#include "SegmentationModel.h"

inline void NormaliseField(Math::Matrix& fprob, Math::Matrix& bprob)
{
    for(uint32_t i=0; i<fprob.Cols(); ++i)
    for(uint32_t j=0; j<fprob.Rows(); ++j)
    {
        double sum = bprob(j,i) + fprob(j,i);
        if( sum > 0 )
        {
            bprob(j,i) /= sum;
            fprob(j,i) /= sum;
        }
        else
        {
            bprob(j,i)=0.5;
            fprob(j,i)=0.5;
        }
    }
}

SegmentationModel::SegmentationModel(std::vector< ImageInfo >& iminfos, const Mesh& mesh,
                                     const std::vector<Math::Matrix>& ProjectionMatrices) :
    m_imageinfos(iminfos),
    m_3dmesh(mesh),
    m_ProjectionMatrices(ProjectionMatrices),
    m_ForegroundColours(3,0),
    m_BackgroundColours(3,0),
    m_pointProbability(mesh.GetNVerts()),
    m_NSamples(10000),
    m_msg("SegmentationModel"),
    m_depthModelsComputed(false),
    m_fprob(iminfos.size()),
    m_bprob(iminfos.size())
{

}

void SegmentationModel::CreateData( )
{
    m_msg.print() << "Create data from sketches" << std::endl;

    std::vector<double> fgarr, bgarr;

    for( const ImageInfo& im : m_imageinfos )
    {
        if(im.m_ForegroundColours.size() > 0)
            fgarr.insert( fgarr.end(), im.m_ForegroundColours.begin(), im.m_ForegroundColours.end() );

        if(im.m_BackgroundColours.size() > 0)
            bgarr.insert( bgarr.end(), im.m_BackgroundColours.begin(), im.m_BackgroundColours.end() );
    }

    Math::Matrix fgcols( 3, fgarr.size()/3, fgarr, false );
    Math::Matrix bgcols( 3, bgarr.size()/3, bgarr, false );

    m_msg.print() << "Sampling colours" << std::endl;

    Math::Generator g;
    m_ForegroundColours = g.UniquelySubSampleColumns( fgcols, m_NSamples );
    m_BackgroundColours = g.UniquelySubSampleColumns( bgcols, m_NSamples );
}

void SegmentationModel::RunSegmentation(double variance, double alpha)
{
    if(m_ForegroundColours.Cols() == 0)
    {
        m_msg.printerr("SegmentationModel::RunSegmentation:: There are no foreground colours");
        return;
    }

    if(m_BackgroundColours.Cols() ==0 )
    {
        m_msg.printerr("SegmentationModel::RunSegmentation:: There are no background colours");
        return;
    }

    m_msg.print() << "Clustering the colour data" << std::endl;
    m_msg.print(1) << "Foreground" << std::endl;
    Math::Matrix Lambda = m_ForegroundColours.mean(1);
    Math::Matrix R = (variance*m_ForegroundColours*m_ForegroundColours.transpose()).inv();
    uint32_t Beta = m_ForegroundColours.Rows() + 2;
    Sample::IGMMHyperparameters hypesfg = {Lambda, R, alpha, Beta};
    Sample::IGMM igmmfg( m_ForegroundColours, hypesfg );
    igmmfg.Estimate(1e-13, 100);
    igmmfg.ConditionBestSigmas();

    m_msg.print(1) << "Background" << std::endl;
    Lambda = m_BackgroundColours.mean(1);
    R = (variance*m_BackgroundColours*m_BackgroundColours.transpose()).inv();
    Beta = m_BackgroundColours.Rows() + 2;
    Sample::IGMMHyperparameters hypesbg = {Lambda, R, alpha, Beta};
    Sample::IGMM igmmbg( m_BackgroundColours, hypesbg );
    igmmbg.Estimate(1e-13, 100);
    igmmbg.ConditionBestSigmas();


    m_msg.print() << "Running segmentations" << std::endl;
    m_fprob.resize(m_imageinfos.size());
    m_bprob.resize(m_imageinfos.size());
    auto fpit = m_fprob.begin();
    auto bpit = m_bprob.begin();
    uint32_t iindex=0;
    for(auto& iminfo : m_imageinfos)
    {
        const uint32_t M = iminfo.m_image.height();
        const uint32_t N = iminfo.m_image.width();

        *fpit = Math::Matrix(M, N);
        *bpit = Math::Matrix(M, N);
        m_msg.print() << "Camera : " << iindex << std::endl;
        m_msg.print(1) << "Evaluating IGMM" << std::endl;
        Math::Matrix imvec(3,M*N);
        uint32_t index=0;
        for(uint32_t i=0; i<N; ++i)
        for(uint32_t j=0; j<M; ++j)
        {
            QColor q(iminfo.m_image.pixel(i,j));
            imvec(0,index) = q.redF();
            imvec(1,index) = q.greenF();
            imvec(2,index) = q.blueF();
            ++index;
        }
        fpit->GetArr() = igmmfg.EvaluateLikelihood( imvec );
        bpit->GetArr() = igmmbg.EvaluateLikelihood( imvec );

        NormaliseField(*fpit,*bpit);

        //fpit->Save("ForegroundIgmm");
        //bpit->Save("BackgroundIgmm");

        if(m_depthModelsComputed)
        {
            if(m_depthModel.find(iindex) != m_depthModel.end())
            {
                auto& dm = m_depthModel[iindex];
                m_msg.print(1) << "Applying known depth models" << std::endl;
                //m_msg.printerr("This model is wrong for multiple views!");
                for(uint32_t i=0; i<N; ++i)
                for(uint32_t j=0; j<M; ++j)
                {
                    double depth = iminfo.m_depthMapValues[ iminfo.getFileHeight() * ( (i*iminfo.getFileWidth()) / iminfo.m_image.width()) + j*iminfo.getFileHeight() / iminfo.m_image.height() ];
                    if(depth > 0)
                    {
                        (*fpit)(j,i) *= std::exp(- (depth-dm.m_foregroundDepth[0])*(depth-dm.m_foregroundDepth[0]) / (2*dm.m_foregroundDepth[1]*dm.m_foregroundDepth[1])) / dm.m_foregroundDepth[1];
                        (*bpit)(j,i) *= std::exp(- (depth-dm.m_backgroundDepth[0])*(depth-dm.m_backgroundDepth[0]) / (2*dm.m_backgroundDepth[1]*dm.m_backgroundDepth[1])) / dm.m_backgroundDepth[1];
                        (*fpit)(j,i) *= dm.m_foregroundUndefinedProbability[0];
                        (*bpit)(j,i) *= dm.m_backgroundUndefinedProbability[0];
                    }
                    else
                    {
                        (*fpit)(j,i) *= dm.m_foregroundUndefinedProbability[1];
                        (*bpit)(j,i) *= dm.m_backgroundUndefinedProbability[1];
                    }
                }
            }
            else
            {
                m_msg.print(1) << "Applying infered depth models" << std::endl;
                //Math::Matrix dm(M,N);
                //dm.SetZero();
                for(uint32_t i=0; i<N; ++i)
                for(uint32_t j=0; j<M; ++j)
                {
                    uint32_t imindex = iminfo.getFileHeight() * ( (i*iminfo.getFileWidth()) / iminfo.m_image.width()) + j*iminfo.getFileHeight() / iminfo.m_image.height();
                    int32_t vertindex = iminfo.m_indexValues[ imindex ];
                    if(vertindex > 0)
                    {
                        auto& fgbg = m_pointProbability[vertindex];
                        (*fpit)(j,i) *= fgbg[0] / ( fgbg[0] + fgbg[1] );
                        (*bpit)(j,i) *= fgbg[1] / ( fgbg[0] + fgbg[1] );
                        (*fpit)(j,i) *= m_fgprobabilityUndefined[0];
                        (*bpit)(j,i) *= m_bgprobabilityUndefined[0];
                        //dm(j,i) = 1;
                    }
                    else
                    {
                        (*fpit)(j,i) *= m_fgprobabilityUndefined[1];
                        (*bpit)(j,i) *= m_bgprobabilityUndefined[1];
                    }
                }
                //dm.Save("dm");

            }
            NormaliseField(*fpit,*bpit);

            //fpit->Save("ForegroundDepth");
            //bpit->Save("BackgroundDepth");
        }
        else
        {
            m_msg.print() << "No depth models applied" << std::endl;
        }

        ++fpit;
        ++bpit;
        ++iindex;
    }

    m_msg.print() << "Evaluating MRFs" << std::endl;
    EvaluateMRF( m_msg);
}

void SegmentationModel::ComputePointProbabilities()
{
    m_msg.print() << "Computing point probabilities" << std::endl;
    auto& verts = m_3dmesh.GetVertices();
    m_pointProbability.resize(m_3dmesh.GetNVerts());
    for(uint32_t i=0; i<m_3dmesh.GetNVerts(); ++i)
    {
        m_pointProbability[i][0] = 0.0;
        m_pointProbability[i][1] = 0.0;
        uint32_t count = 0;

        for(uint32_t j=0; j<m_ProjectionMatrices.size(); ++j)
        {
            if(m_depthModel.find(j) != m_depthModel.end())
            {
                auto& P = m_ProjectionMatrices[j];
                auto& dm = m_depthModel[j];
                double depth = P(2,0)*verts[3*i] + P(2,1)*verts[3*i+1] + P(2,2)*verts[3*i+2] + P(2,3);
                m_pointProbability[i][0] += std::exp(- (depth-dm.m_foregroundDepth[0])*(depth-dm.m_foregroundDepth[0]) / (2*dm.m_foregroundDepth[1]*dm.m_foregroundDepth[1])) / dm.m_foregroundDepth[1];
                m_pointProbability[i][1] += std::exp(- (depth-dm.m_backgroundDepth[0])*(depth-dm.m_backgroundDepth[0]) / (2*dm.m_backgroundDepth[1]*dm.m_backgroundDepth[1])) / dm.m_backgroundDepth[1];

                ++count;
            }
        }

        m_pointProbability[i][0] /= count;
        m_pointProbability[i][1] /= count;
    }
}

void SegmentationModel::ComputeDepthModels()
{
    m_msg.print() << "Computing depth models" << std::endl;

    uint32_t numberofdepthmodels=0;

    m_fgprobabilityUndefined = {0,0};
    m_bgprobabilityUndefined = {0,0};

    for( uint32_t i=0; i<m_imageinfos.size(); ++i )
    {
        auto& iminfo = m_imageinfos[i];
        if( iminfo.m_BackgroundDepths.size() > 0 && iminfo.m_ForegroundDepths.size() > 0 )
        {
            auto& dm = m_depthModel[i];
            m_msg.print(1) << "Camera " << i << std::endl;
            m_msg.print(1) << "Computing background models" << std::endl;
            {
                dm.m_backgroundDepth[0] = 0.0; // The mean
                uint32_t count=0;
                for( auto d : iminfo.m_BackgroundDepths)
                {
                    if(d > 0)
                    {
                        dm.m_backgroundDepth[0] += d;
                        ++count;
                        ++dm.m_backgroundUndefinedProbability[0];
                    }
                    else
                    {
                        ++dm.m_backgroundUndefinedProbability[1];
                    }
                }
                dm.m_backgroundDepth[0] /= count;
                dm.m_backgroundDepth[1] = 0.0;
                for( auto d : iminfo.m_BackgroundDepths)
                {
                    if(d>0)
                    {
                        dm.m_backgroundDepth[1]+=(d - dm.m_backgroundDepth[0])*(d - dm.m_backgroundDepth[0]);
                    }
                }
                dm.m_backgroundDepth[1]/=count;
                double sum = dm.m_backgroundUndefinedProbability[0] + dm.m_backgroundUndefinedProbability[1];
                dm.m_backgroundUndefinedProbability[0] /= sum;
                dm.m_backgroundUndefinedProbability[1] /= sum;

                m_bgprobabilityUndefined[0] += dm.m_backgroundUndefinedProbability[0];
                m_bgprobabilityUndefined[1] += dm.m_backgroundUndefinedProbability[1];

                m_msg.print(1) << "Mean: " << dm.m_backgroundDepth[0] << ", " <<
                                  "Var: " <<  dm.m_backgroundDepth[1] << "\n" <<
                                  "Prob def:" << dm.m_backgroundUndefinedProbability[0] << ", "<<
                                  "Prob undef: " << dm.m_backgroundUndefinedProbability[1] << std::endl;
            }

            m_msg.print(1) << "Computing foreground models" << std::endl;
            {
                dm.m_foregroundDepth[0] = 0.0; // The mean
                uint32_t count=0;
                for( auto d : iminfo.m_ForegroundDepths)
                {
                    if(d > 0)
                    {
                        dm.m_foregroundDepth[0] += d;
                        ++count;
                        ++dm.m_foregroundUndefinedProbability[0];
                    }
                    else
                    {
                        ++dm.m_foregroundUndefinedProbability[1];
                    }
                }
                dm.m_foregroundDepth[0] /= count;
                dm.m_foregroundDepth[1] = 0.0;
                for( auto d : iminfo.m_ForegroundDepths)
                {
                    if(d > 0)
                    {
                        dm.m_foregroundDepth[1]+=(d - dm.m_foregroundDepth[0])*(d - dm.m_backgroundDepth[0]);
                    }
                }
                dm.m_foregroundDepth[1]/=count;
                double sum = dm.m_foregroundUndefinedProbability[0] + dm.m_foregroundUndefinedProbability[1];
                dm.m_foregroundUndefinedProbability[0] /= sum;
                dm.m_foregroundUndefinedProbability[1] /= sum;

                m_fgprobabilityUndefined[0] += dm.m_foregroundUndefinedProbability[0];
                m_fgprobabilityUndefined[1] += dm.m_foregroundUndefinedProbability[1];

                m_msg.print(1) << "Mean: " << dm.m_foregroundDepth[0] << ", " <<
                                  "Var: " <<  dm.m_foregroundDepth[1] << "\n" <<
                                  "Prob def:" << dm.m_foregroundUndefinedProbability[0] << ", "<<
                                  "Prob undef: " << dm.m_foregroundUndefinedProbability[1] << std::endl;
            }

            ++numberofdepthmodels;
            m_depthModelsComputed = true;
        }
    }

    m_fgprobabilityUndefined[0] /= numberofdepthmodels;
    m_fgprobabilityUndefined[1] /= numberofdepthmodels;
    m_bgprobabilityUndefined[0] /= numberofdepthmodels;
    m_bgprobabilityUndefined[1] /= numberofdepthmodels;

    if(!m_depthModelsComputed)
    {
        m_msg.printerr("Warning:: SegmentationModel::ComputeDepthModels :: No depths found in any drawings");
    }
}

inline void EvaluateFieldAtPointList( Math::Matrix& fprob, Math::Matrix& bprob, const std::vector<PenPoint>& points)
{
    const uint32_t M = fprob.Rows();
    const uint32_t N = fprob.Cols();
    for( const PenPoint& pp : points )
    {
        uint32_t radius = pp.m_pen.width()/2;
        uint32_t x = pp.m_point.x();
        uint32_t y = pp.m_point.y();
        for(int32_t i = x-radius; i < x+radius; ++i)
        for(int32_t j = y-radius; j < y+radius; ++j)
        if( i < N && j < M && i > 0 && j > 0)
        if( (i - x)*(i - x) + (j - y)*(j - y) < radius*radius )
        {
            fprob(j,i) = 1 - std::numeric_limits<double>::epsilon();
            bprob(j,i) = std::numeric_limits<double>::epsilon();
        }
    }
}

void SegmentationModel::EvaluateMRF(Stream::Message& msg)
{
    auto fgit = m_fprob.begin();
    auto bgit = m_bprob.begin();
    uint32_t camind = 0;
    for(auto& iminfo : m_imageinfos)
    {
        const uint32_t M = iminfo.m_image.height();
        const uint32_t N = iminfo.m_image.width();
        //Math::Matrix fprob(M, N);
        //Math::Matrix bprob(M, N);

        msg.print(1) << "Creating fields. Camera : " << camind++ << std::endl;
        msg.print(2) << "Adding hard constraints" << std::endl;
        EvaluateFieldAtPointList( *fgit, *bgit, iminfo.m_ForegroundPoints);
        EvaluateFieldAtPointList( *bgit, *fgit, iminfo.m_BackgroundPoints);

        NormaliseField(*fgit, *bgit);

        switch (m_params.m_type) {
        case SegmentationModelParameters::GraphCut :
        {
            msg.print(1) << "Maximising" << std::endl;
            MaxFlow::HyperParameters hypes {
                {m_params.m_Alpha[0], m_params.m_Alpha[1]}, // Alpha the higher the value the more ambiguous it is about the class - choose in relation to N
                m_params.m_NSample,  // N. The number of samples.
                m_params.m_edgeWeight,
                m_params.m_terminalWeight,
                msg // The message interface
            };

            auto assigments = MaxFlow::EvaluateFGBGFromProbabilities(*fgit, *bgit, hypes );

            QColor qwhite( 255,255,255 ), qblack(0,0,0);
            for( uint32_t i=0; i<M; ++i)
            for( uint32_t j=0; j<N; ++j)
            {
                if(assigments(i,j) > 0.5)
                {
                    iminfo.m_segmentation.setPixelColor( QPoint(j,i ), qblack);
                }
                else
                {
                    iminfo.m_segmentation.setPixelColor( QPoint(j,i ), qwhite);
                }
            }
        }
            break;

        case SegmentationModelParameters::ICM :
        {/*
            std::vector<double> probstack;
            probstack.insert(probstack.end(), fgit->GetArr().begin(), fgit->GetArr().end());
            probstack.insert(probstack.end(), bgit->GetArr().begin(), bgit->GetArr().end());

            msg.print(1) << "Maximising MRF" << std::endl;
            auto assigments = MRF::EvaluateFieldFromProbabilityImage( probstack.data(), 2, N, M, MRF::IteratedConditionalModes, 100 );

            QColor qwhite( 255,255,255 ), qblack(0,0,0);
            for( uint32_t i=0; i<M; ++i)
            for( uint32_t j=0; j<N; ++j)
            {
                if(assigments[ M*j + i ] > 0.5)
                {
                    iminfo.m_segmentation.setPixelColor( QPoint(j,i ), qwhite);
                }
                else
                {
                    iminfo.m_segmentation.setPixelColor( QPoint(j,i ), qblack);
                }
            }*/
            msg.printerr("SegmentationModel::EvaluateMRF:: ICM is not supported. It has been removed because it is slow.");
        }
            break;
        default:
            msg.printerr("SegmentationModel::EvaluateMRF:: Option not supported");
            break;
        }

        ++fgit;
        ++bgit;
    }

    msg.print() << "Done." << std::endl;
}
