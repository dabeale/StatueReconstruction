
#include "Sample.h"

namespace Sample
{
    Generator::Generator()
    {
        uint32_t seed = std::chrono::system_clock::now().time_since_epoch().count();
        m_e = std::default_random_engine(seed);
    }

    std::vector<uint32_t> Generator::RandomPermutation( const uint32_t N )
    {
        std::uniform_int_distribution<uint32_t> uniint(0,N-1);
        std::vector<uint32_t> out(N);
        for(uint32_t i=0;i<N;++i) out[i] = i;

        for(uint32_t i=0;i<N;++i)
            std::swap(out[i], out[uniint(m_e)]);
        return out;
    }

    std::vector<uint32_t> Generator::SampleUniformIntegerUI(const uint32_t min, const uint32_t max, const uint32_t NSamples)
    {
        std::uniform_int_distribution<uint32_t> uniint(min,max);
        std::vector<uint32_t> samples(NSamples);
        for(auto& s : samples)
        {
            s = uniint(m_e);
        }
        return samples;
    }

    std::vector<int32_t> Generator::SampleUniformInteger(const int32_t min, const int32_t max, const uint32_t NSamples)
    {
        std::uniform_int_distribution<int32_t> uniint(min,max);
        std::vector<int32_t> samples(NSamples);
        for(auto& s : samples)
        {
            s = uniint(m_e);
        }
        return samples;
    }

    std::vector<uint32_t> Generator::SampleDiscrete1D( const std::vector<double>& dist, const uint32_t N )
    {
        std::discrete_distribution<int> discrete(dist.begin(), dist.end());
        std::vector<uint32_t> samples(N);

        for( auto& sample : samples )
        {
            sample = discrete(m_e);
        }
        return samples;
    }

    std::vector<uint32_t> Generator::SampleDiscrete1D( const Cu::Matrix& distribution, const uint32_t N )
    {
        std::vector<float> dist(distribution.cols());
        uint32_t index=0;
        for( auto& d : dist )
        {
            d = distribution(0, index++);
        }

        std::discrete_distribution<int> discrete(dist.begin(), dist.end());
        std::vector<uint32_t> samples(N);

        for( auto& sample : samples )
        {
            sample = discrete(m_e);
        }
        return samples;
    }


    std::vector<std::pair<uint32_t, uint32_t>> Generator::SampleDiscrete2D( const Cu::Matrix& distribution, const uint32_t N )
    {
        // Sample the prior
        Cu::Matrix prior = distribution.colwise().sum();
        std::vector<uint32_t> pss = SampleDiscrete1D( prior, N );

        std::vector<std::pair<uint32_t, uint32_t>> returnVals(N);
        auto it = returnVals.begin();
        for( auto p : pss )
        {
            int x = *(SampleDiscrete1D( distribution.col( p ).transpose(), 1).begin());
            int y = *(SampleDiscrete1D( distribution.row( x ), 1).begin());
            it->first = x;
            it->second = y;
            it++;
        }
        return returnVals;
    }

    std::vector<Cu::Matrix> Generator::SampleWishart( const Cu::Matrix& Sigma, const uint32_t df, const uint32_t N )
    {
        std::vector<Cu::Matrix> output(N, Cu::Matrix(Sigma.rows(), Sigma.cols()));
        Cu::Matrix Mu(Sigma.rows(),1), nsamples(Sigma.rows(), df);
        Mu.setZero();

        for( uint32_t k=0; k<N; ++k )
        {
            nsamples = SampleMultivariateNormal( Mu, Sigma, df );
            output[k] = nsamples * nsamples.transpose();
        }

        return output;
    }

    std::vector< Cu::Matrix > Generator::Sample3DEdgeResampling( const std::vector<std::vector< std::pair<double, double> > >& edgeMeans,
                                                                 const std::vector< std::vector<double> >& edgeDeviations,
                                                                 const std::vector<Cu::Matrix>& P,
                                                                 const uint32_t K,
                                                                 const double deviation)
    {
        std::vector<Cu::Matrix> retvals(K, Cu::Matrix(3,1));

        const uint32_t NImages = P.size();
        std::uniform_int_distribution<int32_t> uniint(0,NImages-1);

        const uint32_t IterationCount = 0; // The number of gibbs iterations (this can be low).

        // First sample each curve
        std::vector<std::vector< Cu::Matrix > >
                curveSamples(NImages, std::vector< Cu::Matrix >(K, Cu::Matrix(3,1))),
                correspondences(NImages, std::vector< Cu::Matrix >(K, Cu::Matrix(3,1)));

        for(uint32_t tauval,s,e=0; e<NImages; ++e)
        {
            std::uniform_int_distribution<int32_t> pint(0,edgeMeans[e].size()-1);
            for(s=0; s<K; ++s)
            {
                tauval = pint(m_e);
                std::normal_distribution<double> ndist(0.0, edgeDeviations[e][tauval]);
                curveSamples[e][s](0)  = ndist(m_e) + edgeMeans[e][tauval].first;
                curveSamples[e][s](1) = ndist(m_e) + edgeMeans[e][tauval].second;
                curveSamples[e][s](2) = 1;
            }
        }

        // Get the samples from the first image
        uint32_t sim1, sim2;

        sim1 = uniint(m_e); // get a start point
        for(uint32_t s=0; s<K; ++s)
        {
            correspondences[sim1][s](0) = curveSamples[sim1][s](0);
            correspondences[sim1][s](1) = curveSamples[sim1][s](1);
            correspondences[sim1][s](2) = 1;
        }

        Cu::Matrix cv(1,1);
        std::vector<double> probs(K);
        std::vector<Cu::Matrix> Ps(2);
        Cu::Matrix F(3,3);

        for( sim2=0; sim2<NImages; ++sim2 )
        {
            if(sim2==sim1) continue;
            Ps[0] = P[sim1]; Ps[1] = P[sim2];
            F = Triangulate::FFromP( P[sim1], P[sim2] );

            for(uint32_t ind,s,t=0; t<K;++t)
            {
                for(auto& p : probs) p=1.0; // Clear the probabilities
                for(s=0; s<K; ++s)
                {
                    cv = curveSamples[sim1][t].transpose()*F*curveSamples[sim2][s];
                    probs[s] = Funcs::normal_pdf( cv(0), 0.0, deviation );
                }
                double sum=0.0;
                for(auto p : probs) sum+=p;
                for(auto& p : probs) p/=sum;

                ind = SampleDiscrete1D( probs, 1)[0];
                correspondences[sim1][t] = curveSamples[sim2][ind]; // Some people sample some Gaussian noise here.
            }
        }

        /*
         * A gibbs loop to refine the samples (this is optional)
         */
        for(uint32_t its=0; its<IterationCount; ++its)
        {
            sim1 = uniint(m_e);

            for(uint32_t ind,s,t=0; t<K;++t)
            {
                for(auto& p : probs) p=1.0; // Clear the probabilities
                for( sim2=0; sim2 < NImages; ++sim2 )
                {
                    if(sim2==sim1) continue;
                    Ps[0] = P[sim1]; Ps[1] = P[sim2];
                    F = Triangulate::FFromP( P[sim1], P[sim2] ); // This is quite an intensive computation ( perhaps precompute )
                    for(s=0; s<K; ++s)
                    {
                        cv = correspondences[sim1][t].transpose()*F*correspondences[sim2][s];
                        probs[s] *= Funcs::normal_pdf( cv(0), 0.0, deviation );
                    }
                }

                double sum=0.0;
                for(auto p : probs) sum+=p;
                for(auto& p : probs) p/=sum;

                ind = SampleDiscrete1D( probs, 1)[0];
                //correspondences[sim1][t] = curveSamples[sim2][ind]; // Some people sample some Gaussian noise here.
            }
        }

        std::vector<Cu::Matrix> Xs(NImages);
        for(uint32_t t,s=0; s<K; ++s)
        {
            for(t=0; t<NImages; ++t)
            {
                Xs[t] = correspondences[t][s];
            }
            retvals[s] = Triangulate::Triangulate(P, Xs);
        }

        return retvals;
    }

    std::vector<std::vector<double>> Generator::Sample3DEdge(const std::vector<Cu::Matrix>& edges,
                                                  const std::vector<Cu::Matrix>& P,
                                                  const uint32_t K , const AlternateCameraMethod acm)
    {
        std::vector<std::vector<double>> retvals(K, std::vector<double>(3,0.0));
        std::uniform_int_distribution<int32_t> uniint(0,P.size()-1);
        uint32_t sim,sim2; // the start edge
        Cu::Matrix F(3,3);
        Cu::Matrix EpiEdge(edges[0].rows(), edges[0].cols());

        std::vector<Cu::Matrix> NormalEdges(edges.size());
        auto eit=edges.begin();
        for(auto& e : NormalEdges)
        {
            e = *(eit++);
            e /= e.sum();
        }
        Cu::Matrix x1(3,1),x2(3,1),x3(3,1);
        x1.setOnes(3,1); x2.setOnes(3,1);

        Cu::Matrix prod(EpiEdge.rows(), EpiEdge.cols());
        std::vector<Cu::Matrix> Ps(2);
        std::vector<Cu::Matrix> Xs(2);

        std::vector<Cu::Matrix> outwardNormals;
        if(acm == Orthogonal)
        {
            outwardNormals = ComputeOutwardNormals(P);
        }

        double o1,o2;
        for(uint32_t i,j,k=0; k<K; ++k)
        {
            sim = uniint(m_e);
            sim2 = sim;

            switch(acm)
            {
            case Random:
                while(sim2 == sim) sim2 = uniint(m_e); // Choose two random cameras
                break;
            case Orthogonal:
                double dot = outwardNormals[sim].cwiseProduct(outwardNormals[sim2]).sum();
                for(uint32_t l=0; l<P.size(); ++l)
                {
                    double dd= outwardNormals[sim].cwiseProduct(outwardNormals[l]).sum();
                    if(dd < dot)
                    {
                        dot = dd;
                        sim2 = l;
                    }
                }
            }

            Ps[0] = P[sim]; Ps[1] = P[sim2];

            const Cu::Matrix& Edge1 = NormalEdges[sim];
            const Cu::Matrix& Edge2 = NormalEdges[sim2];

            F = Triangulate::FFromP( P[sim], P[sim2] );
            auto x = SampleDiscrete2D( Edge1, 1);

            //Cu::print_matrix_octave(Edge1, "Edge1");
            o1 = x1(0) = x[0].first;
            o2 = x1(1) = x[0].second;
            //Cu::print_matrix_octave(x1, "x1");

            EpiEdge.setZero(EpiEdge.rows(), EpiEdge.cols());

            x3 = (x1.transpose()*F).transpose();
            for(i=0;i<EpiEdge.rows();++i)
            {
                j = static_cast<uint32_t>(std::round((-x3(0)*i - x3(2))/x3(1)));
                if(j < EpiEdge.cols())
                {
                    EpiEdge(i,j)= 1.0;
                }
            }

            //Cu::print_matrix_octave(Edge2, "Edge2");
            //Cu::print_matrix_octave(EpiEdge, "EpiEdge");
/*
            for(i=0;i<EpiEdge.rows();++i)
            {
                for(j=0;j<EpiEdge.cols();++j)
                {
                    ferr = F(0,0)*o1*i + F(1,0)*o2*i +
                           F(0,1)*o1*j + F(1,1)*o2*j +
                           F(0,2)*o1 + F(2,0)*i +
                           F(1,2)*o2 + F(2,1)*j + F(2,2);
                    if(ferr < 1e-1)
                    {
                        EpiEdge(i,j) = Funcs::normal_pdf( ferr, 0.0, 2.0 );
                    }
                }
            }*/

            EpiEdge /= EpiEdge.sum();

            prod = EpiEdge.cwiseProduct(Edge2);
            double psum =  prod.sum();
            if(psum < 1e-8)
            {
                --k;
            }
            else
            {
                prod /= psum;

                //Cu::print_matrix_octave(prod, "prod");

                x = SampleDiscrete2D( prod, 1);
                x2(0) = x[0].first;
                x2(1) = x[0].second;
                //Cu::print_matrix_octave(x2, "x2");
                Xs[0] = x1; Xs[1] = x2;

                auto X = Triangulate::Triangulate(Ps, Xs);

                for(j=0;j<3;++j)
                {
                    retvals[k][j] = X(j);
                }
            }
        }

        return retvals;
    }

    Cu::Matrix Generator::SampleMultivariateNormal( const Cu::Matrix& mu, const Cu::Matrix& Sigma, const uint32_t N )
    {
        // Do cholesky decomposition
        Eigen::LLT<Cu::Matrix> lltofSigma(Sigma);
        Cu::Matrix llt = lltofSigma.matrixL();
        //Cu::print_matrix(llt, "llt");

        // Find the eigenvalues
        Cu::Matrix eigs = (llt.cwiseProduct(llt)).colwise().sum().cwiseSqrt();
        //Cu::print_matrix(eigs, "eigs");

        // colwise divide
        for(uint32_t i=0; i<llt.rows(); i++)
        {
            llt.row(i) = llt.row(i).cwiseQuotient( eigs.row(0) );
        }

        //Cu::print_matrix(llt, "llt scaled");

        Cu::Matrix output( mu.rows(), N);

        // Sample the normal
        for( uint32_t i=0; i<mu.rows(); i++)
        {
            std::normal_distribution<float> distribution ( 0.0, eigs(0,i));
            for( uint32_t k=0; k<N; k++ )
            {
                output(i,k) = distribution(m_e);
            }
        }

        output = (llt*output);
        output.colwise() += (Cu::Vector)mu.col(0);
        //print_matrix(output, "data");

        return output;
    }
}
