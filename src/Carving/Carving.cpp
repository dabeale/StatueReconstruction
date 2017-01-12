#include "Carving.h"

namespace Carve
{
Box GetBoxFromCloud( const std::vector< Cu::Matrix >& P,
                     const std::vector< Image >& Images,
                     Stream::Message msg,
                     double fractionPadding,
                     const std::vector<double >&Samples )
{

    const uint32_t rows = Images[0].rows();
    const uint32_t cols = Images[0].cols();
    const uint32_t N = Samples.size()/3;
    std::vector< const double* > Samples2;
    Samples2.reserve(N);

    for( uint32_t l=0; l< N; ++l)
    {
        const double* s = &Samples[3*l];

        double ave=0.0;
        uint32_t Nin=0;
        for(uint32_t k=0; k<P.size(); ++k)
        {
            double x = P[k](0,0)*s[0] + P[k](0,1)*s[1] + P[k](0,2)*s[2] + P[k](0,3);
            double y = P[k](1,0)*s[0] + P[k](1,1)*s[1] + P[k](1,2)*s[2] + P[k](1,3);
            double scale = P[k](2,0)*s[0] + P[k](2,1)*s[1] + P[k](2,2)*s[2] + P[k](2,3);
            y /= scale;
            x /= scale;
            if(static_cast<uint32_t>(x) < cols &&
               static_cast<uint32_t>(y) < rows )
            {
                ave += static_cast<double>(Images[k].Ptr()[ 3*(rows*static_cast<uint32_t>(x) + static_cast<uint32_t>(y)) ]);
                ++Nin;
            }
        }
        if(Nin > 2 && (ave / Nin) - 0.9 > 0)
        {
            Samples2.push_back( s );
        }
    }

    /*
    Mesh msh(Samples2);
    msh.Write( "/Users/dabeale/Workspace/Libs/ELib/scripts/Carving/Test.ply" );
    */

    msg.print(2) << "Finalising data" << std::endl;

    Box ret;

    for(const auto& s : Samples2)
    {
        for( uint8_t k=0; k<3; ++k )
        {
            ret.center[k]+=s[k];
            if(s[k] > ret.max[k]) ret.max[k] = s[k];
            if(s[k] < ret.min[k]) ret.min[k] = s[k];
        }
    }

    for( uint8_t k=0; k<3; ++k )
    {
        ret.center[k] /= Samples2.size();
        ret.max[k] += fractionPadding*(ret.max[k] - ret.min[k]);
        ret.min[k] -= fractionPadding*(ret.max[k] - ret.min[k]); // Add some padding
    }

    return ret;
}

    Box GetBoxFromCloud( const std::vector< Cu::Matrix >& P,
                         const std::vector< Image >& Images,
                         Stream::Message msg,
                         double fractionPadding,
                         const std::vector<std::vector<double> >&Samples )
    {

        const uint32_t rows = Images[0].rows();
        const uint32_t cols = Images[0].cols();

        std::vector< std::vector<double> > Samples2;
        Samples2.reserve(Samples.size());

        for( uint32_t l=0; l< Samples.size(); ++l)
        {
            const auto& s = Samples[l];

            double ave=0.0;
            uint32_t Nin=0;
            for(uint32_t k=0; k<P.size(); ++k)
            {
                double x = P[k](0,0)*s[0] + P[k](0,1)*s[1] + P[k](0,2)*s[2] + P[k](0,3);
                double y = P[k](1,0)*s[0] + P[k](1,1)*s[1] + P[k](1,2)*s[2] + P[k](1,3);
                double scale = P[k](2,0)*s[0] + P[k](2,1)*s[1] + P[k](2,2)*s[2] + P[k](2,3);
                y /= scale;
                x /= scale;
                if(static_cast<uint32_t>(x) < cols &&
                   static_cast<uint32_t>(y) < rows )
                {
                    ave += static_cast<double>(Images[k].Ptr()[ 3*(rows*static_cast<uint32_t>(x) + static_cast<uint32_t>(y)) ]);
                    ++Nin;
                }
            }
            if(Nin > 2 && (ave / Nin) - 0.9 > 0)
            {
                Samples2.push_back( Samples[l] );
            }
        }

        /*
        Mesh msh(Samples2);
        msh.Write( "/Users/dabeale/Workspace/Libs/ELib/scripts/Carving/Test.ply" );
        */

        msg.print(2) << "Finalising data" << std::endl;

        Box ret;

        for(const auto& s : Samples2)
        {
            for( uint8_t k=0; k<3; ++k )
            {
                ret.center[k]+=s[k];
                if(s[k] > ret.max[k]) ret.max[k] = s[k];
                if(s[k] < ret.min[k]) ret.min[k] = s[k];
            }
        }

        for( uint8_t k=0; k<3; ++k )
        {
            ret.center[k] /= Samples2.size();
            ret.max[k] += fractionPadding*(ret.max[k] - ret.min[k]);
            ret.min[k] -= fractionPadding*(ret.max[k] - ret.min[k]); // Add some padding
        }

        return ret;
    }

    Mesh FindMesh(const std::vector< Cu::Matrix >& P, const std::vector< Image >& Images,
                   const Box& box, Stream::Message msg, const uint32_t density)
    {

        const uint32_t N = P.size();
        const uint32_t rows = Images[0].rows();
        const uint32_t cols = Images[0].cols();

        auto outwardNormals = ComputeOutwardNormals( P );

        auto func = [=](const double* val)
        {
            uint32_t i,j;
            double x[3];
            double ave = 0.0;
            uint32_t Nin=0;

            bool angleCriteria = false; // true of there is an angle of roughly 90 degrees
            uint32_t firstincam =0;
            for(uint32_t k=0; k<N; ++k)
            {
                x[0] = P[k](0,0)*val[0] + P[k](0,1)*val[1] + P[k](0,2)*val[2] + P[k](0,3);
                x[1] = P[k](1,0)*val[0] + P[k](1,1)*val[1] + P[k](1,2)*val[2] + P[k](1,3);
                x[2] = P[k](2,0)*val[0] + P[k](2,1)*val[1] + P[k](2,2)*val[2] + P[k](2,3);
                x[0] /= x[2];
                x[1] /= x[2];

                if(x[0] < cols && x[0] >= 0 && x[1] < rows && x[1] >= 0)
                {
                    i = static_cast<uint32_t> (std::floor(x[0]));
                    j = static_cast<uint32_t> (std::floor(x[1]));
                    ave += static_cast<double>(Images[k].Ptr()[3*(rows*i + j)]);

                    if(Nin > 0)
                    {
                        double angle = std::acos(outwardNormals[k](0)*outwardNormals[firstincam](0) +
                                       outwardNormals[k](1)*outwardNormals[firstincam](1) +
                                       outwardNormals[k](2)*outwardNormals[firstincam](2));
                        if( std::abs(angle - M_PI/2) < M_PI / 8)
                        {
                            angleCriteria = true;
                        }
                    }
                    else
                    {
                        firstincam = k;
                    }

                    ++Nin;
                }
            }

            return (angleCriteria && Nin > P.size()-3) ? (ave / Nin)-0.9 : -0.9;
        };

        /**
         * Create a pair of blobs and save to ply file
         */
        const double iso = 0.01;

        msg.print(2) << "Running marching cubes" << std::endl;
        double max[3], min[3];
        for(uint8_t k=0; k<3;++k)
        {
            max[k] = box.max[k];
            min[k] = box.min[k];
        }
        MC::MarchingCubes mc(iso, density, max, min);
        mc.March(func);

        auto v = mc.GetPoints();
        auto f = mc.GetFaces();

        std::vector< PointTypes::Pointd3 > vo(v.size());
        std::vector< PointTypes::Pointi3 > fo(f.size());

#pragma omp parallel for
        for(uint32_t i=0; i<v.size(); ++i)
        for(uint8_t j=0; j<3; ++j)
            vo[i].p[j] = v[i].p[j];

#pragma omp parallel for
        for(uint32_t i=0; i<f.size(); ++i)
        for(uint8_t j=0; j<3; ++j)
            fo[i].p[j] = f[i].p[j];

        Mesh ret(vo, fo);

        return ret;
    }

    void RemovePointsOutsideHull( Mesh& mesh,
                                 const std::vector< Cu::Matrix >& P, const std::vector< Image >& Images,
                                 Stream::Message msg)
    {
        const uint32_t N = P.size();
        const uint32_t rows = Images[0].rows();
        const uint32_t cols = Images[0].cols();

        std::vector<double> newVertices(0);
        std::vector<double> newColours(0);
        std::vector<double> newNormals(0);
        newVertices.reserve(mesh.GetVertices().size());
        newColours.reserve(mesh.GetVertices().size());
        newNormals.reserve(mesh.GetVertices().size());

        const double* dat = mesh.GetVertices().data();
        const double* datn = mesh.GetNormals().data();
        const double* datc = mesh.GetColours().data();

        for(uint32_t step=0; step<mesh.GetNVerts(); ++step)
        {
            uint32_t i,j;
            double x[3];
            double ave = 0.0;
            uint32_t Nin=0;
            for(uint32_t k=0; k<N; ++k)
            {
                x[0] = P[k](0,0)*dat[step*3] + P[k](0,1)*dat[step*3 + 1] + P[k](0,2)*dat[step*3 + 2] + P[k](0,3);
                x[1] = P[k](1,0)*dat[step*3] + P[k](1,1)*dat[step*3 + 1] + P[k](1,2)*dat[step*3 + 2] + P[k](1,3);
                x[2] = P[k](2,0)*dat[step*3] + P[k](2,1)*dat[step*3 + 1] + P[k](2,2)*dat[step*3 + 2] + P[k](2,3);
                x[0] /= x[2];
                x[1] /= x[2];

                if(x[0] < cols && x[0] >= 0 && x[1] < rows && x[1] >= 0)
                {
                    i = static_cast<uint32_t> (std::floor(x[0]));
                    j = static_cast<uint32_t> (std::floor(x[1]));
                    ave += static_cast<double>(Images[k].Ptr()[3*(rows*i + j)]);
                    ++Nin;
                }
            }


            if( Nin > 2  && (ave / Nin)-0.9 > 0)
            {
                newVertices.push_back( dat[step*3]  );
                newVertices.push_back( dat[step*3 + 1]  );
                newVertices.push_back( dat[step*3 + 2]  );

                if(mesh.GetColours().size() > 0)
                {
                    newColours.push_back(datc[step*3]);
                    newColours.push_back(datc[step*3 + 1]);
                    newColours.push_back(datc[step*3 + 2]);
                }

                if(mesh.GetNormals().size() > 0)
                {
                    newNormals.push_back(datn[step*3]);
                    newNormals.push_back(datn[step*3 + 1]);
                    newNormals.push_back(datn[step*3 + 2]);
                }
            }
        }
        mesh.GetColours() = newColours;
        mesh.GetNormals() = newNormals;
        mesh.GetVertices() = newVertices;
        mesh.SetNVertices( newVertices.size()/3 );
    }

    Mesh Blend( const Mesh& pointcloud, const NearestNeighbours* nnPointCloud,
                const Mesh& visualhull, const NearestNeighbours* nnVHull,
                const Box& box,
                const uint32_t K, const double radius, Stream::Message msg)
    {
        auto func = [=](const double* val)
        {
            auto PCNs = nnPointCloud->Search( val, K );
            auto VHNs = nnVHull->Search( val, K );

            double pcnave = FindNormalisedProbability( PCNs, radius);
            double vhave = FindNormalisedProbability( VHNs, radius);

            double pb, nb, pcloudassigment=0.0, vhullassignment=0.0;
            for(uint32_t k=0; k<K; ++k)
            {
                auto& kppc = PCNs[k];
                auto& vhpc = VHNs[k];

                // Find the visual hull implicit (hoppe) function
                pb = 0.0;
                nb = 0.0;
                for(uint32_t j=0; j<3; ++j)
                {
                    pb += (pointcloud.GetVertices()[3*kppc.index + j] - val[j])*
                            pointcloud.GetNormals()[3*kppc.index + j];
                    nb += (visualhull.GetVertices()[3*vhpc.index + j] - val[j])*
                            visualhull.GetNormals()[3*vhpc.index + j];
                }

                pcloudassigment += pb*PCNs[k].distance;
                vhullassignment += nb*VHNs[k].distance;
            }

            if( pcloudassigment <  vhullassignment)
            {
                pcnave = pdf( pcnave, radius);
                vhave = pdf( vhave, radius);
                double mx = std::max(pcnave, vhave);
                pcnave = std::exp(pcnave - mx);
                vhave = std::exp(vhave - mx);
                pcnave /= pcnave + vhave;
                vhave /= vhave + pcnave;

                return pcnave*pcloudassigment + vhave*vhullassignment;
            }
            else
            {
                return vhullassignment;
            }
        };

        const double iso = 0.01;
        const uint32_t nd = 100;

        msg.print(2) << "Running marching cubes" << std::endl;
        double max[3], min[3];
        for(uint8_t k=0; k<3;++k)
        {
            max[k] = box.max[k];
            min[k] = box.min[k];
        }

        MC::MarchingCubes mc(iso, nd, max, min);
        mc.March(func);

        auto v = mc.GetPoints();
        auto f = mc.GetFaces();

        std::vector< PointTypes::Pointd3 > vo(v.size());
        std::vector< PointTypes::Pointi3 > fo(f.size());

#pragma omp parallel for
        for(uint32_t i=0; i<v.size(); ++i)
        for(uint8_t j=0; j<3; ++j)
            vo[i].p[j] = v[i].p[j];

#pragma omp parallel for
        for(uint32_t i=0; i<f.size(); ++i)
        for(uint8_t j=0; j<3; ++j)
            fo[i].p[j] = f[i].p[j];

        Mesh ret(vo, fo);

        return ret;
    }

    Mesh Hoppe(const Mesh& pointcloud, const NearestNeighbours *nnPointCloud, const Box &box,
               const uint32_t K, const double radius, Stream::Message msg)
    {
        auto func = [=](const double* val)
        {
            auto PCNs = nnPointCloud->Search( val, K );

            FindNormalisedProbability( PCNs, radius);

            double pb, pcloudassigment=0.0;
            for(uint32_t k=0; k<K; ++k)
            {
                auto& kppc = PCNs[k];

                // Find the visual hull implicit (hoppe) function
                pb = 0.0;
                for(uint32_t j=0; j<3; ++j)
                {
                    pb += (pointcloud.GetVertices()[3*kppc.index + j] - val[j])*
                            pointcloud.GetNormals()[3*kppc.index + j];
                }

                pcloudassigment += pb*PCNs[k].distance;
            }

            return pcloudassigment;
        };

        const double iso = 0.01;
        const uint32_t nd = 100;

        msg.print(2) << "Running marching cubes" << std::endl;
        double max[3], min[3];
        for(uint8_t k=0; k<3;++k)
        {
            max[k] = box.max[k];
            min[k] = box.min[k];
        }

        MC::MarchingCubes mc(iso, nd, max, min);
        mc.March(func);

        auto v = mc.GetPoints();
        auto f = mc.GetFaces();

        std::vector< PointTypes::Pointd3 > vo(v.size());
        std::vector< PointTypes::Pointi3 > fo(f.size());

#pragma omp parallel for
        for(uint32_t i=0; i<v.size(); ++i)
        for(uint8_t j=0; j<3; ++j)
            vo[i].p[j] = v[i].p[j];

#pragma omp parallel for
        for(uint32_t i=0; i<f.size(); ++i)
        for(uint8_t j=0; j<3; ++j)
            fo[i].p[j] = f[i].p[j];

        Mesh ret(vo, fo);

        return ret;
    }
}
