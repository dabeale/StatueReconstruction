#ifndef SAMPLER_H
#define SAMPLER_H

#include "Matrix.h"
#include <Eigen/Cholesky>
#include "Triangulate.h"

#include <random>
#include <chrono>

namespace Sample
{

    class Generator
    {
    public:
        Generator( );

        std::vector<int32_t> SampleUniformInteger(const int32_t min, const int32_t max, const uint32_t NSamples);
        std::vector<uint32_t> SampleUniformIntegerUI(const uint32_t min, const uint32_t max, const uint32_t NSamples);
        std::vector<uint32_t> SampleDiscrete1D( const Cu::Matrix& distribution, const uint32_t N );
        std::vector<uint32_t> SampleDiscrete1D(const std::vector<double> &dist, const uint32_t N );
        std::vector<uint32_t> RandomPermutation( const uint32_t N );
        std::vector<std::pair<uint32_t, uint32_t> > SampleDiscrete2D( const Cu::Matrix& distribution, const uint32_t N );

        Cu::Matrix SampleMultivariateNormal( const Cu::Matrix& mu, const Cu::Matrix& Sigma, const uint32_t N );

        /**
         * @brief SampleWishart
         * Sample from a wishart distribution. This implementation is taken from the octave source code, originally written by Nir Krakauer.
         * It is an adaptation for c++.
         * @param Sigma
         * @param df
         * @param N
         * @return
         */
        std::vector<Cu::Matrix> SampleWishart( const Cu::Matrix& Sigma, const uint32_t df, const uint32_t N );

        enum AlternateCameraMethod
        {
            Random, // Choose a random camera
            Orthogonal // Choose an orthogonal camera (slower)
        };

        /**
         * @brief Sample3DEdge
         * This method samples a 3D edge from a collection of curves by first casting the curve mixture
         * to a discrete distribution.
         * Care must be taken that the projection is in the correct order. i.e.
         * if x = PX, x(1) must be the rows of the image, and x(2) must be the columns.
         * This is not consistent between reconstruction libraries. Either swap the rows
         * in the correspoinding projection matrix or transpose the images to ensure that
         * it works.
         * @param edges The 2D edge maps
         * @param P The projection matrices
         * @param K The number of points to sample
         * @return
         */
        std::vector<std::vector<double> > Sample3DEdge( const std::vector<Cu::Matrix>& edges,
                                                       const std::vector<Cu::Matrix>& P,
                                                       const uint32_t K,
                                                       const AlternateCameraMethod acm = Random);

        /**
         * @brief Sample3DEdgeResampling
         * This method does the same as above but uses resampling on each of the conditionals.
         * It avoids the artifacts due to quantisation and also reduces memory requirements.
         * @param edgeMeans The curve means
         * @param edgeDeviations The curve deviations
         * @param P The projection matrices
         * @param K The number of samlpes to take
         * @param deviation The deviation in the fundeamental matrices
         * @return
         */
        std::vector<Cu::Matrix> Sample3DEdgeResampling(const std::vector<std::vector<std::pair<double, double> > > &edgeMeans,
                                                                  const std::vector<std::vector<double> > &edgeDeviations,
                                                                  const std::vector<Cu::Matrix>& P,
                                                                  const uint32_t K, const double deviation);
    protected:
        std::default_random_engine m_e;
        uint32_t seed;
    };

    inline std::vector<Cu::Matrix> ComputeOutwardNormals(const std::vector<Cu::Matrix>& P)
    {
        std::vector<Cu::Matrix> outwardNormals = std::vector<Cu::Matrix>(P.size(), Eigen::MatrixXd::Zero(3,1));
        Cu::Matrix no(3,1),n(3,1);
        no(0) = 0;
        no(1) = 0;
        no(2) = 1.0;
        for(uint32_t k=0; k<P.size(); ++k)
        {
            Cu::Matrix preRQ = P[k].block(0,0,3,3);
            std::swap(preRQ(0,0), preRQ(2,0));
            std::swap(preRQ(0,1), preRQ(2,1));
            std::swap(preRQ(0,2), preRQ(2,2));

            Eigen::HouseholderQR<Cu::Matrix> qr(preRQ.transpose());
            Cu::Matrix Q = qr.householderQ().transpose();
            std::swap(Q(0,0), Q(2,0));
            std::swap(Q(0,1), Q(2,1));
            std::swap(Q(0,2), Q(2,2));

            n = Q.transpose()*no;
            outwardNormals[k] = n;
        }
        return outwardNormals;
    }
}

#endif // SAMPLER_H
