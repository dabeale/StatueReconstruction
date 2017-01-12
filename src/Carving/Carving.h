#ifndef CARVING_H
#define CARVING_H

/* ----------------------------------------------------------------------
 * Copyright (C) 2016 Daniel Beale. All rights reserved.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 * ---------------------------------------------------------------------- */

#include <vector>
#include <functional>
#include <limits>

#include "Image.h"
#include "Matrix.h"
#include "Sample.h"
#include "Mesh.h"
#include "Message.h"

#include "marchingcubes.h"

namespace Carve
{
    class Box
    {
    public:
        Box() :
            center(3,0.0),
            max(3,-std::numeric_limits<double>::max()),
            min(3,std::numeric_limits<double>::max()) {}
        std::vector<double> center;
        std::vector<double> max;
        std::vector<double> min;
    }; ///< A box struct for the marching cubes initialisation

    /**
     * @brief ComputeInitialBox
     * Compute the initial box for the marching cubes
     * @param P
     * @param Images
     * @return
     */
    Box GetBoxFromCloud( const std::vector< Cu::Matrix >& P,
                         const std::vector< Image >& Images,
                         Stream::Message msg,
                         double fractionPadding,
                         const std::vector<std::vector<double> >&Samples );
    Box GetBoxFromCloud( const std::vector< Cu::Matrix >& P,
                         const std::vector< Image >& Images,
                         Stream::Message msg,
                         double fractionPadding,
                         const std::vector< double >&Samples );


    Mesh FindMesh(const std::vector< Cu::Matrix >& P, const std::vector< Image >& Images,
                   const Box& box, Stream::Message msg = Stream::Message("FindMesh") , const uint32_t density=200);

    void RemovePointsOutsideHull(Mesh &mesh,
                                 const std::vector< Cu::Matrix >& P, const std::vector< Image >& Images,
                                 Stream::Message msg = Stream::Message("RemovePointsOutsideHull")); ///< Remove vertices outside the visual hull

    /**
     * @brief Blend
     * This function will blend the point cloud with the visual hull
     * as an alternative to optimisation. It requires the normals to be computed
     * on both of the mesh objects.
     * @param pointcloud
     * @param visualhull
     * @return
     */
    Mesh Blend(const Mesh& pointcloud, const NearestNeighbours *nnPointCloud,
               const Mesh& visualhull, const NearestNeighbours *nnVHull, const Box &box,
               const uint32_t K=20, const double radius=10, Stream::Message msg = Stream::Message("Blend"));

    Mesh Hoppe(const Mesh& pointcloud, const NearestNeighbours *nnPointCloud, const Box &box,
               const uint32_t K=20, const double radius=10, Stream::Message msg = Stream::Message("Blend"));
}

inline double pdf( double x, const double sigma)
{
    x = (x*x) / (sigma*sigma);
    return -0.5 * x;
}

/**
 * @brief FindNormalisedProbability
 *
 * @param X
 * @param sigma
 * @return  void
 */
inline double FindNormalisedProbability( std::vector<KeyPair>& X, const double sigma)
{
    double ave=0.0, sum=0.0, max=-std::numeric_limits<double>::max();
    uint32_t k=0;

    for(k=0; k<X.size(); ++k)
    {
        ave += X[k].distance;
        X[k].distance = pdf(X[k].distance, sigma);
        if(X[k].distance > max) max = X[k].distance;
    }
    ave /= X.size();

    for(k=0; k<X.size(); ++k)
    {
        X[k].distance = std::exp( X[k].distance - max );
        sum += X[k].distance;
    }

    for(k=0; k<X.size(); ++k)
    {
        X[k].distance /= sum;
    }

    return ave;
}

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

#endif
