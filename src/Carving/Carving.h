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

/**
 * \brief A collection of functions to perform space carving.
 * The functions in this namespace are designed to take a calibrated
 * collection of cameras, a point cloud and some segmentations and generate a mesh which represents the
 * object of interest.
 *
 * \details There are three principle methods in order to achieve this.
 * \li The Visual Hull
 * \li Poisson reconstruction
 * \li Hoppe reconstruction
 *
 * \section visualhull The visual hull
 * The visual hull is a computation of the surface of a volume through space carving. The three dimensional
 * space enclosing the object is discretised in to a regular collection of points. Points which are not inside
 * the object are rejected, by projecting them on to the segmentations. More details on the original visual hulll
 * algorithm can be found in the paper,
 *
 *   "The Visual Hull Concept for Sillhouette Based Image Understanding"  - Aldo Laurentini, PAMI 1994.
 *
 * The visual hull algorithm used in this library is slightly different, since it makes use of the Marching Cubes
 * implementation. The Marching Cubes algorithm takes a functional defined on 3D space, which is positive on the inside of
 * an object, and negative on the outside; and returns a watertight mesh representing the boundary. In order to compute the
 * visual hull given a collection of calibrated segmentations the functional projects the given 3D point on to each, and
 * statistically combines then to find the probability that the 3D point is inside the object. If that probability
 * is high enough it returns a one, if not it will return a minus one.
 *
 * The probabilistic formulation means that the segmentations can be 'fuzzy' or probabilistic. Generally speaking, a segmentation
 * image will be white on the inside of an object and black on the outside. Regions where the segmentation is grey identify areas
 * which are unknown, such as on the boundary. Computing the visual hull in this way leads to smoother surfaces when the segmentations
 * are low quality.
 *
 * \section poisson Poisson reconstruction
 * The Poisson reconstruction is a method for reconstructing a surface from a point cloud. In order to comptute the surface in this
 * fashion, one first has to compute the surface normals for each point. A function is then fit to the points which is negative
 * behind the normals, zero on the surface, and positive in front of the normals. The Posson reconstruction method makes this fit
 * by solving the Poisson equation. More details can be found in the paper,
 *
 * "Poisson Surface Reconstruction" - Michael Kazhdan et al. Eurographics. 2006.
 *
 * The principle problem in this computation is computing an oriented set of normals. If the point cloud contains holes, for example,
 * there may be no suitable heuristic to decide how to make them consistent without object recognition. The Poisson reconstruction
 * can be computed by transferring the normals from the visual hull on to the point cloud. There is no implementation in this code, but it can be used
 * to transfer the normals for Poisson reconstruction in a third party library, such as MeshLab.
 *
 * \section hoppe Hoppe reconstruction
 * Hoppe surface reconstruction is an earlier example of implicit surface reconstruction. More information can be found in the following
 * paper,
 *
 *  "Surface Reconstruction from Unorganised Points" - Hughes Hoppe et al. ACM 1992.
 *
 * The method forulates a signed distance function between a point in 3D space, and the surface using the point cloud normals. The
 * signed distance function is beneficial because it can easily be put in to a probabilistic framework if the point cloud
 * exhibits a large amount of noise. In a probabilistic framework, it can also be combined with another functor such as the visual hull functor described
 * above. The combined Hoppe / Visual Hull method has been implemented in the Blend function of this library.
 *
 */
namespace Carve
{
    /**
     * @brief A box class for the marching cubes initialisation.
     * This class contains the extents of a 3D box for the initialisation of the Marching Cubes algorihtm.
     */
    class Box
    {
    public:
        Box() :
            center(3,0.0),
            max(3,-std::numeric_limits<double>::max()),
            min(3,std::numeric_limits<double>::max()) {}
        std::vector<double> center; ///< The center of the box
        std::vector<double> max; ///< The maximum of each dimension
        std::vector<double> min; ///< The minimum of each dimension
    };

    /**
     * @brief ComputeInitialBox
     * Compute the initial box for the marching cubes algortihm
     * @param P A vector of projection matrices
     * @param Images A vector of calibrated images ( in the same order as P )
     * @param msg An interface to write messages to stdout.
     * @param fractionPadding A fraction of the sample extends to add as padding.
     * @param Samples a vector of 3D vectors. Samples on the inside of the object for initialisation.
     * @return A box which completely encloses the object.
     */
    Box GetBoxFromCloud( const std::vector< Cu::Matrix >& P,
                         const std::vector< Image >& Images,
                         Stream::Message msg,
                         double fractionPadding,
                         const std::vector<std::vector<double> >&Samples );

    /**
     * @brief ComputeInitialBox
     * Compute the initial box for the marching cubes algortihm
     * @param P A vector of projection matrices
     * @param Images A vector of calibrated images ( in the same order as P )
     * @param msg An interface to write messages to stdout.
     * @param fractionPadding A fraction of the sample extends to add as padding.
     * @param Samples A column major vector of 3D points. Samples on the inside of the object for initialisation.
     * @return A box which completely encloses the object.
     */
    Box GetBoxFromCloud( const std::vector< Cu::Matrix >& P,
                         const std::vector< Image >& Images,
                         Stream::Message msg,
                         double fractionPadding,
                         const std::vector< double >&Samples );

    /**
     * @brief FindMesh - Compute the visual hull
     * Find the visual hull of an object given a collection of calibrated images. The function passes a
     * functor to a marching cubes algorithm.
     * @param P A vector of camera matrices
     * @param Images A vector of calibrated images (in the same order as P)
     * @param box The initial box for the marching cubes algorithm.
     * @param msg An interface to write messages to stdout.
     * @param density The number of divisions of each dimension for marching cubes. i.e. 200 initialises marching cubes with 200^3 boxes.
     * @return The mesh representing the boundary of the object (or visual hull).
     */
    Mesh FindMesh(const std::vector< Cu::Matrix >& P, const std::vector< Image >& Images,
                   const Box& box, Stream::Message msg = Stream::Message("FindMesh") , const uint32_t density=200);

    /**
     * @brief RemovePointsOutsideHull - Remove points outside the visual hull
     * This function will take a large and noisy point cloud and remove all of the points outside of
     * the visual hull.
     * @param mesh A reference to the point cloud (or mesh)
     * @param P A vector of projection matrices
     * @param Images A vector of calibrated images (in the same order as P)
     * @param msg An interface to write messages to stdout.
     */
    void RemovePointsOutsideHull(Mesh &mesh,
                                 const std::vector< Cu::Matrix >& P, const std::vector< Image >& Images,
                                 Stream::Message msg = Stream::Message("RemovePointsOutsideHull")); ///< Remove vertices outside the visual hull



    /**
     * @brief Hoppe - Compute the Hoppe surface reconstruction algorithm
     * This function computes the surface of a point cloud on which the surface normals have been computed. It
     * is an implementation of the paper "Surface Reconstruction from Unorganised Points" - Hughes Hoppe et al. ACM 1992.
     * @param pointcloud The point cloud (it must have the normals computed).
     * @param nnPointCloud A pointer to a nearest neighbours object of pointcloud.
     * @param box The initial box for the marching cubes algorithm.
     * @param K The number of neighbours to use when computing the signed distance function
     * @param radius The the radius parameter is the assumed std deviation of the error between the 3D point and the surface.
     * @param msg An interface to write messages to stdout.
     * @return A mesh representing the surface of the object.
     */
    Mesh Hoppe(const Mesh& pointcloud, const NearestNeighbours *nnPointCloud, const Box &box,
               const uint32_t K=20, const double radius=10, Stream::Message msg = Stream::Message("Blend"));

    /**
     * @brief Blend - A blend between Hoppe surface reconstruction and visual hull.
     * This function will blend the Hoppe surface reconstruction with the visual hull
     * as an alternative to optimisation. It requires the normals to be computed
     * on both of the mesh objects.
     * @param pointcloud The point cloud
     * @param nnPointCloud The nearest neighbour object for pointcloud
     * @param visualhull The visual hull (computed using FindMesh)
     * @param nnVhull The nearest neighbour object for visualhull
     * @param box  The initial box for the marching cubes algorithm.
     * @param K The number of neighbours to use when computing the signed distance function
     * @param radius The the radius parameter is the assumed std deviation of the error between the 3D point and the surface.
     * @param msg An interface to write messages to stdout.
     * @return A mesh representing the surface of the object.
     */
    Mesh Blend(const Mesh& pointcloud, const NearestNeighbours *nnPointCloud,
               const Mesh& visualhull, const NearestNeighbours *nnVHull, const Box &box,
               const uint32_t K=20, const double radius=10, Stream::Message msg = Stream::Message("Blend"));


    /**
     * @brief pdf The exponent of the zero mean Gaussian distribution
     * @param x The value
     * @param sigma The std deviation
     * @return
     */
    inline double pdf( double x, const double sigma)
    {
        x = (x*x) / (sigma*sigma);
        return -0.5 * x;
    }

    /**
     * @brief FindNormalisedProbability
     * Find the joint probability of the distance element of each keypair is zero
     * given that they are normally distributed with zero mean and std. deviation zero.
     * Replace the distance element with the probability and normalise using the
     * exp-sum-log calculation, which avoids numerical underflow.
     * @param X A vector of KeyPairs
     * @param sigma The std deviation of the distances
     * @return The average probability
     */
    inline double FindNormalisedProbability( std::vector<KD::KeyPair>& X, const double sigma)
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

    /**
     * @brief ComputeOutwardNormals
     * Compute the outward facing normals from a collection of cameras. i.e. the direction
     * that the camera is pointing.
     * @param P A vector of camera matrices
     * @return A vector of normals
     */
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



#endif
