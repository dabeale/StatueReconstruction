#ifndef ZBUFFER_H
#define ZBUFFER_H

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
#include <list>
#include <algorithm>
#include <stdint.h>
#include <iostream>
#include <cassert>
#include <cmath>

#include "Matrix.h"

/**
 * \brief A collection of functions, structs and classes for computing the ZBuffer of a mesh, or point cloud, from a particular view.
 */
namespace Buffer
{
    typedef std::vector<double> BVector; ///< \brief Vector for matrix algebra

    /**
     * @brief The ProjectionMatrix class.
     * An imlementation of a homogeneous 3x4 matrix. This imlpementation
     * is distinct from other matrix classes because the overloaded multiplication operator
     * ensures that the final element of the returned vector is one.
     *
     * The implementation is compatible with the auto differentiation scheme used in
     * ceres, and so it can be used for mesh optimisation.
     *
     */
    class ProjectionMatrix
    {
    public:
        /**
         * @brief Construct an empty projection matrix.
         * The matrix is 3x4 with all elements set to zero.
         */
        ProjectionMatrix();

        /**
         * @brief Construct a 3x4 projection matrix filled with a particular value
         * @param val The value to fill the matrix with
         */
        ProjectionMatrix(double val);

        /**
         * @brief Construct a projection matrix and fill it with the values in the input vector.
         * The matrix is filled in a column major fashion.
         * @param arr A vector of size 12 to fill the matrix with.
         */
        ProjectionMatrix(const BVector& arr );


        BVector operator *(const BVector& X);          ///< \brief Overload the multiplication operator for a matrix with a BVector
        BVector operator *(const BVector& X) const;    ///< \brief Overload the const multiplication operator for a matrix with a BVector

        /**
         *  \brief operator *.
         * Overload the multiplication operator for a matrix with a pointer to a type T. The method
         * is templated so that is can accept standard types and also the auto differentiation
         * types used in ceres.
         */
        template<typename T>
        std::array<T,3> operator *(const T* const X)
        {
            std::array<T,3> x;
            x[0] = m_arr[0]*X[0] +m_arr[3]*X[1] +m_arr[6]*X[2] +m_arr[9]*X[3];
            x[1] = m_arr[1]*X[0] +m_arr[4]*X[1] +m_arr[7]*X[2] +m_arr[10]*X[3];
            x[2]=  m_arr[2]*X[0] +m_arr[5]*X[1] +m_arr[8]*X[2] +m_arr[11]*X[3];
            x[0] /= x[2];
            x[1] /= x[2];
            x[2] = T(1);
            return x;
        }

        /**
         *  \brief operator *.
         *  Overload the const multiplication operator for a matrix with a pointer to a type T. The method
         *  is templated so that is can accept standard types and also the auto differentiation
         *  types used in ceres.
         */
        template<typename T>
        std::array<T,3> operator *(const T* const X) const
        {
            std::array<T,3> x;
            x[0] = m_arr[0]*X[0] +m_arr[3]*X[1] +m_arr[6]*X[2] +m_arr[9]*X[3];
            x[1] = m_arr[1]*X[0] +m_arr[4]*X[1] +m_arr[7]*X[2] +m_arr[10]*X[3];
            x[2] = m_arr[2]*X[0] +m_arr[5]*X[1] +m_arr[8]*X[2] +m_arr[11]*X[3];
            x[0] /= x[2];
            x[1] /= x[2];
            x[2] = T(1);
            return x;
        }


        /**
         * @brief ComputeLocation.
         * Return the location of the camera in 3D space. This method assumes that the camera is canonical,
         * in the sense that the first 3x3 submatrix is unitary.
         * @return A BVector of size 3, containing the location of the focus of the camera.
         */
        BVector ComputeLocation() const;

        /**
         * @brief data.
         * Return a reference to the internal vector used to store the data. It is of size 12 and column major.
         * @return A vector reference.
         */
        inline std::vector<double>& data(){ return m_arr;}

        /**
         * @brief at.
         * Return the (i,j)th element of the matrix
         * @param i The row
         * @param j The column
         * @return The value at (i,j)
         */
        inline double at( const uint32_t i, const uint32_t j) const
        {
            return m_arr[j*3 + i];
        }

        /**
         * @brief at.
         * Return a refernce to the (i,j)th element of the matrix
         * @param i The row
         * @param j The column
         * @return The value at (i,j)
         */
        inline double& at( const uint32_t i, const uint32_t j)
        {
            return m_arr[j*3 + i];
        }
    private:
        BVector m_arr;  ///< A vector of length 12, in column major format, containing the entries of the matrix.

        /**
         * @brief at.
         * Return the (i,j)th element of the matrix
         * @param i The row
         * @param j The column
         * @return The value at (i,j)
         */
        double el( const uint32_t i, const uint32_t j) const;

        /**
         * @brief at.
         * Return a refernce to the (i,j)th element of the matrix
         * @param i The row
         * @param j The column
         * @return The value at (i,j)
         */
        double& el( const uint32_t i, const uint32_t j);
    };

    /**
     * @brief The ZBuffer class.
     * A Z buffer for a collection of vertices and faces. There is a mesh implementation which will efficiently compute the
     * depths at every pixel and also a point cloud implementation which models each point as a blob before projecting it
     * on to the target image.
     *
     * The implementation is designed to be fast, but does not use any external libraries except octave for its
     * matrices. There is no use of the GPU or of multiple processors. A faster implementation
     * exists in OpenGL which uses the GPU and can do this computation in real time, compared to
     * about half a second for an image in the order of 10M pixels, and a mesh with 1M faces.
     */
    class ZBuffer
    {
    public:
        /**
         * @brief The Mesh struct.
         * A container for the vertices and faces of a mesh. This is the simplest possible representation
         * for a mesh.
         */
        struct Mesh
        {
            const double* vertices; ///< A pointer to the column major vertex array
            const uint32_t* faces;  ///< A pointer to the column major face array
            const uint32_t Nv;      ///< The number of vertices
            const uint32_t Nf;      ///< The number of faces
        };

        /**
         * @brief Construct an empty ZBuffer
         */
        ZBuffer();

        /**
         * @brief ZBuffer
         * @param mesh
         * @param height
         * @param width
         * @param P
         */
        ZBuffer(const Mesh mesh,
                const uint32_t height,
                const uint32_t width,
                const ProjectionMatrix& P);
        ~ZBuffer(); ///< Empty destructor

        /**
         * @brief ComputePointCloudDepthBuffer.
         * Compute the depth buffer based on the vertices. This method models each vertex as a blob and fills the
         * depth buffer with the closest depths to the camera.
         * @param blobsize The size, in pixels, of each blob.
         */
        void ComputePointCloudDepthBuffer(  const int32_t blobsize = 5 );

        /**
         * @brief ComputePointCloudVisibilities.
         * Copmute the mesh  visibilities based on the point cloud. This method requires ComputePointCloudDepthBuffer() to be called first,
         * and computed the per vertex visibilities on the point cloud for the given camera.
         * @param blobsize
         */
        void ComputePointCloudVisibilities( const int32_t blobsize = 5 );

        /**
         * @brief ComputeDepthBuffer.
         * Compute the depth buffer, based on vertices and faces. This method will fill an image with the depths
         * of the nearest side of the mesh to the camera.
         */
        void ComputeDepthBuffer();

        /**
         * @brief ComputeVisibilities.
         * Compute the vertex visibilities for a triangle mesh. This method requires ComputeDepthBuffer() to be called first,
         * and computes the per vertex visibilities on the mesh of the given camera.
         */
        void ComputeVisibilities();

        const std::vector<double>& ReturnDepthBuffer() const;                   ///< \brief Return a reference to the current depth buffer
        const std::vector<int32_t>& ReturnIndexBuffer() const {return m_ib;}    ///< \brief Return a reference to the index buffer. It is the same size as the depth buffer and contains the vertex that is visible at each pixel.
        const std::vector<bool>& ReturnVertexVisibilities() const;              ///< \brief Return a reference to the vertex visibilities. A bool array the same size as the number of vertices.
        const std::vector<double>& ReturnVertexDepth();                         ///< \brief Compute and return the current vertex depths for camera @param P. An array the same size as the number of vertices, indicating the depths for each vertex.

    private:
        const double* m_vertices;            ///< A pointer to the vertex array
        const uint32_t* m_faces;             ///< A pointer to the face index array
        uint32_t m_Nv;                       ///< The number of vertices
        uint32_t m_Nf;                       ///< The number of faces
        uint32_t m_height;                   ///< The height of the buffer
        uint32_t m_width;                    ///< The width of the buffer

        ProjectionMatrix m_P;                ///< The projection matrix

        std::vector<double> m_db;            ///< The depth buffer. A column major matrix m_height*m_width containing the depths for the camera view (it contains a -1 where the depth is undefined).
        std::vector<int32_t> m_ib;           ///< The index buffer. A vector of indicies in to the vertex list. The value -1 if is has an undefined depth.
        std::vector<double> m_vertexDepth;   ///< The depth of each vertex. A vector the same size as the number of vertices.
        std::vector<double> m_faceDepth;     ///< The depth of each face. A vector the same size as the number of faces.
        std::vector<bool> m_visibilities;    ///< True if the vertex is visible. A vector the same size as the number of vertices.


        /**
         * @brief GetDepthBufferValue.
         * Get the value of the depth buffer at a given index
         * @param a A vector of length two containing the location on the image (depth buffer)
         * @return The depth
         */
        inline double GetDepthBufferValue( const BVector& a )
        {
            return m_db[m_height*static_cast<uint32_t>(a[0])+static_cast<uint32_t>(a[1])];
        }

        /**
         * @brief GetDistanceFromTriangle.
         * A private method to compute the distance of an arbitrary point from a triangle. To save time
         * this method simply computes the distance from the center of the triangle. All of the vectors
         * are assumed to be in 3D space.
         * @param Center The point from which to compute the distance
         * @param a The first point on the triangle
         * @param b The second point on the triangle
         * @param c Thie third point on the triangle
         * @return The distance to the triangle
         */
        inline double GetDistanceFromTriangle( const BVector& Center, const BVector& a, const BVector& b, const BVector& c )
        {
            double distance = 0.0;
            for (int k=0;k<3;k++)
            {
                distance += (Center[k] - (a[k] + b[k] + c[k])/3.0) *
                            (Center[k] - (a[k] + b[k] + c[k])/3.0);
            }
            return std::sqrt(distance);
        }

        /**
         * @brief GetDistanceFromPoint.
         * A private method to compute the distance to a point.
         * @param Center The point from which to compute the distance
         * @param a The point to compute the distance to
         * @return The Euclidean distance between Center and a.
         */
        inline double GetDistanceFromPoint( const BVector& Center, const BVector& a )
        {
            double distance = 0.0;
            for (int k=0;k<3;k++)
            {
                distance += (Center[k] - a[k]) *
                            (Center[k] - a[k]);
            }
            return std::sqrt(distance);
        }
    };

    /**
     * @brief InsideTriangle.
     * Check whether a point (x) lies within a triangle (a,b,c), on a plane.
     * @param a The first vertex on the triangle
     * @param b The second vertex on the triangle
     * @param c The third vertex on the triangle
     * @param x The point to test
     */
    inline bool InsideTriangle( const BVector& a, const BVector&b, const BVector& c, const BVector& x );

    /**
     * @brief SetVertices.
     * Fill a 3 vectors with the vertices on the given face in a mesh.
     * @param vertices A pointer to the array of vertices 3xN column major
     * @param faces A pointer to the array of faces 3xK column major
     * @param a The first vertex on the face
     * @param b The second vertex on the face
     * @param c The third vertex on the face
     * @param i The face index to fill the vectors with
     */
    inline void SetVertices(const double* vertices, const double* faces, BVector &a, BVector& b, BVector& c, const int& i);

    /**
     * @brief FindCameraVisibilities.
     * Given a collection of camera matrices and a mesh find out which vertices are visible from which cameras. This method
     * constructs a collection of depth buffer objects and computes the visibilities for each camera. A vector of lists is returned,
     * one list for each vertex, indicating which cameras are visible from each vertex.
     * @param mesh An input mesh
     * @param height The height of the camera view
     * @param width The width of the camera view
     * @param Ps A vector of projection metrices.
     * @return A vector of lists, of size N, where N is the number of vertices. The lists contain indexes in to the projection matrix vector.
     */
    std::vector< std::list<uint32_t> > FindCameraVisibilities( const ZBuffer::Mesh mesh, uint32_t height, uint32_t width, const std::vector< ProjectionMatrix >& Ps  );
}
#endif
