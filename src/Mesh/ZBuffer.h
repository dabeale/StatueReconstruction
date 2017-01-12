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

namespace Buffer
{
    typedef std::vector<double> BVector; ///< Vector for matrix algebra

    class ProjectionMatrix
    {
    public:
        ProjectionMatrix();
        ProjectionMatrix(double val);
        ProjectionMatrix(const BVector& arr );

        BVector operator *(const BVector& X);
        BVector operator *(const BVector& X) const;

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


        BVector ComputeLocation() const; ///< Get the location

        inline std::vector<double>& data(){ return m_arr;}

        inline double at( const uint32_t i, const uint32_t j) const
        {
            return m_arr[j*3 + i];
        }

        inline double& at( const uint32_t i, const uint32_t j)
        {
            return m_arr[j*3 + i];
        }
    private:
        BVector m_arr;

        double el( const uint32_t i, const uint32_t j) const;
        double& el( const uint32_t i, const uint32_t j);
    };

    /**
     * @brief The ZBuffer class
     * This class represents a z buffer for a collection of vertices and faces. The implementation
     * is supposed to be fast, but does not use any external libraries except octave for its
     * matrices. There is no use of the GPU or of multiple processors. A faster implementation
     * exists in OpenGL which uses the GPU and can do this computation in real time, compared to
     * about half a second for an image in the order of 10M pixels, and a mesh with 1M faces.
     */
    class ZBuffer
    {
    public:
        struct Mesh
        {
            const double* vertices;
            const uint32_t* faces;
            const uint32_t Nv;
            const uint32_t Nf;
        };

        /**
         * @brief ZBuffer
         * Contructor for the ZBuffer class
         * @param vertices A pointer to the vertex list (column major with 3 rows)
         * @param faces A pointer to the face list (column major 3 rows)
         * @param Nv The number of vertices
         * @param Nf The number fo faces
         * @param height The height of the image
         * @param width the width of the image
         */
        ZBuffer();
        ZBuffer(const Mesh mesh,
                const uint32_t height,
                const uint32_t width,
                const ProjectionMatrix& P);
        ~ZBuffer(); ///< Empty destructor

        void ComputePointCloudDepthBuffer(const int32_t blobsize = 5); ///< Compute the depth buffer based on the vertices
        void ComputePointCloudVisibilities( const int32_t blobsize = 5 ); ///< Copmute the visibilities based on the point cloud (requires ComputePointCloudDepthBuffer() to be called first)


        void ComputeDepthBuffer();           ///< Compute the depth buffer, based on vertices and faces
        void ComputeVisibilities();          ///< Compute the vertex visibilities for a triangle mesh (requires ComputeDepthBuffer() to be called)

        const std::vector<double>& ReturnDepthBuffer() const;      ///< Return the current depth buffer
        const std::vector<int32_t>& ReturnIndexBuffer() const {return m_ib;} ///< Return the index buffer
        const std::vector<bool>& ReturnVertexVisibilities() const; ///<
        const std::vector<double>& ReturnVertexDepth(); ///< compute and return the current vertex depths for camera @param P

    private:
        const double* m_vertices; ///< A pointer to the vertex array
        const uint32_t* m_faces;    ///< A pointer to the face index array
        uint32_t m_Nv;            ///< the number of vertices
        uint32_t m_Nf;            ///< The number of faces
        uint32_t m_height;        ///< The height of the buffer
        uint32_t m_width;         ///< The width of the buffer

        ProjectionMatrix m_P; ///< The projection matrix

        std::vector<double> m_db;              ///< The depth buffer
        std::vector<int32_t> m_ib;            ///< The index buffer. The the point has a depth it contains the index of the vertex (-1 otherwise)
        std::vector<double> m_vertexDepth;     ///< The depth of each vertex
        std::vector<double> m_faceDepth;       ///< The depth of each face
        std::vector<bool> m_visibilities;           ///< True if the vertex is visible


        inline double GetDepthBufferValue( const BVector& a )
        {
            return m_db[m_height*static_cast<uint32_t>(a[0])+static_cast<uint32_t>(a[1])];
        }

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
     * @brief InsideTriangle
     * Check whether a point (x) lies within a triangle (a,b,c)
     * This is not used in favour of more efficient inline code
     */
    inline bool InsideTriangle( const BVector& a, const BVector&b, const BVector& c, const BVector& x );
    inline void SetVertices(const double* vertices, const double* faces, BVector &a, BVector& b, BVector& c, const int& i);

    std::vector< std::list<uint32_t> > FindCameraVisibilities( const ZBuffer::Mesh mesh, uint32_t height, uint32_t width, const std::vector< ProjectionMatrix >& Ps  );
}
#endif
