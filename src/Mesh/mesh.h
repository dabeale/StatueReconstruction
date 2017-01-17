#ifndef MESH_H
#define MESH_H

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

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <set>

#include <cassert>

#include "rply.h"
#include "NearestNeighbours.h"

namespace PointTypes
{
    typedef struct
    {
        double p[3];
    } Pointd3;

    typedef struct
    {
        uint32_t p[3];
    } Pointi3;
}

/**
 * @brief The Mesh class
 * This class represents a triangle mesh. It is a simple, non-extensible, implementation; designed
 * to provide a thread-safe interface to mesh algorithms. The class has the potential to store
 * the vertices and faces, and also the vertex colours and normals. It is possible to efficiently compute the
 * nearest neighbours on the mesh also.
 */
class Mesh
{
public:
    /**
     * @brief The Primitive enum
     * An enumeration of primitives which can be used to
     * initialise the mesh object. Currently only the unit sphere is implemented.
     */
    enum Primitive
    {
        UnitSphere
    };

    /**
     * @brief Mesh
     * Construct an empty mesh.
     */
    Mesh();

    /**
     * @brief Mesh
     * Construct a mesh with the given primitive type
     * @param primitive The primitive
     */
    Mesh( const Primitive primitive );

    /**
     * @brief Mesh
     * Construct a mesh from a file. Currently only the stanford ply format is
     * supported.
     * @param filename
     */
    Mesh( const std::string& filename );

    /**
     * @brief Mesh
     * Construct a mesh from a vector of vertices. The faces, normals and colours are all
     * set to be empty in this constructor.
     * @param vertices A vector of vectors
     */
    Mesh(const std::vector< std::vector<double> >& vertices);

    /**
     * @brief Mesh
     * Construct a mesh from a vector of vertices and a vector of faces. The faces must be integer
     * indexes in to the vertex vector.
     * @param vertices A vector of double precision arrays
     * @param faces A vector of integral arrays
     */
    Mesh( const std::vector< PointTypes::Pointd3 >& vertices, const std::vector< PointTypes::Pointi3 >& faces);

    /**
     * @brief Mesh
     * Construct a mesh from a vector of vertices, faces, per vertex normals and per vertex colours
     * @param vertices A column major matrix of vertices of size 3N, where N is the number of vertices
     * @param faces A column major matrix of faces os size 3K where K is the number of faces
     * @param normals   A column major matrix of per vertex normals of size 3N
     * @param colours A column major matrix or per vertex colours of size 3N
     */
    Mesh( std::vector< double > vertices, std::vector< uint32_t > faces,
          std::vector< double > normals, std::vector< double > colours );

    /**
     * @brief Load
     * Load a ply from file. Currently only the stanford ply format is supported.
     * @param filename The path to the file.
     */
    void Load( const std::string& filename );

    /**
     * @brief Write
     * Write the mesh to disk. Currently only the stanford ply format is supoprted.
     * @param filename the path to the file.
     */
    void Write( const std::string& filename );

    inline uint32_t GetNVerts() const {return m_NVerts;}        ///< \brief Get the number of vertices
    inline uint32_t GetNFaces() const {return m_NFaces;}        ///< \brief Get the number of faces

    std::vector<double> GetVertex( const uint32_t i) const;     ///< \brief Return a vector of size 3 corresponding to the vertex at index i
    std::vector<uint32_t> GetFace( const uint32_t i) const;     ///< \brief Return a vector of size 3 corresponding to the face at index i
    std::vector<double> GetColour( const uint32_t i) const;     ///< \brief Return a vector of size 3 corresponding to the colour at index i
    std::vector<double> GetNormal( const uint32_t i) const;     ///< \brief Return a vector of size 3 corresponding to the normal at index i

    /**
     * @brief ComputeNormalsFromPointsAndFaces
     * Compute the per vertex normals by iterating through the faces and adding the face normal to a
     * per vertex list. The final vertex normal is computed as an average of all the neighbouring
     * face normals. The face normals are computed using the cross product.
     */
    void ComputeNormalsFromPointsAndFaces();

    /**
     * @brief ApplyNormals
     * Apply the normals from an input mesh to 'this'. The method requires a nearest neighbours object
     * to find the nearest vertex in terms of Euclidean distance to apply the normal from.
     *
     * NOTE: There is no thresholding applied in this method, it will apply the normal no matter how far
     * away the corresponding normal on the input mesh is.
     * @param mesh A mesh from whish to apply the normals
     * @param nn The pre-computed NearestNeighbour object for mesh
     */
    void ApplyNormals( const Mesh& mesh, const NearestNeighbours& nn );

    /**
     * @brief ComputeAverageEdgeLength
     * Return the average edge length
     * @return The average edge length
     */
    double ComputeAverageEdgeLength() const;

    /**
     * @brief ComputeMean
     * Compute the mean vertex
     * @return A vector of length 3, the mean vertex.
     */
    std::vector<double> ComputeMean() const;

    /**
     * @brief ComputeClosestPoint
     * Compute the closest vertex to a given point. This algorithm search all of the vertices, it is not
     * based on a nearest neighbours algorithm.
     * @param pt A vector of length three, the test point.
     * @return A vector of length three, the closest vertex on the mesh.
     */
    std::vector<double> ComputeClosestPoint( const std::vector<double>& pt ) const;

    /**
     * @brief ComputeAverageDistance
     * Compute the average distance from the given test point, to all of the points in the mesh.
     * @param pt A vector of length 3, the test point.
     * @return The average distance
     */
    double ComputeAverageDistance( const std::vector<double>& pt ) const ;

    /**
     * @brief FillNRing
     * Compute the NRing for each of the vertices in the mesh, that is, all of the neighbouring
     * vertices which are connected by 1 edge.
     */
    void FillNRing();

    /**
     * @brief GetNRing
     * Get the NRing around a given vertex; all of the vertex indices which are connected to
     * the given vertex by an edge. This method requires FillNRing to be called beforehand.
     * @param i Return the NRing for this vertex
     * @return The NRing
     */
    std::set<uint32_t> GetNRing( const uint32_t i);

    /**
     * @brief TranslateVertices
     * Perform a affine transformation on each of the vertices. Supposing that
     * v is an arbitrary vertex and A and b are the parameters for the affine  transformation,
     * compute a new vertex w = Av + b.
     * @param A A vector of length 9 representing a 3x3 matrix, in column major form.
     * @param b A vector of length 3.
     */
    void TranslateVertices(const std::vector<double> &A, const std::vector<double> &b );

    inline std::vector<double>& GetVertices(){ return m_vertices; }                 ///< \brief Return a reference to the vector of vertices
    inline std::vector<uint32_t>& GetFaces(){ return m_faces; }                     ///< \brief Return a reference to the vector of faces
    inline std::vector<double>& GetNormals(){ return m_normals; }                   ///< \brief Return a reference to the vector of normals
    inline std::vector<double>& GetColours(){ return m_colours; }                   ///< \brief Return a reference to the vector of colours
    inline const std::vector<double>& GetVertices() const{ return m_vertices; }     ///< \brief Return a const reference to the vector of vertices
    inline const std::vector<uint32_t>& GetFaces() const{ return m_faces; }         ///< \brief Return a const reference to the vector of faces
    inline const std::vector<double>& GetNormals() const{ return m_normals; }       ///< \brief Return a const reference to the vector of normals
    inline const std::vector<double>& GetColours() const{ return m_colours; }       ///< \brief Return a const reference to the vector of colours

    inline void SetNVertices( const uint32_t nverts ){ m_NVerts = nverts; }         ///< \brief Set the number of vertices to be equal to NVerts

private:
    std::vector< double > m_vertices;                   ///< A vector of vertices of size 3N, where N is the number of vertices (column major).
    std::vector< uint32_t > m_faces;                    ///< A vector of faces of size 3K, where K is the number of faces (column major).
    std::vector< double > m_normals;                    ///< A vector of per vertex normals of size 3N, where N is the number of vertices (column major).
    std::vector< double > m_colours;                    ///< A vector of per vertex colours of size 3N, where N is the number of vertices (column major).

    uint32_t m_NVerts;                                  ///< The number of vertices
    uint32_t m_NFaces;                                  ///< The number of faces


    std::vector< std::set<uint32_t> > m_NRing;          ///< A vector of sets of indexes, representing the NRing of each vertex. This vector is of size m_NVerts.

    /**
     * @brief ComputeTangents
     * Compute a pair of tangent vectors for a given face. This pair of vectors form a basis for
     * the tangent space at the face, but are arbitrarily chosen from the set of all possible bases.
     * The tangent vectors are othonormal.
     * @param n1 The first tangent vector
     * @param n2 The second tangent vector.
     * @param face The face at which to compute the tangents.
     */
    void ComputeTangents(double *n1, double *n2, const uint32_t face);
};

#endif // MESH_H
