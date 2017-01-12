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


class Mesh
{
public:
    enum Primitive
    {
        UnitSphere
    }; ///< A collection of primitive types

    Mesh();
    Mesh( const Primitive primitive );
    Mesh( const std::string& filename );
    Mesh(const std::vector< std::vector<double> >& vertices);
    Mesh( const std::vector< PointTypes::Pointd3 >& vertices, const std::vector< PointTypes::Pointi3 >& faces);
    Mesh( std::vector< double > vertices, std::vector< uint32_t > faces,
          std::vector< double > normals, std::vector< double > colours );

    void Load( const std::string& filename );
    void Write( const std::string& filename );

    inline uint32_t GetNVerts() const {return m_NVerts;} ///< Get the number of vertices
    inline uint32_t GetNFaces() const {return m_NFaces;} ///< Get the number of faces

    std::vector<double> GetVertex( const uint32_t i) const;
    std::vector<uint32_t> GetFace( const uint32_t i) const;
    std::vector<double> GetColour( const uint32_t i) const;
    std::vector<double> GetNormal( const uint32_t i) const;

    void ComputeNormalsFromPointsAndFaces();
    void ApplyNormals( const Mesh& mesh, const NearestNeighbours& nn );
    double ComputeAverageEdgeLength() const;
    std::vector<double> ComputeMean() const;
    std::vector<double> ComputeClosestPoint( const std::vector<double>& pt ) const;
    double ComputeAverageDistance( const std::vector<double>& pt ) const ;

    void FillNRing(); ///< Fill the map from each vertex to its neighbours
    std::set<uint32_t> GetNRing( const uint32_t i);

    void TranslateVertices(const std::vector<double> &A, const std::vector<double> &b );

    inline std::vector<double>& GetVertices(){ return m_vertices; }
    inline std::vector<uint32_t>& GetFaces(){ return m_faces; }
    inline std::vector<double>& GetNormals(){ return m_normals; }
    inline std::vector<double>& GetColours(){ return m_colours; }
    inline const std::vector<double>& GetVertices() const{ return m_vertices; }
    inline const std::vector<uint32_t>& GetFaces() const{ return m_faces; }
    inline const std::vector<double>& GetNormals() const{ return m_normals; }
    inline const std::vector<double>& GetColours() const{ return m_colours; }

    inline void SetNVertices( const uint32_t nverts ){ m_NVerts = nverts; }

private:
    std::vector< double > m_vertices;
    std::vector< uint32_t > m_faces;
    std::vector< double > m_normals;
    std::vector< double > m_colours;

    uint32_t m_NVerts;
    uint32_t m_NFaces;


    std::vector< std::set<uint32_t> > m_NRing;

    void ComputeTangents(double *n1, double *n2, const uint32_t face);
};

#endif // MESH_H
