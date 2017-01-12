#ifndef MARCHINGCUBES_H
#define MARCHINGCUBES_H

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

#if defined(_OPENMP)
    #include <libiomp/omp.h>
#endif

#include <numeric>
#include <complex>
#include <iostream>
#include <vector>

//#include "mlists.cpp"

extern int wq[12];
extern int ww[12];
extern int pw[12];
extern int edgeTable[256];
extern int triTable[256][16];

#ifdef RENDER_USING_GLUT
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#endif

#include <vector>
#include <functional>
#include <map>
#include <fstream>
#include <string>

namespace MC
{
    typedef struct
    {
        double p[3];
    } Pointd3;

    typedef struct
    {
        uint32_t p[3];
    } Pointi3;

    typedef struct {
       Pointd3 p[8];         ///< The cube of points
       Pointd3 n[8];         ///< The cube of normals
       double val[8];        ///< The value of the points in the grid
    } GridCell;

    typedef struct {
       Pointd3 p[3];         ///< Vertices
    } Triangle;

    Pointd3 VertexInterp(const double &isolevel, const Pointd3 &p1, const Pointd3 &p2, const double &valp1, const double &valp2);

    /**
     * @brief The MarchingCubes class
     * This is an implementation of the marching cubes algorithm. It is a c++
     * rehash of the original c code by Paul Bourke which can be found here
     * http://astronomy.swin.edu.au/pbourke/modelling/polygonise/
     * There are some open mp directives to speed the process up.
     */
    class MarchingCubes
    {
    public:
        MarchingCubes( const double isoval , const double npoints, const double max[3], const double min[3]);
        ~MarchingCubes();

        const std::vector<Pointd3>& GetPoints() const;  ///< Return the point list
        const std::vector<Pointi3>& GetFaces() const;   ///< Return the face list
        void SetRenderToScreen(bool render);
        void EnableBoolFunc( const std::function<bool(const double*)>& fun );
            ///< Render the implicit surface if true

        void March(const std::function<double(const double*)> &fun); ///< Compute the marching cubes algorithm
        void MarchKD(const std::function<double(const double*)>& fun);
        void WriteToPly( const std::string& filename ); ///< Write the current mesh to file
    private:
        const double   m_isoval;      ///< The iso value
        const uint32_t m_npoints;     ///< The number of cubes per dimension

        double m_max[3]; ///< The maximum grid values
        double m_min[3]; ///< The minumum grid values

        std::vector<Pointd3> m_pointList;
        std::vector<Pointi3> m_faceList;

        int cubeindex, vi, j, i;
        int m_map[12]; ///< Temporary storage
        bool m_renderToScreen;

        double m_step[3];

        std::function<bool(const double*)> m_boolFunc; ///< The boolean function
        bool m_doboolfunc;

        void ComputeExtremalValues(); ///< Compute extreme values on the grid
        void PolygoniseCube(
                const GridCell &g,
                std::vector<Pointd3>::iterator& vertit,
                std::vector<Pointi3>::iterator &faceit,
                uint32_t& vertindex);
        void PolygoniseCube(const GridCell &g, std::vector<Pointd3> &pointList, std::vector<Pointi3> &faceList, uint32_t cindex);

        void KD( const GridCell& gc, uint32_t dim,
                 std::vector<Pointd3>::iterator& vertit,
                 std::vector<Pointi3>::iterator& faceit,
                 uint32_t& vertindex,
                 const std::function<double(const double*)>& fun);
    };

    void WritePly( std::string filename,
                   const std::vector<Pointd3>& vertices,
                   const std::vector<Pointi3>& faces);
}

#endif // MARCHINGCUBES_H
