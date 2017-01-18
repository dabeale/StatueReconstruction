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

/**
 * \brief A collection of classes, structs and functions for hte Marching Cubes algorithm.
 */
namespace MC
{
    /**
     * @brief A double precision 3 dimensional vector
     */
    struct Pointd3
    {
        double p[3];
    };

    /**
     * @brief An integral 3 dimensional vector
     */
    struct Pointi3
    {
        uint32_t p[3];
    };

    /**
     * @brief The GridCell struct.
     * A structure representing a cube in the
     * marchine cubes algorithm.
     */
    struct GridCell
    {
       Pointd3 p[8];         ///< The cube of points
       Pointd3 n[8];         ///< The cube of normals
       double val[8];        ///< The value of the points in the grid
    } ;

    /**
     * @brief A structure which represents a triangle, or 3 vertices.
     */
    struct Triangle
    {
       Pointd3 p[3];         ///< Vertices
    } ;

    /**
     * @brief VertexInterp.
     * Interpolate between the vertices according to the function values at each of the vertices
     * @param isolevel The isolevel, or threshold.
     * @param p1 The first vertex
     * @param p2 The second vertex
     * @param valp1 The function value at the first vertex
     * @param valp2 The function value at the second vertex
     * @return
     */
    Pointd3 VertexInterp(const double &isolevel, const Pointd3 &p1, const Pointd3 &p2, const double &valp1, const double &valp2);

    /**
     * @brief The MarchingCubes class.
     * This is an implementation of the marching cubes algorithm. It separates the 3D space in
     * to a collection of cubes and then turns each cube in to a surface segment.
     *
     * The algorithm is described in detail in the paper:
     *   "Marching cubes: A high resolution 3D surface construction algorithm" - Lorensen and Cline, ACM Siggraph, 1987.
     *
     * For a different, perhaps simpler, Marching Cubes implementation
     * take a look at http://astronomy.swin.edu.au/pbourke/modelling/polygonise/
     */
    class MarchingCubes
    {
    public:
        /**
         * @brief onstruct a marching cubes object.
         * @param isoval The threshold to apply to the input function
         * @param npoints The number of cubes per dimension
         * @param max An array specifying the maximum for each dimension
         * @param min An array specifying the minimum for each dimension
         */
        MarchingCubes( const double isoval , const double npoints, const double max[3], const double min[3]);
        ~MarchingCubes();

        const std::vector<Pointd3>& GetPoints() const;  ///< \brief Return the point list
        const std::vector<Pointi3>& GetFaces() const;   ///< \brief Return the face list

        /**
         * @brief SetRenderToScreen.
         * Set the render to screen flag. If this is set to true, then the algorithm will
         * copy data into open gl inline.
         * @param render True if the algorithm should render the result.
         */
        void SetRenderToScreen(bool render);

        /**
         * @brief EnableBoolFunc.
         * The bool functional returns true if a cube should even be considered for surface construction.
         * It provides a way of quickly rejecting cubes if there is another heuristic for deciding whether
         * the main function if crosses the ise value in the cube.
         * @param fun A functional which takes a pointer to a double precision array and returns a bool
         */
        void EnableBoolFunc( const std::function<bool(const double*)>& fun );

        /**
         * @brief March.
         * Compute the surface by checking each cube in the volume. This method has openmp directives
         * to parallelise the procedure.
         * @param fun The functional to evaluate at each point. It takes a pointer to a double precision array and returns a double.
         */
        void March(const std::function<double(const double*)> &fun);

        /**
         * @brief MarchKD.
         * A recursive method for the Marching cubes. It performs a coarse to fine strategy starting from
         * a maximum cube size. If all corners of the cube are below, or above, the isovalue, the cube is ignored;
         * otherwise the cube is divided and rechecked down to a minimum cube size.
         *
         * This method has openmp directives to parallelise the procedure.
         * @param fun The functional to evaluate at each point. It takes a pointer to a double precision array and returns a double.
         */
        void MarchKD(const std::function<double(const double*)>& fun);

        /**
         * @brief Write the current mesh to file
         * @param filename The math to the output file
         */
        void WriteToPly( const std::string& filename );
    private:
        const double   m_isoval;      ///< The iso value
        const uint32_t m_npoints;     ///< The number of cubes per dimension

        double m_max[3];              ///< The maximum grid values
        double m_min[3];              ///< The minumum grid values

        std::vector<Pointd3> m_pointList; ///< A vector of points on the surface
        std::vector<Pointi3> m_faceList;  ///< A vector of faces

        int cubeindex, vi, j, i;    ///< Temporary storage
        int m_map[12];              ///< Temporary storage
        bool m_renderToScreen;      ///< If true copy the vertices and faces in to opengl

        double m_step[3];           ///< The size of each cube

        std::function<bool(const double*)> m_boolFunc; ///< The boolean function
        bool m_doboolfunc;                             ///< If true, evaluate the bool func.

        /**
         * @brief PolygoniseCube.
         * Find the surface inside the cube. Append the vertices and faces to the list, and update
         * the vertex index accordingly. This method assumes that the vertices and faces are preallocated
         * @param g The current grid cell
         * @param vertit A reference to the vertex iterator
         * @param faceit A reference to the face iterator
         * @param vertindex The current vertex index
         */
        void PolygoniseCube(  const GridCell &g, std::vector<Pointd3>::iterator& vertit, std::vector<Pointi3>::iterator &faceit, uint32_t& vertindex);

        /**
         * @brief PolygoniseCube.
         * Find the surface inside the cube. Append the vertices and faces to the list, and update
         * the vertex index accordingly. This method does not assume that the vertex and face vectors are preallocated,
         * but they can be reserved.
         * @param g The current grid cell
         * @param vertit A reference to the vertex vector
         * @param faceit A reference to the face vector
         * @param cindex The current vertex index
         */
        void PolygoniseCube(const GridCell &g, std::vector<Pointd3> &pointList, std::vector<Pointi3> &faceList, uint32_t cindex);

        /**
         * @brief KD.
         * The method for the recursive KD implementation of marching cubes
         * @param gc the current grid cell.
         * @param dim the current dimension
         * @param vertit A reference to the vertex iterator
         * @param faceit A reference to the face iterator
         * @param vertindex A reference to the current vertex index
         * @param fun The functional to evaluate at each point. It takes a pointer to a double precision array and returns a double.
         */
        void KD( const GridCell& gc, uint32_t dim,
                 std::vector<Pointd3>::iterator& vertit,
                 std::vector<Pointi3>::iterator& faceit,
                 uint32_t& vertindex,
                 const std::function<double(const double*)>& fun);
    };

    /**
     * @brief Write a mesh to the stanford ply format.
     * @param filename The path to the output file
     * @param vertices A vector of vertices
     * @param faces A vector of triangle faces.
     */
    void WritePly( std::string filename,
                   const std::vector<Pointd3>& vertices,
                   const std::vector<Pointi3>& faces);
}

#endif // MARCHINGCUBES_H
