#ifndef NEARESTNEIGHBOURS_H
#define NEARESTNEIGHBOURS_H


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
#include "tree.h"
#include <vector>
#include <set>
#include <cmath>

class NearestNeighbours
{
public:
    enum Algorithm
    {
        RandomTrees,
        KDTree,
        RandomisedMorton,
        KDMap,
        Direct
    };

    NearestNeighbours(const Algorithm& algorithm,
                       const double* data,
                       const uint32_t M,
                       const uint32_t N,
                       const uint32_t Permutations = 6,
                       const bool verbose = false );
    ~NearestNeighbours();

    std::set<uint32_t> SearchRadius(const double* point, const double radius, const uint32_t sparsity=1) const;
    std::set<uint32_t> NearestToPointCloud(const double* pts, const uint32_t N, const double thresh) const;

    std::vector<KeyPair> Search(const double *point, const uint32_t K ) const;
    //std::pair< std::vector<double>, std::vector<uint32_t> >
      //  GetNearestNeighbours(const double *point, const uint32_t K ) const;
private:
    struct DataStore
    {
        KD::RandomisedTrees* rtree;
        KD::Tree* tree;
    } m_memstore;

    NearestNeighbours( const NearestNeighbours& nn ); ///< Do not allow copying
    void operator=(const NearestNeighbours& nn); ///< Do not allow assignment

    Algorithm m_algorithm; ///< The type of algorithm to use
    const double* m_data;    ///< A pointer to the data (column major)

    const uint32_t m_M; ///< The dimension
    const uint32_t m_N; ///< The number of points
    const uint32_t m_perms; ///< The number of random variations of the data (randomised nns)

    const uint32_t m_verbose; ///< If true print output

    void PrintMessage(const std::string message , const uint32_t level=0) const;
};

#endif // NEARESTNEIGHBOURS_H
