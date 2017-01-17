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

/**
 * @brief The NearestNeighbours class
 * The nearest neighbour class was set up as an interface for a number of different nearest neighbour algorithms.
 * In this implementation, only the KD and Randomised KD tree are used, the other methods will not work.
 */
class NearestNeighbours
{
public:
    /**
     * @brief The Algorithm enum
     * An enumeration of possible nearest neighbour algorithms.
     */
    enum Algorithm
    {
        RandomTrees
    };

    /**
     * @brief NearestNeighbours
     * A constructor for the nearest neighbour algorithm
     * @param algorithm The nearest neighbour algorithm to use
     * @param data A pointer to the column wise data matrix
     * @param M The number of rows
     * @param N The number of columns
     * @param Permutations The number of trees to use
     * @param verbose If true, produce verbose output
     */
    NearestNeighbours(const Algorithm& algorithm,
                       const double* data,
                       const uint32_t M,
                       const uint32_t N,
                       const uint32_t Permutations = 6,
                       const bool verbose = false );
    ~NearestNeighbours(); ///< An empty destructor

    /**
     * @brief Search
     * Search for the nearest neighbours to a given point.
     * @param point A pointer to the M dimensional point
     * @param K The number of neighbours to find
     * @return A vector of KeyPairs, length K, sorted by the distance to the point.
     */
    std::vector<KD::KeyPair> Search(const double *point, const uint32_t K ) const;
private:

    /**
     * @brief The DataStore struct
     * A pointer to the KD::RandomisedTrees object
     */
    struct DataStore
    {
        KD::RandomisedTrees* rtree; ///< A pointer to the random trees object
    } m_memstore;

    NearestNeighbours( const NearestNeighbours& nn );   ///< \brief A private copy constructor to prevent copying
    void operator=(const NearestNeighbours& nn);        ///< \brief A private assignment operator to prevent copyingD

    Algorithm m_algorithm;                              ///< The type of algorithm to use. This can only be RandomisedTrees
    const double* m_data;                               ///< A pointer to the data ( M*N column major)

    const uint32_t m_M;                                 ///< The dimension or number of rows in the data matrix
    const uint32_t m_N;                                 ///< The number of points
    const uint32_t m_perms;                             ///< The number of random variations of the data (randomised nns). The number of trees.

    const uint32_t m_verbose;                           ///< If true print output

    /**
     * @brief PrintMessage
     * Print a message to stdout preceded by arrows of length specified by the input variable 'level'.
     * For example
     * \code
     *     PrintMessage("Hello, World!", 2)
     * \endcode
     * produces the output
     *  --> Hello, World!
     * @param message The message to be printed
     * @param level The level at which to print it, or the length of the arrow.
     */
    void PrintMessage(const std::string message , const uint32_t level=0) const;
};

#endif // NEARESTNEIGHBOURS_H
