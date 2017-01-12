#ifndef TREE_H
#define TREE_H

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
#include <stdint.h>
#include <map>
#include <vector>
#include <vector>
#include <numeric>
#include <iostream>
#include <random>
#include <algorithm>
#include <chrono>
#include "node.h"

namespace KD
{
    /**
     * @brief The Tree class
     * The KD::Tree is a 'no bells' implementation of the KD tree. It is not
     * randomised, only initialised on the first dimension. It does, therefore,
     * suffer the limitations. However, this implementation is fast and written
     * entirely in stl.
     */
    class Tree
    {
    public:
        Tree(); ///< Empty tree
        Tree(const double* vertices, const uint32_t dims, const uint32_t N, const bool randomize=false, uint32_t maxDepth=-1);
        ~Tree();

        Tree(const Tree& tree); ///< Copy constructor for the tree structure

        /**
         * @brief Search
         * Search the KD tree for K near neighbours
         * @param point The dims dimensional point to find the nearest neighbours of
         * @param K The number of near neighbours
         * @return A map of each node that has been checked (euclidean distance, index) pair
         * the K nearest neighbours are the first K entries in the map
         */
        std::vector<KeyPair> Search(const double* point, uint32_t K ) const;
        void Search( const double* point, uint32_t K, std::vector<KeyPair>& mp ) const;

        /**
         * @brief Histogram
         * Create a histogram for the depth limitied tree from input data.
         * @param inputdata A pointer to the column major, contiguous, data matrix
         * @param Npts The number of points in the input data
         * @return A map histogram
         */
        std::map<uint32_t, uint32_t> Histogram(const double* inputdata, const uint32_t Npts) const;

        uint32_t GetMaxDepth() const;
        uint32_t GetDataNumberOfPoints() const;

    private:
        std::shared_ptr<uint32_t> m_indecies; ///< A pointer to indexes in to the data matrix
        std::shared_ptr<Data> m_data;         ///< The data
        NodeData m_nd;        ///< The first node data
        std::shared_ptr<Node> m_headNode;     ///< The head node (points to all of the data).

        std::shared_ptr<uint32_t> m_ref;      ///< Number of references to the indecies data
    };

    class RandomisedTrees
    {
    public:
        RandomisedTrees(const double* vertices, const uint32_t dims, const uint32_t N, const uint32_t NumberOfTrees);
        std::vector<KeyPair> Search(const double* point, const uint32_t K ) const;

    private:
        std::vector< std::vector<double> > m_treedata;
        std::vector< Tree > m_trees;
        std::vector< std::vector<uint32_t> > m_permutes;
    };
}
#endif // TREE_H
