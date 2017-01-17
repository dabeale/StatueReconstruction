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
     * This class is an implementation of the KD Tree. It is the simplest possible
     * implementation there is no balancing or statistics involved, making it easy
     * to modify. It was first designed so that a functional could be passed down the tree.
     */
    class Tree
    {
    public:
        /**
         * @brief Tree
         * Construct an empty tree. Implementing this constructor allows stl containers to be used.
         */
        Tree();

        /**
         * @brief Tree
         * Construct a tree with some data. The tree can be chosen as random, for the randomised KD tree implementation.
         * Setting this flag to be true randomises the first search dimension of the tree.
         * @param vertices A pointer to the (column major) data matrix
         * @param dims The number of dimensions of the matrix
         * @param N The number of columns
         * @param randomize If true, randomise the first search dimension.
         * @param maxDepth The maximum depth of the tree
         */
        Tree(const double* vertices, const uint32_t dims, const uint32_t N, const bool randomize=false, uint32_t maxDepth=-1);

        /**
          * \brief ~Tree
          * Destroy the tree
          */
        ~Tree();

        /**
         * @brief Tree
         * Copy constructor for the tree structure. The latest implementation
         * uses smart pointers for each of the nodes and so the copy constructor is trivial. Care should
         * be taken though since it is still copies of pointers.
         * @param tree
         */
        Tree(const Tree& tree);

        /**
         * @brief Search
         * Search the KD tree for K near neighbours. This implementation recurses to the bottom
         * of the tree, since data is stored at each level and it is possible that the lower levels contain a closer point.
         * This means that the algorithm may pass more than K nodes, which are all appended to the KeyPair vector. Specifying a
         * value for K ensures that, at least, K points are returned. It is up to the programmer to truncate the vector if necessary.
         *
         * NOTE: This method will recurse the nodes and sort the returned vector, meaning it -should- be used to find the nearest neighbours.
         * @param point The dims dimensional point to find the nearest neighbours of
         * @param K The minimum number of near neighbours to return
         * @return A vector of KeyPairs.
         */
        std::vector<KeyPair> Search(const double* point, uint32_t K ) const;

        /**
         * @brief Search
         * This search method simply passes a vector of KeyPairs to the head node. The vector is
         * filled with at least K KeyPairs. The vector is not truncated or sorted.
         * @param point A pointer to a point with the same number of dimensions as the input data
         * @param K The number of neighbours to find
         * @param mp A reference to a vector of KeyPairs.
         */
        void Search( const double* point, uint32_t K, std::vector<KeyPair>& mp ) const;

        /**
         * @brief Histogram
         * Create a histogram for the depth limitied tree from input data.
         * @param inputdata A pointer to the column major, contiguous, data matrix
         * @param Npts The number of points in the input data
         * @return A map histogram
         */
        std::map<uint32_t, uint32_t> Histogram(const double* inputdata, const uint32_t Npts) const;

        uint32_t GetMaxDepth() const;              ///< \brief Return the maximum depth of the tree (set on initialisation)
        uint32_t GetDataNumberOfPoints() const;    ///< \brief Get the number of datapoints (set on initialisation).

    private:
        std::shared_ptr<uint32_t> m_indecies; ///< A pointer to indexes in to the data matrix
        std::shared_ptr<Data> m_data;         ///< The data
        NodeData m_nd;                          ///< The first node data
        std::shared_ptr<Node> m_headNode;     ///< The head node (points to all of the data).

        std::shared_ptr<uint32_t> m_ref;      ///< Number of references to the indecies data
    };

    /**
     * @brief The RandomisedTrees class
     * This class contains a collection of randomly initialised KD trees.
     *
     * In the KD tree implementation
     * each dimension is considered independently, and so points which are close can end up on different branches of the
     * tree. This makes a KD tree an approximate nearest neighbour algorithm. Randomly selecting the
     * intial search dimension and permuting the input data on multiple trees, and then searching each of them
     * lowers the probability of getting an incorrect selection of neighbours.
     */
    class RandomisedTrees
    {
    public:
        /**
         * @brief RandomisedTrees
         * Construct the randomised trees. The constructor will create an array of indexes
         * for the permutations and then create a KD tree for each of them
         * @param vertices A pointer to the data matrix
         * @param dims The number of dimensions
         * @param N The number of columns
         * @param NumberOfTrees The number of trees
         */
        RandomisedTrees(const double* vertices, const uint32_t dims, const uint32_t N, const uint32_t NumberOfTrees);

        /**
         * @brief Search
         * Search the randomised trees. This method searches each tree and then sorts the
         * vector of results based on their distance from the intput point. The result is truncated
         * so that only K KeyPairs are returned.
         * @param point The point to find the nearest neighbours of
         * @param K The number of neighbours
         * @return A vector of K KeyPairs
         */
        std::vector<KeyPair> Search(const double* point, const uint32_t K ) const;

    private:
        std::vector< std::vector<double> > m_treedata;    ///< The permuted data matrices
        std::vector< Tree > m_trees;                      ///< The randomised KD trees
        std::vector< std::vector<uint32_t> > m_permutes;  ///< The indexes from the m_treedatas in to the original data matrix
    };
}
#endif // TREE_H
