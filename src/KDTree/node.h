#ifndef NODE_H
#define NODE_H

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

#include <stdint.h>
#include <map>
#include <vector>
#include <memory>
#include "KeyPair.h"

/**
 *\brief A collection of classes, structs and functions to create a KD tree from data.
 *
 * Please see https://en.wikipedia.org/wiki/K-d_tree for more information about the implementation. A randomised
 * KD tree is also implemented in this namespace.
 */
namespace KD
{

/**
 * @brief print_array.
 * Print a double array to stdout, which contains a column major matrix.
 * @param array A pointer to the array
 * @param dims The number of dimensions
 * @param N The number of columns
 */
void print_array(const double *array, const uint32_t dims, const uint32_t N );

/**
 * @brief print_array.
 * Print an integer array to stdout, which contains a column major matrix.
 * @param array A pointer to the array
 * @param dims The number of dimensions
 * @param N The number of columns
 */
void print_array(const uint32_t *array, const uint32_t dims, const uint32_t N );

/**
 * @brief The NodeData struct POD.
 * All of the data required for the node. A copy of this struct is
 * contained within each node.
 */
struct NodeData
{
    uint32_t m_index;            ///< Index for the current node
    uint32_t m_currentDimension; ///< The current dimension, or level in the heirarchy
    uint32_t m_level;            ///< The current level of the node
    uint32_t m_min;              ///< The minimum index
    uint32_t m_max;              ///< The maximum index
    uint32_t m_maximumDepth;     ///< The maximum depth for the data
};

/**
 * @brief The Data struct.
 * A storage space for the data. This struct is
 * unique, and pointers to it are passed down the trees.
 */
struct Data
{
    const double* m_data;  ///< The data
    uint32_t* m_indecies;  ///< The indexes
    const uint32_t m_dims; ///< The number of dimensions
    const uint32_t m_N;    ///< the number of points
};

/**
 * @brief The Node class.
 * This is an abstract class for a node. It contains a pointer to another
 * node on the left and right, which are null if the node is either empty or a leaf.
 * To instantiate the entire heirarchy instantiate a single data node with a
 * reference to the data.
 */
class Node
{
public:
    /**
     * @brief The Type enum.
     * An enumeration of the possible types for a node (Leaf or Data).
     * The empty type was included to begin with, but became unecessary since smart pointers can be null.
     */
    enum Type
    {
        LeafType,
        EmptyType, // Not used in this implementation
        DataNodeType
    };

    /**
     * @brief Instantiate a single node
     * @param type The node type
     * @param data A reference to the data struct
     * @param nd A node data object.
     */
    Node(Type type, Data &data, const NodeData nd);

    /**
     * @brief ~Node.
     * The virtual destructor for the node
     */
    virtual ~Node(){}

    /**
     * @brief Node.
     * Copy this node and all of its subnodes. Calling this function
     * effectively copies the entire tree. The tree is stored as using the
     * shared pointer paradigm.
     * @param node
     */
    Node(const Node &node);

    /**
     * @brief Search.
     * Search the node for a particular point. Starting from the head, recurse down the
     * branches until a leaf is reached. The method is virtual since only a DataNode has
     * pointers to lower level nodes, the Leaf Nodes cannot search left or right.
     *
     * @param DistanceIndex A map containing every index that has been visited
     * @param point The dims dimensional point to search for
     * @param K The number of branches to search. If it is greater than zero the algorithm
     * will search K alternative branches from the bottom of the tree. To find one nearest neighbour
     * the two lowermost Leaf Nodes should be checked.
     */
    virtual void Search( std::vector<KeyPair>& DistanceIndex, const double* point, uint32_t& K ) const;

    Type GetType() const;                        ///< \brief Get the node type
    std::pair<std::shared_ptr<Node>, std::shared_ptr<Node> > GetLeftRight() const; ///< \brief Get left and right subnode pointers
    void setHeadNode(bool isheadnode);           ///< \brief Set the 'is head node' property.

protected:
    Type m_type;   ///< The type of node
    Data& m_data;  ///< A reference to the data matrices.
    NodeData m_nd; ///< The node data

    /**
     * @brief DistanceToPoint.
     * Return the distance from the node data element to the point
     * @param point A pointer to the point which should be compared
     * @return The Euclidean distance
     */
    double DistanceToPoint( const double* point ) const;

    /**
     * @brief NodeFactory.
     * Create a new node based on the node data. This method will select the type of node
     * depending on how much data is left, and either return a leaf or data node.
     * @param data A reference to the data store
     * @param nd The nodedata for level n of the tree
     * @param left If true, create the left hand node, else, create the right hand node
     * @return A node at level n+1
     */
    std::shared_ptr<Node> NodeFactory(Data &data, const NodeData &nd, bool left) const;

    uint32_t GetIndeciesIndex() const;   ///< \brief Get the nodes index in to the indecies array
    uint32_t GetDataColIndex() const;    ///< \brief Get the column index of the data point in the data array
    uint32_t GetDataLinearIndex() const; ///< \brief Get the data index

    std::shared_ptr<Node> m_left;  ///< A smart pointer to the left node.
    std::shared_ptr<Node> m_right; ///< A smart pointer the right node.

    bool m_isheadNode; ///< true if this is the head node.
};

/**
 * @brief The Leaf class.
 * The lowermost node in the tree. The pointers to left and right sub nodes
 * are null.
 */
class Leaf : public Node
{
public:
    /**
     * @brief Create a new leaf node
     * @param data A reference to the data struct
     * @param nd The node data
     */
    Leaf( Data& data, const NodeData nd);

    /**
     * @brief Leaf.
     * A copy constructor for the leaf node from a Node
     * @param node
     */
    Leaf( const Node& node );

    /**
     * @brief Leaf.
     * A copy constructor for the leaf node from a Leaf
     * @param node
     */
    Leaf( const Leaf& node );

    /**
      * @brief ~Leaf.
      * The destructor for the LeafNode.
      */
    ~Leaf(){}

    /**
     * @brief Search.
     * Search the leaf. This search method is trivial since we are at the bottom of the
     * tree. The recursion stops and so we add the index and distance to the DistanceIndex vector
     * and decrement K, if it is greater than zero.
     * @param DistanceIndex
     * @param point
     * @param K
     */
    void Search( std::vector<KeyPair>& DistanceIndex, const double* point, uint32_t& K) const;
};

/**
 * @brief The DataNode class.
 * This class represents a Node on the tree which is not a leaf. It, therefore, contains
 * pointers to a left and right node, both of which could contain a leaf node, data node,
 * or be null. The search method recurses down the tree.
 */
class DataNode : public Node
{
public:
    /**
     * @brief DataNode.
     * Construct a data node from a reference to the data store and a NodeData object.
     * This constructor is used when the tree is first created.
     * @param data A reference to the data store
     * @param nd The node data
     */
    DataNode( Data& data, const NodeData nd);

    /**
     * @brief DataNode.
     * Create a data node from a Node.
     * @param node
     */
    DataNode( const Node& node );

    /**
     * @brief DataNode.
     * Create a data node from a DataNode.
     * @param node
     */
    DataNode( const DataNode& node);
    ~DataNode();

    /**
     * @brief Search.
     * Search the tree for the data pointed to by point. The DistanceIndex vector is filled as
     * the tree is parsed. It may not contain unique values, since there is data on each of the
     * nodes, not just the leaves.
     * @param DistanceIndex
     * @param point
     * @param K
     */
    void Search( std::vector<KeyPair>& DistanceIndex, const double* point, uint32_t& K ) const;
};

/**
 * @brief The Empty class.
 * The empty class was originally designed so that the KD tree could be stored on the
 * stack. The latest implementation uses smart pointers, which can be null, and so it is not necessary any more.
 */
class Empty : public Node
{
public:
    Empty(Data &data, const NodeData nd);
    Empty(Data &data);
    ~Empty();
    void Search( std::vector<KeyPair>& DistanceIndex, const double* point, uint32_t& K ) const;
};

}

#endif // NODE_H
