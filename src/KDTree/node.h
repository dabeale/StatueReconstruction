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

namespace KD
{
void print_array(const double *array, const uint32_t dims, const uint32_t N );
void print_array(const uint32_t *array, const uint32_t dims, const uint32_t N );

/**
 * @brief The NodeData struct POD
 * All of the data required for the node
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

struct Data
{
    const double* m_data;  ///< The data
    uint32_t* m_indecies;  ///< The indexes
    const uint32_t m_dims; ///< The number of dimensions
    const uint32_t m_N;    ///< the number of points
};

/**
 * @brief The Node class
 * This is an abstract class for a node. It contains a pointer to another
 * node on the left and right, which are null if the node is either empty or a leaf.
 * To instantiate the entire heirarchy instantiate a single data node with a
 * reference to the data.
 */
class Node
{
public:
    enum Type
    {
        LeafType,
        EmptyType, // Not used in this implementation
        DataNodeType
    };

    Node(Type type, Data &data, const NodeData nd); ///< Instantiate a node
    virtual ~Node(){} ///< Destructor for the node (virtual)
    Node(const Node &node); ///< Copy node and all subnodes.

    /**
     * @brief Search
     * Search the node for a particular point. Starting from the head, recurse down the
     * branches until a leaf is reached. The method is virtual since only a DataNode has
     * pointers to lower level nodes, the Leaf Nodes cannot search left or right.
     * @param DistanceIndex A map containing every index that has been visited
     * @param point The dims dimensional point to search for
     * @param K The number of branches to search. If it is greater than zero the algorithm
     * will search K alternative branches from the bottom of the tree. To find one nearest neighbour
     * the two lowermost Leaf Nodes should be checked.
     */
    virtual void Search( std::vector<KeyPair>& DistanceIndex, const double* point, uint32_t& K ) const;

    Type GetType() const;                        ///< Get the node type
    std::pair<std::shared_ptr<Node>, std::shared_ptr<Node> > GetLeftRight() const; ///< Get left and right subnode pointers
    void setHeadNode(bool isheadnode); ///< Set the 'is head node' property.

protected:
    Type m_type;   ///< The type of node
    Data& m_data;  ///< A reference to the data matrices.
    NodeData m_nd; ///< The node data

    double DistanceToPoint( const double* point ) const;                 ///< Return the distance from the node data element to the point
    std::shared_ptr<Node> NodeFactory(Data &data, const NodeData &nd, bool left) const; ///< Create a new node based on the node data - method is used by the constructor

    uint32_t GetIndeciesIndex() const;   ///< Get the nodes index in to the indecies array
    uint32_t GetDataColIndex() const;    ///< Get the column index of the data point in the data array
    uint32_t GetDataLinearIndex() const; ///< Get the data index

    std::shared_ptr<Node> m_left;  ///< Contain the left node.
    std::shared_ptr<Node> m_right; ///< Contain the right node.

    bool m_isheadNode; ///< true if this is the head node.
};

/**
 * @brief The Leaf class
 * The lowermost node in the tree. The pointers to left and right sub nodes
 * are null.
 */
class Leaf : public Node
{
public:
    Leaf( Data& data, const NodeData nd);
    Leaf( const Node& node );
    Leaf( const Leaf& node );
    ~Leaf(){}
    void Search( std::vector<KeyPair>& DistanceIndex, const double* point, uint32_t& K) const;
};

/**
 * @brief The DataNode class
 * A DataNode contains pointers to the left and right sub nodes.
 */
class DataNode : public Node
{
public:
    DataNode( Data& data, const NodeData nd);
    DataNode( const Node& node );     ///< Create a data node from super Node. If it is a leaf it
    DataNode( const DataNode& node);
    ~DataNode();

    /**
     * @brief Search
     * Search the tree for the data pointed to by point. The DistanceIndex vector is filled as
     * the tree is parsed. It may not contain unique values, since there is data on each of the
     * nodes, not just the leaves.
     * @param DistanceIndex
     * @param point
     * @param K
     */
    void Search( std::vector<KeyPair>& DistanceIndex, const double* point, uint32_t& K ) const;
};

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
