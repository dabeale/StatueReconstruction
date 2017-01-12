#include "node.h"
#include <cmath>
#include <iostream>
#include <algorithm>


namespace KD
{

void print_array( const double* array, const uint32_t dims, const uint32_t N )
{
    for(uint32_t j=0; j<dims; j++)
    {
        for(uint32_t i=0; i<N; i++)
        {
            std::cout << array[ i*dims + j ] << ",";
        }
        std::cout << std::endl;
    }
}

void print_array(const uint32_t *array, const uint32_t dims, const uint32_t N )
{
    for(uint32_t j=0; j<dims; j++)
    {
        for(uint32_t i=0; i<N; i++)
        {
            std::cout << array[ i*dims + j ] << ",";
        }
        std::cout << std::endl;
    }
}

Node::Node(Type type,
      Data& data,
      const NodeData nd) :
         m_type(type),
         m_data(data),
         m_nd(nd),
         m_left(NULL),
         m_right(NULL),
         m_isheadNode(false)
{

}

Node::Node( const Node& node) :
    m_type(node.m_type),
    m_data(node.m_data),
    m_nd(node.m_nd),
    m_left(NULL),
    m_right(NULL),
    m_isheadNode(node.m_isheadNode)
{

}

inline Node* CreateNode(Node* node)
{
    switch( node->GetType() )
    {
        case Node::DataNodeType :
            return new DataNode(*node);
        case Node::LeafType :
            return new Leaf(*node);
        case Node::EmptyType :
            return new Node(*node);
    }
}

DataNode::DataNode( const Node& node ) :
    Node(node)
{
    auto nodes(node.GetLeftRight());
    if(node.GetType()==Node::DataNodeType) // This is correct since node should be a DataNode
    {
        m_left = std::shared_ptr<Node>(CreateNode(nodes.first.get()));
        m_right = std::shared_ptr<Node>(CreateNode(nodes.second.get()));
    }
    else // Bad - but no error thrown
    {
        std::cerr << "DataNode: Data node cannot be created from a leaf" << std::endl;
        m_left = std::shared_ptr<Node>(new Leaf(m_data, m_nd));
        m_right = std::shared_ptr<Node>(new Leaf(m_data, m_nd));
    }
}

DataNode::DataNode( const DataNode& node) :
    Node(node)
{
    auto nodes(node.GetLeftRight());
    m_left = std::shared_ptr<Node>(CreateNode(nodes.first.get()));
    m_right = std::shared_ptr<Node>(CreateNode(nodes.second.get()));
}


DataNode::~DataNode()
{
}


void Node::Search( std::vector<KeyPair>& DistanceIndex, const double* point, uint32_t& K ) const
{
   (void) DistanceIndex;
   (void) point;
   (void) K;
}

std::pair<std::shared_ptr<Node>,std::shared_ptr<Node>> Node::GetLeftRight() const
{
    return std::make_pair(m_left, m_right);
}

Node::Type Node::GetType() const
{
    return m_type;
}

void Node::setHeadNode(bool isheadnode)
{
    m_isheadNode = isheadnode;
}

double Node::DistanceToPoint( const double* point ) const
{
    double dist=0;
    double pt, dt;
    uint32_t dim = m_data.m_dims;
    for(int i=0; dim > 0 ; dim--)
    {
        i = GetDataColIndex()*m_data.m_dims+dim-1;
        pt = point[dim-1];
        dt = m_data.m_data[i];
        dist += (pt-dt)*(pt-dt);
    }
    return std::sqrt(dist);
}

uint32_t Node::GetIndeciesIndex() const
{
    return m_nd.m_index;
}

uint32_t Node::GetDataColIndex() const
{
    return m_data.m_indecies[GetIndeciesIndex()];
}

uint32_t Node::GetDataLinearIndex() const
{
    return m_data.m_dims*GetDataColIndex() + m_nd.m_currentDimension;
}

Empty::Empty(Data &data, const NodeData nd) :
        Node(Node::LeafType, data, nd)
{

}

Empty::Empty( Data & data) :
    Node(Node::EmptyType, data, NodeData())
{

}

Empty::~Empty()
{

}

Leaf::Leaf(Data &data, const NodeData nd) :
        Node(Node::LeafType, data, nd)
{
    // Sort the indexes based on the data and current dimension
    std::sort( static_cast<uint32_t*>(data.m_indecies + nd.m_min),
               static_cast<uint32_t*>(data.m_indecies + nd.m_max),
               [&](const uint32_t& a, const uint32_t& b){
                     return data.m_data[a*data.m_dims + nd.m_currentDimension] <
                            data.m_data[b*data.m_dims + nd.m_currentDimension];
               } );
}

Leaf::Leaf( const Node& node ) :
    Node(node)
{

}

Leaf::Leaf(const Leaf &node) :
    Node(node)
{
}


void Empty::Search(std::vector<KeyPair>& DistanceIndex, const double* point, uint32_t& K) const
{
    (void) DistanceIndex;
    (void) point;
    (void) K;
}

void Leaf::Search( std::vector<KeyPair>& DistanceIndex, const double* point, uint32_t& K) const
{
   DistanceIndex.push_back( KeyPair(DistanceToPoint(point), GetDataColIndex()) );
   if( K > 0 ) --K;
}

std::shared_ptr<Node> Node::NodeFactory(Data &data, const NodeData &nd, bool left) const
{
   // Create child nodes
   uint32_t newindex, newmax, newmin;

   if(left)
   {
       newmax = nd.m_index;
       newmin = nd.m_min;
       newindex = (newmin+newmax)/2;

       NodeData ndi({
           newindex,
           (nd.m_currentDimension+1)%data.m_dims,
           nd.m_level+1,
           newmin,
           newmax,
           nd.m_maximumDepth});

       if( newindex == newmax )
       {
           return std::shared_ptr<Node>(new Empty(data, ndi));
       }

       if( newmax - newmin > 1 && nd.m_level + 1 < nd.m_maximumDepth)
       {
           return std::shared_ptr<Node>(new DataNode(data, ndi));
       }
       else
       {
           newindex = newmin;
           return std::shared_ptr<Node>(new Leaf(data, ndi));
       }
   }
   else
   {
       newmax = nd.m_max;
       newmin = nd.m_index+1;
       newindex = (newmin+newmax)/2;

       NodeData ndi({
           newindex,
           (nd.m_currentDimension+1)%data.m_dims,
           nd.m_level+1,
           newmin,
           newmax,
           nd.m_maximumDepth});

       if( newindex == newmax )
       {
           return std::shared_ptr<Node>(new Empty(data, ndi));
       }

       if(newmax - newmin > 1 && nd.m_level + 1 < nd.m_maximumDepth )
       {
           return std::shared_ptr<Node>(new DataNode(data, ndi));
       }
       else
       {
           newindex = newmin;
           return std::shared_ptr<Node>(new Leaf(data, ndi));
       }
   }
}

DataNode::DataNode(Data &data,
        const NodeData nd) :
            Node(Node::DataNodeType, data, nd)
{
    // Sort the indexes based on the data and current dimension
    std::sort( static_cast<uint32_t*>(data.m_indecies + nd.m_min),
               static_cast<uint32_t*>(data.m_indecies + nd.m_max),
               [&](const uint32_t& a, const uint32_t& b){
                     return data.m_data[a*data.m_dims + nd.m_currentDimension] <
                            data.m_data[b*data.m_dims + nd.m_currentDimension];
               } );
/*
    for(uint32_t i=nd.m_min; i<nd.m_max; ++i)
    {
        std::cout << data.m_indecies[i] << ",";
    }
    std::cout << std::endl;*/

    m_left = Node::NodeFactory(m_data, m_nd, true);
    m_right = Node::NodeFactory(m_data, m_nd, false);
}

void DataNode::Search( std::vector<KeyPair>& DistanceIndex, const double* point, uint32_t& K ) const
{
    if(!m_isheadNode)
    {
        double dist = DistanceToPoint(point);
        uint32_t colindex = GetDataColIndex();
        DistanceIndex.push_back(KeyPair(dist, colindex));
    }

   bool goLeft = point[m_nd.m_currentDimension] < m_data.m_data[ GetDataLinearIndex() ];

   goLeft ? m_left->Search( DistanceIndex, point, K ) : m_right->Search( DistanceIndex, point, K );

   if( K > 0 )
   {
       --K;
       if(goLeft)
       {
           m_right->Search( DistanceIndex, point, K );
       }
       else
       {
            m_left->Search( DistanceIndex, point, K );
       }
   }
}
}
