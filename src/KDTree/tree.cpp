#include "tree.h"

namespace KD
{
    inline std::vector< uint32_t > RandomPermutation( const uint32_t N )
    {
        unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
        std::default_random_engine generator(seed1);
        std::uniform_int_distribution<uint32_t> distribution(0,N-1);

        std::vector< uint32_t > out(N);
        int i=0;
        for(auto& o : out)
        {
            o = i++;
        }
        for(auto& o : out)
        {
            std::swap(o, out[distribution(generator)]);
        }
        return out;
    }

    RandomisedTrees::RandomisedTrees(const double* vertices, const uint32_t dims, const uint32_t N, const uint32_t NumberOfTrees)
        : m_treedata(NumberOfTrees, std::vector<double>(dims*N)),
          m_trees(NumberOfTrees),
          m_permutes(NumberOfTrees)
    {

#pragma omp parallel for
        for( uint32_t j=0; j<NumberOfTrees; ++j)
        {
            m_permutes[j] = RandomPermutation(N);
        }

#pragma omp parallel for
        for(uint32_t j=0; j<NumberOfTrees;++j)
        {
            for(uint32_t pindex,k,i=0; i<N; i++)
            {
                pindex = m_permutes[j][i];
                for(k=0; k<dims; k++)
                {
                    m_treedata[j][i*dims+k] = vertices[pindex*dims+k];
                }
            }
            m_trees[j]=Tree(&m_treedata[j].front(), dims, N, true);
        }
    }

    std::vector<KeyPair> RandomisedTrees::Search(const double* point, const uint32_t K ) const
    {
        std::vector<KeyPair> out(K*m_trees.size());

#pragma omp parallel for
        for(uint32_t t =0; t<m_trees.size(); t++)
        {
            std::vector<KeyPair> conc;
            conc.reserve(4*K);

            m_trees[t].Search( point, K, conc);
            std::sort(conc.begin(), conc.end(), [&](const KeyPair& a, const KeyPair& b){return a.distance < b.distance;});
            for(uint32_t it=0;it<conc.size();++it)
            {
                conc[it].index = m_permutes[t][conc[it].index];
            }

            if(conc.size() > K)
            {
                std::swap_ranges(conc.begin(), conc.begin()+K, out.begin()+K*t);
            }
            else
            {
                std::swap_ranges(conc.begin(), conc.end(), out.begin()+K*t);
            }
            ///out.insert(out.begin()+K*t, conc.begin(), conc.begin()+K);
        }

        // Make the solution unique
        std::sort(out.begin(), out.end(), [&](const KeyPair& a, const KeyPair& b){return a.index < b.index;});
        std::stable_sort(out.begin(), out.end(), [&](const KeyPair& a, const KeyPair& b){return a.distance < b.distance;});
        std::vector<KeyPair> outunique(K);
        uint32_t ouback=std::numeric_limits<uint32_t>::max();
        uint32_t k=0;
        for(const auto& o : out)
        {
            if(o.index != ouback)
            {
                outunique[k++] = o;
                ouback = o.index;
                if(k>=K) break;
            }
        }
        return outunique;
    }

    Data DataFactory( const double* vertices, uint32_t* indecies, const uint32_t dims, const uint32_t N )
    {
       for(uint32_t i=0;i<N;i++)
       {
           indecies[i] = i;
       };

       return {vertices, indecies, dims, N};
    }

    Tree::Tree() :
        m_indecies(NULL),
        m_data(new Data({NULL, NULL, 0,0})),
        m_nd({0,0,0,0,0,static_cast<uint32_t>(-1)}),
        m_headNode(new Node(Node::LeafType, *m_data, m_nd))
    {

    }


    Tree::Tree(const double* vertices, const uint32_t dims, const uint32_t N, const bool randomize, uint32_t maxDepth) :
       m_indecies(new uint32_t[N]),
       m_data(new Data(DataFactory(vertices, m_indecies.get() , dims, N))),
       m_nd({(N)/2, 0, 0, 0, N-1, maxDepth}),
       m_headNode(NULL),
       m_ref(new uint32_t)
    {
        *m_ref = 1;
        if(randomize)
        {
            unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
            std::default_random_engine generator(seed1);
            std::uniform_int_distribution<int> distribution(0,dims-1);
            m_nd.m_currentDimension = distribution(generator);
        }
        m_headNode = std::shared_ptr<Node>(new DataNode(*m_data, m_nd));
        m_headNode->setHeadNode(true);
    }

    Tree::Tree(const Tree &tree):
              m_data(tree.m_data),
              m_nd(tree.m_nd),
              m_headNode(new DataNode(*tree.m_headNode)),
              m_ref(tree.m_ref)
    {
        (*m_ref)++; // The data matrix and indecies are ref counted since they do not change.
    }

    Tree::~Tree()
    {
    }


    void Tree::Search(const double* point, uint32_t K, std::vector<KeyPair> &mp ) const
    {
       m_headNode->Search( mp, point, K);
    }

    std::vector<KeyPair> Tree::Search( const double* point, uint32_t K ) const
    {
       std::vector<KeyPair> mp;
       mp.reserve(2*K);
       Search( point, K, mp);
       if(mp.size() > K)
       {
           std::partial_sort(mp.begin(), mp.begin()+K, mp.end(),
                             [&](const KeyPair& a, const KeyPair& b)
           {
               return a.distance < b.distance;
           });
           mp.resize(K);
       }
       else
       {
           std::sort(mp.begin(), mp.end(),
                     [&](const KeyPair& a, const KeyPair& b)
           {
               return a.distance < b.distance;
           });
       }
       return mp;
    }

    uint32_t Tree::GetMaxDepth() const
    {
        return m_nd.m_maximumDepth;
    }

    uint32_t Tree::GetDataNumberOfPoints() const
    {
        return m_nd.m_index*2;
    }

    std::map<uint32_t, uint32_t> Tree::Histogram(const double* inputdata, const uint32_t Npts) const
    {
        std::vector<KeyPair> mpout;
        std::map<uint32_t, uint32_t> histog;
        for(uint32_t i=0;i<Npts;++i)
        {
            mpout.clear();
            mpout.reserve(1);
            Search(&inputdata[i*m_data->m_dims], 1, mpout);
            histog[mpout.begin()->index]++;
        }
        return histog;
    }
}
