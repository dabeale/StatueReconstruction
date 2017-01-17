
#include "NearestNeighbours.h"

#include <iostream>

void NearestNeighbours::PrintMessage( const std::string message, const uint32_t level ) const
{
    if(m_verbose)
    {
        for(uint32_t k=0; k<level;k++)
        {
            std::cout << "--";
        }
        if(level>0)
        {
            std::cout << ">";
        }
        std::cout << " " << message << std::endl;
    }
}

NearestNeighbours::NearestNeighbours(const Algorithm& algorithm , const double *data,
                                     const uint32_t M,
                                     const uint32_t N,
                                     const uint32_t Permutations,
                                     const bool verbose) :
    m_algorithm(algorithm),
    m_data(data),
    m_M(M),
    m_N(N),
    m_perms(Permutations),
    m_verbose(verbose)
{
    switch(m_algorithm)
    {
        case RandomTrees:
        {
            PrintMessage("Algorithm: Random Trees");
            m_memstore.rtree = new KD::RandomisedTrees(data, M, N, m_perms);
        } break;
        default:
        {
            std::cerr << "NearestNeighbours: No such implementataion" << std::endl;
        }
    }
}

std::vector<KD::KeyPair> NearestNeighbours::Search(const double *point, const uint32_t K ) const
{
    switch(m_algorithm)
    {
        case RandomTrees:
        {
            auto out = m_memstore.rtree->Search(point, K);
            return out;
        } break;
    default:
        std::cout << "NearestNeighbours : No such algorithm" << std::endl;
        return std::vector<KD::KeyPair>();
    }
}

NearestNeighbours::~NearestNeighbours()
{
    switch(m_algorithm)
    {
        case RandomTrees:
        {
            delete m_memstore.rtree;
        } break;
    default :
        break;
    }
}
