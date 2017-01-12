
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
        case KDTree:
        case RandomisedMorton:
        case KDMap :
        case Direct:
            std::cerr << "NearestNeighbours: This implementation has been removed" << std::endl;
        break;
        default:
        {
            std::cerr << "NearestNeighbours: No such implementataion" << std::endl;
        }
    }
}

std::set<uint32_t> NearestNeighbours::SearchRadius(const double* point, const double radius, const uint32_t sparsity) const
{
    (void) point;
    (void) radius;
    (void) sparsity;

    if(m_algorithm==RandomisedMorton)
    {
        std::cerr << "Currently out of action" << std::endl;
        return std::set<uint32_t>();
        //return m_memstore.m_morton->SearchRadius(point, radius, sparsity);
    }
    else
    {
        std::cerr << "Algorthm not available" << std::endl;
        return std::set<uint32_t>();
    }
}


std::set<uint32_t> NearestNeighbours::NearestToPointCloud(const double* pts, const uint32_t N, const double thresh ) const
{
    (void) pts;
    (void) N;
    (void) thresh;

    switch(m_algorithm)
    {
        case RandomisedMorton:
        {
            std::cout << "Currently out of action " << std::endl;
            return std::set<uint32_t>();
        /*
            std::set<uint32_t> out;
#pragma omp parallel for
            for(uint32_t k=0; k<N;++k)
            {
                std::set<uint32_t> set = m_memstore.m_morton->SearchRadius(&pts[m_M*k], thresh);
                out.insert(set.begin(), set.end());
            }
            return out;*/
        } break;
    default:
        std::cout << "NearestNeighbours: No such algorithm " << std::endl;
        return std::set<uint32_t>();
    }
}

std::vector<KeyPair> NearestNeighbours::Search(const double *point, const uint32_t K ) const
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
    }
}
/*
std::pair<std::vector<double>, std::vector<uint32_t> > NearestNeighbours::GetNearestNeighbours( const double* point, const uint32_t K ) const
{
    std::vector<uint32_t> ret(K,0);
    std::vector<double> dist(K,0.0);

    std::map<double, uint32_t> out(Search(point, K));

    auto oit = out.begin();
    auto dit = dist.begin();

    uint32_t i=0;
    for(auto& k : ret)
    {
        if(++i > K) break;
        if(dit != dist.end())
            *(dit++) = oit->first;
        else
            *(dit++) = NAN;
        if(oit != out.end())
            k = (oit++)->second;
        else
            k = NAN;
    }
    return std::make_pair(dist,ret);
}*/

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
