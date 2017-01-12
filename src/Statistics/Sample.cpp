
#include "Sample.h"

namespace Math
{
    Generator::Generator()
    {
        uint32_t seed = std::chrono::system_clock::now().time_since_epoch().count();
        m_e = std::default_random_engine(seed);
    }

    std::vector<uint32_t> Generator::SampleUniformIntegerUI(const uint32_t min, const uint32_t max, const uint32_t NSamples)
    {
        std::uniform_int_distribution<uint32_t> uniint(min,max);
        std::vector<uint32_t> samples(NSamples);
        for(auto& s : samples)
        {
            s = uniint(m_e);
        }
        return samples;
    }

    std::vector<int32_t> Generator::SampleUniformInteger(const int32_t min, const int32_t max, const uint32_t NSamples)
    {
        std::uniform_int_distribution<int32_t> uniint(min,max);
        std::vector<int32_t> samples(NSamples);
        for(auto& s : samples)
        {
            s = uniint(m_e);
        }
        return samples;
    }

    std::vector<uint32_t> Generator::SampleDiscrete1D( const std::vector<double>& dist, const uint32_t N )
    {
        std::discrete_distribution<int> discrete(dist.begin(), dist.end());
        std::vector<uint32_t> samples(N);

        for( auto& sample : samples )
        {
            sample = discrete(m_e);
        }
        return samples;
    }

    std::vector<uint32_t> Generator::SampleDiscrete1D( const Matrix& distribution, const uint32_t N )
    {
        /**
        std::vector<float> dist(distribution.cols());
        uint32_t index=0;
        for( auto& d : dist )
        {
            d = distribution(0, index++);
        }**/

        std::discrete_distribution<int32_t> discrete(distribution.GetArr().begin(), distribution.GetArr().end());
        std::vector<uint32_t> samples(N);

        for( auto& sample : samples )
        {
            sample = discrete(m_e);
        }
        return samples;
    }


    std::vector<std::pair<uint32_t, uint32_t>> Generator::SampleDiscrete2D( const Matrix& distribution, const uint32_t N )
    {
        // Sample the prior
        Matrix prior = distribution.Sum(0);
        std::vector<uint32_t> pss = SampleDiscrete1D( prior, N );

        std::vector<std::pair<uint32_t, uint32_t>> returnVals(N);
        auto it = returnVals.begin();
        for( auto p : pss )
        {
            int x = *(SampleDiscrete1D( distribution.Column( p ).transpose(), 1).begin());
            int y = *(SampleDiscrete1D( distribution.Row( x ), 1).begin());
            it->first = x;
            it->second = y;
            it++;
        }
        return returnVals;
    }

    Matrix Generator::SampleMultivariateNormal( const Matrix& mu, const Matrix& Sigma, const uint32_t N )
    {
        // Do cholesky decomposition
        Matrix llt = Sigma.chol();

        Matrix output( mu.Rows(), N);
        std::normal_distribution<double> distribution ( 0.0, 1.0);
        for( uint32_t i=0; i<output.numel(); i++)
        {
            output(i) = distribution(m_e);
        }
        output = (llt.transpose()*output);
        output.AddVectorToColumns( mu );

        return output;
    }

    Matrix Generator::SampleGMM( const Matrix& Mu, const std::vector<Matrix>& Sigma,
                                 const Matrix& Pi, const uint32_t N )
    {
        const std::vector<uint32_t> zs = SampleDiscrete1D( Pi, N );

        Matrix returnValues(Mu.Rows(), N);
        uint32_t index=0;
        for( auto z : zs )
        {
            returnValues.SetColumn(index++, SampleMultivariateNormal( Mu.Column(z), Sigma.at(z), 1 ));
        }
        return returnValues;
    }

    Matrix Generator::UniformlyRandomRotation(const uint32_t M)
    {
        std::normal_distribution<double> distribution ( 0.0, 1.0);
        Matrix norm(M,M);
        for(uint32_t k=0; k<norm.numel(); ++k)
        {
            norm(k) = distribution(m_e);
        }
        Matrix Q = norm.QR().first;
        if(Q.det()<0)
        {
            for(uint32_t k=0;k<M;++k)
            {
                Q(k,0) = -Q(k,0);
            }
        }
        return Q;
    }

    std::vector<uint32_t> Generator::SamplePermutation( const uint32_t N )
    {
        std::vector<uint32_t> ret(N);
        std::uniform_int_distribution<uint32_t> uniint(0,N-1);
        uint32_t index=0;
        for(auto& r : ret)
        {
            r = index++;
        }
        for(uint32_t k=0; k<N; ++k)
        {
            std::swap(ret[k], ret[uniint(m_e)]);
        }
        return ret;
    }

    Matrix Generator::UniquelySubSampleColumns( Matrix data, uint32_t N, double sigma )
    {
        if(N > data.Cols())
        {
            std::cerr << "Warning::UniquelySubSampleColumns:: N is too large" << std::endl;
            return data;
        }

        auto perm = SamplePermutation( data.Cols() );
        Matrix out(data.Rows(), N);

        std::normal_distribution<double> distribution ( 0.0, sigma);

        for(uint32_t k=0; k<N; ++k)
        for(uint32_t i=0; i<data.Rows();++i)
            out(i,k) = data(i,perm[k]) + distribution(m_e);

        return out;
    }

    std::vector<Math::Matrix> Generator::SampleWishart( const Math::Matrix& Sigma, const uint32_t df, const uint32_t N )
    {
       std::vector<Matrix> output(N, Matrix(Sigma.Rows(), Sigma.Cols()));
       Matrix Mu(Sigma.Rows(),1), nsamples(Sigma.Rows(), df);
       Mu.SetZero();

       for( uint32_t k=0; k<N; ++k )
       {
           nsamples = SampleMultivariateNormal( Mu, Sigma, df );
           output[k] = nsamples * nsamples.transpose();
       }

       return output;
    }
}
