#include "Distributions.h"

namespace Distribution
{
void NormaliseLogarithmicMatrix(double* R, double *sum, const uint32_t K, const uint32_t N)
{
    uint32_t i,j;
    double *maxes = new double[N];
    double apmin = std::numeric_limits<double>::min();
    for(i=0;i<N;++i)
    {
        maxes[i] = apmin;
        for(j=0;j<K;++j)
        {
            maxes[i]= R[K*i + j] > maxes[i] ? R[K*i + j] : maxes[i];
        }
        sum[i]=0.0;
    }

    for(i=0;i<N;++i)
    for(j=0;j<K;++j)
    {
        R[K*i + j] = std::exp( R[K*i + j] - maxes[i] );
        sum[i] += R[K*i + j];
    }

    for(i=0;i<N;++i)
    for(j=0;j<K;++j)
    {
        if(sum[i]>0)
            R[K*i + j] /= sum[i] ;
    }

    delete [] maxes;
}

void distance2(const double* data1, const double* data2, double* D, const uint32_t N1, const uint32_t N2, const uint32_t K)
{
    uint32_t k,j,i;
    for(i=0;i<N1; ++i)
    {
        for(j=0;j<N2;++j)
        {
            D[N1*j + i]=0.0;
            for(k=0;k<K;++k)
            {
                D[N1*j + i]+=(data1[K*i + k] - data2[K*j+k])*(data1[K*i + k] - data2[K*j+k]);
            }
        }
    }
}

std::vector< double > MvnPdf( const Cu::Matrix& data, const Cu::Matrix& Mu, const Cu::Matrix& Sigma  )
{
    if(data.rows() != Mu.rows() || data.rows() != Sigma.rows() || Sigma.rows() != Sigma.cols())
    {
        std::cerr << "Distribution::MvnPdf:: The dimensions are incorrect" << std::endl;
        assert(0);
    }

    const uint32_t K = data.rows();
    const uint32_t N = data.cols();

    Eigen::LLT<Cu::Matrix> llt(Sigma);
    Cu::Matrix L = llt.matrixL();

    double proddiag=1.0;
    for(uint32_t k=0; k<K; ++k)
        proddiag *= L(k,k);

    Cu::Matrix reorient(K, 1);

    std::vector<double> p(N, 0);
    if(proddiag == 0)
    {
        std::cerr << "Distribution::MvnPdf:: warning - Cholesky returned a zero diagonal" << std::endl;
        return p;
    }

    double coeff = std::pow(2*M_PI, static_cast<double>(K)/2.0);
    coeff = 1.0/(coeff*proddiag);

    double sum;
    for(uint32_t k=0; k<N; ++k)
    {
        reorient = L.fullPivHouseholderQr().solve(data.block(0, k, K, 1) - Mu);
        sum=0.0;
        for(uint32_t i=0; i<K; ++i)
            sum += reorient(i)*reorient(i);

        p[k] = coeff * std::exp( -sum/2.0 );
    }

    return p;
}

std::vector< double > StudentT3(const Cu::Matrix& data, const Cu::Matrix& Sigma, const Cu::Matrix& Mu, const uint32_t df )
{
    const uint32_t N = data.cols();
    const uint32_t K = data.rows();

    if( Sigma.rows() != K || Sigma.cols() != K || Mu.rows() != K || Mu.cols() != 1 )
    {
        std::cerr << "Distribution::StudentT3:: the dimensions are incorrect" <<std::endl;
        assert(0);
    }

    Eigen::LLT<Cu::Matrix> llt(Sigma);
    Cu::Matrix L = llt.matrixL();
    double logsum=0.0;
    for(uint32_t i=0; i<K; ++i)
        logsum += std::log(L(i,i));

    double coef = std::lgamma( (K+df)/2.0 ) - (std::lgamma(df/2.0) + logsum + (K/2.0)*std::log(df*M_PI));

    std::vector<double> pout(N, 0.0);
    Cu::Matrix reorient(K,1);
    double sum;

    for(uint32_t k=0; k<N; ++k)
    {
        reorient = L.lu().solve(data.block(0,k,K,1) - Mu);
        sum = 0.0;
        for(uint32_t i=0; i<K; ++i)
            sum += reorient(i)*reorient(i);
        sum /= df;
        pout[k] = std::exp( coef - ((df+K)/2.0)*std::log(1.0 + sum) );
    }

    return pout;
}


std::vector< double > LatentGaussianWishart(const Cu::Matrix& data, const uint32_t df, const double k, const Cu::Matrix& Mu, const Cu::Matrix& Precision)
{
    const uint32_t K = data.rows();

    if(Mu.rows() != K || Precision.rows() != K || Precision.cols() != K || Mu.cols() != 1)
    {
        std::cerr << "Distribution::LatentGaussianWishart:: Dimensions do not match" << std::endl;
        assert(0);
    }

    Cu::Matrix Sigma = (static_cast<double>(k*(df-K+1)) / (k + 1)) * Precision;
    Sigma = Sigma.inverse();
    return StudentT3(data, Sigma, Mu, df-K+1);
}
}
