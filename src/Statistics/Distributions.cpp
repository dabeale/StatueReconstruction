#include "Distributions.h"

namespace Distribution
{
void NormaliseLogarithmicMatrix(double* R, double *sum, const uint32_t K, const uint32_t N)
{
    uint32_t i,j;
    double *maxes = new double[N];
    double apmin = -std::numeric_limits<double>::max();
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

Math::Matrix NormaliseLogarithmicMatrix(const Math::Matrix& R, std::vector<double>& sum)
{

    uint32_t i,k;
    uint32_t N = R.Cols();
    uint32_t K = R.Rows();

    Math::Matrix out(K,N);

    std::vector<double> maxes(N);

    double apmin = -std::numeric_limits<double>::max();
    for(i=0;i<N;++i)
    {
        maxes[i] = apmin;
        for(k=0;k<K;++k)
        {
            if(R(k,i) > maxes[i])
                maxes[i]= R(k,i);
        }
        sum[i]=0.0;
    }

    for(i=0;i<N;++i)
    for(k=0;k<K;++k)
    {
        out(k,i) = std::exp( R(k,i) - maxes[i] );
        sum[i] += out(k,i);
    }

    for(i=0;i<N;++i)
    for(k=0;k<K;++k)
    {
        if(sum[i]>0)
            out(k,i) /= sum[i] ;
    }

    return out;
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

std::vector< double > LogMvnPdf( const Math::Matrix& data, const Math::Matrix& Mu, const Math::Matrix& Sigma  )
{
    if(data.Rows() != Mu.Rows() || data.Rows() != Sigma.Rows() || Sigma.Rows() != Sigma.Cols())
    {
        std::cerr << "Distribution::LogMvnPdf:: The dimensions are incorrect" << std::endl;
        assert(0);
    }

    const uint32_t K = data.Rows();
    const uint32_t N = data.Cols();

    Math::Matrix L = Sigma.chol().transpose();

    double proddiag=1.0;
    for(uint32_t k=0; k<K; ++k)
        proddiag *= L(k,k);


    std::vector<double> p(N, 0);
    if(proddiag == 0)
    {
        std::cerr << "Distribution::LogMvnPdf:: warning - Cholesky returned a zero diagonal" << std::endl;
        return p;
    }

    const double coeff = -(static_cast<double>(K)/2.0)*std::log(2*M_PI) - std::log(proddiag);

#pragma omp parallel for
    for(uint32_t k=0; k<N; ++k)
    {
        Math::Matrix reorient = L.solve(data.Column(k) - Mu);
        double sum=0.0;
        for(uint32_t i=0; i<K; ++i)
            sum += reorient(i)*reorient(i);

        p[k] = coeff - (sum/2.0);
    }

    return p;
}

std::vector< double > MvnPdf( const Math::Matrix& data, const Math::Matrix& Mu, const Math::Matrix& Sigma  )
{
    if(data.Rows() != Mu.Rows() || data.Rows() != Sigma.Rows() || Sigma.Rows() != Sigma.Cols())
    {
        std::cerr << "Distribution::MvnPdf:: The dimensions are incorrect" << std::endl;
        assert(0);
    }

    const uint32_t K = data.Rows();
    const uint32_t N = data.Cols();

    Math::Matrix L = Sigma.chol().transpose();

    double proddiag=1.0;
    for(uint32_t k=0; k<K; ++k)
        proddiag *= L(k,k);

    std::vector<double> p(N, 0);

    if(Sigma.det() < 1e-16)
    {
        std::cerr << "Distribution::MvnPdf:: warning - Sigma is ill conditioned" << std::endl;
        return p;
    }

    double coeff = std::pow(2*M_PI, static_cast<double>(K)/2.0);
    coeff = 1.0/(coeff*proddiag);

#pragma omp parallel for
    for(uint32_t k=0; k<N; ++k)
    {
        Math::Matrix reorient = L.solve(data.Column(k) - Mu);
        double sum=0.0;
        for(uint32_t i=0; i<K; ++i)
            sum += reorient(i)*reorient(i);

        p[k] = coeff * std::exp( -sum/2.0 );
    }

    return p;
}

std::vector< double > StudentT3(const Math::Matrix& data, const Math::Matrix& Sigma, const Math::Matrix& Mu, const uint32_t df )
{
    const uint32_t N = data.Cols();
    const uint32_t K = data.Rows();

    if( Sigma.Rows() != K || Sigma.Cols() != K || Mu.Rows() != K || Mu.Cols() != 1 )
    {
        std::cerr << "Distribution::StudentT3:: the dimensions are incorrect" <<std::endl;
        assert(0);
    }

    Math::Matrix L = Sigma.chol().transpose();
    double logsum=0.0;
    for(uint32_t i=0; i<K; ++i)
        logsum += std::log(L(i,i));

    double coef = std::lgamma( (K+df)/2.0 ) - (std::lgamma(df/2.0) + logsum + (K/2.0)*std::log(df*M_PI));

    std::vector<double> pout(N, 0.0);


#pragma omp parallel for
    for(uint32_t k=0; k<N; ++k)
    {
        Math::Matrix cc = data.Column(k) - Mu;
        Math::Matrix reorient = L.solve(cc);
        double sum = 0.0;
        for(uint32_t i=0; i<K; ++i)
            sum += reorient(i)*reorient(i);
        sum /= df;
        pout[k] = std::exp( coef - ((df+K)/2.0)*std::log(1.0 + sum) );
    }

    return pout;
}


std::vector< double > LatentGaussianWishart(const Math::Matrix& data, const uint32_t df, const double k, const Math::Matrix& Mu, const Math::Matrix& Precision)
{
    const uint32_t K = data.Rows();

    if(Mu.Rows() != K || Precision.Rows() != K || Precision.Cols() != K || Mu.Cols() != 1)
    {
        std::cerr << "Distribution::LatentGaussianWishart:: Dimensions do not match" << std::endl;
        assert(0);
    }

    Math::Matrix Sigma = (static_cast<double>(k*(df-K+1)) / (k + 1)) * Precision;
    Sigma = Sigma.inv();
    return StudentT3(data, Sigma, Mu, df-K+1);
}
}
