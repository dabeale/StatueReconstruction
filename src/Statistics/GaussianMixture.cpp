
#include "GaussianMixture.h"

namespace Math
{
GaussianMixture::GaussianMixture(Matrix mu, std::vector<Matrix> sigma, Matrix pi) :
    m_mu(mu), m_sigma(sigma), m_prec(m_sigma.size()), m_dets(m_sigma.size()), m_pi(pi), m_g(), m_M(m_mu.Rows()), m_K(m_mu.Cols())
{
    assert(m_K == m_sigma.size() && m_pi.Cols() == m_K);
    auto pit = m_prec.begin();
    auto dit = m_dets.begin();
    for( auto s : m_sigma)
    {
        assert(s.Rows() == m_M);
        *(pit++) = s.inv();
        *(dit++) = s.det();
    }
}

Matrix GaussianMixture::DrawSamples(const uint32_t N)
{
    return m_g.SampleGMM(m_mu, m_sigma, m_pi, N);
}

// This method is taken from the paper
//
Matrix GaussianMixture::FindModes()
{
    const uint32_t NSamples = m_K*10;

    // Draw N samples
    Matrix samples = DrawSamples(NSamples);

    // Find the most likely
    std::vector<double> lik = Evaluate(samples);
    std::vector<std::pair<uint32_t,double>> indexedLik(lik.size());
    for(uint32_t k=0;k<lik.size();++k)
    {
        indexedLik[k] = std::make_pair(k,lik[k]);
    }


    std::sort(indexedLik.begin(), indexedLik.end(),
              [&](const std::pair<uint32_t,double>& l, const std::pair<uint32_t,double>& r)
    {
        return l.second < r.second;
    });

    Matrix Modes(m_M,m_K);
    for(uint32_t k=0; k<m_K; k++)
    {
        Modes.SetColumn(k, samples.Column(indexedLik[k].first));
    }

    const uint32_t MaxIter = 50;
    const double tol = 1e10;

    // Iterate each of the modes using fixed point iterations until convergence
#pragma omp parallel for
    for(uint32_t k=0; k<m_K; ++k)
    {
        double err=std::numeric_limits<double>::max();
        Matrix prevval = Matrix(m_M, 1, std::numeric_limits<double>::max());
        uint32_t iter=0;
        Matrix currentx = Modes.Column(k);
        Matrix H(m_M,m_M), f(m_M,1);

        while(err > tol && iter++ < MaxIter)
        {
            H.SetZero();
            f.SetZero();
            for(uint32_t l=0;l<m_K;++l)
            {
                H+=ComponentProbability(currentx, l)*m_prec[l];
                f+=ComponentProbability(currentx, l)*m_prec[l]*m_mu.Column(l);
            }

            currentx = H.solve(f);
            err = (currentx - prevval).Norm();
            prevval = currentx;
        }
        Modes.SetColumn(k,currentx);
    }
    return Modes;
}

double GaussianMixture::ComponentProbability(const Matrix& x, const uint32_t i) const
{
    const double coef = (1.0 / std::sqrt(std::pow(2*M_PI,m_M)));
    Matrix r(m_M,1);
    Matrix c=m_mu.Column(i);
    Multiply( m_prec[i], ( x - c), r );
    double m=r.dot(x - c);
    return coef * (1.0/std::sqrt(m_dets[i])) * std::exp( -0.5*m );
}

std::vector<double> GaussianMixture::Evaluate(const Matrix& x) const
{
    std::vector<double> ret(x.Cols());
    const double coef = (1.0 / std::sqrt(std::pow(2*M_PI,m_M)));
    double innerprod, m;
    uint32_t i,k;
    Matrix r(m_M,1), c(m_M,1);
    for(i=0;i<x.Cols();++i)
    {
        innerprod=0;
        c = x.Column(i);
        for(k=0;k<m_K;++k)
        {
            Multiply( m_prec[k], ( c - m_mu.Column(k)), r );
            m=r.dot(c - m_mu.Column(k));
            innerprod += m_pi(k)* coef * (1.0/std::sqrt(m_dets[k])) * std::exp( -0.5*m );
        }
        ret[i] = innerprod;
    }
    return ret;
}
}
