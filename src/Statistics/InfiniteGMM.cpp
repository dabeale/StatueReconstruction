#include "InfiniteGMM.h"

namespace Math
{
    IGMM::IGMM( const Math::Matrix& data, const DataStore& datastore, const IGMMHyperparameters& hypes ) :
        m_hypes(hypes),
        m_data(data),
        m_M(data.Rows()),
        m_N(data.Cols()),
        m_P( 1, std::vector< double >(m_N, 1.0) ),
        m_PSum(m_N, 0.0),
        m_ProbAny( m_N, 0),
        m_BadClass( 1, false),
        m_current(datastore),
        m_best(datastore),
        m_K(1),
        m_error(std::numeric_limits<double>::infinity()),
        m_loglikelihood(std::numeric_limits<double>::infinity()),
        m_previousLikelihood(std::numeric_limits<double>::infinity()),
        m_bestLikelidood( -std::numeric_limits<double>::infinity())
    {

    }

    IGMM::IGMM(const Math::Matrix& data, const IGMMHyperparameters& hypes) :
        m_hypes(hypes),
        m_data(data),
        m_M(data.Rows()),
        m_N(data.Cols()),
        m_P( 1, std::vector< double >(m_N, 1.0) ),
        m_PSum(m_N, 0.0),
        m_ProbAny( m_N, 0),
        m_BadClass( 1, false),
        m_current({
            std::vector< uint32_t >(m_N, 0), // Classes
            std::vector< Math::Matrix >(1, Math::Matrix(m_M, 1)),
            std::vector< Math::Matrix >(1, Math::Matrix(m_M, m_M)),
            std::vector< uint32_t >(1, m_N)
        }),
        m_K(1),
        m_error(std::numeric_limits<double>::infinity()),
        m_loglikelihood(std::numeric_limits<double>::infinity()),
        m_previousLikelihood(std::numeric_limits<double>::infinity()),
        m_bestLikelidood( -std::numeric_limits<double>::infinity())
    {
        CalculateParametersFromExpectedLikelihood();

        // Fix the initial mu - the mean may not be near any datapoints.
        for(uint32_t i=0; i<m_M; ++i)
            m_current.m_Mu[0](i) = m_data(i,0);

        CalculateProbabilityOfAnyClass();
    }

    std::vector<uint32_t> IGMM::Classify( const Math::Matrix& data ) const
    {
        if( m_best.m_Ns.size() == 0  )
        {
            std::cerr << "IGMM::Classify:: Estimate returned no parameters, or it was not called" <<std::endl;
            assert(false);
        }

        std::vector<uint32_t> classes(data.Cols(), 0);
        std::vector< std::vector<double> > MVNprobs(m_best.m_Ns.size());

#pragma omp parallel for
        for(uint32_t m=0;m<m_best.m_Ns.size(); ++m)
        {
            MVNprobs[m] = Math::MvnPdf(data, m_best.m_Mu[m], m_best.m_Sigma[m]);
            for(auto& mv : MVNprobs[m])
            {
                mv *= static_cast<double>(m_best.m_Ns[m]) / m_data.Cols();
            }
        }

#pragma omp parallel for
        for(uint32_t k=0; k<data.Cols(); ++k)
        {
            uint32_t bestindex=0;
            double max= -std::numeric_limits<double>::max();
            double ma;
            for(uint32_t m=0;m<m_best.m_Ns.size(); ++m)
            {
                ma = MVNprobs[k][m] ;
                if(ma > max)
                {
                    max = ma;
                    bestindex = m;
                }
            }
            classes[k] = bestindex;
        }

        return classes;
    }

    std::vector<double> IGMM::EvaluateLikelihood(const Math::Matrix& data) const
    {
        if( m_best.m_Ns.size() == 0  )
        {
            std::cerr << "IGMM::Classify:: Estimate returned no parameters, or it was not called" <<std::endl;
            assert(false);
        }

        std::vector<double> probs(data.Cols(), 0);
        std::vector< std::vector<double> > MVNprobs(m_best.m_Ns.size());

        DataStore besttmp = m_best;

#pragma omp parallel for
        for(uint32_t m=0;m<besttmp.m_Ns.size(); ++m)
        {
            MVNprobs[m] = Math::MvnPdf(data, besttmp.m_Mu[m], besttmp.m_Sigma[m]);
            for(auto& mv : MVNprobs[m])
            {
                mv *= static_cast<double>(besttmp.m_Ns[m]) / m_data.Cols();
            }
        }

#pragma omp parallel for
        for(uint32_t k=0; k<data.Cols(); ++k)
        {
            double sum = 0.0;
            for(uint32_t m=0;m<besttmp.m_Ns.size(); ++m)
            {
                sum += MVNprobs[m][k];
            }
            probs[k] = sum;
        }

        return probs;
    }

    Math::Matrix IGMM::FindAssignmentLikelihoods(const Math::Matrix& data) const
    {
        Math::Matrix AssignmentProbs( m_best.m_Ns.size(), data.Cols(), 0.0 );

#pragma omp parallel for
        for(uint32_t m=0;m<m_best.m_Ns.size(); ++m)
        {
            std::vector<double> MVNprobs = Math::MvnPdf(data, m_best.m_Mu[m], m_best.m_Sigma[m]);
            for(uint32_t i=0; i<MVNprobs.size(); ++i)
            {
                AssignmentProbs(m,i) = MVNprobs[i]*(static_cast<double>(m_best.m_Ns[m]) / m_data.Cols());
            }
        }

        return AssignmentProbs;
    }

    Math::Matrix IGMM::FindNormalisedAssignmentLikelihoods(const Math::Matrix& data) const
    {
        Math::Matrix AssignmentProbs( m_best.m_Ns.size(), data.Cols(), 0.0 );

#pragma omp parallel for
        for(uint32_t m=0;m<m_best.m_Ns.size(); ++m)
        {
            std::vector<double> MVNprobs = Math::LogMvnPdf(data, m_best.m_Mu[m], m_best.m_Sigma[m]);
            for(uint32_t i=0; i<MVNprobs.size(); ++i)
            {
                AssignmentProbs(m,i) = MVNprobs[i] + std::log(static_cast<double>(m_best.m_Ns[m]) / m_data.Cols());
            }
        }

        std::vector<double> sums(AssignmentProbs.Cols(),0.0);
        auto nassignments = Math::NormaliseLogarithmicMatrix( AssignmentProbs, sums);

        return nassignments;
    }

    void IGMM::Estimate( const double tol, const uint32_t MaxIter )
    {
        uint32_t currentIter = 0;
        while( currentIter < MaxIter )
        {
            ComputeExpectedLikelihood();
            CalculateLogLikelihoodAndErrors();
            NormaliseProbabilities();

            if( SampleClassProbabilities(false) )
            {
                ++m_K;
            }

            CalculateParametersFromExpectedLikelihood();

            RemoveEmptyComponents();

            SaveBest();

            SampleParameters();

            ++currentIter;
            std::cout << "Iteration : " << currentIter << " / " << MaxIter << ", Best Likelihood : " << m_bestLikelidood << ", NComps : " << m_K << "  \r";
            std::cout.flush();
        }
        std::cout << std::endl;
    }

    void IGMM::RemoveEmptyComponents()
    {
        int32_t badindex=0;

        for(auto bs : m_BadClass)
            if(bs) ++badindex;

        for(uint32_t k=m_K; k>0; --k)
        if(m_BadClass[k-1])
        {
            m_current.m_Mu.erase( m_current.m_Mu.begin() + k - 1 );
            m_current.m_Sigma.erase( m_current.m_Sigma.begin() + k - 1 );
            m_current.m_Ns.erase( m_current.m_Ns.begin() + k - 1 );
            m_P.erase( m_P.begin() + k - 1 );
            m_K--;
        }
        else
        {
            if(badindex > 0)
            {
                for(auto& label : m_current.m_Classes)
                if(label == k-1)
                {
                    label -= badindex;
                }
            }
            --badindex;
        }
    }

    void IGMM::SampleParameters()
    {
        const double Kappa = 1;

#pragma omp parallel for
        for(uint32_t k=0; k<m_K; ++k)
        {
            // Find out if Sigma is cholesky
            if(m_hypes.Beta + m_current.m_Ns[k] > m_M) // Make sure we do not sample a invalid covariance
            {
                Math::Matrix Mun, Tn;
                double vn, kappan;
                ComputeParameters(Mun, vn, Tn, kappan, k);

                m_current.m_Sigma[k] = m_gen.SampleWishart( Tn, vn, 1)[0].inv() / kappan;
                m_current.m_Mu[k] = m_gen.SampleMultivariateNormal(Mun, m_current.m_Sigma[k], 1);
            }

            if(m_hypes.Beta + m_current.m_Ns[k] <= m_M || !isSPD(m_current.m_Sigma[k]))
            { // if we do, only sample the prior
                m_current.m_Mu[k] = (Kappa*m_hypes.Lambda + m_current.m_Ns[k]*m_current.m_Mu[k]) / (Kappa+m_current.m_Ns[k]);
                m_current.m_Sigma[k] = m_hypes.R.inv();
            }
        }
    }

    void IGMM::ComputeParameters( Math::Matrix& Mun, double& vn, Math::Matrix& Tn, double& kappan, const uint32_t k )
    {
        const uint32_t Kappa = 1;
        Mun = (m_hypes.Lambda*Kappa + m_current.m_Ns[k]*m_current.m_Mu[k]) / (Kappa*m_current.m_Ns[k]);
        kappan = Kappa + m_current.m_Ns[k];
        vn = m_hypes.Beta + m_current.m_Ns[k];
        Tn = m_hypes.R.inv() + m_current.m_Ns[k]*m_current.m_Sigma[k] +
             (static_cast<double>(Kappa*m_current.m_Ns[k]) / (m_current.m_Ns[k] + Kappa)) *
             (m_hypes.Lambda - m_current.m_Mu[k])*(m_hypes.Lambda - m_current.m_Mu[k]).transpose();
        Tn = Tn.inv();
    }

    bool IGMM::SampleClassProbabilities(const bool maxdebug)
    {

        bool newClass = false;
#pragma omp parallel for
        for(uint32_t i=0; i<m_N; ++i)
        {
            std::vector<double> distro(m_P.size());
            for(uint32_t k=0; k<m_P.size(); ++k)
            {
                distro[k] = m_P[k][i];
            }

            if(maxdebug)
                m_current.m_Classes[i] = std::max_element(distro.begin(), distro.end()) - distro.begin();
            else
                m_current.m_Classes[i] = m_gen.SampleDiscrete1D(distro, 1)[0];

#pragma omp critical
            if(!newClass && m_current.m_Classes[i] == m_K)
            {
                newClass=true;
            }
        }
        return newClass;
    }

    void IGMM::SaveBest()
    {
        if(m_loglikelihood > m_bestLikelidood)
        {
            m_best = m_current;
            m_bestLikelidood = m_loglikelihood;
            if(m_best.m_Mu.size() != m_best.m_Sigma.size() ||
                    m_best.m_Mu.size() != m_best.m_Ns.size() ||
                    m_best.m_Sigma.size() != m_best.m_Ns.size())
            {
                std::cerr<<"Error::IGMM::SaveBest() The sizes do not match" <<std::endl;
                assert(false);
            }
        }
    }

    void IGMM::NormaliseProbabilities()
    {
        for(auto& p : m_PSum)
            p=0.0;

        const uint32_t K = m_P.size();

#pragma omp parallel for
        for(uint32_t k=0; k<K; ++k)
        for(uint32_t i=0; i<m_N; ++i)
        {
            m_PSum[i] += m_P[k][i];
        }

#pragma omp parallel for
        for(uint32_t k=0; k<K; ++k)
        for(uint32_t i=0; i<m_N; ++i)
        {
            if(m_PSum[i] < 1e-13)
                m_P[k][i] = (k < K-1) ? 0.0 : 1.0;
            else
                m_P[k][i] /= m_PSum[i];
        }
    }

    void IGMM::CalculateLogLikelihoodAndErrors()
    {
        m_previousLikelihood = m_loglikelihood;
        m_loglikelihood = 0.0;

#pragma omp parallel for
        for( uint32_t i=0; i<m_N; ++i)
        {
            double sum = 0.0;
            for( uint32_t k=0; k<m_P.size(); ++k )
            {
                sum += m_P[k][i];
            }

#pragma omp critical
            if(sum > 1e-13)
            {
                m_loglikelihood += std::log(sum);
            }
        }
        m_loglikelihood /= m_N;

        m_error = std::abs( m_loglikelihood - m_previousLikelihood );
    }

    void IGMM::CalculateProbabilityOfAnyClass()
    {
        const double Kappa = 1;
        Math::Matrix BT = (1.0 / (Kappa + m_N))*(m_hypes.R.inv() + m_N*m_current.m_Sigma[0] +
                (static_cast<double>(Kappa*m_N)/(Kappa + m_N))*(m_hypes.Lambda - m_current.m_Mu[0])*(m_hypes.Lambda - m_current.m_Mu[0]).transpose());

        Math::Matrix StudentSigma = BT*(static_cast<double>(Kappa + m_N + 1 )/(Kappa*(m_hypes.Beta+m_N-m_M+1)));
        Math::Matrix StudentMu = (Kappa*m_hypes.Lambda + m_N*m_current.m_Mu[0]) / (Kappa*(m_hypes.Beta + m_N - m_M + 1));
        m_ProbAny = Math::StudentT3( m_data, StudentSigma, StudentMu, Kappa + m_N - m_M + 1 );
        for(auto& p : m_ProbAny)
        {
            p *= m_hypes.Alpha / (m_N - 1 + m_hypes.Alpha);
        }
    }

    void IGMM::ComputeExpectedLikelihood()
    {
        m_P.resize(m_K+1);

#pragma omp parallel for
        for(uint32_t k=0; k<m_K; ++k)
        {
            m_P[k] = Math::MvnPdf( m_data, m_current.m_Mu[k], m_current.m_Sigma[k]);
            for(auto& p : m_P[k])
            {
                p *= m_current.m_Ns[k] / (m_N - 1 + m_hypes.Alpha);
            }
        }
        m_P[m_K] = m_ProbAny;
    }

    void IGMM::CalculateParametersFromExpectedLikelihood()
    {
        m_current.m_Mu.resize(m_K, Math::Matrix(m_M,1));
        m_current.m_Sigma.resize(m_K, Math::Matrix(m_M,m_M));
        m_current.m_Ns.resize(m_K, 0);
        m_BadClass.resize(m_K, false);

#pragma omp parallel for
        for(uint32_t k=0; k<m_K; ++k)
        {
            m_current.m_Mu[k].SetZero();
            m_current.m_Sigma[k].SetZero();
            m_BadClass[k] = false;
        }

#pragma omp parallel for
        for(uint32_t k=0; k<m_K; ++k)
        {
            uint32_t index=0;
            for(uint32_t j=0; j<m_N; ++j)
            {
                if(m_current.m_Classes[j] == k)
                {
                    for(uint32_t i=0; i<m_M; ++i)
                        m_current.m_Mu[k](i) += m_data(i,j);
                    ++index;
                }
            }
            m_current.m_Ns[k] = index;

            if(index > 0)
            {
                m_current.m_Mu[k] /= index;

                if(index > 1)
                {
                    for(uint32_t j=0; j<m_N; ++j)
                    if(m_current.m_Classes[j] == k)
                    for(uint32_t i=0; i<m_M; ++i)
                    for(uint32_t l=0; l<m_M; ++l)
                        m_current.m_Sigma[k](i,l) += (m_data(i,j) - m_current.m_Mu[k](i)) *
                                                     (m_data(l,j) - m_current.m_Mu[k](l));

                    m_current.m_Sigma[k] /= index - 1;
                }
                else
                {
                    m_current.m_Sigma[k].SetIdentity();
                }
            }
            else
            {
                m_current.m_Mu[k].SetZero();
                m_current.m_Sigma[k].SetZero();
                m_BadClass[k] = true;
            }
        }
    }
}
