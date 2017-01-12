#ifndef INFINITEGMM_H
#define INFINITEGMM_H

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

#include <limits>
#include <vector>

#include "Matrix.h"
#include "Distributions.h"
#include "Sample.h"

namespace Sample
{
    struct IGMMHyperparameters
    {
        Math::Matrix Lambda;
        Math::Matrix R;
        double Alpha;
        uint32_t Beta;
    };

    struct DataStore
    {
        std::vector< uint32_t > m_Classes; ///< The classes
        std::vector< Math::Matrix > m_Mu;
        std::vector< Math::Matrix > m_Sigma;
        std::vector< uint32_t > m_Ns; ///< The number of elements per class

        inline void RemoveIllConditionedCovariances()
        {
            for(uint32_t k= m_Sigma.size()-1; k<m_Sigma.size(); ++k)
            {
                Math::Matrix L = m_Sigma[k].chol().transpose();

                double proddiag=1.0;
                for(uint32_t k=0; k<m_Sigma[k].Cols(); ++k)
                    proddiag *= L(k,k);

                if(proddiag == 0)
                {
                    m_Mu.erase(m_Mu.begin() + k);
                    m_Sigma.erase(m_Sigma.begin() + k);
                    if(m_Ns[k]>0)
                    {
                        std::cerr << "Warning::DataStore::RemoveIllConditionedCovariances() - There are values assigned to the ill-conditioned covariance " << std::endl;
                    }
                    m_Ns.erase(m_Ns.begin() + k);
                }
            }
        }
    };

    /**
     * @brief The IGMM class
     * This class is an imlpementation of the Infinite Gaussian Mixture Model.
     * It does not work because of the limitations of the Eigen::Matrix, but it
     * has been tested against an octave implementation, and appears to be working numerically.
     * A segmentation fault may happen when the mixture is estimated.
     * Rather that trying to fix the Eigen problem I have migrated some code in to
     * FLib which is for hierachical statistics and computer vision. It has a bespoke
     * Matrix library which will not suffer the same problems, at the expense of
     * potentially higher complexity.
     */
    class IGMM
    {
    public:
        IGMM( const Math::Matrix& data, const IGMMHyperparameters& hypes );
        IGMM( const Math::Matrix& data, const DataStore& datastore, const IGMMHyperparameters& hypes );

        void Estimate( const double tol=1e-6, const uint32_t MaxIter=500 );

        inline std::vector< Math::Matrix > GetMu() const { return m_best.m_Mu; }
        inline std::vector< Math::Matrix > GetSigma() const { return m_best.m_Sigma; }
        inline std::vector< uint32_t > GetNs() const { return m_best.m_Ns;}
        inline std::vector< uint32_t > GetClasses() const { return m_best.m_Classes; }

        /**
         * @brief Classify
         * classify a data matrix in to each of the Gaussian components
         * @param data
         * @return
         */
        std::vector<uint32_t> Classify( const Math::Matrix& data ) const;

        /**
         * @brief EvaluateLikelihood
         * Evaluate the likelhood of a collection of data opints against the GMM
         * @return
         */
        std::vector<double> EvaluateLikelihood(const Math::Matrix& data) const;

        /**
         * @brief FindAssignmentLikelihoods
         * Find the probabilities that each data point bleongs to a component
         * @return
         */
        Math::Matrix FindAssignmentLikelihoods(const Math::Matrix& data) const;
        Math::Matrix FindNormalisedAssignmentLikelihoods(const Math::Matrix& data) const;

        inline void ConditionBestSigmas()
        {
            Math::Matrix eye(m_M, m_M);
            eye.SetIdentity();

            for(  Math::Matrix& S : m_best.m_Sigma )
            {
                if(S.det() < 1e-6)
                {
                    double mineig = S.mineig();
                    if(mineig < 0)
                        S = S - mineig*eye;
                    S = S + 1e-6*eye;
                }
            }
        }

    private:
        const IGMMHyperparameters& m_hypes; ///< The hyperparameters

        const Math::Matrix& m_data; ///< The data (MxN) column major.
        const uint32_t m_M; ///< The dimensionality
        const uint32_t m_N; ///< The number of data points

        std::vector< std::vector< double > > m_P; ///< The membership probabilities
        std::vector< double > m_PSum; ///< The probability sums (temporary)
        std::vector< double > m_ProbAny; ///< The probability that the data belongs to any class
        std::vector< bool > m_BadClass; ///< True if a class is corrupt

        DataStore m_current; ///< The most recent parameters
        DataStore m_best; ///< The parameters with the best log likelihood

        uint32_t m_K; ///< The number of classes

        double m_error; ///< The log likelihood error
        double m_loglikelihood; ///< The log-likelihood
        double m_previousLikelihood; ///< The error on the previous step
        double m_bestLikelidood;

        Math::Generator m_gen;

        void CalculateParametersFromExpectedLikelihood(); ///< Calculate the means and covariances given the expected likelihood
        void CalculateProbabilityOfAnyClass(); ///< Calculate the pribability that the data belongs to any class (should only be called once)
        void ComputeExpectedLikelihood(); ///< Compute the expected likelihood given the parameters
        void CalculateLogLikelihoodAndErrors(); ///< Compute the log likelihood and errors
        void NormaliseProbabilities();
        void SaveBest(); ///< Saves the current parameters if they are the best log -likelihood
        bool SampleClassProbabilities(const bool maxdebug=false); ///< Sample the class probabilities returns true if there is a new class
        void ComputeParameters(Math::Matrix& Mun, double &vn, Math::Matrix& Tn, double& kappan, const uint32_t k  );
        void SampleParameters();
        void RemoveEmptyComponents();

        inline bool isSPD( const Math::Matrix& A, const double tol=1e-6 ) const
        {
            if(A.Rows() != A.Cols()) return false;
            auto usv = A.svd();
            auto& S = std::get<1>(usv);
            double minsing = std::numeric_limits<double>::max();

            for(uint32_t i=0; i<A.Rows(); ++i)
                if(S(i,i) < minsing)
                    minsing = S(i,i);
            return minsing > tol;
        }
    };
}

#endif
