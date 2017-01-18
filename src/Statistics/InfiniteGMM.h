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

namespace Math
{
    /**
     * @brief The IGMMHyperparameters struct.
     * This structure contains all of the priors for the IGMM.
     */
    struct IGMMHyperparameters
    {
        Math::Matrix Lambda; ///< An Mx1 matrix, center prior
        Math::Matrix R;      ///< An MxM matrix, precision prior
        double Alpha;        ///< The Dirichlet 'new class probability' prior
        uint32_t Beta;       ///< The degree of freedom prior
    };

    /**
     * @brief The DataStore struct.
     * All of the information required for the Gaussian mixture model
     */
    struct DataStore
    {
        std::vector< uint32_t > m_Classes;      ///< The class assignments, length N
        std::vector< Math::Matrix > m_Mu;       ///< A vector of means, length K
        std::vector< Math::Matrix > m_Sigma;    ///< A vector of covariances, length K
        std::vector< uint32_t > m_Ns;           ///< The number of elements per class, length K

        /**
         * @brief RemoveIllConditionedCovariances.
         * Remove all of the covariances from the parameter set which are not
         * symmetric positive definite, and also remove their associated means and
         * number of elements per class. This method will update the m_Classes variable
         * accordingly.
         */
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
     * @brief The IGMM class.
     * This class is for extimating the parameters of a Gaussian mixture model by Gibbs
     * sampling a Dirichlet process. The implementation is not fully documented, but more details
     * on the method can be found in the following papers,
     *   \li "The Infinite Gaussian Mixture Model" - Carl Rasmussen, NIPS 1999.
     *   \li "Markov Chain Sampling Methods for Dirichlet Process Mixture Models" - Radford Neal, Technical Report, Unverisity of Toronto 1998.
     *   \li "Conjugate Bayesian analysis of the Gaussian distribution" - Kevin P. Murphy, 2007.
     * The method not only estimates the parameters of the Gaussian mixture, but also the number of components, given some suitable
     * hyperparameters. The class does not estimate the hyperparameters.
     */
    class IGMM
    {
    public:
        /**
         * @brief IGMM.
         * Construct the Infinite Mixture with a data matrix (MxN) and the hyper parameters.
         * @param data An MxN matrix of observations one datapoint per column
         * @param hypes The hyper parameters.
         */
        IGMM( const Math::Matrix& data, const IGMMHyperparameters& hypes );

        /**
         * @brief IGMM.
         * Construct the Infinite Mixture with a data matrix, the hyperparameters, and also the parameters.
         * @param data  An MxN matrix of observations one datapoint per column
         * @param datastore The parameters
         * @param hypes The hyperparameters
         */
        IGMM( const Math::Matrix& data, const DataStore& datastore, const IGMMHyperparameters& hypes );

        /**
         * @brief Estimate.
         * Estimate the parameters of the distribution. This method will take a collections of samples,
         * up to MaxIter, and return the parameters which achieve the best log-likelihood.
         * @param tol The tolerance - not used in this implementation.
         * @param MaxIter The maximum number of iterations, ot number of samples.
         */
        void Estimate( const double tol=1e-6, const uint32_t MaxIter=500 );

        inline std::vector< Math::Matrix > GetMu() const { return m_best.m_Mu; }            ///< \brief Return the means
        inline std::vector< Math::Matrix > GetSigma() const { return m_best.m_Sigma; }      ///< \brief Return the covariances
        inline std::vector< uint32_t > GetNs() const { return m_best.m_Ns;}                 ///< \brief Return the number of datapoints per class
        inline std::vector< uint32_t > GetClasses() const { return m_best.m_Classes; }      ///< \brief Return the class assignments

        /**
         * @brief Classify.
         * Classify a data matrix in to each of the Gaussian components.
         * @param data A MxL data matrix, each column an observation.
         * @return A vector of length L of class assignments
         */
        std::vector<uint32_t> Classify( const Math::Matrix& data ) const;

        /**
         * @brief EvaluateLikelihood.
         * Evaluate the likelhood of a collection of data opints against the GMM
         * @param data An MxL data matrix, each column an observation
         * @return A vector of length L, containing the likelihoods given the distribution.
         */
        std::vector<double> EvaluateLikelihood(const Math::Matrix& data) const;

        /**
         * @brief FindAssignmentLikelihoods
         * Find the probabilities that each data point belongs to each of the components.
         * @param data An MxL data matrix, each column an observation
         * @return A KxL matrix, where K is the number of classes, containing the likelihood that the lth datapoint belongs to the kth class.
         */
        Math::Matrix FindAssignmentLikelihoods(const Math::Matrix& data) const;

        /**
         * @brief FindAssignmentLikelihoods.
         * Find the probabilities that each data point belongs to each of the components. This method copmutes the
         * log-Multivariate normal distribution and then uses ExpSumLog to normalise the matrix. It is completely robust to numberical
         * underflow, meaning it will never return inf or nan in the returned matrix.
         * @param data An MxL data matrix, each column an observation
         * @return A KxL matrix, where K is the number of classes, containing the likelihood that the lth datapoint belongs to the kth class. Summing over the rows will produce a row vector of ones.
         */
        Math::Matrix FindNormalisedAssignmentLikelihoods(const Math::Matrix& data) const;

        /**
         * @brief ConditionBestSigmas.
         * If any of the covariances in the best estimate are singular, condition them
         * by adding the identity. There are several methods for getting round the problem
         * of ill conditioned matrices in the computation, this is only one of them. There
         * are more details in the paper,
         *  "Markov Chain Sampling Methods for Dirichlet Process Mixture Models" - Radford Neal, Technical Report, Unverisity of Toronto 1998.
         */
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
        const IGMMHyperparameters& m_hypes;             ///< The hyperparameters

        const Math::Matrix& m_data;                     ///< The data matrix (MxN) column major.
        const uint32_t m_M;                             ///< The dimensionality
        const uint32_t m_N;                             ///< The number of data points

        std::vector< std::vector< double > > m_P;       ///< The membership probabilities i.e. the probabilities that a datapoint belongs to the kth class.
        std::vector< double > m_PSum;                   ///< The probability sums (temporary)
        std::vector< double > m_ProbAny;                ///< The probability that the data belongs to any class
        std::vector< bool > m_BadClass;                 ///< True if a class is corrupt (ill conditioned, or containing no points)

        DataStore m_current;                            ///< The most recent parameters
        DataStore m_best;                               ///< The parameters with the best log likelihood

        uint32_t m_K;                                   ///< The current number of classes

        double m_error;                                 ///< The log likelihood error
        double m_loglikelihood;                         ///< The log-likelihood
        double m_previousLikelihood;                    ///< The error on the previous step
        double m_bestLikelidood;                        ///< The best likelihood

        Math::Generator m_gen;                          ///< A random sample generator

        void CalculateParametersFromExpectedLikelihood();         ///< \brief Calculate the means and covariances given the expected likelihood
        void CalculateProbabilityOfAnyClass();                    ///< \brief Calculate the pribability that the data belongs to any class (should only be called once)
        void ComputeExpectedLikelihood();                         ///< \brief Compute the expected likelihood given the parameters
        void CalculateLogLikelihoodAndErrors();                   ///< \brief Compute the log likelihood and errors
        void NormaliseProbabilities();                            ///< \brief Normalise the membership probabilities
        void SaveBest();                                          ///< \brief Saves the current parameters if they are the best log -likelihood
        bool SampleClassProbabilities(const bool maxdebug=false); ///< \brief Sample the class probabilities returns true if there is a new class
        void SampleParameters();                                  ///< \brief Sample the parameters given the data
        void RemoveEmptyComponents();                             ///< \brief Remove any components which do not contain a datapoint

        /**
         * @brief ComputeParameters.
         * Compute the hyper-parameters given the estimated parameters and the data.
         * @param Mun A Mx1 matrix, the center hyperparameter
         * @param vn The degrees of freedom
         * @param Tn An MxM matrix, the precision hyperparameter
         * @param kappan The scale hyper-parameter
         * @param k The Gaussian component
         */
        void ComputeParameters(Math::Matrix& Mun, double &vn, Math::Matrix& Tn, double& kappan, const uint32_t k  );


        /**
         * @brief isSPD.
         * Return true if the given matrix is symmetric positive definite.
         * @param A A square test matrix
         * @param tol The tolerance
         * @return True if the matrix is spd.
         */
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
