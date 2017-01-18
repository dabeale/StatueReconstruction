#ifndef SAMPLER_H
#define SAMPLER_H

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

#include "Matrix.h"
#include <random>
#include <chrono>

namespace Math
{
    /**
     * @brief The Generator class.
     * A class used to generate samples from a collection of probability distributions. This class is
     * necessary because it generates pseudo random number and so keeps an update of how many samples
     * have been taken since the seed.
     */
    class Generator
    {
    public:
        /**
         * @brief Generator.
         * Construct an empty random number generator.
         */
        Generator( );

        /**
         * @brief SampleUniformInteger.
         * Sample a uniform signed 32bit integer with the given parameters
         * @param min The minimum number
         * @param max The maximum number
         * @param NSamples The number of samples to draw
         * @return A vector of length NSamples
         */
        std::vector<int32_t> SampleUniformInteger(const int32_t min, const int32_t max, const uint32_t NSamples);

        /**
         * @brief SampleUniformInteger.
         * Sample a uniform unsigned 32bit integer with the given parameters
         * @param min The minimum number
         * @param max The maximum number
         * @param NSamples The number of samples to draw
         * @return A vector of length NSamples
         */
        std::vector<uint32_t> SampleUniformIntegerUI(const uint32_t min, const uint32_t max, const uint32_t NSamples);


        /**
         * @brief SampleDiscrete1D.
         * Sample from a discrete distribution.
         * @param distribution A 1xK matrix
         * @param N The number of samples to draw
         * @return A vector of unsigned 32bit integers - the indecies
         */
        std::vector<uint32_t> SampleDiscrete1D( const Matrix& distribution, const uint32_t N );

        /**
         * @brief SampleDiscrete1D.
         * Sample from a discrete distribution.
         * @param distribution A vector of length K
         * @param N The number of samples to draw
         * @return A vector of unsigned 32bit integers - the indecies
         */
        std::vector<uint32_t> SampleDiscrete1D(const std::vector<double> &dist, const uint32_t N );

        /**
         * @brief SampleDiscrete2D.
         * Sample from a discrete 2D distribution, using Gibbs sampling.
         * @param distribution A KxL matrix
         * @param N The number of samples to draw
         * @return A pair of 32bit unsigned integers, the 2D indecies.
         */
        std::vector<std::pair<uint32_t, uint32_t> > SampleDiscrete2D( const Matrix& distribution, const uint32_t N );


        /**
         * @brief SamplePermutation.
         * Sample a permutation vector of length N.
         * @param N The length of the vector.
         * @return A vector of length N.
         */
        std::vector<uint32_t> SamplePermutation( const uint32_t N );

        /**
         * @brief SampleMultivariateNormal.
         * Draw samples from a multivariate normal distribution.
         * @param mu A Mx1 matrix containing the mean
         * @param Sigma An MxM matrix containing the covariance
         * @param N The number of samples to draw
         * @return A MxN matrix of samples, each column containing a different sample.
         */
        Matrix SampleMultivariateNormal( const Matrix& mu, const Matrix& Sigma, const uint32_t N );

        /**
         * @brief SampleGMM.
         * Sample from a Gaussian Mixture Model.
         * @param Mu An MxK matrix of means, each column containing a separate mean.
         * @param Sigma A vector of length K, each element a MxM covariance matrix
         * @param Pi a 1xK matrix of coefficients. This vector must sum to one.
         * @param N The number of samples to take
         * @return A MxN matrix of floating point samples, each column contains a different sample.
         */
        Matrix SampleGMM(const Matrix &Mu, const std::vector<Matrix> &Sigma, const Matrix &Pi, const uint32_t N );

        /**
        * @brief SampleWishart.
        * Sample from a Wishart distribution.
        * @param Sigma The covariance prior
        * @param df The degrees of freedom
        * @param N The number of samples to take
        * @return A vector of length N, containing MxM matrices.
        */
        std::vector<Matrix> SampleWishart( const Matrix& Sigma, const uint32_t df, const uint32_t N );

        /**
         * @brief UniformlyRandomRotation.
         * Sample a random rotation matrix which is uniformly distributed.
         * @param M The dimensionality of the matrix
         * @return An MxM matrix
         */
        Matrix UniformlyRandomRotation(const uint32_t M);

        /**
         * @brief SubSampleColumns.
         * Sample the colums of a matrix. This sampler avoids repetition.
         * @param data An MxK matrix of samples
         * @param N The number of columns to sample (N <= K)
         * @return an MxN matrix with each column a random column of data.
         */
        Matrix UniquelySubSampleColumns(Matrix data, uint32_t N , double sigma=0.01);
    protected:
        std::default_random_engine m_e; ///< The stl pseudo random number generator
        uint32_t seed;                  ///< The seed
    };
}

#endif // SAMPLER_H
