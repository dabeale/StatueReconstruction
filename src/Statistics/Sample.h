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

    class Generator
    {
    public:
        Generator( );

        /**
         * @brief Generator::SampleUniformInteger
         * Sample a uniform integer with the given parameters
         * @param min The minimum number
         * @param max The maximum number
         * @param NSamples The number of samples to draw
         * @return
         */
        std::vector<int32_t> SampleUniformInteger(const int32_t min, const int32_t max, const uint32_t NSamples);
        std::vector<uint32_t> SampleUniformIntegerUI(const uint32_t min, const uint32_t max, const uint32_t NSamples);


        /**
         * @brief SampleDiscrete1D
         * Sample from a discrete distribution
         * @param distribution A 1xK matrix
         * @param N The number of samples to draw
         * @return A vector of uints - the indecies
         */
        std::vector<uint32_t> SampleDiscrete1D( const Matrix& distribution, const uint32_t N );
        std::vector<uint32_t> SampleDiscrete1D(const std::vector<double> &dist, const uint32_t N );
        /**
         * @brief SampleDiscrete2D
         * Sample from a discrete 2D distribution. The process is Gibbs Sampling.
         * @param distribution A KxL matrix
         * @param N The number of samples to draw
         * @return A pair of uints, the 2D indecies
         */
        std::vector<std::pair<uint32_t, uint32_t> > SampleDiscrete2D( const Matrix& distribution, const uint32_t N );


        std::vector<uint32_t> SamplePermutation( const uint32_t N );

        /**
         * @brief SampleMultivariateNormal
         * Draw samples from a multivariate normal distribution.
         * @param mu
         * @param Sigma
         * @param N
         * @return
         */
        Matrix SampleMultivariateNormal( const Matrix& mu, const Matrix& Sigma, const uint32_t N );

        /**
         * @brief SampleGMM
         * Sample from the given GMM
         * @param gmm The GMM
         * @param N The number of samples
         * @return A matrix of floating point samples
         */
        Matrix SampleGMM(const Matrix &Mu, const std::vector<Matrix> &Sigma, const Matrix &Pi, const uint32_t N );

        /**
        * @brief SampleWishart
        * Sample from a wishart distribution. This implementation is taken from the octave source code, originally written by Nir Krakauer.
        * It is an adaptation for c++.
        * @param Sigma
        * @param df
        * @param N
        * @return
        */
        std::vector<Matrix> SampleWishart( const Matrix& Sigma, const uint32_t df, const uint32_t N );

        Matrix UniformlyRandomRotation(const uint32_t M);

        /**
         * @brief SubSampleColumns
         * Sample the colums of a matrix. This sampler avoids repetition.
         * @param data MxK
         * @param N The number of columns to sample (N <= K)
         * @return
         */
        Matrix UniquelySubSampleColumns(Matrix data, uint32_t N , double sigma=0.01);
    protected:
        std::default_random_engine m_e;
        uint32_t seed;
    };
}

#endif // SAMPLER_H
