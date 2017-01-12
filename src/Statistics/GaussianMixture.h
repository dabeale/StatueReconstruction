#ifndef GAUSSIANMIXTURE_CPP
#define GAUSSIANMIXTURE_CPP

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

#if defined(_OPENMP)
    #include <libiomp/omp.h>
#endif

#include "Matrix.h"
#include "Sample.h"

#include <cassert>
#include <functional>
#include <algorithm>

namespace Math
{
    class GaussianMixture
    {
    public:
        GaussianMixture(Matrix mu, std::vector<Matrix> sigma, Matrix pi);

        std::vector<double> Evaluate(const Matrix& x) const;

        Matrix DrawSamples(const uint32_t N);
        Matrix FindModes();

    private:
        Matrix m_mu; ///< A matrix containing hte center of each gaussian in a column
        std::vector<Matrix> m_sigma; ///< An array of covariances
        std::vector<Matrix> m_prec; ///< Array of precisions, computed on construction
        std::vector<double> m_dets; ///< Array of determinants, computed on construction
        Matrix m_pi; ///< Array of priors

        Generator m_g; ///< Sampling generator

        uint32_t m_M; ///< The dimension
        uint32_t m_K; ///< The number of components

        double ComponentProbability(const Matrix& x, const uint32_t i) const;
    };
}

#endif // GAUSSIANMIXTURE_CPP
