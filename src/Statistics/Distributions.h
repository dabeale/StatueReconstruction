#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H

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

#include <stdint.h>
#include <cmath>
#include <numeric>
#include <limits>
#include <vector>

#include "Matrix.h"

namespace Distribution
{
    /**
     * @brief NormaliseLogarithmicMatrix
     * Normalise a logarithmic matrix while avoiding numerical underflow.
     * @param R - The KxN column major data matrix
     * @param sum - The rowwise sums of the data matrix
     * @param K - The number of rows
     * @param N - The number of columns
     */
    void NormaliseLogarithmicMatrix(double* R, double *sum, const uint32_t K, const uint32_t N);
    Math::Matrix NormaliseLogarithmicMatrix(const Math::Matrix& R, std::vector<double>& sum);

    /**
     * @brief distance2
     * Find the matrix of squared distances between two data matrices
     * @param data1 - KxN1 column major
     * @param data2 - KxN2 column major
     * @param D - N1xN2 squared distance matrix
     * @param N1 - The number of columns in data1
     * @param N2 - The number of columns of data2
     * @param K - The number of rows ( dimensionality )
     */
    void distance2(const double* data1, const double* data2, double* D,  const uint32_t N1, const uint32_t N2, const uint32_t K); ///< The distance 2 function

    std::vector< double > LogMvnPdf( const Math::Matrix& data, const Math::Matrix& Mu, const Math::Matrix& Sigma  );
    std::vector< double > MvnPdf( const Math::Matrix& data, const Math::Matrix& Mu, const Math::Matrix& Sigma  );
    std::vector< double > StudentT3(const Math::Matrix& data, const Math::Matrix& Sigma, const Math::Matrix& Mu, const uint32_t df=3 );
    std::vector< double > LatentGaussianWishart(const Math::Matrix& data, const uint32_t df, const double k,
                                                const Math::Matrix& Mu, const Math::Matrix& Precision);
}

#endif
