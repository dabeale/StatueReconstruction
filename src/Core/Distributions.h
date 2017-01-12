#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H

#include <stdint.h>
#include <cmath>
#include <numeric>
#include <limits>

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

    std::vector< double > MvnPdf( const Cu::Matrix& data, const Cu::Matrix& Mu, const Cu::Matrix& Sigma  );
    std::vector< double > StudentT3(const Cu::Matrix& data, const Cu::Matrix& Sigma, const Cu::Matrix& Mu, const uint32_t df=3 );
    std::vector< double > LatentGaussianWishart(const Cu::Matrix& data, const uint32_t df, const double k,
                                                const Cu::Matrix& Mu, const Cu::Matrix& Precision);
}

#endif
