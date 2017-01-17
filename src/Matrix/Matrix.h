#ifndef MATHMATRIX_H
#define MATHMATRIX_H

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

#include <vector>
#include <tuple>
#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>
#include <limits>
#include <initializer_list>

#include "Eigen/Dense"
#include "Eigen/LU"
#include "Eigen/Cholesky"
#include <Eigen/Sparse>
#include <Eigen/SPQRSupport>
#include <Eigen/Eigenvalues>

/**
 * \brief This namespace is an interface into the Eigen library. It is required for some of the
 * functions and classes and really only serves to avoid naming conflicts. Most of the matrices use
 * the Math::Matrix implementations, which are slower, but less prone to threading and memory management
 * problems.
 */
namespace Cu
{
    typedef Eigen::MatrixXd Matrix;                 ///< \brief A dynamic double precision Eigen matrix
    typedef Eigen::VectorXd Vector;                 ///< \brief A dynamic double precision Eigen vector
    typedef Eigen::MatrixXcd ComplexMatrix;         ///< \brief A dynamic complex Eigen matrix
    typedef Eigen::VectorXcd ComplexVector;         ///< \brief A dynamic complex Eigen vector
    typedef Eigen::SparseMatrix<double> SPMatrix;   ///< \brief A dynamic complex sparse Eigen matrix

    typedef std::pair< std::vector<double> , std::vector<uint32_t> > DataDims;  ///< An abstraction of an N tensor

    /**
     * @brief read_matrix_octave
     * Read a matrix stored in the octave format. This implementation currently only reads text files.
     * @param filename The path to the file
     * @param name A reference to the name of the matrix. This value is filled from the name stored in the file
     * @return A DataDims object containing the tensor.
     */
    DataDims read_matrix_octave(const std::string &filename, std::string& name);

    /**
     * @brief read_matrix_bundler
     * Read a matrix stored in the Bundler format. A bundler matrix is necessarily
     * a camera matrix (i.e. 3x4 homogeneous). The first line should be CONTOUR, for example,
     *
     * > CONTOUR
     * > 1 2 3 4
     * > 5 6 7 8
     * > 9 10 11 12
     *
     * @param filename
     * @return
     */
    Matrix read_matrix_bundler( const std::string& filename );

    /**
     * @brief ConvertToEigen
     * Map a double precision array on to a Eigen matrix. This imlpements the Eigen map
     * functionality and so it is as fast as a mem copy.
     * @param in
     * @param M
     * @param N
     * @return
     */
    Cu::Matrix ConvertToEigen(const double *in, const uint32_t M, const uint32_t N );
}


/**
 * \brief A collection of classes and functions for performing mathematical routines. This includes matrix algorithms and
 * statistical models.
 */
namespace Math
{
    class Matrix;
    static void Multiply(const Matrix& A, const Matrix& B, Matrix &C); ///< A fast multiplication routine (no checking)
    static void Add(const Matrix& A, const Matrix& B, Matrix &C); ///< Fast add (no checking)
    static void Subtract(const Matrix& A, const Matrix& B, Matrix &C); ///< Fast subtract (no checking)

    /**
     * @brief The Matrix class
     * This class implements a collection of algorithms for matrix operations. All of the methods are the most basic
     * and are not heavily optimised for speed. The matrix factorisation algorithms are taken from a collection of papers,
     * but also the Rosetta Code project website (http://rosettacode.org/wiki/Rosetta_Code).
     *
     * There are no dependencies on any other matrix, or linear algebra library, making the code completely thread safe.
     * It is slower than LAPACK and not parallelised.
     */
    class Matrix
    {
    public:
        friend void Multiply(const Matrix& A, const Matrix& B, Matrix &C); ///< A fast multiplication routine (no checking)
        friend void Add(const Matrix& A, const Matrix& B, Matrix &C); ///< Fast add (no checking)
        friend void Subtract(const Matrix& A, const Matrix& B, Matrix &C); ///< Fast subtract (no checking)

        /**
         * @brief Matrix
         * Construct an empty matrix
         */
        Matrix();

        /**
         * @brief Matrix
         * Construct an MxN matrix filled with zeros.
         * @param M The number of rows
         * @param N The number of columns
         */
        Matrix(const uint32_t M, const uint32_t N);

        /**
         * @brief Matrix
         * Construct an MxN matrix and fill it with a single value
         * @param M The number of rows
         * @param N The number of columns
         * @param val The value with which to fill the matrix
         */
        Matrix(const uint32_t M, const uint32_t N, const double val);

        /**
         * @brief Matrix
         * Construct an MxN matrix and fill it with the contents of an array. The array stored in
         * this matrix implementation is, by default, column major. The input array is assumed to
         * be row major by default since it allows a matrix to be constructed using an initializer_list,
         * as follows,
         * \code
         *  Matrix Mat( 3, 4, { 1, 2, 3, 4,
         *                      5, 6, 7, 8,
         *                      9, 10, 11, 12} );
         * \endcode
         *
         * If one wants to construct the matrix with a column major array, the default parameter rowmajor
         * should be set to false.
         * @param M The number of rows
         * @param N The number of columns
         * @param array The array to be copied
         * @param rowmajor True if the input array is row major.
         */
        Matrix(const uint32_t M, const uint32_t N, const std::vector<double>& array, bool rowmajor=true);

        /**
         * @brief Matrix
         * Construct a matrix using an initializer_list. This will assume that the
         * matrix is a column vector, and allows the matrix to be constructed as follows,
         * \code
         *  Matrix Mat {1, 2, 3, 4};
         * \endcode
         * @param array
         */
        Matrix(std::initializer_list<double> array);

        uint32_t Rows() const;                                      ///< \brief Return the number of rows
        uint32_t Cols() const;                                      ///< \brief Return the number of columns
        uint32_t numel() const;                                     ///< \brief Return the number of elements

        double& operator () (uint32_t row, uint32_t col);           ///< \brief Return a reference to the (row, col) element of the matrix
        double  operator () (uint32_t row, uint32_t col) const;     ///< \brief Return the (row,col) element of the matrix
        double& operator () (uint32_t elem);                        ///< \brief Return a reference to the  elem index of the matrix (column major)
        double  operator () (uint32_t elem) const;                  ///< \brief Return the elem index of the matrix (column major)

        Matrix& operator += (const Matrix& m);                      ///< \brief Overload the += operator for matrices
        Matrix& operator -= (const Matrix& m);                      ///< \brief Overload the -= operator for matrices
        Matrix& operator *= (const Matrix& m);                      ///< \brief Overload the *= operator for matrices
        Matrix& operator *= (const double& c);                      ///< \brief Overload the *= operator for a matrix with a double
        Matrix& operator /= (const double& c);                      ///< \brief Overload the /= operator for a matrix with a double
        Matrix& operator ^= (const double& pow);                    ///< \brief Overload the ^= operator for a matrix with a double. This will compute the power of each element to pow.
        Matrix& operator ^= (const Matrix& m);                      ///< \brief Overload the ^= operator for a matrix with a matrix. This will take the ith element of 'this' to the powr of the ith element of m.

        Matrix transpose() const;                                   ///< \brief Return the transpose of 'this'.

        friend std::ostream& operator<<(std::ostream& s, const Matrix& f); ///< \brief Overload the << operator for matrices on to the output stream, for printing to stdout.

        friend Matrix operator +(const Matrix& A, const Matrix& B); ///< \brief Overload the + operator for matrices
        friend Matrix operator -(const Matrix& A, const Matrix& B); ///< \brief Overload the - operator for matrices
        friend Matrix operator *(const Matrix& A, const Matrix& B); ///< \brief Overload the * operator for matrices
        friend Matrix operator +(const Matrix& A, const double B);  ///< \brief Overload the + operator for a matrix with a double
        friend Matrix operator -(const Matrix& A, const double B);  ///< \brief Overload the - operator for a matrix with a double
        friend Matrix operator *(const Matrix& A, const double B);  ///< \brief Overload the * operator for a matrix with a double
        friend Matrix operator /(const Matrix& A, const double B);  ///< \brief Overload the / operator for a matrix with a double
        friend Matrix operator +(const double B, const Matrix& A);  ///< \brief Overload the + operator for a double with a matrix
        friend Matrix operator -(const double B, const Matrix& A);  ///< \brief Overload the - operator for a double with a matrix
        friend Matrix operator *(const double B, const Matrix& A);  ///< \brief Overload the * operator for a double with a matrix
        friend Matrix operator /(const double B, const Matrix& A);  ///< \brief Overload the / operator for a double with a matrix
        friend Matrix operator -(const Matrix& A);

        /**
         * @brief det
         * Compute the determinant of a square matrix. It is a recursive algorithm
         * which computes the matrix of cofactors and muliplies each element by a power of -1.
         * @return The determinant.
         */
        double det() const;

        /**
         * @brief inv
         * Compute the inverse of the matrix usign Gauss-Jordan elimination.
         * @return The inverse
         */
        Matrix inv() const;

        /**
         * @brief solve
         * Solve the system of equations Ax = b using LU decomposition.
         * @param b
         * @return
         */
        Matrix solve( const Matrix& b ) const;

        /**
         * @brief mean
         * Compute the column or rowwise mean of the matrix. For example, calling
         * mean(0) returns a row vector containing the mean across the rows. It follows the
         * same semantics as the associated matlab function.
         * @param dim The dimension to take the mean across
         * @return The mean
         */
        Matrix mean( const uint32_t dim=0 ) const;

        /**
         * @brief cov
         * Compute the covariance across the specified dimension. For example, calling
         * mean(0) will return a square matrix with the same number of columns as 'this', containing
         * the covarainces across each of the rows.
         * @param dim The dimension to take the covaraince matrix across
         * @return The covariance matrix
         */
        Matrix cov( const uint32_t dim=0 ) const;

        /**
         * @brief LU
         * Return the LU decomposition of the matrix. This implementation is taken from the Rosetta Code project.
         * It will return three matrices L, U, and P, so that PA = LU; where A is 'this'.
         * @return A tuple containing L, U, and P; the lower triangular, upper triangular and permutation matrix, respectively.
         */
        std::tuple<Matrix, Matrix, Matrix> LU() const;

        /**
         * @brief svd
         * Return the singular value decomposition of 'this'. It will return three matrices U, S, V,
         * where U and V are left and right orthogonal transformation matrices and S is a diagonal matrix of
         * singular values.
         *
         * This routine is adapted from svdecomp.c in XLISP-STAT 2.1 which is
         * code from Numerical Recipes adapted by Luke Tierney and David Betz.
         * @return A tuple containing U, S and V, respectively.
         */
        std::tuple<Matrix,Matrix,Matrix> svd() const; ///< singular value decomposition for a square matrix

        void SetIdentity();                                             ///< \brief Set this matrix to the identity
        void SetZero();                                                 ///< \brief Set this matrix to be zero.
        void SwapRow( const uint32_t i, const uint32_t j);              ///< \brief Swap rows i and j

        Matrix Column(const uint32_t i) const;                          ///< \brief Return the ith column of the matrix
        void SetColumn(const uint32_t i, const Matrix& a);              ///< \brief Set the ith column to be 'a'
        Matrix Row(const uint32_t i) const;                             ///< \brief Return the ith row of the matrix
        void SetRow(const uint32_t i, const Matrix& a);                 ///< \brief Set the ith row to be 'a'

        Matrix Block(const uint32_t imin, const uint32_t jmin,
                     const uint32_t imax, const uint32_t jmax) const;   ///< \brief Return a submatrix with rows between imin and imax (inclusive), and columns between jmin and jmax (inclusive).

        Matrix Sum( const uint32_t dimension) const;                    ///< \brief Compute the sum across the 'dimension'.

        Matrix SolveZero() const;                                       ///< \brief Find the solution to Ax=0 using singular value decomposition

        Matrix pinv() const;                                            ///< \brief Return the moore penrose pseudo inverse.
        double dot(const Matrix& c);                                    ///< \brief Take the dot product of this matrix with 'c'

        /**
         * @brief QR
         * Return a pair of matrices Q and R such that A = QR, where A is 'this'.
         * @return A pair of matrices, Q and R, respectively.
         */
        std::pair<Matrix, Matrix> QR() const; ///< Return the QR decomposition

        /**
         * @brief eig
         * Compute the eigenvalue decomposition of a square matrix; matrices U and S such that
         * A = USU^{-1}.
         *
         * This implementation uses QR decomposition iteratively.
         * @return A pair of matrices containing
         */
        std::pair<Matrix, Matrix> eig() const;
        double mineig() const;                          ///< \brief Return the minimum eigenvalue

        double Norm() const;                            ///< \brief Return the Frobenius norm of the matrix

        const double *Data() const;                     ///< \brief Return a pointer to the first element of the matrix

        /**
         * @brief Jacobi
         * Compute the Eigen decompoosition for a symmetric matrix.
         * This code is an adaptation of the numerical recipes version
         * Originally written by David Squire (DMS), David.Squire@infotech.monash.edu.au
         * @return Matrices U and S, respectively.
         */
        std::pair<Matrix, Matrix> Jacobi() const;

        /**
         * @brief chol
         * Compute the cholesky decomposition for a square positive semidefinite matrix; a matrix L
         * such that A = L'L, where A is 'this'.
         *
         *
         * The original implementation was taken from the Rosetta Code project but it was
         * adapted for for positive semidefinite matrices according to the papers,
         *  'Singular Values using Cholesky Decomposition' - Krishnamoorthy and Kocagoez
         *  'Analysis of the Cholesky Decomposition of a Semi-definite Matrix' - N. J. Higham
         * @return
         */
        Matrix chol() const;

        const std::vector<double>& GetArr() const;      ///< \brief Return a constant reference to the column major vector of values
        std::vector<double>& GetArr();                  ///< \brief Return a reference to the column major vector of values

        void AddVectorToColumns(const Matrix& mat);     ///< \brief Add a vector matrix to each of the columns

        void Save( const std::string& filename) const;  ///< \brief Save the matrix to filenam in the octave format.
    private:
        std::vector<double> m_arr;                  ///< The internal vector of values (column major)
        uint32_t m_M;                               ///< The number of rows
        uint32_t m_N;                               ///< The number of columns

        /**
         * @brief pivot
         * Internal pivoting for Gauss-Jordan elimination.
         * @param row The row to pivot
         * @return
         */
        inline int pivot (uint32_t row);

        /**
         * @brief det
         * The internal determinant method. Used for recursion.
         * @param a A square matrix
         * @param n The number of rows
         * @return
         */
        double det(const Matrix& a, uint32_t n) const;
    };
}
#endif // MATHMATRIX_H
