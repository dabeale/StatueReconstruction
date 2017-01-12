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

namespace Cu
{ // The namespace is necessary to avoid conflicts with another library, so it is also quite short
    typedef Eigen::MatrixXd Matrix;
    typedef Eigen::VectorXd Vector;
    typedef Eigen::MatrixXcd ComplexMatrix;
    typedef Eigen::VectorXcd ComplexVector;
    typedef Eigen::SparseMatrix<double> SPMatrix;
/*
    void print_matrix(const Cu::Matrix&, const std::string&);
    void print_matrix_octave(const Cu::Matrix&, const std::string&);
    void print_matrix_octave(const std::vector<std::vector<double>>& A, const std::string& name);
    void print_matrix_octave(const double* data, const uint32_t M, const uint32_t N, const std::string& name, bool binary=false);
    void print_matrix_octave(const uint32_t* data, const uint32_t M, const uint32_t N, const std::string& name, bool binary=false);
    void print_matrix_octave(const uint32_t* data, const uint32_t M, const uint32_t N, const uint32_t K, const std::string& name);
    void print_matrix_octave(const double* data, const uint32_t M, const uint32_t N, const uint32_t K, const std::string& name, bool binary=false);
    void transpose(std::vector<double>& A, uint32_t M, uint32_t N);
*/
    void print_matrix_octave(const int32_t* data, const uint32_t M, const uint32_t N, const std::string& name, bool binary=false);


    typedef std::pair< std::vector<double> , std::vector<uint32_t> > DataDims;
    DataDims read_matrix_octave(const std::string &filename, std::string& name);
    Matrix read_matrix_bundler( const std::string& filename );

    Cu::Matrix ConvertToEigen(const double *in, const uint32_t M, const uint32_t N );
}

namespace Math
{
    class Matrix;
    static void Multiply(const Matrix& A, const Matrix& B, Matrix &C); ///< A fast multiplication routine (no checking)
    static void Add(const Matrix& A, const Matrix& B, Matrix &C); ///< Fast add (no checking)
    static void Subtract(const Matrix& A, const Matrix& B, Matrix &C); ///< Fast subtract (no checking)

    class Matrix
    {
    public:
        friend void Multiply(const Matrix& A, const Matrix& B, Matrix &C); ///< A fast multiplication routine (no checking)
        friend void Add(const Matrix& A, const Matrix& B, Matrix &C); ///< Fast add (no checking)
        friend void Subtract(const Matrix& A, const Matrix& B, Matrix &C); ///< Fast subtract (no checking)

        Matrix();
        Matrix(const uint32_t M, const uint32_t N);
        Matrix(const uint32_t M, const uint32_t N, const double val);
        Matrix(const uint32_t M, const uint32_t N, const std::vector<double>& array, bool rowmajor=true);
        Matrix(std::initializer_list<double> array);

        uint32_t Rows() const; ///< Return the number of rows
        uint32_t Cols() const; ///< Return the number of columns
        uint32_t numel() const; ///< Return the number of elements

        double& operator () (uint32_t row, uint32_t col);
        double  operator () (uint32_t row, uint32_t col) const;
        double& operator () (uint32_t elem);
        double  operator () (uint32_t elem) const;

        Matrix& operator += (const Matrix& m);
        Matrix& operator -= (const Matrix& m);
        Matrix& operator *= (const Matrix& m);
        Matrix& operator *= (const double& c);
        Matrix& operator /= (const double& c);
        Matrix& operator ^= (const double& pow);
        Matrix& operator ^= (const Matrix& m);

        Matrix transpose() const;

        friend std::ostream& operator<<(std::ostream& s, const Matrix& f);

        friend Matrix operator +(const Matrix& A, const Matrix& B);
        friend Matrix operator -(const Matrix& A, const Matrix& B);
        friend Matrix operator *(const Matrix& A, const Matrix& B);
        friend Matrix operator +(const Matrix& A, const double B);
        friend Matrix operator -(const Matrix& A, const double B);
        friend Matrix operator *(const Matrix& A, const double B);
        friend Matrix operator /(const Matrix& A, const double B);
        friend Matrix operator +(const double B, const Matrix& A);
        friend Matrix operator -(const double B, const Matrix& A);
        friend Matrix operator *(const double B, const Matrix& A);
        friend Matrix operator /(const double B, const Matrix& A);
        friend Matrix operator -(const Matrix& A);

        double det() const;
        double det(const Matrix& a, uint32_t n) const;
        Matrix inv() const;
        Matrix solve( const Matrix& b ) const;

        Matrix mean( const uint32_t dim=0 ) const;
        Matrix cov( const uint32_t dim=0 ) const;

        std::tuple<Matrix, Matrix, Matrix> LU() const; ///< @returns LUP LU decomposition of a square matrix
        std::tuple<Matrix,Matrix,Matrix> svd() const; ///< singular value decomposition for a square matrix

        void SetIdentity(); ///< Set this matrix to the identity
        void SetZero();
        void SwapRow( const uint32_t i, const uint32_t j);

        Matrix Column(const uint32_t i) const;
        void SetColumn(const uint32_t i, const Matrix& a);
        Matrix Row(const uint32_t i) const;
        void SetRow(const uint32_t i, const Matrix& a);
        Matrix Block(const uint32_t imin, const uint32_t jmin,
                     const uint32_t imax, const uint32_t jmax) const;

        Matrix Sum( const uint32_t dimension) const;

        Matrix SolveZero() const; ///< Find the solution to Ax=0

        Matrix pinv() const; ///< Return the moore penrose pseudo inverse

        std::pair<Matrix, Matrix> QR() const; ///< Return the QR decomposition

        // This may need speeding up.
        // There is a mem copy on every matrix operation.
        // This would be improved if there were fast element wise operations
        std::pair<Matrix, Matrix> eig() const; ///< return the eigenvalue decomposition
        double mineig() const; ///< Return the minimum eigenvalue

        double Norm() const;

        const double *Data() const;
        std::pair<Matrix, Matrix> Jacobi() const; ///< Eigendecomposition of a symmetric matrix
        Matrix chol() const; ///< Cholesky decomposition, works for positive semidefinite

        const std::vector<double>& GetArr() const;
        std::vector<double>& GetArr();

        void AddVectorToColumns(const Matrix& mat); ///< Add a vector matrix to each of the columns

        void Save( const std::string& filename) const;
        double dot(const Matrix& c);
    private:
        std::vector<double> m_arr;
        uint32_t m_M;
        uint32_t m_N;

        inline int pivot (uint32_t row);
    };
}
#endif // MATHMATRIX_H
