#include "Matrix.h"

#define SIGN(a, b) ((b) >= 0.0 ? std::fabs(a) : -std::fabs(a))
#define MIN(x,y) ( (x) < (y) ? (x) : (y) )
#define MAX(x,y) ((x)>(y)?(x):(y))

namespace Cu
{
Cu::Matrix ConvertToEigen(const double *in, const uint32_t M, const uint32_t N )
{
    Cu::Matrix out(M,N);
    out = Eigen::Map<Cu::Matrix>(const_cast<double *>(in), M, N);
    return out;
}

void print_matrix_octave(const int32_t* data, const uint32_t M, const uint32_t N, const std::string& name, bool binary)
{
    uint32_t k=name.find_last_of('/');
    std::ofstream ostream(name);
    //ostream.open(name);
    ostream << " # Created by CLib" << std::endl;
    ostream << " # name: " << name.substr(k+1,name.size()-k-1) << std::endl;
    ostream << " # type: matrix " << std::endl;
    ostream << " # rows : " << M << std::endl;
    ostream << " # columns : " << N << std::endl;
    for( uint32_t j=0; j<M; j++ )
    {
        for( uint32_t i=0; i<N; i++ )
        {
            //ostream.write((char *)&data[M*i + j],sizeof(uint32_t));
            ostream << " " << data[M*i+j];
        }
        ostream << std::endl;
    }
     ostream.close();
}

DataDims read_matrix_octave(const std::string &filename, std::string& name)
{
    std::string line;
    std::ifstream file(filename);
    std::vector<uint32_t> dims;
    std::vector<double> data;
    bool lastwashash=true;
    uint32_t linenum=0;
    if (file.is_open())
    {
        while ( std::getline (file,line) )
        {
            if(line.size() > 0)
            {
                line.shrink_to_fit();
                if(line.find('#') !=std::string::npos)
                {
                    if(line.find("name") != std::string::npos )
                    {
                        uint32_t l = line.find_last_of(':');
                        name = line.substr(l+2, line.length() - l -1);
                    }
                    else if(line.find("type") !=std::string::npos)
                    {
                        uint32_t l = line.find_last_of(':');
                        std::string type =line.substr(l+2, line.length() - l -1);
                        if(type.find("matrix") == std::string::npos)
                        {
                            std::cerr << "read_matrix_octave : the file is not a matrix" << std::endl;
                            return DataDims();
                        }
                    }
                    else if(line.find("rows") !=std::string::npos)
                    {
                        if(dims.size()==0)
                        {
                            dims.resize(2);
                        }
                        uint32_t l = line.find_last_of(':');
                        std::string rows = line.substr(l+2, line.length() - l -1);
                        dims[0] = atoi(rows.c_str());
                    }
                    else if(line.find("cols") !=std::string::npos ||
                            line.find("columns") !=std::string::npos)
                    {
                        if(dims.size()==0)
                        {
                            dims.resize(2);
                        }
                        uint32_t l = line.find_last_of(':');
                        std::string cols = line.substr(l+2, line.length() - l -1);
                        dims[1] = atoi(cols.c_str());
                    }
                    else if(line.find("ndims") !=std::string::npos)
                    {
                        uint32_t l = line.find_last_of(':');
                        std::string dimsize = line.substr(l+2, line.length() - l -1);
                        dims.resize(atoi(dimsize.c_str()));
                    }
                    lastwashash =true;
                }
                else
                {
                    std::stringstream ss(line);

                    if(lastwashash)
                    {
                        if( dims.size()!=2)
                        {
                            uint32_t in, index=0;
                            while(!ss.eof() && index < dims.size())
                            {
                                ss >> in;
                                dims[index] = in;
                                ++index;
                            }
                        }
                        uint32_t in=1;
                        for(auto d : dims)
                        {
                            in*=d;
                        }
                        data.resize(in);
                        lastwashash = false;
                        if( dims.size()!=2) continue;
                    }

                    if(!lastwashash || dims.size()==2)
                    {
                        double in;
                        uint32_t index=0,i,j;
                        while(!ss.eof())
                        {
                            ss >> in;
                            i = linenum;
                            j = index;
                            data[ j*dims[0] + i] = in;
                            ++index;
                        }

                        ++linenum;
                    }
                }
            }
        }
        file.close();

        return std::make_pair(data, dims);
    }
    else
    {
        std::cerr << "read_matrix_octave : Unable to open file" << std::endl;
        return DataDims();
    }
}

Matrix read_matrix_bundler( const std::string& filename )
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Matrix::read_matrix_bundler :: Cannot open the file." << std::endl;
        return Matrix(3,4);
    }

    std::string line;
    std::getline (file,line);
    if(line.find("CONTOUR")==std::string::npos)
    {
        std::cerr << "Matrix::read_matrix_bundler :: The file is not a bundler file." << std::endl;
        return Matrix(3,4);
    }

    Matrix P(3,4);
    uint32_t l,k;
    for(l=0;l<3;++l)
    {
        std::getline (file,line);
        std::stringstream ss(line);
        for(k=0;k<4;++k)
            ss >> P(l,k);
    }
    return P;
}
}


namespace Math
{
    void Multiply(const Matrix& A, const Matrix& B, Matrix &C)
    {
        uint32_t r,c,k;
        for(r=0;r<A.Rows();++r)
        {
            for(c=0;c<B.Cols();++c)
            {
                C.m_arr[C.m_M*c + r] = 0;
                for(k=0;k<B.Rows();++k)
                {
                    C.m_arr[C.m_M*c + r] += A.Data()[A.m_M*k + r]*B.Data()[B.m_M*c + k];
                }
            }
        }
    }

    void Add(const Matrix& A, const Matrix& B, Matrix &C)
    {
        for(uint32_t k=0;k<A.numel();++k)
        {
            C.m_arr[k] = A.Data()[k]+B.Data()[k];
        }
    }

    void Subtract(const Matrix& A, const Matrix& B, Matrix &C)
    {
        for(uint32_t k=0;k<A.numel();++k)
        {
            C.m_arr[k] = A.Data()[k]-B.Data()[k];
        }
    }

    const double* Matrix::Data() const
    {
        return static_cast<const double*>(m_arr.data());
    }

    Matrix::Matrix() :
        m_M(0), m_N(0)
    {

    }

    Matrix::Matrix(const uint32_t M, const uint32_t N) :
        m_arr(M*N), m_M(M), m_N(N)
    {

    }

    Matrix::Matrix(const uint32_t M, const uint32_t N, const double val) :
        m_arr(M*N, val), m_M(M), m_N(N)
    {

    }

    Matrix::Matrix(const uint32_t M, const uint32_t N, const std::vector<double>& array, bool rowmajor) :
        m_arr(array), m_M(rowmajor?N:M), m_N(rowmajor?M:N)
    {
        if(m_M*m_N != m_arr.size())
        {
            std::cerr << "Matrix::Matrix (uin32_t, uint32_t, vector<double>) :: sizes do not match" <<std::endl;
            assert(false);
        }
        if(rowmajor)
            *this = this->transpose();
    }

    Matrix::Matrix(std::initializer_list<double> array) :
        m_arr(array), m_M(array.size()), m_N(1)
    {

    }

    Matrix Matrix::mean( const uint32_t dim ) const
    {
        const double N = (dim==0)?Rows():Cols();
        return (1.0/N)*this->Sum(dim);
    }

    Matrix Matrix::cov( const uint32_t dim ) const
    {
        const double N = (dim==0)?Rows():Cols();
        Matrix Cov = (1.0/N)*((dim==0)?this->transpose()*(*this) : (*this)*(this->transpose()) );
        return Cov;
    }

    const std::vector<double>& Matrix::GetArr() const
    {
        return m_arr;
    }

    std::vector<double>& Matrix::GetArr()
    {
        return m_arr;
    }

    void Matrix::AddVectorToColumns(const Matrix& mat)
    {
        uint32_t i,j;
        for(i=0; i<m_N; ++i)
        {
            for(j=0; j<m_M; ++j)
            {
                m_arr[m_M*i + j] += mat.GetArr()[j];
            }
        }
    }

    void Matrix::SetColumn(const uint32_t i, const Matrix& a)
    {
        for(uint32_t k=0; k<m_M;++k)
        {
            m_arr[m_M*i + k] = a.GetArr()[k];
        }
    }

    void Matrix::SetRow(const uint32_t i, const Matrix& a)
    {
        for(uint32_t k=0; k<m_N;++k)
        {
            m_arr[m_M*k + i] = a.GetArr()[k];
        }
    }

    Matrix Matrix::Sum( const uint32_t dimension) const
    {
        switch (dimension)
        {
            case 0:
            {
                Matrix out(1, m_N);
                uint32_t i,j;
                for(i=0;i<m_N;i++)
                {
                    for(j=0;j<m_M;j++)
                    {
                        out(i)+=m_arr[m_M*i + j];
                    }
                }
                return out;
            } break;
            case 1:
            {
                Matrix out(m_M, 1);
                uint32_t i,j;
                for(j=0;j<m_M;j++)
                {
                    for(i=0;i<m_N;i++)
                    {
                        out(j) += m_arr[m_M*i + j];
                    }
                }
                return out;
            } break;
            default:
                assert(false);
        }
    }

    void Matrix::Save( const std::string& filename) const
    {
        std::ofstream ostream(filename, std::ofstream::out);
        ostream << " # Created by Matrix " << std::endl;
        ostream << " # name: " << filename << std::endl;
        ostream << " # type: matrix " << std::endl;
        ostream << " # rows : " << m_M << std::endl;
        ostream << " # columns : " << m_N << std::endl;
        uint32_t i,j;
        for( i=0; i<m_M; i++ )
        {
            for( j=0; j<m_N; j++ )
            {
                ostream << " " << m_arr[m_M*j +i];
            }
            ostream << std::endl;
        }
    }

    uint32_t Matrix::Rows() const
    {
        return m_M;
    }

    uint32_t Matrix::Cols() const
    {
        return m_N;
    }

    uint32_t Matrix::numel() const
    {
        return m_M*m_N;
    }

    double& Matrix::operator () (uint32_t row, uint32_t col)
    {
        assert( row < m_M && col < m_N );
        return m_arr[col*m_M + row];
    }

    double Matrix::operator () (uint32_t row, uint32_t col) const
    {
        assert( row < m_M && col < m_N );
        return m_arr.at(col*m_M + row);
    }

    double& Matrix::operator () (uint32_t elem)
    {
        assert( elem < m_M*m_N );
        return m_arr[elem];
    }

    double  Matrix::operator () (uint32_t elem) const
    {
        assert( elem < m_M*m_N );
        return m_arr.at(elem);
    }

    Matrix& Matrix::operator += (const Matrix& m)
    {
        assert( m_N == m.Cols() && m_M == m.Rows() );
        for(uint32_t i=0;i<m_M*m_N;++i)
        {
            m_arr[i] += m(i);
        }
        return *this;
    }

    Matrix& Matrix::operator -= (const Matrix& m)
    {
        assert( m_N == m.Cols() && m_M == m.Rows() );
        for(uint32_t i=0;i<m_M*m_N;++i)
        {
            m_arr[i] -= m(i);
        }
        return *this;
    }

    Matrix& Matrix::operator *= (const Matrix& m)
    {
        assert( m_N == m.Rows());
        uint32_t i,j,k;
        std::vector<double> temp(m_M*m.Cols());

        for(i=0;i<m_M;++i)
        {
            for(j=0;j<m.Cols();++j)
            {
                for(k=0;k<m_N;++k)
                {
                    temp[j*m.Cols() + i] += m_arr[k*m_M + i]*m(k,j);
                }
            }
        }
        m_N = m.Cols();
        m_arr = temp;
        return *this;
    }

    Matrix& Matrix::operator *= (const double& c)
    {
        for(uint32_t i=0;i<m_M*m_N;++i)
        {
            m_arr[i]*=c;
        }
        return *this;
    }

    Matrix& Matrix::operator /= (const double& c)
    {
        for(uint32_t i=0;i<m_M*m_N;++i)
        {
            m_arr[i]/=c;
        }
        return *this;
    }

    Matrix& Matrix::operator ^= (const double& pow)
    {
        for(uint32_t i=0;i<m_M*m_N;++i)
        {
            m_arr[i] = std::pow(m_arr[i], pow);
        }
        return *this;
    }

    Matrix& Matrix::operator ^= (const Matrix& m)
    {
        assert(m_N = m.Cols() && m_M == m.Rows());
        for(uint32_t i=0;i<m_M*m_N;++i)
        {
            m_arr[i] = std::pow(m_arr[i], m(i));
        }
        return *this;
    }

    Matrix Matrix::transpose() const
    {
        Matrix out(m_N, m_M);
        uint32_t i,j;
        for(i=0;i<m_M;++i)
        {
            for(j=0;j<m_N;++j)
            {
                out(j,i) = m_arr[j*m_M + i];
            }
        }
        return out;
    }

    std::ostream& operator<<(std::ostream& s, const Matrix& m)
    {
        uint32_t i,j;
        s << std::endl;
        for(i=0;i<m.Rows();++i)
        {
            for(j=0;j<m.Cols();++j)
            {
                s << m(i,j);
                if(j<m.Cols()-1)
                {
                    s << ",";
                }
            }
            s << std::endl;
        }
        s << std::endl;
        return s;
    }

    Matrix operator +(const Matrix& A, const Matrix& B)
    {
        assert( A.m_M == B.m_M && A.m_N == B.m_N );
        Matrix C(A.m_M, A.m_N);
        uint32_t i,j;
        for(i=0;i<A.m_M;++i)
        {
            for(j=0;j<A.m_N;++j)
            {
                C.m_arr[j*A.m_M + i] = A.m_arr[j*A.m_M + i] + B.m_arr[j*A.m_M + i];
            }
        }
        return C;
    }


    Matrix operator -(const Matrix& A)
    {
        return (0.0-A);
    }

    Matrix operator -(const Matrix& A, const Matrix& B)
    {
        assert( A.m_M == B.m_M && A.m_N == B.m_N );
        Matrix C(A.m_M, A.m_N);
        uint32_t i,j;
        for(i=0;i<A.m_M;++i)
        {
            for(j=0;j<A.m_N;++j)
            {
                C.m_arr[j*A.m_M + i] = A.m_arr[j*A.m_M + i] - B.m_arr[j*A.m_M + i];
            }
        }
        return C;
    }

    Matrix operator *(const Matrix& A, const Matrix& B)
    {
        assert( A.m_N == B.m_M );
        Matrix C(A.m_M, B.m_N);
        uint32_t i,j,k;
        for(i=0;i<A.m_M;++i)
        {
            for(j=0;j<B.m_N;++j)
            {
                for(k=0;k<A.m_N;++k)
                {
                    C.m_arr[j*C.m_M + i] += A.m_arr[k*A.m_M + i]*B.m_arr[j*B.m_M + k];
                }
            }
        }
        return C;
    }

    Matrix operator *(const Matrix& A, const double B)
    {
        Matrix C(A.m_M, A.m_N);
        for(uint32_t i=0;i<C.numel(); ++i)
        {
            C(i) = A.Data()[i]*B;
        }
        return C;
    }

    Matrix operator /(const Matrix& A, const double B)
    {
        Matrix C(A.m_M, A.m_N);
        for(uint32_t i=0;i<C.numel(); ++i)
        {
            C(i) = A.Data()[i]/B;
        }
        return C;
    }

    Matrix operator /(const double B, const Matrix& A )
    {
        Matrix C(A.m_M, A.m_N);
        for(uint32_t i=0;i<C.numel(); ++i)
        {
            C(i) = B/A.Data()[i];
        }
        return C;
    }

    Matrix operator +(const Matrix& A, const double B)
    {
        Matrix C(A.m_M, A.m_N);
        for(uint32_t i=0;i<C.numel(); ++i)
        {
            C(i) = A.Data()[i]+B;
        }
        return C;
    }

    Matrix operator -(const Matrix& A, const double B)
    {
        Matrix C(A.m_M, A.m_N);
        for(uint32_t i=0;i<C.numel(); ++i)
        {
            C(i) = A.Data()[i]-B;
        }
        return C;
    }

    Matrix operator *(const double B,const Matrix& A)
    {
        return A*B;
    }

    Matrix operator +(const double B,const Matrix& A)
    {
        return A+B;
    }

    Matrix operator -(const double B, const Matrix& A)
    {
        Matrix C(A.m_M, A.m_N);
        for(uint32_t i=0;i<C.numel(); ++i)
        {
            C(i) = B-A.Data()[i];
        }
        return C;
    }

    void Matrix::SwapRow( const uint32_t i, const uint32_t j)
    {
        for(uint32_t k=0;k<m_N;k++)
        {
            std::swap(m_arr[m_M*k + i], m_arr[m_M*k + j]);
        }
    }

    inline int Matrix::pivot (uint32_t row)
    {
        uint32_t k = row;
        double amax,temp;

        amax = -1;
        for (uint32_t i=row; i < m_M; i++)
        {
            if ( (temp = std::fabs( m_arr[ i*m_M+row] )) > amax && temp != 0.0)
            {
                amax = temp;
                k = i;
            }
        }
        if (m_arr[ k*m_M+row] == 0.0)
        {
            return -1;
        }
        if (k != row)
        {
            SwapRow(k,row);
            return k;
        }
        return 0;
    }

    /*
    double Matrix::det() const
    {
        size_t i,j,k;
        double piv,detVal = 1.0;

        assert( m_M == m_N );

        Matrix temp(*this);

        for (k=0; k < m_M; k++)
        {
            int indx = temp.pivot(k);
            if (indx == -1)
            {
                return 0;
            }
            if (indx != 0)
            {
                detVal = - detVal;
            }
            detVal = detVal * temp(k,k);
            for (i=k+1; i < m_M; i++)
            {
                piv = temp(i,k) / temp(k,k);
                for (j=k+1; j < m_M; j++)
                {
                    temp(i,j) -= piv * temp(k,j);
                }
            }
        }
        return detVal;
    }
    */

    // This code was taken from paulborke.net
    double Matrix::det() const
    {
        assert(m_N == m_N);
        return det(*this, m_N);
    }
    double Matrix::det(const Matrix& a, uint32_t n) const
    {
        uint32_t i,j,j1,j2;
        double d = 0;

        if (n < 1)
        { /* Error */
            assert(false);
        } else if (n == 1)
        { /* Shouldn't get used */
          d = a(0,0);
        } else if (n == 2)
        {
          d = a(0,0) * a(1,1) - a(1,0) * a(0,1);
        } else
        {
          d = 0;
          for (j1=0;j1<n;j1++)
          {
             Matrix m(n-1, n-1);
             for (i=1;i<n;i++)
             {
                j2 = 0;
                for (j=0;j<n;j++)
                {
                   if (j == j1)
                   {
                      continue;
                   }
                   m(i-1,j2) = a(i,j);
                   j2++;
                }
             }
             d += std::pow(-1.0,1.0+j1+1.0) * a(0,j1) * det(m,n-1);
          }
        }
        return(d);
    }

    Matrix Matrix::inv() const
    {
        size_t i,j,k;
        double a1,a2;

        assert(m_M == m_N);

        Matrix temp(m_M,m_N);
        Matrix ti = *this;

        temp.SetIdentity();
        for (k=0; k < m_M; k++)
        {
            int indx = ti.pivot(k);
            if (indx == -1)
            {
                std::cout << " Warning: Matrix is singular " << std::endl;
                return *this;
            }

            if (indx != 0)
            {
                temp.SwapRow(k, indx);
            }
            a1 = ti(k,k);
            for (j=0; j <m_M; j++)
            {
                ti(k,j) /=a1;
                temp(k,j) /= a1;
            }
            for (i=0; i < m_M; i++)
            {
                if (i != k)
                {
                    a2 = ti(i,k);
                    for (j=0; j < m_M; j++)
                    {
                        ti(i,j) -= a2 * ti(k,j);
                        temp(i,j) -= a2 * temp(k,j);
                    }
                }
            }
        }
        return temp;
    }

    // from http://rosettacode.org/wiki/LU_decomposition
    std::tuple<Matrix, Matrix, Matrix> Matrix::LU() const
    {
        assert(m_M == m_N);
        Matrix U(m_M, m_N), L(m_M, m_N), P(m_M,m_N);
        P.SetIdentity();
        L.SetIdentity();
        uint32_t i,j,k;
        double s;
        for(i=0;i<m_M;++i)
        {
            uint32_t max_j=i;
            for(j=i;j<m_N;++j)
            {
                if (std::fabs(m_arr[m_M*i +j] ) > std::fabs(m_arr[m_M*i + max_j]))
                {
                    max_j = j;
                }
            }
            if(max_j != i)
            {
                P.SwapRow(i,max_j);
            }
        }

        Matrix Aprime = P*(*this);
        for(i=0; i<m_M;++i )
        {
            for(j=0; j<m_M;++j)
            {
                if(j<=i)
                {
                    s=0.0;
                    for(k=0;k<j;++k)
                    {
                        s+=L(j,k)*U(k,i);
                    }
                    U(j,i) = Aprime(j,i) - s;
                }
                if(j>=i)
                {
                    s=0.0;
                    for(k=0;k<i;++k)
                    {
                        s+=L(j,k)*U(k,i);
                    }

                    if(std::abs(U(i,i)) > 0)
                    {
                        L(j,i) = (Aprime(j,i)-s) / U(i,i);
                    }
                    else
                    {
                        L(j,i) = 0;
                    }
                }
            }
        }
        return std::make_tuple(L, U, P);
    }

    Matrix Matrix::solve( const Matrix& b ) const
    {
        assert(m_M == m_N && b.Rows()==m_N && b.Cols()==1);
        Matrix z(m_N, 1), x(m_N,1);
        auto lu = LU();
        auto L = std::get<0>(lu);
        auto U = std::get<1>(lu);
        auto P = std::get<2>(lu);
        auto bp = P*b;

        double sum;
        uint32_t i,p;
        for(i=0;i<m_N;i++)
        {
            sum=0;
            for(p=0;p<i;p++)
            {
                sum+=L(i,p)*z(p);
            }
            if(std::abs(L(i,i)) > 0)
            {
                z(i)=(bp(i)-sum)/L(i,i);
            }
            else
            {
                std::cerr << "Warning:: solve :: There was a zero divisor" << std::endl;
                z(i) = 0;
            }
        }

        for(i=m_M-1;i!=static_cast<uint32_t>(-1);i--)
        {
            sum=0;
            for(p=m_N-1;p!=i;p--)
            {
                sum+=U(i,p)*x(p);
            }
            if(std::abs(U(i,i)) > 0)
            {
                x(i)=(z(i)-sum)/U(i,i);
            }
            else
            {
                std::cerr << "Warning:: solve :: There was a zero divisor" << std::endl;
                x(i) = 0;
            }
        }

        return x;
    }

    inline double PYTHAG(double a, double b)
    {
        double at = std::fabs(a), bt = std::fabs(b), ct, result;

        if (at > bt)       { ct = bt / at; result = at * sqrt(1.0 + ct * ct); }
        else if (bt > 0.0) { ct = at / bt; result = bt * sqrt(1.0 + ct * ct); }
        else result = 0.0;
        return(result);
    }


    /*
     * svdcomp - SVD decomposition routine.
     * Takes an mxn matrix a and decomposes it into udv, where u,v are
     * left and right orthogonal transformation matrices, and d is a
     * diagonal matrix of singular values.
     *
     * This routine is adapted from svdecomp.c in XLISP-STAT 2.1 which is
     * code from Numerical Recipes adapted by Luke Tierney and David Betz.
     *
     * Input to dsvd is as follows:
     *   a = mxn matrix to be decomposed, gets overwritten with u
     *   m = row dimension of a
     *   n = column dimension of a
     *   w = returns the vector of singular values of a
     *   v = returns the right orthogonal transformation matrix
    */
    std::tuple<Matrix,Matrix,Matrix> Matrix::svd() const
    {
        bool dotranspose;
        uint32_t M, N;
        Matrix U;
        if(m_M < m_N)
        {
            dotranspose = true;
            M = m_N;
            N = m_M;
            U = this->transpose();
        }
        else
        {
            dotranspose = false;
            M = m_M;
            N = m_N;
            U = *this;
        }


        uint32_t flag, i, its, j, jj, k, l, nm;
        double c, f, h, s, x, y, z;
        double anorm = 0.0, g = 0.0, scale = 0.0;
        double *rv1 = new double[N];

        Matrix V( N, N);
        Matrix S( N, N);

         /* Householder reduction to bidiagonal form */
        for (i = 0; i <  N; i++)
        {
            /* left-hand reduction */
            l = i + 1;
            rv1[i] = scale * g;
            g = s = scale = 0.0;
            if (i <  M)
            {
                for (k = i; k <  M; k++)
                    scale += std::fabs(U(k,i));
                if (scale)
                {
                    for (k = i; k <  M; k++)
                    {
                        U(k,i) = U(k,i)/scale;
                        s += U(k,i) * U(k,i);
                    }
                    f = U(i,i);
                    g = -SIGN(sqrt(s), f);
                    h = f * g - s;
                    U(i,i) = (f - g);
                    if (i !=  N - 1)
                    {
                        for (j = l; j <  N; j++)
                        {
                            for (s = 0.0, k = i; k <  M; k++)
                                s += U(k,i) * U(k,j);
                            f = s / h;
                            for (k = i; k <  M; k++)
                                U(k,j) += (f * U(k,i));
                        }
                    }
                    for (k = i; k <  M; k++)
                        U(k,i) = U(k,i)*scale;
                }
            }
            S(i,i) = (scale * g);

            /* right-hand reduction */
            g = s = scale = 0.0;
            if (i <  M && i !=  N - 1)
            {
                for (k = l; k <  N; k++)
                    scale += std::fabs(U(i,k));
                if (scale)
                {
                    for (k = l; k <  N; k++)
                    {
                        U(i,k) = U(i,k)/scale;
                        s += U(i,k) * U(i,k);
                    }
                    f = U(i,l);
                    g = -SIGN(sqrt(s), f);
                    h = f * g - s;
                    U(i,l) = (f - g);
                    for (k = l; k <  N; k++)
                        rv1[k] = U(i,k) / h;
                    if (i !=  M - 1)
                    {
                        for (j = l; j <  M; j++)
                        {
                            for (s = 0.0, k = l; k <  N; k++)
                                s += U(j,k) * U(i,k);
                            for (k = l; k <  N; k++)
                                U(j,k) += (s * rv1[k]);
                        }
                    }
                    for (k = l; k <  N; k++)
                        U(i,k) = U(i,k)*scale;
                }
            }
            anorm = MAX(anorm, (std::fabs(S(i,i)) + fabs(rv1[i])));
        }

        /* accumulate the right-hand transformation */
        for (i =  N - 1; i != static_cast<uint32_t>(-1); i--)
        {
            if (i <  N - 1)
            {
                if (g)
                {
                    for (j = l; j <  N; j++)
                        V(j,i) = ((U(i,j) / U(i,l)) / g);
                        /* double division to avoid underflow */
                    for (j = l; j <  N; j++)
                    {
                        for (s = 0.0, k = l; k <  N; k++)
                            s += (U(i,k) * V(k,j));
                        for (k = l; k <  N; k++)
                            V(k,j) += (s * V(k,i));
                    }
                }
                for (j = l; j <  N; j++)
                    V(i,j) = V(j,i) = 0.0;
            }
            V(i,i) = 1.0;
            g = rv1[i];
            l = i;
        }

        /* accumulate the left-hand transformation */
        for (i =  N - 1; i != static_cast<uint32_t>(-1); i--)
        {
            l = i + 1;
            g = S(i,i);
            if (i <  N - 1)
                for (j = l; j <  N; j++)
                    U(i,j) = 0.0;
            if (g)
            {
                g = 1.0 / g;
                if (i !=  N - 1)
                {
                    for (j = l; j <  N; j++)
                    {
                        for (s = 0.0, k = l; k <  M; k++)
                            s += (U(k,i) * U(k,j));
                        f = (s / U(i,i)) * g;
                        for (k = i; k <  M; k++)
                            U(k,j) += (f * U(k,i));
                    }
                }
                for (j = i; j <  M; j++)
                    U(j,i) = U(j,i)*g;
            }
            else
            {
                for (j = i; j <  M; j++)
                    U(j,i) = 0.0;
            }
            ++U(i,i);
        }

        /* diagonalize the bidiagonal form */
        for (k =  N - 1; k != static_cast<uint32_t>(-1); k--)
        {                             /* loop over singular values */
            for (its = 0; its < 30; its++)
            {                         /* loop over allowed iterations */
                flag = 1;
                for (l = k; l != static_cast<uint32_t>(-1); l--)
                {                     /* test for splitting */
                    nm = l - 1;
                    if (std::fabs(rv1[l]) + anorm == anorm)
                    {
                        flag = 0;
                        break;
                    }
                    if (std::fabs(S(nm,nm)) + anorm == anorm)
                        break;
                }
                if (flag)
                {
                    c = 0.0;
                    s = 1.0;
                    for (i = l; i <= k; i++)
                    {
                        f = s * rv1[i];
                        if (std::fabs(f) + anorm != anorm)
                        {
                            g = S(i,i);
                            h = PYTHAG(f, g);
                            S(i,i) = h;
                            h = 1.0 / h;
                            c = g * h;
                            s = (- f * h);
                            for (j = 0; j <  M; j++)
                            {
                                y = U(j,nm);
                                z = U(j,i);
                                U(j,nm) = (y * c + z * s);
                                U(j,i) = (z * c - y * s);
                            }
                        }
                    }
                }
                z = S(k,k);
                if (l == k)
                {                  /* convergence */
                    if (z < 0.0)
                    {              /* make singular value nonnegative */
                        S(k,k) = (-z);
                        for (j = 0; j <  N; j++)
                            V(j,k) = (-V(j,k));
                    }
                    break;
                }

                assert(its <= 30);

                /* shift from bottom 2 x 2 minor */
                x = S(l,l);
                nm = k - 1;
                y = S(nm,nm);
                g = rv1[nm];
                h = rv1[k];
                f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
                g = PYTHAG(f, 1.0);
                f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;

                /* next QR transformation */
                c = s = 1.0;
                for (j = l; j <= nm; j++)
                {
                    i = j + 1;
                    g = rv1[i];
                    y = S(i,i);
                    h = s * g;
                    g = c * g;
                    z = PYTHAG(f, h);
                    rv1[j] = z;
                    c = f / z;
                    s = h / z;
                    f = x * c + g * s;
                    g = g * c - x * s;
                    h = y * s;
                    y = y * c;
                    for (jj = 0; jj <  N; jj++)
                    {
                        x = V(jj,j);
                        z = V(jj,i);
                        V(jj,j) = (x * c + z * s);
                        V(jj,i) = (z * c - x * s);
                    }
                    z = PYTHAG(f, h);
                    S(j,j) = z;
                    if (z)
                    {
                        z = 1.0 / z;
                        c = f * z;
                        s = h * z;
                    }
                    f = (c * g) + (s * y);
                    x = (c * y) - (s * g);
                    for (jj = 0; jj <  M; jj++)
                    {
                        y = U(jj,j);
                        z = U(jj,i);
                        U(jj,j) = (y * c + z * s);
                        U(jj,i) = (z * c - y * s);
                    }
                }
                rv1[l] = 0.0;
                rv1[k] = f;
                S(k,k) = x;
            }
        }
        delete [] rv1;
        if(!dotranspose)
        {
            return std::tie(U,S,V);
        }
        else
        {
            return std::tie(V,S,U);
        }
    }

    Matrix Matrix::Column(const uint32_t i) const
    {
        assert(i<m_N);
        Matrix ret(m_M, 1);
        for(uint32_t k=0;k<m_M;++k)
        {
            ret(k) = m_arr[m_M*i+k];
        }
        return ret;
    }

    Matrix Matrix::Row(const uint32_t i) const
    {
        assert(i<m_M);
        Matrix ret(1,m_N);
        for(uint32_t k=0; k<m_N; ++k)
        {
            ret(k) = m_arr[m_M*k+i];
        }
        return ret;
    }

    Matrix Matrix::Block(const uint32_t imin, const uint32_t jmin,
                 const uint32_t imax, const uint32_t jmax) const
    {
        assert( imin<m_M && imax < m_M && jmin < m_N && jmax < m_N );
        Matrix ret(imax-imin+1,jmax-jmin+1);
        uint32_t i,j;
        for(i=imin; i<=imax; i++)
        {
            for(j=jmin; j<=jmax; j++)
            {
                ret(i-imin,j-jmin) = m_arr[m_M*j + i];
            }
        }
        return ret;
    }

    Matrix Matrix::SolveZero() const
    {
        if(m_M < m_N)
        {
            Matrix temp = (*this).transpose()*(*this);
            Matrix v = std::get<2>(temp.svd());
            return v.Column(v.Rows()-1);
        }
        else
        {
            Matrix v = std::get<2>(svd());
            return v.Column(v.Rows()-1);
        }
    }

    Matrix Matrix::pinv() const
    {
        if( m_M <= m_N )
        {
            return (*this).transpose()*(((*this)*((*this).transpose())).inv());
        }
        else
        {
            return (((*this).transpose())*(*this)).inv()*((*this).transpose());
        }
    }

    double Matrix::Norm() const
    {
        double sum=0.0;
        for(auto s : m_arr)
        {
            sum+=(s*s);
        }
        return std::sqrt(sum);
    }

    inline Matrix MatrixMinor(const Matrix& x, uint32_t d)
    {
        Matrix m(x.Rows(), x.Cols());
        for (uint32_t i = 0; i < d; i++)
            m(i,i) = 1;
        for (uint32_t i = d; i < x.Rows(); i++)
            for (uint32_t j = d; j < x.Cols(); j++)
                m(i,j) = x(i,j);
        return m;
    }

    void Matrix::SetZero()
    {
        for(auto& s : m_arr)
        {
            s=0.0;
        }
    }

    Matrix vmul(const Matrix& v, uint32_t n)
    {
        Matrix x(n, n);
        for (uint32_t i = 0; i < n; i++)
            for (uint32_t j = 0; j < n; j++)
                x(i,j) = -2 *  v(i) * v(j);
        for (uint32_t i = 0; i < n; i++)
            x(i,i) += 1;

        return x;
    }

    // This code was adapted from rosettacode.net
    std::pair<Matrix, Matrix> Matrix::QR() const
    {
        std::vector<Matrix> q(m_M);
        Matrix z = *this, z1;
        Matrix e(m_M,1), x(m_M,1);
        double a, enorm;
        uint32_t l;
        Matrix Q(m_M,m_M), R(m_M,m_N);

        for (uint32_t k = 0; k < m_N && k < m_M - 1; k++)
        {
            a=0;
            e.SetZero();
            x.SetZero();

            z1 = MatrixMinor(z, k);
            z = z1;

            x = z.Column(k);
            a = x.Norm();
            if (m_arr[m_M*k+k] > 0) a = -a;

            e(k)=1;
            for(l=0;l<m_M;++l)
            {
                e(l) = x(l) + e(l)*a;
            }
            enorm = e.Norm();
            for(l=0;l<m_M;++l)
            {
                e(l) = e(l)/enorm;
            }

            q[k] = vmul(e, m_M);
            z1 = q[k]*z;
            z = z1;
        }
        Q = q[0];
        R = q[0]*(*this);
        for (uint32_t i = 1; i < m_N && i < m_M - 1; i++)
        {
            z1 =q[i]*Q;
            Q = z1;
        }
        z = Q*(*this);
        R = z;
        Q = Q.transpose();
        return std::make_pair(Q,R);
    }

    // This code is an implementation of the wikipedia article
    std::pair<Matrix, Matrix> Matrix::eig() const
    {
        assert( m_M == m_N );
        const double tol = 1e-6;
        const uint32_t MaxIter = 300;
        Matrix prevval(m_M,m_N,std::numeric_limits<double>::max());
        double err = std::numeric_limits<double>::max();
        uint32_t iter=0;
        Matrix A = *this;
        std::pair<Matrix, Matrix> qr;
        Matrix U(m_N,m_N);
        Matrix subtemp(m_M,m_N);
        U.SetIdentity();
        while(err > tol && iter < MaxIter)
        {
            qr = A.QR();
            U*=qr.first; // This involves creating new memory
            prevval = A;
            Multiply(qr.second,qr.first,A);
            Subtract(A,prevval,subtemp);
            err = (subtemp).Norm();
            iter++;
        }
        Matrix S(m_M,m_M);
        for(uint32_t k=0;k<m_M;++k)
        {
            S(k,k) = A(k,k);
        }
        return std::make_pair(U, S);
    }

    double Matrix::mineig() const
    {
        auto us = eig();
        double mineig= std::numeric_limits<double>::max();
        for(uint32_t k=0;k<m_M;++k)
        {
            if(us.second(k,k)<mineig)
                mineig = us.second(k,k);
        }
        return mineig;
    }

    #define ROTATE(A,i,j,k,l) g=A(i,j);h=A(k,l);A(i,j)=g-s*(h+g*tau);\
        A(k,l)=h+s*(g-h*tau);

    // This code is an adaptation of the numerical recipes version
    // Originall written by
    //  David Squire (DMS), David.Squire@infotech.monash.edu.au
    // but adapted for the matrix format
    //
    std::pair<Matrix, Matrix> Matrix::Jacobi() const
    {
        assert(m_M == m_N);
        /* Computes all eigenvalues and eigenvectors of a real symmetric matrix
           a[0..n-1][0..n-1]. On output, elements of a above the diagonal are
           destroyed. d[0..n-1] returns the eigenvalues of a. v[0..n-1][0..n-1] is
           a matrix whose columns contain, on output, the normalized eigenvectors
           of a. nrot returns the number of Jacobi rotations that were required.
        */
        Matrix A(*this);
        Matrix U(m_M, m_N), S(m_M, m_N);
        U.SetIdentity();
        std::vector<double> b(m_M), z(m_M);
        int nrot;
        int32_t j, iq, ip, i;
        double tresh, theta, tau, t, sm, s, h, g, c;

        /* initialise b and d to the diagonal of a */
        for (ip = 0; ip < static_cast<int32_t>(m_M); ip++) {
            b[ip] = S(ip,ip) = A(ip,ip);
            z[ip] = 0.0;
        }
        nrot = 0;
        for (i = 1; i <= 50; i++) {
            /* sum off-diagonal elements */
            /* DMS: surely this is just the lower triangle? */
            sm = 0.0;
            for (ip = 0; ip < static_cast<int32_t>(m_N)-1; ip++) {
                for (iq = ip+1; iq < static_cast<int32_t>(m_N); iq++)
                    sm +=  std::fabs(A(ip,iq));
            }
            if (sm == 0.0)
            {
                return std::make_pair(U,S);
            }
            if (i < 4)
            {
                tresh = 0.2*sm/(m_N*m_M);	/* On the first three sweeps... */
            }
            else
            {
                tresh = 0.0;			/* ...thereafter */
            }
            for (ip = 0; ip < static_cast<int32_t>(m_N)-1; ip++)
            {
                for (iq = ip+1; iq < static_cast<int32_t>(m_N); iq++)
                {
                    g = 100.0*std::fabs(A(ip,iq));
                    /* After four sweeps, skip the rotation if the off-diagonal
                       element is small */
                    if (i > 4 && std::fabs(S(ip,ip))+g == std::fabs(S(ip,ip))
                        && std::fabs(S(iq,iq))+g == std::fabs(S(iq,iq)))
                    {
                        A(ip,iq) = 0.0;
                    }
                    else if (std::fabs(A(ip,iq)) > tresh)
                    {
                        h = S(iq,iq)-S(ip,ip);
                        if (std::fabs(h)+g == std::fabs(h))
                        {
                            t = (A(ip,iq))/h;
                        }
                        else {
                            theta = 0.5*h/(A(ip,iq));
                            t = 1.0/(std::fabs(theta)+std::sqrt(1.0+theta*theta));
                            if (theta < 0.0)
                                t = -t;
                        }
                        c = 1.0/sqrt(1+t*t);
                        s = t*c;
                        tau = s/(1.0+c);
                        h = t*A(ip,iq);
                        z[ip] -=  h;
                        z[iq] +=  h;
                        S(ip,ip) -=  h;
                        S(iq,iq) +=  h;
                        A(ip,iq) = 0.0;
                        /* Case of rotations 0 <= j < p */
                        for (j = 0; j <= ip-1; j++)
                        {
                            ROTATE(A, j, ip, j, iq)
                        }
                        /* Case of rotations p < j < q */
                        for (j = ip+1; j <= iq-1; j++)
                        {
                            ROTATE(A, ip, j, j, iq)
                        }
                        /* Case of rotations q < j < n */
                        for (j = iq+1; j < static_cast<int32_t>(m_N); j++)
                        {
                            ROTATE(A, ip, j, iq, j)
                        }
                        for (j = 0; j < static_cast<int32_t>(m_N); j++)
                        {
                            ROTATE(U, j, ip, j, iq)
                        }
                        ++nrot;
                    }
                }
            }
            for (ip = 0; ip < static_cast<int32_t>(m_N); ip++) {
                b[ip] +=  z[ip];
                A(ip,ip) = b[ip];
                z[ip] = 0.0;
            }
        }
        std::cerr << "Too many iterations in routine JACOBI" << std::endl;
        return std::make_pair(U,S);
    }

    #undef ROTATE

    // This is an adaptation of the rosettacode.net implementation.
    // The adaptation for positive semidefinite matrices is take from
    // 'Singular Values using Cholesky Decomposition' - Krishnamoorthy and Kocagoez
    // 'Analysis of the Cholesky Decomposition of a Semi-definite Matrix' - N. J. Higham
    Matrix Matrix::chol() const
    {
        assert(m_N == m_M);
        Matrix L(m_N, m_N);
        uint32_t i,j,k;
        double s,ii;
        for (i = 0; i < m_N; i++)
        {
            for (j = 0; j < (i+1); j++)
            {
                s = 0;
                for (k = 0; k < j; k++)
                {
                    s += L(i * m_N + k) * L(j * m_N + k);
                }
                ii = std::sqrt(m_arr[i * m_N + i] - s);
                L(i * m_N + j) = (i == j) ?
                               ii :
                               (1.0 / L(j * m_N + j) * (m_arr[i * m_N + j] - s));

            }
        }

        return L;
    }

    double Matrix::dot(const Matrix& c)
    {
        double s=0;
        assert( c.numel() == m_arr.size() );
        for(uint32_t k=0;k<c.numel();++k)
        {
            s+=m_arr[k]*c(k);
        }
        return s;
    }

    void Matrix::SetIdentity()
    {
        uint32_t i,j;
        for(i=0;i<m_M;++i)
        {
            for(j=0;j<m_N;++j)
            {
                m_arr[m_M*j+i] = (i==j) ? 1.0 : 0.0;
            }
        }
    }
}
