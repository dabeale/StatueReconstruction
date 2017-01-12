#include "Utils.h"

namespace Utils
{
void print_matrix_octave(const Math::Matrix& A, const std::string& name)
{
    std::ofstream ostream(name);
    //ostream.open(name);
    uint32_t k=name.find_last_of('/');
    ostream << " # Created by FLib" << std::endl;
    ostream << " # name: " << name.substr(k+1,name.size()-k-1) << std::endl;
    ostream << " # type: matrix " << std::endl;
    ostream << " # rows : " << A.Rows() << std::endl;
    ostream << " # columns : " << A.Cols() << std::endl;
    for( uint32_t i=0; i<A.Rows(); i++ )
    {
        for( uint32_t j=0; j<A.Cols(); j++ )
        {
            ostream << " " << A(i,j);
        }
        ostream << std::endl;
    }
     ostream.close();
}

void print_matrix_octave(const std::vector< Math::Matrix >& A, const std::string& name)
{
    if(A.size()==0)
    {
        std::cerr <<"Utils::print_matrix_octave:: A must have non-zero size" <<std::endl;
        assert(false);
    }
    std::ofstream ostream(name);
    //ostream.open(name);
    uint32_t k=name.find_last_of('/');
    ostream << " # Created by FLib" << std::endl;
    ostream << " # name: " << name.substr(k+1,name.size()-k-1) << std::endl;
    ostream << " # type: matrix " << std::endl;
    ostream << " # ndims : 3" << std::endl;
    ostream << " " << A[0].Rows() << " " << A[0].Cols() << " " << A.size() << std::endl;

    for(auto & a : A)
    {
        if( a.Rows() != A[0].Rows() || a.Cols() != A[0].Cols())
        {
            std::cerr << "Utils::print_matrix_octave:: All matrices should have the same dimensions" << std::endl;
            assert(false);
        }
        for( uint32_t i=0; i<a.Rows(); i++ )
        {
            for( uint32_t j=0; j<a.Cols(); j++ )
            {
                ostream << " " << a(i,j) << "\n";
            }
        }
    }
    ostream << std::endl;
    ostream.close();
}

void print_matrix_octave(const std::vector<std::vector<double> > &A, const std::string& name)
{
    std::ofstream ostream(name);
    //ostream.open(name);

    uint32_t k=name.find_last_of('/');
    ostream << " # Created by FLib" << std::endl;
    ostream << " # name: " << name.substr(k+1,name.size()-k-1)  << std::endl;
    ostream << " # type: matrix " << std::endl;
    ostream << " # rows : " << A[0].size() << std::endl;
    ostream << " # columns : " << A.size() << std::endl;
    for( uint32_t i=0; i<A.size(); i++ )
    {
        for( uint32_t j=0; j<A[0].size(); j++ )
        {
            ostream << " " << A[i][j];
        }
        ostream << std::endl;
    }
    ostream.close();
}

void print_matrix_octave(const double* data, const uint32_t M, const uint32_t N, const std::string& name, bool binary)
{
    uint32_t k=name.find_last_of('/');
    if(binary)
    {
        std::ofstream ostream(name, std::ofstream::binary);
        ostream << "Octave-1-L";
        unsigned char tt = 0;
        ostream.write(reinterpret_cast<char *>(&tt), 1);
        //ostream.open(name);
        std::string nm = name.substr(k+1,name.size()-k-1);
        int32_t nlength = nm.length();
        ostream.write(reinterpret_cast<char *>(&nlength), 4);
        ostream << nm;

        std::string doc("");
        int32_t dlength = doc.length();
        ostream.write(reinterpret_cast<char *>(&dlength), 4);
        ostream << doc;

        unsigned char mark_as_global;
        mark_as_global = false;
        ostream.write (reinterpret_cast<char *> (&mark_as_global), 1);

        unsigned char tmp = 255;
        ostream.write (reinterpret_cast<char *> (&tmp), 1);

        std::string typ("matrix");
        int32_t len = typ.length();
        ostream.write (reinterpret_cast<char *> (&len), 4);
        const char *btmp = typ.data ();
        ostream.write (btmp, len);

        int32_t NDims = -2;
        ostream.write (reinterpret_cast<char *> (&NDims), 4);

        int32_t tmpint = M;
        ostream.write (reinterpret_cast<char *> (&tmpint), 4);
        tmpint = N;
        ostream.write (reinterpret_cast<char *> (&tmpint), 4);

        char tmp_type = 7; //This means double type (it is an enumeration in octave)
        ostream.write (&tmp_type, 1);

        for( uint32_t j=0; j<M*N; j++ )
        {
            double d = data[j];
            ostream.write(reinterpret_cast<char *>(&d),sizeof(double));
        }
        ostream.close();
    }
    else
    {
        std::ofstream ostream(name);
        //ostream.open(name);
        ostream << " # Created by FLib" << std::endl;
        ostream << " # name: " <<  name.substr(k+1,name.size()-k-1)  << std::endl;
        ostream << " # type: matrix " << std::endl;
        ostream << " # rows : " << M << std::endl;
        ostream << " # columns : " << N << std::endl;
        for( uint32_t j=0; j<M; j++ )
        {
            for( uint32_t i=0; i<N; i++ )
            {
                //ostream.write((char *)&data[M*i + j],sizeof(double));
                ostream << " " << data[M*i+j];
            }
            ostream << std::endl;
        }
        ostream.close();
    }
}


void print_matrix_octave(const uint32_t* data, const uint32_t M, const uint32_t N, const std::string& name, bool binary)
{
    uint32_t k=name.find_last_of('/');
    std::ofstream ostream(name);
    //ostream.open(name);
    ostream << " # Created by FLib" << std::endl;
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

void print_matrix_octave(const uint32_t* data, const uint32_t M, const uint32_t N, const uint32_t K, const std::string& name)
{
    uint32_t k=name.find_last_of('/');
        std::ofstream ostream;
        ostream.open(name);
        ostream << " # Created by FLib" << std::endl;
        ostream << " # name: " << name.substr(k+1,name.size()-k-1)  << std::endl;
        ostream << " # type: matrix "   << std::endl;
        ostream << " # ndims: 3 "       << std::endl;
        ostream << " " << M << " " << N << " " << K << std::endl;

        for(uint32_t k=0; k < M*N*K; ++k)
        {
            ostream << " " << data[k];
        }
        ostream << std::endl;
        ostream.close();
}

void print_matrix_octave(const double* data, const uint32_t M, const uint32_t N, const uint32_t K, const std::string& name, bool binary)
{
    uint32_t k=name.find_last_of('/');
    if(binary)
    {
        std::ofstream ostream(name, std::ofstream::binary);
        ostream << "Octave-1-L";
        unsigned char tt = 0;
        ostream.write(reinterpret_cast<char *>(&tt), 1);
        //ostream.open(name);
        std::string nm = name.substr(k+1,name.size()-k-1);
        int32_t nlength = nm.length();
        ostream.write(reinterpret_cast<char *>(&nlength), 4);
        ostream << nm;

        std::string doc("");
        int32_t dlength = doc.length();
        ostream.write(reinterpret_cast<char *>(&dlength), 4);
        ostream << doc;

        unsigned char mark_as_global;
        mark_as_global = false;
        ostream.write (reinterpret_cast<char *> (&mark_as_global), 1);

        unsigned char tmp = 255;
        ostream.write (reinterpret_cast<char *> (&tmp), 1);

        std::string typ("matrix");
        int32_t len = typ.length();
        ostream.write (reinterpret_cast<char *> (&len), 4);
        const char *btmp = typ.data ();
        ostream.write (btmp, len);

        int32_t NDims = -3;
        ostream.write (reinterpret_cast<char *> (&NDims), 4);

        int32_t tmpint = M;
        ostream.write (reinterpret_cast<char *> (&tmpint), 4);
        tmpint = N;
        ostream.write (reinterpret_cast<char *> (&tmpint), 4);
        tmpint = K;
        ostream.write (reinterpret_cast<char *> (&tmpint), 4);

        char tmp_type = 7; //This means double type (it is an enumeration in octave)
        ostream.write (&tmp_type, 1);

        for( uint32_t j=0; j<M*N*K; j++ )
        {
            double d = data[j];
            ostream.write(reinterpret_cast<char *>(&d),sizeof(double));
        }
        ostream.close();
    }
    else
    {
        std::ofstream ostream;
        ostream.open(name);
        ostream << " # Created by FLib" << std::endl;
        ostream << " # name: " << name.substr(k+1,name.size()-k-1)  << std::endl;
        ostream << " # type: matrix "   << std::endl;
        ostream << " # ndims: 3 "       << std::endl;
        ostream << " " << M << " " << N << " " << K << std::endl;

        for(uint32_t k=0; k < M*N*K; ++k)
        {
            ostream << " " << data[k];
        }
        ostream << std::endl;
        ostream.close();
    }
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

Math::Matrix read_matrix_bundler( const std::string& filename )
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Matrix::read_matrix_bundler :: Cannot open the file." << std::endl;
        return Math::Matrix(3,4);
    }

    std::string line;
    std::getline (file,line);
    if(line.find("CONTOUR")==std::string::npos)
    {
        std::cerr << "Matrix::read_matrix_bundler :: The file is not a bundler file." << std::endl;
        return Math::Matrix(3,4);
    }

    Math::Matrix P(3,4);
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

#ifdef DEBUG
    void print_matrix(const Math::Matrix& P, const std::string& name)
    {
        std::cout << name << std::endl;
        for( int i=0; i<P.Rows(); i++)
        {
            std::cout << " ";
            for( int j=0; j<P.Cols(); j++)
            {
                std::cout << P(i,j) << " ";
            }
            std::cout << std::endl;
        }
    }
#else
    void print_matrix(const Math::Matrix&, const std::string&)
    {
    }
#endif
}
