#include "Matrix.h"

Cu::Matrix ConvertToEigen(const double *in, const uint32_t M, const uint32_t N )
{
    Cu::Matrix out(M,N);
    out = Eigen::Map<Cu::Matrix>(const_cast<double *>(in), M, N);
    return out;
}

namespace Cu
{
void print_matrix_octave(const Cu::Matrix& A, const std::string& name)
{
    std::ofstream ostream(name);
    //ostream.open(name);
    uint32_t k=name.find_last_of('/');
    ostream << " # Created by CLib" << std::endl;
    ostream << " # name: " << name.substr(k+1,name.size()-k-1) << std::endl;
    ostream << " # type: matrix " << std::endl;
    ostream << " # rows : " << A.rows() << std::endl;
    ostream << " # columns : " << A.cols() << std::endl;
    for( uint32_t i=0; i<A.rows(); i++ )
    {
        for( uint32_t j=0; j<A.cols(); j++ )
        {
            ostream << " " << A(i,j);
        }
        ostream << std::endl;
    }
     ostream.close();
}

void print_matrix_octave(const std::vector<std::vector<double> > &A, const std::string& name)
{
    std::ofstream ostream(name);
    //ostream.open(name);

    uint32_t k=name.find_last_of('/');
    ostream << " # Created by CLib" << std::endl;
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
        ostream << " # Created by CLib" << std::endl;
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

void print_matrix_octave(const uint32_t* data, const uint32_t M, const uint32_t N, const uint32_t K, const std::string& name)
{
    uint32_t k=name.find_last_of('/');
        std::ofstream ostream;
        ostream.open(name);
        ostream << " # Created by CLib" << std::endl;
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
        ostream << " # Created by CLib" << std::endl;
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

void transpose(std::vector<double>& A, uint32_t M, uint32_t N)
{
    std::vector<double> B(A);
#pragma omp parallel for
    for(uint32_t i=0;i<M;++i)
    for(uint32_t j=0;j<N;++j)
    {
        A[i*N + j] = B[j*M + i];
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

#ifdef DEBUG
    void print_matrix(const Cu::Matrix& P, const std::string& name)
    {
        std::cout << name << std::endl;
        for( int i=0; i<P.rows(); i++)
        {
            std::cout << " ";
            for( int j=0; j<P.cols(); j++)
            {
                std::cout << P(i,j) << " ";
            }
            std::cout << std::endl;
        }
    }
#else
    void print_matrix(const Cu::Matrix&, const std::string&)
    {
    }
#endif
}

namespace MatrixFuncs
{
    void HomogeneousNormalise( Cu::Matrix& in )
    {
        Cu::Matrix bottomRow = in.bottomRows(1);
        Cu::Matrix brep = bottomRow.replicate(in.rows(),1);
        in = in.cwiseQuotient(brep);
    }

    std::vector<uint32_t> find( const Cu::Matrix& in, const std::function<bool(float)>& func  )
    {
        std::vector<uint32_t> out(in.cols()*in.rows());
        uint32_t ind=0;
        for(uint32_t i=0; i<in.cols()*in.rows(); i++)
        {
            if(func(in.data()[i]))
            {
                out[ind++]=i;
            }
        }
        out.resize(ind);
        return out;
    }

    std::vector<uint32_t> unique( uint32_t* vec, uint32_t M )
    {
        std::vector<uint32_t> vector(vec, vec+M);
        std::sort(vector.begin(), vector.end());
        std::vector<uint32_t> out(M);
        out[0]=vector[0];
        uint32_t l,k;
        for(l=0,k=1; k<M; ++k)
        {
            if(out[l] != vector[k])
            {
                out[l++] = vector[k];
            }
        }
        out.resize(l);
        return out;
    }

    std::pair<Cu::Matrix, Cu::Matrix> EigenValues( const Cu::Matrix& A  )
    {
        Eigen::EigenSolver<Cu::Matrix> es(A);
        return std::make_pair(es.eigenvectors().real() / es.eigenvectors().real().determinant(), es.eigenvalues().real());
    }

    std::pair<Cu::ComplexMatrix, Cu::ComplexMatrix> EigenValuesLPCK( const Cu::Matrix& A  )
    {
        int32_t n = A.rows ();
        if(n != A.cols())
        {
            std::cerr << "EigenValuesLPCK :: The matrix must be square" << std::endl;
            assert(0);
        }

        int32_t info = 0;

        std::vector<double> atmp(A.data(), A.data()+A.rows()*A.cols());
        double *tmp_data = atmp.data();

        std::vector<double> wr(n);
        double *pwr = wr.data ();
        std::vector<double> wi (n);
        double *pwi = wi.data ();

        const bool calc_ev = true; // calculate eigenvectors

        int32_t tnvr = calc_ev ? n : 0;
        Cu::Matrix vr (tnvr, tnvr);
        double *pvr = vr.data ();

        int32_t lwork = -1;
        double dummy_work;
        double *dummy = 0;
        int32_t idummy = 1;

        const char cN = 'N';
        const char vN = 'V';

        DGEEV( &cN,&vN,n, tmp_data, n, pwr, pwi, dummy,idummy, pvr, n, &dummy_work, lwork, info, 1,1);

        Cu::ComplexMatrix lambda(n,1);
        Cu::ComplexMatrix v(n, n);

        if (info == 0)
        {
            lwork = static_cast<int32_t> (dummy_work);
            std::vector<double> work (lwork);
            double *pwork = work.data ();
            DGEEV(&cN,&vN,n, tmp_data, n, pwr, pwi, dummy,idummy, pvr, n, pwork, lwork, info, 1, 1);

            if (info < 0)
            {
                std::cerr << "EigenValuesLPCK :: unrecoverable error in dgeev" << std::endl;
                assert(0);
            }

            if (info > 0)
            {
                std::cerr << "EigenValuesLPCK :: dgeev failed to converge" << std::endl;
                assert(0);
            }



            // This loop ignores the comoplex parts of the eigen values and eigen vectors
            for (int32_t j = 0; j < n; j++)
            {
                if (wi[j] == 0.0)
                {
                    lambda (j) = std::complex<double>(wr[j], 0.0);
                    for (int32_t i = 0; i < n; i++)
                        v (i, j) = std::complex<double>(vr (i, j),0.0);
                }
                else
                {
                    if (j+1 >= n)
                    {
                        std::cerr << "EigenValuesLPCK :: internal error" << std::endl;
                        assert(0);
                    }

                    lambda(j) = std::complex<double> (wr[j], wi[j]);
                    lambda(j+1) = std::complex<double> (wr[j+1], wi[j+1]);

                    for (int32_t i = 0; i < n; i++)
                    {
                        double real_part = vr (i, j);
                        double imag_part = vr (i, j+1);
                        v (i, j) = std::complex<double> (real_part, imag_part);
                        v (i, j+1) = std::complex<double> (real_part, -imag_part);
                    }
                    j++;
                }
            }
        }

        return std::make_pair(v,lambda);
    }

    std::pair<Cu::ComplexMatrix, Cu::ComplexMatrix> EigenValuesARPK( const Cu::SPMatrix& A, int32_t K, const bool takehighest )
    {
        uint32_t N = A.rows();

        const uint32_t maxit = 500; // The maximum number of iterations
        uint32_t mode = 1; // The mode

        Eigen::SPQR<Cu::SPMatrix > *solver = NULL;
        if(takehighest)
        {
            mode=1;
        }
        else
        {
            // We need to compute the inverse for the mode 3 driver
            solver = new Eigen::SPQR<Cu::SPMatrix >(A);
            if(solver->info()!=Eigen::Success)
            {
                std::cerr << " EigenValuesARPK :: There was an error in QR decomposition " << std::endl;
                assert(0);
            }
            mode=3;
        }

        const double tol = std::numeric_limits<double>::epsilon();
        const double sigmar =0.0;
        const double sigmai = 0.0;

        int32_t p = K*2+1;
        if(p<20) p = 20;
        if(p>N-1) p = N-1;

        char bmat = 'I';

        std::vector<double> resid(N, 0.0);
        for(auto& r: resid) r = static_cast<double>(rand()%100)/100;

        std::vector<int32_t> ip(11,0);
        ip[0] = 1; // ishift
        ip[2] = maxit;
        ip[3] = 1;
        ip[6] = mode;
        int32_t* iparam = ip.data();

        std::vector<int32_t> iptr(14,0.0);
        int32_t* ipntr = iptr.data();

        int32_t ido=0;
        //int32_t iter=0;
        int32_t lwork = 3 * p * (p + 2);

        Cu::Vector workd(3*N + 1);
        std::vector<double> v(N*(p+1)), workl(lwork+1);
        double *vptr=v.data();
        double *worklptr = workl.data();
        double *workdptr = workd.data();

        int32_t info=0;

        double* presid = resid.data();

        std::string typ;
        if(takehighest)
        {
            typ = "LM"; // largest
        }
        else
        {
            typ = "LM"; // smallest
        }

        do
        {
            dnaupd_ (ido, &bmat, N, typ.c_str(), K, tol, presid, p, vptr, N, iparam,
                        ipntr, workdptr, worklptr, lwork, info, 1,2);

            if(ido == -1 || ido == 1 || ido == 2)
            {
                switch(ido)
                {
                    case -1:
                    {
                        if(takehighest)
                        {
                            workd.segment(iptr[1] - 1, N) = A * workd.segment(iptr[0] - 1, N);
                        }
                        else
                        {
                            workd.segment(iptr[1] - 1, N) = solver->solve( workd.segment(iptr[0] - 1, N));
                        }
                    } break;
                    case 1 :
                    {
                        std::copy(workdptr + iptr[0] -1, workdptr + iptr[0] -1 +N, workdptr + iptr[2] -1);
                        if(takehighest)
                        {
                            workd.segment(iptr[1] - 1, N) = A * workd.segment(iptr[2] - 1, N);
                        }
                        else
                        {
                            workd.segment(iptr[1] - 1, N) = solver->solve( workd.segment(iptr[2] - 1, N));
                        }
                    } break;
                    case 2 :
                    {
                        std::copy(workdptr + iptr[0] -1, workdptr + iptr[0] -1 +N, workdptr + iptr[1] -1);
                    } break;
                }
            }
            else
            {
                if(info < 0)
                {
                    std::cerr << "EigenValuesARPK:: Error " << info << " from DNAUPD" << std::endl;
                    assert(0);
                }
                break;
            }
        } while(true);

        int32_t info2;

        std::vector<int32_t> s(p);
        int32_t *sptr=s.data();

        Cu::Matrix evec2(N, K+1);
        evec2.setZero();
        Cu::ComplexMatrix eval(K+1,1);
        eval.setZero();
        double* evecptr = evec2.data();
        std::complex<double>* evalptr = eval.data();
        Cu::ComplexMatrix evec(N,K);
        evec.setZero();

        std::string Achr = "A";
        bool rvec=true;

        std::vector<double> dr(K+1,0.0), di(K+1,0.0), workev(3*p);
        double* drptr = dr.data();
        double* diptr = di.data();
        double* workevptr = workev.data();

        dneupd_(rvec, Achr.c_str(), sptr, drptr, diptr, evecptr, N, sigmar, sigmai, workevptr,
               &bmat, N, typ.c_str (), K, tol, presid, p, vptr, N, iparam,
               ipntr, workdptr, worklptr, lwork, info2, 1,1,2);


        if(info2==0)
        {
            int32_t jj=0;
            for(int32_t i=0; i<K+1; ++i)
            {
                if(drptr[i]==0.0 && diptr[i]==0.0 && jj==0)
                    ++jj;
                else
                    evalptr[i-jj] = std::complex<double>( drptr[i], diptr[i] );
            }
            if(jj==0 && !rvec)
            {
                for(int32_t i=0; i<K; ++i)
                    evalptr[i] = evalptr[i+1];
            }
            for(int32_t i=0; i<K/2; ++i)
            {
                std::complex<double> dtmp = evalptr[i];
                evalptr[i] = evalptr[K - i - 1];
                evalptr[K - i - 1] = dtmp;
            }

            eval.conservativeResize(K,1);
            evalptr = eval.data();

            if(rvec)
            {
                std::vector<double> dtmp(N);
                for(int32_t i=0; i< K/2; ++i)
                {
                    int32_t off1 = i*N;
                    int32_t off2 = (K - i - 1) * N;
                    if(off1 == off2)
                        continue;

                    for(int32_t j=0; j<N; ++j)
                        dtmp[j] = evecptr[off1 + j];

                    for(int32_t j=0; j<N; ++j)
                        evecptr[off1 + j] = evecptr[off2 + j];

                    for(int32_t j=0; j<N; ++j)
                        evecptr[off2 + j] = dtmp[j];
                }

                evec.conservativeResize(N,K);
                int32_t i=0;
                while(i<K)
                {
                    int32_t off1 = i*N;
                    int32_t off2 = (i+1)*N;
                    if(std::imag(evalptr[i])==0)
                    {
                        for(int32_t j=0; j<N; ++j)
                            evec(j,i) = std::complex<double>(evecptr[j+off1], 0.0);
                        ++i;
                    }
                    else
                    {
                        for(int32_t j=0; j<N; ++j)
                        {
                            evec(j,i) = std::complex<double>(evecptr[j+off1], evecptr[j+off2]);
                            if(i < K-1)
                                evec(j,i+1) = std::complex<double>(evecptr[j+off1], -evecptr[j+off2]);
                        }
                        i+=2;
                    }
                }
            }
        }
        else
        {
            std::cerr << " EigenValuesARPK:: error " << info2 << " in DNEUPD" << std::endl;
            assert(0);
        }

        if(!takehighest)
        {
            delete solver;
        }

        return std::make_pair(evec, eval);
    }

    std::pair<std::vector<Cu::BlasVector>, std::vector<double>> EigenValues( const Cu::SparseMatrix& A, bool takehighest, uint32_t N )
    {
        uint32_t M = A.size1();

        typedef ietl::vectorspace<Cu::BlasVector> Vecspace;
        typedef boost::lagged_fibonacci607 Gen;

        Vecspace vec(M);
        Gen mygen;
        ietl::lanczos<Cu::SparseMatrix,Vecspace> lanczos(A,vec);

        // Creation of an iteration object:
        int max_iter = 200;
        double rel_tol = 500*std::numeric_limits<double>::epsilon();
        double abs_tol = std::pow(std::numeric_limits<double>::epsilon(),2./3);

        std::vector<double> eigen;
        std::vector<double> err;
        //std::vector<int> multiplicity;

        try
        {
            if( takehighest )
            {
                ietl::lanczos_iteration_nhighest<double> iter(max_iter, N, rel_tol, abs_tol);
                lanczos.calculate_eigenvalues(iter,mygen);
            }
            else
            {
                ietl::lanczos_iteration_nlowest<double> iter(max_iter, N, rel_tol, abs_tol);
                lanczos.calculate_eigenvalues(iter,mygen);
            }
        }
        catch (std::runtime_error& e)
        {
            std::cerr <<  "Error:: Eigenvalues::Sparse:: " << e.what() << std::endl;
        }

        eigen = lanczos.eigenvalues();
        err = lanczos.errors();
        auto multiplicity = lanczos.multiplicities();


        std::vector<double>::iterator start = eigen.begin();
        std::vector<double>::iterator end = eigen.begin()+N;
        std::vector<Cu::BlasVector> eigenvectors; // for storing the eigen vectors.
        ietl::Info<double> info; // (m1, m2, ma, eigenvalue, residualm, status).

        try
        {
            lanczos.eigenvectors( start, end, std::back_inserter(eigenvectors), info, mygen );
        }
        catch (std::runtime_error& e)
        {
            std::cerr <<  "Error:: Eigenvalues::Sparse:: " << e.what() << std::endl;
        }

        for(int i = 0; i < info.size(); i++)
        {
            if(info.residual(i) > 1e-10)
            {
                std::cerr << " Warning:: Eigenvalues::Sparse The residual for vector "<<i << "is " << info.residual(i) <<std::endl;
            }
        }

        return std::make_pair(eigenvectors, eigen);
    }

    /*
    void EigenValuesSS(const Cu::SPMatrix& A )
    {
        Cu::SPMatrix currentA = A;

        Eigen::SparseQR<Cu::SPMatrix, Eigen::COLAMDOrdering<int> > solver(A);

        auto Q = solver.matrixQ();
        auto R = solver.matrixR();
        Cu::Matrix AA;

        while(true)
        {
            Eigen::SparseQR<Cu::SPMatrix, Eigen::COLAMDOrdering<int> > solver2(Q*AA);
        }


    }*/

     Cu::Matrix diag(const Cu::Matrix& x)
     {
         const uint32_t N = x.rows()*x.cols();
         Cu::Matrix ret(N,N);
         ret.setZero();
         for(uint32_t k=0;k<N; ++k)
         {
             ret(k,k) = x(k);
         }
         return ret;
     }

     Cu::Matrix Cos(const Cu::Matrix& x)
     {
         Cu::Matrix ret(x.rows(),x.cols());
         for(uint32_t k=0; k<ret.rows()*ret.cols(); ++k)
         {
             ret(k) = std::cos(x(k));
         }
         return ret;
     }

     Cu::Matrix Sin(const Cu::Matrix& x)
     {
         Cu::Matrix ret(x.rows(),x.cols());
         for(uint32_t k=0; k<ret.rows()*ret.cols(); ++k)
         {
             ret(k) = std::sin(x(k));
         }
         return ret;
     }

     Cu::Matrix ones(const uint32_t i, const uint32_t j)
     {
         Cu::Matrix ret(i,j);
         ret.setOnes();
         return ret;
     }

     Cu::Matrix tril(const uint32_t i, const uint32_t j)
     {
         Cu::Matrix ret(i,j);
         ret.setZero();
         for(uint32_t k,l=0; l<i; ++l)
         for(k=l; k<j; ++k)
         {
            ret(k,l) = 1;
         }
         return ret;
     }

     Cu::Matrix triu(const uint32_t i, const uint32_t j)
     {
         Cu::Matrix ret(i,j);
         ret.setZero();
         for(uint32_t k,l=0; l<i; ++l)
         for(k=l; k<j; ++k)
         {
            ret(l,k) = 1;
         }
         return ret;
     }
}

namespace MeshFuncs
{
    template<typename T>
    inline void EraseFromList( std::vector<std::vector<T>>& list, const std::vector<bool>& indicators )
    {
        uint32_t revfit =  list.size()-1;
        uint32_t remit = indicators.size()-1;
        for(; true; --revfit)
        {
            if(indicators[remit])
            {
                list.erase( list.begin()+revfit );
            }
            if(revfit==0) break;
            --remit;
        }
    }

    void RemoveVerticesFaces( VerticesFaces& vf, std::vector<uint32_t>& VerticesToRemove )
    {
        std::vector<bool> removev(vf.first.size(),false);
        std::vector<bool> removef(vf.second.size(),false);
        std::vector<uint32_t> newFaceNames(vf.first.size(),0);

        //VerticesToRemove.sort();
        std::sort(VerticesToRemove.begin(), VerticesToRemove.end());

        for(auto vtr : VerticesToRemove)
        {
            removev[vtr] = true;
        }

        uint32_t i=0, index=0;
        for(auto& nfn : newFaceNames)
        {
            if(!removev[i++])
                nfn = index++;
        }

        auto bf = removef.begin();
        for(auto& f : vf.second)
        {
            for(auto& fi : f)
            {
                if(removev[fi])
                {
                    *bf = true;
                    break;
                }
            }
            ++bf;
        }

        // Remove faces
        EraseFromList<uint32_t>(vf.second,removef);

        for(auto& f1 : vf.second)
            for(auto& f2 : f1)
              f2 = newFaceNames[f2];

        // Remove vertices
        EraseFromList<double>(vf.first,removev);
    }

    void WritePly(const VerticesFaces& vf, const std::string filename)
    {
        std::ofstream f(filename, std::ofstream::out);
        f << "ply" << std::endl;
        f << "format ascii 1.0" << std::endl;
        f << "comment CLib::Core generated" << std::endl;
        f << "element vertex " << vf.first.size() << std::endl;
        f << "property float x" << std::endl;
        f << "property float y" << std::endl;
        f << "property float z" << std::endl;
        if(vf.second.size() > 0)
        {
            f << "element face " << vf.second.size() << std::endl;
            f << "property list uchar int vertex_indices" << std::endl;
        }
        f << "end_header" << std::endl;

        for(auto& v : vf.first)
        {
            f << v[0] << " " << v[1] << " " << v[2] << std::endl;
        }

        if(vf.second.size() > 0)
        {
            for(auto& fa : vf.second)
            {
                f << 3 << " " << fa[0] << " " << fa[1] << " " << fa[2] << std::endl;
            }
        }
        f.close();
    }

     void WritePly(const VerticesFaces& vf, const std::vector< std::vector<double>>& vertexcolours, const std::string filename)
     {
         std::ofstream f(filename, std::ofstream::out);
         f << "ply" << std::endl;
         f << "format ascii 1.0" << std::endl;
         f << "comment CLib::Core generated" << std::endl;
         f << "element vertex " << vf.first.size() << std::endl;
         f << "property float x" << std::endl;
         f << "property float y" << std::endl;
         f << "property float z" << std::endl;
         f << "property uchar red" << std::endl;
         f << "property uchar green" << std::endl;
         f << "property uchar blue" << std::endl;

         if(vf.second.size() > 0)
         {
             f << "element face " << vf.second.size() << std::endl;
             f << "property list uchar int vertex_indices" << std::endl;
         }
         f << "end_header" << std::endl;
         auto vcit = vertexcolours.begin();
         for(auto& v : vf.first)
         {
             f << v[0] << " " << v[1] << " " << v[2] << " ";

             for(auto& vc : *vcit)
                 f << vc << " ";

             f << std::endl;
             ++vcit;
         }
         if(vf.second.size() > 0)
         {
             for(auto& fa : vf.second)
             {
                 f << 3 << " " << fa[0] << " " << fa[1] << " " << fa[2] << std::endl;
             }
         }
         f.close();
     }

     inline double ComputeDistance(const std::vector<double>& vertices, const uint32_t M, const uint32_t i, const uint32_t j)
     {
         double d=0.0;
         for(uint32_t k=0; k<M; ++k)
         {
             d += ( vertices[ M*i +k ] - vertices[M*j + k] ) * ( vertices[ M*i +k ] - vertices[M*j + k] );
         }
         return std::sqrt(d);
     }

     inline void ComputeDistances( std::vector<double>& distances, const uint32_t i,
                                   const std::vector<double>& vertices, const uint32_t M,
                                   const std::vector<std::vector<uint32_t>>& nbs,
                                   const std::vector<uint32_t>& vertexclass)
     {
         uint32_t ind,k;
         double d;
         for(k=0; k<nbs[i].size(); ++k)
         {
             if(nbs[i].size() > 0)
             {
                 ind = nbs[i][k];
                 d = ComputeDistance( vertices, M, i, ind ) + distances[i];
                 if( vertexclass[ind] < 2 && d < distances[ind] )
                 {
                    distances[ind] = d;
                 }
             }
         }
     }

     inline void InsertElementsToSet( DistanceSet& diset, const uint32_t ind,
                                      const std::vector<double>& distances,
                                      const std::vector<std::vector<uint32_t>>& nbs,
                                      const std::vector<uint32_t>& vertexclass)
     {
        uint32_t k,kind;
        for(k=0; k<nbs[ind].size(); ++k)
        {
            kind = nbs[ind][k];
            if(vertexclass[kind] < 2)
                diset.insert({ kind, distances[ kind ] });
        }
     }

     std::vector< double > Dijkstra(const std::vector<double>& vertices, const uint32_t M,
                                     const std::vector<std::vector<uint32_t>>& nbs,
                                     const uint32_t i, const int32_t j)
     {
        std::vector<uint32_t> vertexclass( vertices.size()/M, 0 ); // 0-unvisited, 1-fringe, 2-visited
        DistanceSet diset( [&]( const DistanceIndex& l, const DistanceIndex& r ){
            if(l.distance == r.distance)
            {
                return l.index < r.index;
            }
            else
            {
                return l.distance < r.distance;
            }
        }); // This should order the elements by their distance but maintain equality based on the index.

        std::vector<double> distances( vertices.size()/M, std::numeric_limits<double>::infinity() );
        distances[i] = 0;

        diset.insert({0,0.0});

        InsertElementsToSet( diset, i, distances, nbs, vertexclass );
        while (diset.size() > 0)
        {
            auto dit = diset.begin();
            vertexclass[dit->index] = 2;

            if(dit->index == j) // If j is a visited node return
                return distances;

            if(dit->index > 0 && dit->index == static_cast<uint32_t>(j))
            {
                break;
            }

            // Compute new distances for fringe vertices
            ComputeDistances( distances, dit->index, vertices, M, nbs, vertexclass);

            // Add the elements to the set.
            InsertElementsToSet( diset, dit->index, distances, nbs, vertexclass);

            // Erase the beginning vertex
            diset.erase( dit );
        }

        return distances;
     }

     std::list<uint32_t> FindShortestPath(const std::vector<std::vector<uint32_t> >& neighbours,
                                           const std::vector<double>& distances,
                                           const uint32_t i, const uint32_t j)
     {
        uint32_t nind,minindex,h,l=j;
        double min;
        std::list<uint32_t> path;
        path.push_back(j);
        while( l!=i )
        {
            min=std::numeric_limits<double>::max();
            minindex=l;
            for( h=0; h<neighbours[l].size(); ++h )
            {
                nind = neighbours[l][h];
                if( distances[nind] < min )
                {
                    min = distances[nind];
                    minindex = nind;
                }
            }

            if( minindex == l )
            {
                std::cerr << "error:: FindShortestPath:: \
                              Minimum index is not the target. \
                              The graph must be undirected." <<std::endl;
                return path;
            }
            else
            {
                path.push_back( minindex );
                l = minindex;
            }
        }
        path.reverse();
        return path;
     }
}


std::set<uint32_t> FindConnectedVertices(const uint32_t vertex , const Faces &faces)
{
    std::set<uint32_t> vertexset;
    for( const auto& f : faces )
    {
        if(f[0]==vertex || f[1] == vertex || f[2] == vertex)
        {
            vertexset.insert(f[0]);
            vertexset.insert(f[1]);
            vertexset.insert(f[2]);
        }
    }
    return vertexset;
}

std::set<uint32_t> FindConnectedVertices(const std::list<uint32_t> &connectedfacelist , const Faces &faces)
{
    std::set<uint32_t> vertexset;
    for(const auto& c : connectedfacelist )
    {
        for(uint8_t k=0; k<3; ++k)
            vertexset.insert(faces[c][k]);
    }
    return vertexset;
}

std::list<uint32_t> FindConnectedFaces( const uint32_t vertex, const Faces& faces )
{
    std::list<uint32_t> faceset;
    uint32_t index=0;
    for( const auto& f : faces )
    {
        if(f[0]==vertex || f[1] == vertex || f[2] == vertex)
        {
            faceset.push_back(index);
        }
        ++index;
    }
    return faceset;
}

std::list<uint32_t> FindQuadrilateral( const uint32_t v1, const uint32_t v2, const std::list<uint32_t> connectedfacelist, const Faces& faces )
{
    std::list<uint32_t> vertices;
    vertices.push_back(v1);
    vertices.push_back(v2);
    for(auto cfl : connectedfacelist)
    {
        const Face& face = faces[cfl];
        if( (face[0] == v1 || face[1] == v1 || face[2] == v1) &&
            (face[0] == v2 || face[1] == v2 || face[2] == v2))
        {
            for(uint8_t k=0; k<3; ++k)
            {
                if(face[k]!=v1 && face[k]!=v2)
                    vertices.push_back(face[k]);
            }
        }
    }
    return vertices;
}


namespace Funcs
{
    double normal_pdf(const double x, const double m, const double s)
    {
        static const double inv_sqrt_2pi = 0.398942280401433;
        double a = (x - m) / s;

        return inv_sqrt_2pi / s * std::exp(-0.5f * a * a);
    }

    int hsv2rgb(const float h, const float s , const float v, float &r, float &g, float &b)
    {
       /*
        * Purpose:
        * Convert HSV values to RGB values
        * All values are in the range [0.0 .. 1.0]
        */
       float S, H, V, F, M, N, K;
       int   I;

       S = s;  /* Saturation */
       H = h;  /* Hue */
       V = v;  /* value or brightness */

       if ( S == 0.0 ) {
          /*
           * Achromatic case, set level of grey
           */
          r = V;
          g = V;
          b = V;
       } else {
          /*
           * Determine levels of primary colours.
           */
          if (H >= 1.0) {
             H = 0.0;
          } else {
             H = H * 6;
          } /* end if */
          I = (int) H;   /* should be in the range 0..5 */
          F = H - I;     /* fractional part */

          M = V * (1 - S);
          N = V * (1 - S * F);
          K = V * (1 - S * (1 - F));

          if (I == 0) { r = V; g = K; b = M; }
          if (I == 1) { r = N; g = V; b = M; }
          if (I == 2) { r = M; g = V; b = K; }
          if (I == 3) { r = M; g = N; b = V; }
          if (I == 4) { r = K; g = M; b = V; }
          if (I == 5) { r = V; g = M; b = N; }
       } /* end if */

       return 0;
    } /* end function hsv2rgb */
}
