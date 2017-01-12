
#include "Rotation.h"

namespace NN
{
    inline double Norm( double* r, uint32_t M)
    {
        double out=0.0;
        for(M--; M>0; M--)
        {
            out += r[M]*r[M];
        }
        return std::sqrt(out + r[0]*r[0]);
    }

    inline double Dot( double* r, double* c, uint32_t M )
    {
        double out=0.0;
        for(M--;M>0;M--)
        {
            out += c[M]*r[M];
        }
        return out + c[0]*r[0];
    }

    void RandomRotation( double* R, const uint32_t M, std::default_random_engine& e,
                         std::uniform_real_distribution<double>& discrete )
    {
        uint32_t i,j,k;
        for(i=0;i<M;++i)
        {
            for(j=0;j<M;++j)
            {
                R[i +M*j]=discrete(e);
            }
        }

        double nu, ku;
        for(i=0;i<M;++i)
        {
            nu = Norm(&R[M*i], M);
            for(k=0;k<M;++k)
            {
                R[M*i+k] /= nu;
            }
            nu = Dot(&R[M*i],&R[M*i], M);
            for(j=i+1;j<M;++j)
            {
                ku = Dot(&R[M*i],&R[M*j], M);
                for(k=0;k<M;++k)
                {
                    R[M*j+k] -= (ku / nu) * R[M*i+k];
                }
            }
        }
    }
}
