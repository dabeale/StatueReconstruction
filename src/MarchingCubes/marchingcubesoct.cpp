
#include "marchingcubes.h"
#include <octave/oct.h>
#include <vector>
#include <functional>

DEFUN_DLD(ImplicitSurface, args, , 
          "[v f] = ImplicitSurface(Data, N, isoval, verbose) \n \
          Input: \n \
             Data     - 1xN The data values (column major) \n \
             N        - The number of points per dimension")
{
    int nargin = args.length ();
    if( nargin > 4 || nargin < 3 )
    {
      print_usage();
    }
    
    int verbose=0;
    if( nargin == 4 )
    {
      verbose = args(3).scalar_value();
    }

    if(verbose==1)
    {
        std::cout << "Marching Cubes" << std::endl;
    }

    const Matrix Data = args(0).array_value();
    const double N = args(1).scalar_value ();
    const double isoval = args(2).scalar_value();
    
    if(verbose == 1)
    {
        std::cout << "N: " << N << std::endl;
    } 
    
    const double min[3] = {0.0,0.0,0.0};
    const double max[3] = {N-1, N-1, N-1};

    if( Data.rows()*Data.cols() != N*N*N )
    {
      std::cerr << "Incorrect number of elements" << std::endl;
      return octave_value(0);
    }
        
    MC::MarchingCubes mc(isoval, N, max, min);
    
    auto fun = [&] ( const double* val ) 
    {
      return Data.data()[ 
                static_cast<uint32_t>(
                  N*N*std::floor(val[2]) + 
                  N*std::floor(val[1]) + 
                  std::floor(val[0])) ];
    };
    
    mc.March(fun);
    
    std::vector< MC::Pointd3 > pts(mc.GetPoints());
    std::vector< MC::Pointi3 > fc(mc.GetFaces());

    Matrix v(pts.size(),3);
    Matrix f(fc.size(),3);
    v.fill(0.0);
    f.fill(0.0);

    
    //#pragma omp parallel for
    for(uint32_t j=0;j<3;j++)
    {
      uint32_t i=0;
      for( const auto& p : pts )
      {
        v(i,j) = p.p[j];
        i++;
      }
    }
    
    //#pragma omp parallel for
    for(uint32_t j=0;j<3;j++)
    {
      uint32_t i=0;
      for(const auto& p : fc)
      {
        f(i,j) = p.p[j];
        i++;
      }
    }
    
    octave_value_list ovl;
    ovl(0) = v;
    ovl(1) = f;
    return ovl;
}
