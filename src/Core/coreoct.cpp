
#include <octave/oct.h>
#include <octave/Cell.h>
#include <vector>

#include <stdio.h>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>

#include "Sample.h"
#include "Matrix.h"

void handler(int sig) {
  void *array[10];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 10);

  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDERR_FILENO);
  exit(1);
}

DEFUN_DLD(SampleEdgesGibbsResampling, args, , "[X] = SampleEdges( EdgeMu, EdgeSigma, P, K, deviation) ")
{
  signal(SIGSEGV, handler);   // install our handler

  int nargin = args.length ();
  if( nargin != 5 )
  {
    print_usage();
  }

  const Cell EdgeMu = args(0).cell_value();
  const Cell EdgeSigma = args(1).cell_value();
  const Cell P = args(2).cell_value();
  const uint32_t K = args(3).scalar_value();
  const double deviation = args(4).scalar_value();
    
  Matrix Temp;
  std::vector< std::vector< std::pair<double, double> > > Means(EdgeMu.numel());
  for(uint32_t s,k=0; k<EdgeMu.numel(); ++k)
  {
    Temp = EdgeMu(k).matrix_value();
    Means[k] = std::vector< std::pair<double, double> >(Temp.rows());
    for(s=0;s<Temp.rows(); ++s)
    {
      Means[k][s].first = Temp(0,s);
      Means[k][s].second = Temp(1,s);
    }
  }
  
  std::vector< std::vector< double > > Covs(EdgeSigma.numel());
  for(uint32_t s,k=0; k<EdgeMu.numel(); ++k)
  {
    Temp = EdgeSigma(k).matrix_value();
    Covs[k] = std::vector< double >(Temp.rows());
    for(s=0;s<Temp.rows(); ++s)
    {
      Covs[k][s] = Temp(s);
    }
  }
  
  std::vector< Cu::Matrix > Ps(P.numel(), Cu::Matrix(3,4));
  for(uint32_t k=0; k<P.numel(); ++k)
  {
    Matrix POTemp = P(k).matrix_value();
    for(uint32_t i,j=0;j<4;++j)
    for(i=0;i<3;++i)
      Ps[k](i,j) = POTemp(i,j);
  }
  
  Sample::Generator g;
  std::cout << "Start resampling... " << std::endl;
  auto X = g.Sample3DEdgeResampling(Means, Covs, Ps, K, deviation);
  std::cout << "Finished resampling. " << std::endl;
  
  Matrix out(3,X.size());
  for(uint32_t k=0;k<X.size();++k)
  {
    out(0,k) = X[k](0);
    out(1,k) = X[k](1);
    out(2,k) = X[k](2);
  }
  
  return octave_value(out);
}

DEFUN_DLD(SampleEdgesGibbs, args, , "[X] = SampleEdges( Edge, P, K ) ")
{
  signal(SIGSEGV, handler);   // install our handler

  int nargin = args.length ();
  if( nargin != 3 )
  {
    print_usage();
  }
  
  std::cout << "This method id not yet ported to octave" << std::endl;
}
