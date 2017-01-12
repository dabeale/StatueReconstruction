
#include "tree.h"
#include <octave/oct.h>
#include <vector>
#include <stdio.h>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>

static std::vector<KD::Tree> trees;

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

DEFUN_DLD(KDTreeCreate, args, , "[KDTree] = KDTreeCreate(Data, Maxdepth, verbose)")
{
    signal(SIGSEGV, handler);   // install our handler
    
    int nargin = args.length ();
    if( nargin > 3 || nargin ==0 )
    {
      print_usage();
      return octave_value(0);
    }
    
    int verbose=0;
    if( nargin == 3 )
    {
      verbose = args(2).scalar_value();
    }
    
    uint32_t Maxdepth = static_cast<uint32_t>(-1);
    if( nargin >= 2)
    {
      Maxdepth =  static_cast<uint32_t>(args(1).scalar_value());
    }

    if(verbose==1)
    {
        std::cout << "Running KDTree" << std::endl;
    }

    const Matrix Data = args(0).array_value ();
    if(verbose == 1)
    {
        std::cout << "N: " << Data.cols() <<
                     ", Dims: " << Data.rows() <<
                     ", MaxDepth: " << Maxdepth << std::endl;
    } 
    trees.push_back(KD::Tree(Data.data(), Data.rows(), Data.cols(), false, Maxdepth));
    return octave_value(trees.size()-1);
}

DEFUN_DLD(KDTreeHistogram, args, , "[centers hist] = KDTreeSearch( KDTree, Points)")
{
    signal(SIGSEGV, handler);   // install our handler
    int nargin = args.length ();
    if( nargin > 2 )
    {
      print_usage();
      return octave_value(0);
    }
    
    const uint32_t kdtree = static_cast<uint32_t>(args(0).scalar_value ());
    const Matrix Points = args(1).array_value();
    const uint32_t N = Points.cols();
    const KD::Tree& kdt = trees[kdtree];
    const uint32_t K = std::pow(2,kdt.GetMaxDepth()+1)-1;
    std::cout << "N:" << N << "  K:" << K << std::endl;
    if(K > kdt.GetDataNumberOfPoints())
    {
      std::cerr << "KDHistogram: The value of K is greater than N/2" << std::endl;
      return octave_value(0);
    }
    
    auto histmap = kdt.Histogram(Points.data(), N);
    
    Matrix histograms(1, K);
    Matrix centers(1, K);
    uint32_t j=0;
    for(auto& q : histmap)
    {
      histograms(0,j) = q.second;
      centers(0,j) = q.first + 1;
      ++j;
    }

    octave_value_list ovl;
    ovl(0) = centers;
    ovl(1) = histograms;
    return ovl;
}

DEFUN_DLD(KDTreeSearch, args, , "[distances, ids] = KDTreeSearch( KDTree, Points, K, verbose )")
{
    signal(SIGSEGV, handler);   // install our handler
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
        std::cout << "Searching KDTree" << std::endl;
    }
    
    const uint32_t kdtree = args(0).scalar_value ();
    const Matrix Points = args(1).array_value();
    const uint32_t K = args(2).scalar_value();
    
    if(verbose == 1)
    {
        std::cout << "kdindex: " << kdtree
                  << ", K: " << K <<
                     ", NPoints:" << Points.cols() << 
                     ", Tree Size" << trees.size() << std::endl;
    }
    
    Matrix OutIndexes(K, Points.cols());
    OutIndexes.fill(0.0);
    Matrix OutDistances(K, Points.cols());
    OutDistances.fill(0.0);

    //#pragma omp parallel for
    for(uint32_t j,k=0; k<Points.cols(); k++)
    {
      const KD::Tree& kdt = trees[kdtree];
      auto res = kdt.Search(&Points.data()[k*Points.rows()], K);
      
      j=0;
      for(auto& r : res)
      {
        OutDistances.fortran_vec()[K*k + j] = r.distance;
        OutIndexes.fortran_vec()[K*k + j] = r.index + 1;
        ++j;
        if(j>=K) break;
      }
    }
    octave_value_list ovl;
    ovl(0) = OutIndexes;
    ovl(1) = OutDistances;
    return ovl;
}

DEFUN_DLD(KDTreeDeleteAll, args, , "KDTreeDeleteAll()")
{
  signal(SIGSEGV, handler);   // install our handler
  trees.clear();
  return octave_value(0);
}
