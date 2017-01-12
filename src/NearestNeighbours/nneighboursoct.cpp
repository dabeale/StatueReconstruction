
#include "NearestNeighbours.h"
#include <octave/oct.h>
#include <vector>
#include <set>
#include <stdio.h>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>

static std::vector<NearestNeighbours::NearestNeighbours*> trees;
static std::vector<std::vector<double>> datas;


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

DEFUN_DLD(NNCreate, args, , "[NN] = NNCreate(Data, type, perms, verbose) \
                                 Data - MxN matrix \
                                 type - 0 : randomised kdtree \
                                        1 : kdtree \
                                        2 : randomised morton \
                                 perms - number of random trees \
                                 verbose")
{
    signal(SIGSEGV, handler);   // install our handler
    int nargin = args.length ();
    if( nargin > 4 || nargin <3 )
    {
      print_usage();
      return octave_value(-1);
    }
    int verbose=0;
    if( nargin == 4 )
    {
      verbose = args(2).scalar_value();
    }

    if(verbose==1)
    {
        std::cout << "Running KDTree" << std::endl;
    }
  
    const uint32_t ti = static_cast<uint32_t>(args(1).scalar_value());
    const uint32_t pi = static_cast<uint32_t>(args(2).scalar_value());
    Matrix Data = args(0).array_value ();
    datas.reserve(10);
    datas.push_back(std::vector<double>(Data.data(), Data.data()+Data.rows()*Data.cols()));
    
    if(verbose == 1)
    {
        std::cout << "N: " << Data.cols() <<
                     ", Dims: " << Data.rows() << std::endl;
    } 
    
    uint32_t i = datas.size()-1;
    
    switch(ti)
    {
      case 0 : 
        trees.push_back( new NearestNeighbours::NearestNeighbours(
          NearestNeighbours::RandomTrees, 
          datas[i].data(), Data.rows(), Data.cols(), pi));
        break;
      case 1 :
        trees.push_back( new NearestNeighbours::NearestNeighbours(
          NearestNeighbours::KDTree, 
          datas[i].data(), Data.rows(), Data.cols(), pi));
        break;
      case 2 :
        trees.push_back( new NearestNeighbours::NearestNeighbours(
          NearestNeighbours::RandomisedMorton, 
          datas[i].data(), Data.rows(), Data.cols(), pi));
      break;
      case 3 :
        trees.push_back( new NearestNeighbours::NearestNeighbours(
          NearestNeighbours::KDMap, 
          datas[i].data(), Data.rows(), Data.cols(), pi));
      break;
      case 4 :
        trees.push_back( new NearestNeighbours::NearestNeighbours(
          NearestNeighbours::Direct, 
          datas[i].data(), Data.rows(), Data.cols(), pi));
      break;
      default:
        std::cerr << "NearestNeighbours: This method is not implemented" << std::endl;
    }
    
    
    return octave_value(i);
}
/*
DEFUN_DLD(NNSearchRadius, args, , "ids = NNSearchRadius( KDTree, Points, radius, verbose )")
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
      return octave_value(-1);
    }

    if(verbose==1)
    {
        std::cout << "Searching KDTree" << std::endl;
    }
    
    const uint32_t kdtree = static_cast<uint32_t>(args(0).scalar_value ());
    const Matrix Points = args(1).array_value();
    const double radius = args(2).scalar_value();
    
    if(verbose == 1)
    {
        std::cout << "kdindex: " << kdtree
                  << ", radius: " << radius <<
                     ", NPoints:" << Points.cols() << 
                     ", Tree Size" << trees.size() << std::endl;
    }
    
   

    auto res = trees[kdtree]->NearestToPointCloud(Points.data(), Points.cols(), radius);
    Matrix OutIndexes(1, res.size());
    OutIndexes.fill(0.0);
    
    uint32_t i=0;
    for(const auto r : res)
    {
      OutIndexes.fortran_vec()[i++] = r + 1;
    }
    
    return octave_value(OutIndexes);
}
*/
DEFUN_DLD(NNSearch, args, , "[distances, ids] = NNSearch( KDTree, Points, K, verbose )")
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
    
    const uint32_t kdtree = static_cast<uint32_t>(args(0).scalar_value ());
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
    
    #pragma omp parallel for
    for(uint32_t k=0; k<Points.cols(); k++)
    {
      auto res = trees[kdtree]->Search(&Points.data()[k*Points.rows()], K);
      for(uint32_t j=0; j<K; j++)
      {
        OutDistances(j,k) = res[j].distance;
        OutIndexes(j,k) = res[j].index + 1;
      }
    }
    octave_value_list ovl;
    ovl(0) = OutIndexes;
    ovl(1) = OutDistances;
    return ovl;
}

DEFUN_DLD(NNDeleteAll, args, , "NNDeleteAll()")
{
  signal(SIGSEGV, handler);   // install our handler
  for(auto& t : trees)
  {
    delete t;
  }
  trees.clear();
  datas.clear();
  return octave_value(0);
}
