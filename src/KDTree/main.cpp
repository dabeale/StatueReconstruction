#include "tree.h"
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <vector>
static std::vector<KD::Tree> kdv;
double* point;
double* data;

void CreateData()
{
    const uint32_t NumberOfPoints = 30;
    const uint32_t dimension = 2;

    data = new double[NumberOfPoints*dimension];
    for(uint32_t i=0; i< NumberOfPoints*dimension; i++)
    {
        data[i]=static_cast<double>(rand()%10000)/10000;
    }

    int n = 20;

    point = new double[dimension];
    point[0] = data[n*dimension] + 0.001;
    point[1] = data[n*dimension+1] + 0.001;

    for(uint32_t k=0; k<NumberOfPoints; k++)
    {
        std::cout << std::sqrt((point[0] - data[k*dimension + 0])*(point[0] - data[k*dimension + 0]) +
                     (point[1] - data[k*dimension + 1])*(point[1] - data[k*dimension + 1])) << ",";
    }
    std::cout << std::endl;

    std::cout << "Point:" << std::endl;
    KD::print_array(point, dimension, 1);

    std::cout << "Vertices: " << std::endl;
    KD::print_array( data, dimension, NumberOfPoints );

    KD::Tree t(data, dimension, NumberOfPoints, false, 1);

    kdv.push_back(t);
}

int main(int argc, char *argv[])
{
    (void) argc;
    (void) argv;
    CreateData();
    const KD::Tree& kdp = kdv[0];
    auto results = kdp.Search(&point[0], 20);
    //test
    for(uint32_t i=0;i<30; i++)
    {
      kdp.Search(&data[i*2], 20);
    }

    uint32_t i=0;
    for(auto tr : results)
    {
        std::cout << tr.distance << "," << tr.index << std::endl;
        if((i++)==3)
        {
            break;
        }
    }

    KD::RandomisedTrees rdt(data, 2, 30, 5);
    auto resultsrandom  = rdt.Search(&point[0], 3);
    std::cout << std::endl << "Random Results" << std::endl;
    i=0;
    for(auto tr : resultsrandom)
    {
        std::cout << tr.distance << "," << tr.index << std::endl;
        if((i++)==3)
        {
            break;
        }
    }

    kdv.clear();
    delete point;
    delete data;
}
