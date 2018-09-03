#ifndef __2DTree_hpp__
#define __2DTree_hpp__

#include <cstdlib>
#include "MatrixData.hpp"
// #include "./Node.hpp"
#include "utils/determineCenterAndRadius.hpp"

class FMM2DTree 
{
private:
    unsigned N;        // Total number of source particles
    unsigned n_levels; // Number of levels in the tree. This needs to be 
                       // determined depending upon the number of particles in the input

public:
    // Constructor for the class
    FMM2DTree(MatrixData &M);
    // Destructor for the class:
    // ~FMM2DTree();

    // Create Tree method:
    void createTree();
};

FMM2DTree::FMM2DTree(MatrixData &M)
{   
    this->N = M.getNumColumns(); // since this is number of sources
}

#endif
