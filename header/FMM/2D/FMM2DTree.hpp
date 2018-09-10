#ifndef __2DTree_hpp__
#define __2DTree_hpp__

#include <cstdlib>
#include "MatrixData.hpp"
#include "utils/getStandardNodes.hpp"
#include "utils/scalePoints.hpp"
#include "MatrixFactorizations/getInterpolation.hpp"

class FMM2DTree 
{
private:
    MatrixData *M_ptr;   // Pointer to the structure describing the underlying problem
    unsigned N;          // Total number of source particles
    unsigned N_nodes;    // Number of nodes used for interpolation
    unsigned rank;       // Rank of the low-rank interaction
    unsigned max_levels; // Number of levels in the tree.
                       
    const array charges;     // Holds the information for the charges of the points
    array standard_nodes_1d; // Standard nodes alloted as given by user choice
    std::string nodes_type;  // Whether Chebyshev, or Legendre or Equispaced

    std::vector<std::vector<FMM2DBox>> tree; // The tree storing all the information.

public:
    // Constructor for the class
    FMM2DTree(const MatrixData &M);
    // Destructor for the class:
    ~FMM2DTree(){};
    
    // Transfer from parent to child:
    void getTransferParentToChild();

    // Create Tree method:
    void createTree();
};

FMM2DTree::FMM2DTree(MatrixData &M, unsigned N_nodes, std::string nodes_type, const array &charges)
{   
    this->M_ptr      = &M;
    this->N          = M.getNumColumns(); // since this is number of sources
    this->N_nodes    = N_nodes;
    this->nodes_type = nodes_type;
    this->rank       = N_nodes * N_nodes
    this->N_levels   = 0; // Initializing to zero
    this->charges    = charges;
    
    // Finding the standard nodes(in [-1, 1]):
    getStandardNodes(N_nodes, nodes_type, this->standard_nodes_1d);
}
