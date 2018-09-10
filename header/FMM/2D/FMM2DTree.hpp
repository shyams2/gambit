#ifndef __2DTree_hpp__
#define __2DTree_hpp__

#include <cstdlib>
#include "MatrixData.hpp"
#include "utils/determineCenterAndRadius.hpp"
#include "utils/getStandardNodes.hpp"
#include "utils/scalePoints.hpp"
#include "MatrixFactorizations/getInterpolation.hpp"
#include "FMM2DBox.hpp"

class FMM2DTree 
{
private:
    unsigned N;          // Total number of source particles
    unsigned N_nodes;    // Number of nodes used for interpolation
    unsigned rank;       // Rank of the low-rank interaction
    unsigned max_levels; // Number of levels in the tree.
                       
    const array charges;     // Holds the information for the charges of the points
    array standard_nodes_1d; // Standard nodes alloted as given by user choice

    std::string nodes_type; // Whether Chebyshev, or Legendre or Equispaced

public:
    // Constructor for the class
    FMM2DTree(MatrixData &M);
    // Destructor for the class:
    ~FMM2DTree(){};
    // Transfer from parent to child:
    void getTransferParentToChild();

    // Create Tree method:
    void createTree();
};

FMM2DTree::FMM2DTree(MatrixData &M, unsigned N_nodes, std::string nodes_type, const array &charges)
{   
    this->N          = M.getNumColumns(); // since this is number of sources
    this->N_nodes    = N_nodes;
    this->nodes_type = nodes_type;
    this->rank       = N_nodes * N_nodes
    this->N_levels   = 0; // Initializing to zero
    this->charges    = charges;
    
    // Finding the standard nodes(in [-1, 1]):
    getStandardNodes(N_nodes, nodes_type, this->standard_nodes_1d);

    // Creating the root-level box:
    FMM2DBox *root = new FMM2DBox;

    root->is_root     = true;
    root->N_level     = 0;       // since this is the root level
    root->N           = this->N; // since number of particles would be the same
    root->inds_in_box = af::range(this->N);

    determineCenterAndRadius(M.getSourceCoords(af::span, 0), root->c_x, root->r_x);
    determineCenterAndRadius(M.getSourceCoords(af::span, 1), root->c_y, root->r_y);
}
