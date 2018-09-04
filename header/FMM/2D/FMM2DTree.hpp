#ifndef __2DTree_hpp__
#define __2DTree_hpp__

#include <cstdlib>
#include "MatrixData.hpp"
#include "utils/determineCenterAndRadius.hpp"
#include "utils/getStandardNodes.hpp"
#include "MatrixFactorizations/getInterpolation.hpp"
#include "FMM2DBox.hpp"

class FMM2DTree 
{
private:
    unsigned N;        // Total number of source particles
    unsigned N_nodes;  // Number of nodes used for interpolation
    unsigned rank;     // Rank of the low-rank interaction
    unsigned N_levels; // Number of levels in the tree. This needs to be 
                       // determined depending upon the number of particles in the input

    const array *charges; // holds the information for the charges of the points
    const array *locations; // holds the information for the location of points
    array standard_nodes; // Standard nodes alloted as given by user choice

    std::string nodes_type; // Whether Chebyshev, or Legendre or Equispaced
    std::vector<std::vector<FMM2DBox>> tree; // used in storing all the boxes in a hierarchical fashion

public:
    // Constructor for the class
    FMM2DTree(MatrixData &M);
    // Destructor for the class:
    ~FMM2DTree(){};
    // Transfer from parent to child:
    void getTransferParentToChild();

    // Create Tree method:
    void createTree();
    // Assign children:
    void assignChildren()
};

FMM2DTree::FMM2DTree(MatrixData &M, unsigned N_nodes, std::string nodes_type, const array &charges)
{   
    this->N         = M.getNumColumns(); // since this is number of sources
    this->N_nodes   = N_nodes;
    this->rank      = N_nodes * N_nodes
    this->N_levels  = 0; 
    this->charges   = &charges;
    this->locations = M.getSourceCoordsPtr();
    
    getStandardNodes(N_nodes, nodes_type, this->standard_nodes);

}

// The following nomenclature is used to describe the child cell number:
// ============================
// ||            |           ||
// ||            |           ||
// ||     4      |      3    ||
// ||            |           ||
// ||            |           ||
// ============================
// ||            |           ||
// ||            |           ||
// ||     1      |     2     ||
// ||            |           ||
// ||            |           ||
// ============================
FMM2DTree::getTransferParentToChild(int N_child, array &L2L)
{
    // Dividing the domain [-1, 1] to [-1, 0] and [0, 1]:
    // Child Nodes which are < 0:
    array child_nodes_1 = 0.5 * (this->standard_nodes - 1);
    // Child Nodes which are > 0:
    array child_nodes_2 = 0.5 * (this->standard_nodes + 1);

    // Getting nodes_x and nodes_y depending upon the quadrant:
    array child_nodes_x, child_nodes_y;

    if(N_child == 1)
    {
        child_nodes_x = af::flat(af::tile(child_nodes_1, 1, this->N_nodes));
        child_nodes_y = af::flat(af::tile(child_nodes_1.T(), this->N_nodes));
    }

    else if(N_child == 2)
    {
        child_nodes_x = af::flat(af::tile(child_nodes_2, 1, this->N_nodes));
        child_nodes_y = af::flat(af::tile(child_nodes_1.T(), this->N_nodes));
    }

    else if(N_child == 3)
    {
        child_nodes_x = af::flat(af::tile(child_nodes_2, 1, this->N_nodes));
        child_nodes_y = af::flat(af::tile(child_nodes_2.T(), this->N_nodes));
    }

    else if(N_child == 4)
    {
        child_nodes_x = af::flat(af::tile(child_nodes_1, 1, this->N_nodes));
        child_nodes_y = af::flat(af::tile(child_nodes_2.T(), this->N_nodes));
    }

    else
    {
        cout << "INVALID!!!" << endl;
        exit(1);
    }

    getL2L2D(child_nodes_x, child_nodes_y, this->standard_nodes, L2L);
}

void assignChildren()
{
    if(Root.N == 0)
    {
        node.isLeaf = True
        node.isEmpty = True
    }
}