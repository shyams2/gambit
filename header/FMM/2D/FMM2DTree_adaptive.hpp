#ifndef __2DTree_hpp__
#define __2DTree_hpp__

#include <cstdlib>
#include "MatrixData.hpp"
#include "FMM2DBox.hpp"
#include "utils/getStandardNodes.hpp"
#include "utils/determineCenterAndRadius.hpp"
#include "utils/scalePoints.hpp"
#include "MatrixFactorizations/getInterpolation.hpp"

using af::matmul;

class FMM2DTree 
{
private:
    
    MatrixData *M_ptr;       // Pointer to the structure describing the underlying problem
    bool user_def_locations; // User defined locations for the points in the domain
    unsigned N;              // Total number of source particles
    unsigned N_nodes;        // Number of nodes used for interpolation
    unsigned rank;           // Rank of the low-rank interaction
    unsigned max_levels;     // Number of levels in the tree.

    array standard_nodes_1d;  // Standard nodes alloted as given by user choice in [-1, 1]
    array standard_nodes;     // Gets the points in 2D
    std::string nodes_type;   // Whether Chebyshev, or Legendre or Equispaced
    const array *charges_ptr; // Pointer to the array that contains the charges
    array potentials;         // Evaluated potentials for all the particles

    double c_x, c_y, r_x, r_y; // Radii and the center coordinates for the domain
    // Flags which determine structure for the underlying kernel:
    bool is_homogeneous, is_log_homogeneous, is_translation_invariant;

    std::vector<size_t> number_of_boxes; // Number of boxes at each level
    // The following functions will only be used if the function is homog / loghomog
    double alpha; // degree of homog

    std::array<af::array, 4> standard_nodes_child;    // Nodes of the child given that parent has standard nodes in [-1, 1]
    std::array<af::array, 4> L2L;                     // L2L operators from the parent to the four children
    std::array<af::array, 4> M2M;                     // M2M operators from the children to the parent
    std::vector<std::array<af::array, 16>> M2L_inner; // M2L operators from the inner well separated clusters
    std::vector<std::array<af::array, 24>> M2L_outer; // M2L operators from the outer well separated clusters
    std::vector<std::vector<FMM2DBox>> tree;          // The tree storing all the information.

    // ==================================================================
    // ========== FUNCTIONS USED IN EVALUATION OF POTENTIAL =============
    // Step 1 of algo mentioned in Page 5 of Fong. et al.(P2M):
    // Gets the charges at the leaf level:
    void assignLeafCharges();

    // Step 2 of algo mentioned in Page 5 of Fong. et al.(M2M):
    // Computes the multipoles for all the boxes using upward sweep:
    void upwardTraveral();

    // Step 3 of algo mentioned in Page 5 of Fong. et al.(M2L):
    // Computing M2Ls for all the boxes:
    // Used when kernel is homogeneous / log-homogeneous
    void evaluateAllM2LHomogeneous();
    // Otherwise:
    void evaluateAllM2L();

    // Step 4 of algo mentioned in Page 5 of Fong. et al.(L2L):
    // Performs the L2L using a downwardsweep:
    void downwardTraversal();

    // Step 5 of algo mentioned in Page 5 of Fong. et al.(L2L):
    // Getting the potentials at the particles using node potentials:
    void evaluateL2P();
    // Performs the leaf level interactions(i.e. neighbor and self interactions):
    void evaluateLeafLevelInteractions();

public:
    // Constructor for the class
    FMM2DTree(MatrixData &M, unsigned N_nodes, std::string nodes_type, unsigned N_levels, double radius);
    // Destructor for the class:
    ~FMM2DTree(){};

    // Gives the box details of the prescribed box and level number:
    void printBoxDetails(unsigned N_level, unsigned N_box);
    // Lists details of all boxes in the tree
    void printTreeDetails();
    // Returns the potential:
    array& getPotential(const array& charges);
    // Checks the potential in FMM(used in validation)
    void checkPotentialInBox(int N_box);
};

// =================== PUBLIC FUNCTIONS =============================
void FMM2DTree::printBoxDetails(unsigned N_level, unsigned N_box)
{
    tree[N_level][N_box].printBoxDetails();
}

void FMM2DTree::printTreeDetails()
{
    for(unsigned N_level = 0; N_level < this->max_levels; N_level++)
    {
        for(unsigned N_box = 0; N_box < pow(4, N_level); N_box++)
        {
            FMM2DTree::printBoxDetails(N_level, N_box);
            cout << "==========================================================================================================================================" << endl;
        }
        cout << endl << endl;
    }
}
