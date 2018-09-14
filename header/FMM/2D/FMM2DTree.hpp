#ifndef __2DTree_hpp__
#define __2DTree_hpp__

#include <cstdlib>
#include "MatrixData.hpp"
#include "FMM2DBox.hpp"
#include "utils/getStandardNodes.hpp"
#include "utils/determineCenterAndRadius.hpp"
#include "MatrixFactorizations/getInterpolation.hpp"

// #include "utils/scalePoints.hpp"

class FMM2DTree 
{
private:
    MatrixData *M_ptr;       // Pointer to the structure describing the underlying problem
    unsigned N;              // Total number of source particles
    unsigned N_nodes;        // Number of nodes used for interpolation
    unsigned rank;           // Rank of the low-rank interaction
    unsigned max_levels;     // Number of levels in the tree.
                       
    const array *charges_ptr; // Holds the information for the charges of the points
    array standard_nodes_1d;  // Standard nodes alloted as given by user choice in [-1, 1]
    array standard_nodes;     // Gets the points in 2D
    std::string nodes_type;   // Whether Chebyshev, or Legendre or Equispaced

    std::array<af::array, 4> standard_nodes_child; 
    std::array<af::array, 8> neighbor_interaction;
    std::array<af::array, 4> L2L;
    std::array<af::array, 4> M2M;
    std::array<af::array, 16> M2L_inner;
    std::array<af::array, 24> M2L_outer;
    array self_interaction;
    std::vector<std::vector<FMM2DBox>> tree; // The tree storing all the information.

public:
    // Constructor for the class
    FMM2DTree(MatrixData &M, unsigned N_nodes, std::string nodes_type, const array &charges);
    // Destructor for the class:
    ~FMM2DTree(){};
  
    // Create Tree method:
    void buildTree();

    // Assigns relations amongst boxes:
    void assignChild0Relations(int N_level, int N_box);
    void assignChild1Relations(int N_level, int N_box);
    void assignChild2Relations(int N_level, int N_box);
    void assignChild3Relations(int N_level, int N_box);
    // Uses above functions to assign to all boxes in tree:
    void assignTreeRelations();

    // Gives the box details of the prescribed box and level number:
    void printBoxDetails(unsigned N_level, unsigned N_box);
    // Lists details of all boxes in the tree
    void printTreeDetails();

    // Transfer from parent to child:
    void getTransferMatrices();
    // Getting M2L operators:
    void getM2L();
    // Getting leaf level operators(neighbor and self interactions):
    void getLeafLevelInteractions();
};

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

FMM2DTree::FMM2DTree(MatrixData &M, unsigned N_nodes, std::string nodes_type, const array &charges)
{   
    this->M_ptr       = &M;
    this->N_nodes     = N_nodes;
    this->nodes_type  = nodes_type;
    this->rank        = N_nodes * N_nodes;
    this->max_levels  = 0; // Initializing to zero. Updated during tree building
    this->charges_ptr = &charges;
    
    // Finding the standard nodes(in [-1, 1]):
    getStandardNodes(N_nodes, nodes_type, this->standard_nodes_1d);

    // Array which stores the tiled version of the nodes:
    array nodes_tiled = af::tile(this->standard_nodes_1d, 1, N_nodes);

    // Getting the standard nodes in 2D:
    this->standard_nodes = af::join(1, 
                                    af::flat(nodes_tiled),
                                    af::flat(af::reorder(nodes_tiled, 1, 0))
                                   );
    this->standard_nodes.eval();

    // Getting the standard nodes for the children:
    array child_nodes;

    // We now proceed to divide the standard domain [-1, 1] to [-1, 0] and [0, 1]:
    // =================== NOTE ===========================
    // The following nomenclature is used to describe the child cell number:
    // ============================
    // ||            |           ||
    // ||            |           ||
    // ||     3      |      2    ||
    // ||            |           ||
    // ||            |           ||
    // ============================
    // ||            |           ||
    // ||            |           ||
    // ||     0      |     1     ||
    // ||            |           ||
    // ||            |           ||
    // ============================

    // For child0:
    child_nodes = 0.5 * (this->standard_nodes - 1);
    child_nodes.eval();
    this->standard_nodes_child[0] = child_nodes;

    // For child1:
    child_nodes = af::join(1,
                           0.5 * (this->standard_nodes(af::span, 0) + 1),
                           0.5 * (this->standard_nodes(af::span, 1) - 1)
                          );
    child_nodes.eval();
    this->standard_nodes_child[1] = child_nodes;

    // For child2:
    child_nodes = 0.5 * (this->standard_nodes + 1);
    child_nodes.eval();
    this->standard_nodes_child[2] = child_nodes;

    // For child3:
    child_nodes = af::join(1,
                           0.5 * (this->standard_nodes(af::span, 0) - 1),
                           0.5 * (this->standard_nodes(af::span, 1) + 1)
                          );
    child_nodes.eval();
    this->standard_nodes_child[3] = child_nodes;

    // Getting Transfer Matrices(that is M2M and L2L):
    cout << "Getting Transfer Matrices..." << endl;
    FMM2DTree::getTransferMatrices();
    cout << "Building Tree..." << endl;
    FMM2DTree::buildTree();
    cout << "Number of Levels in the Tree: " << this->max_levels << endl;
    cout << "Assigning Relations Amongst Boxes in the tree..." << endl;
    FMM2DTree::assignTreeRelations();
    cout << "Getting M2L and neighbor interactions" << endl;
}

void FMM2DTree::getTransferMatrices()
{
    array nodes_x, nodes_y, L2L_array;
    // Initializing L2L array
    L2L_array = array(N_nodes * N_nodes, N_nodes * N_nodes, f64);
    
    // For child0:
    nodes_x = (this->standard_nodes_child[0])(af::span, 0);
    nodes_y = (this->standard_nodes_child[0])(af::span, 1);

    getL2L2D(nodes_x, nodes_y, this->standard_nodes_1d, L2L_array);
    L2L_array.eval();
    this->L2L[0] = L2L_array;

    // For child1:
    nodes_x = (this->standard_nodes_child[1])(af::span, 0);
    nodes_y = (this->standard_nodes_child[1])(af::span, 1);

    getL2L2D(nodes_x, nodes_y, this->standard_nodes_1d, L2L_array);
    L2L_array.eval();
    this->L2L[1] = L2L_array;

    // For child2:
    nodes_x = (this->standard_nodes_child[2])(af::span, 0);
    nodes_y = (this->standard_nodes_child[2])(af::span, 1);

    getL2L2D(nodes_x, nodes_y, this->standard_nodes_1d, L2L_array);
    L2L_array.eval();
    this->L2L[2] = L2L_array;

    // For child3:
    nodes_x = (this->standard_nodes_child[3])(af::span, 0);
    nodes_y = (this->standard_nodes_child[3])(af::span, 1);

    getL2L2D(nodes_x, nodes_y, this->standard_nodes_1d, L2L_array);
    L2L_array.eval();
    this->L2L[3] = L2L_array;

    for(unsigned i = 0; i < 4; i++)
    {
        this->M2M[i] = this->L2L[i].T();
        this->M2M[i].eval();
    }
}

void FMM2DTree::buildTree()
{
    // Box at the root level. This box will contain all the particles
    // under consideration. This is the box that will now be recursively
    // subdivided until the smallest box only contains 4 * N_nodes^2 particles
    FMM2DBox root;

    root.N_box   =  0; // only box on it's level
    root.N_level =  0; // root is always on level 0
    root.parent  = -1; // since it doesn't have a parent

    // ===================================================================
    // Similarly since the root box doesn't have neighbors, inner or outer 
    // boxes, we hold the default value -1 given by the constructor
    // ===================================================================

    #pragma omp parallel for
    for (unsigned i = 0; i < 4; i++) 
    {
        root.children[i] = i;
    }

    // Since the root level would consist of all the boxes:
    root.inds_in_box = af::range(this->M_ptr->getNumCols(), 1, 1, 1, -1, u32);

    // Finding the radius and the centers for root:
    determineCenterAndRadius((*(this->M_ptr->getSourceCoordsPtr()))(af::span, 0), root.c_x, root.r_x);
    determineCenterAndRadius((*(this->M_ptr->getSourceCoordsPtr()))(af::span, 1), root.c_y, root.r_y);

    // Pushing back the data onto the tree variable:
    std::vector<FMM2DBox> root_level;
    root_level.push_back(root);
    tree.push_back(root_level);

    // Using this flag to check whether we have reached leaf level:
    bool reached_leaf = false;
    while(!reached_leaf)
    {
        this->max_levels++;
        std::vector<FMM2DBox> level;
        level.reserve(pow(4, this->max_levels));
        unsigned N_leaf_criterion = 0; // number of cells that satisfy leaf criterion
        for(unsigned i = 0; i < (unsigned) pow(4, this->max_levels); i++) 
        {
            FMM2DBox box;
            box.N_level = this->max_levels;
            box.N_box   = i;
            box.parent  = i / 4;

            FMM2DBox &parent_node = (tree.back())[box.parent];

            // Assigning the new centers and radii:
            box.c_x = parent_node.c_x + (((i - 4 * box.parent) % 2) - 0.5) * parent_node.r_x;
            box.c_y = parent_node.c_y + (((i - 4 * box.parent) / 2) - 0.5) * parent_node.r_y;
            box.r_x = 0.5 * parent_node.r_x;
            box.r_y = 0.5 * parent_node.r_y;

            // Locations local to the parent box:
            array locations_local_x = (*(this->M_ptr->getSourceCoordsPtr()))(parent_node.inds_in_box, 0); 
            array locations_local_y = (*(this->M_ptr->getSourceCoordsPtr()))(parent_node.inds_in_box, 1); 

            // Finding the indices in the children:
            box.inds_in_box =
            parent_node.inds_in_box(af::where(   locations_local_x < box.c_x + box.r_x
                                              && locations_local_x > box.c_x - box.r_x
                                              && locations_local_y < box.c_y + box.r_y
                                              && locations_local_y > box.c_y - box.r_y
                                             )
                                   );

            box.inds_in_box.eval();

            // We are saying that if N < 4 * N_nodes^2, then it's 
            // meets the criterion of being a leaf cell:
            if(box.inds_in_box.elements()<= 4 * this->rank)
            {
                N_leaf_criterion++;
            }

            // Assigning children for the box:
            #pragma omp parallel for
            for(unsigned j = 0; j < 4; j++) 
            {
                box.children[j] = 4 * i + j;
            }

            level.push_back(box);
        }

        tree.push_back(level);
        // If all the cells in the level have N < 4 * rank:
        if(N_leaf_criterion == pow(4, this->max_levels))
            reached_leaf = true;
    }
}

// Assigns the relations for the children all boxes in the tree
void FMM2DTree::assignTreeRelations() 
{
    // DO NOT USE OPENMP HERE!!!
    for(unsigned N_level = 0; N_level < this->max_levels; N_level++) 
    {
        // #pragma omp parallel for
        // Says Invalid control predicate??
        for(unsigned N_box = 0; N_box < pow(4, N_level); N_box++) 
        {
            assignChild0Relations(N_level, N_box);
            assignChild1Relations(N_level, N_box);
            assignChild2Relations(N_level, N_box);
            assignChild3Relations(N_level, N_box);
        }
    }
}

// For diagrams which elucidate the logic in the 
// the following functions refer to https://goo.gl/TffP2B

// Assigns the relations for child0 of a box
void FMM2DTree::assignChild0Relations(int N_level, int N_box) 
{
    // Level number for the child:
    int N_lc = N_level + 1;
    // Box number for the child:
    int N_bc = 4 * N_box;
    // Neighbor number in consideration:
    int N_neighbor;

    // Assign siblings
    tree[N_lc][N_bc].neighbor[3] = N_bc + 1;
    tree[N_lc][N_bc].neighbor[4] = N_bc + 2;
    tree[N_lc][N_bc].neighbor[5] = N_bc + 3;

    // Assign children of parent's zeroth neighbor
    N_neighbor = tree[N_level][N_box].neighbor[0];
    if(N_neighbor != -1) 
    {
        tree[N_lc][N_bc].inner[0]    = tree[N_level][N_neighbor].children[0];
        tree[N_lc][N_bc].inner[1]    = tree[N_level][N_neighbor].children[1];
        tree[N_lc][N_bc].neighbor[0] = tree[N_level][N_neighbor].children[2];
        tree[N_lc][N_bc].inner[15]   = tree[N_level][N_neighbor].children[3];             
    }

    // Assign children of parent's first neighbor
    N_neighbor = tree[N_level][N_box].neighbor[1];
    if(N_neighbor != -1) 
    {
        tree[N_lc][N_bc].inner[2]    = tree[N_level][N_neighbor].children[0];
        tree[N_lc][N_bc].inner[3]    = tree[N_level][N_neighbor].children[1];
        tree[N_lc][N_bc].neighbor[2] = tree[N_level][N_neighbor].children[2];
        tree[N_lc][N_bc].neighbor[1] = tree[N_level][N_neighbor].children[3];
    }

    // Assign children of parent's second neighbor
    N_neighbor = tree[N_level][N_box].neighbor[2];
    if(N_neighbor != -1) 
    {
        tree[N_lc][N_bc].inner[4] = tree[N_level][N_neighbor].children[0];
        tree[N_lc][N_bc].outer[7] = tree[N_level][N_neighbor].children[1];
        tree[N_lc][N_bc].outer[8] = tree[N_level][N_neighbor].children[2];
        tree[N_lc][N_bc].inner[5] = tree[N_level][N_neighbor].children[3];             
    }

    // Assign children of parent's third neighbor
    N_neighbor = tree[N_level][N_box].neighbor[3];
    if(N_neighbor != -1) 
    {
        tree[N_lc][N_bc].inner[6]  = tree[N_level][N_neighbor].children[0];
        tree[N_lc][N_bc].outer[9]  = tree[N_level][N_neighbor].children[1];
        tree[N_lc][N_bc].outer[10] = tree[N_level][N_neighbor].children[2];
        tree[N_lc][N_bc].inner[7]  = tree[N_level][N_neighbor].children[3];
    }

    // Assign children of parent's fourth neighbor
    N_neighbor = tree[N_level][N_box].neighbor[4];
    if(N_neighbor != -1) 
    {
        tree[N_lc][N_bc].inner[8]  = tree[N_level][N_neighbor].children[0];
        tree[N_lc][N_bc].outer[11] = tree[N_level][N_neighbor].children[1];
        tree[N_lc][N_bc].outer[12] = tree[N_level][N_neighbor].children[2];
        tree[N_lc][N_bc].outer[13] = tree[N_level][N_neighbor].children[3];             
    }

    //  Assign children of parent's fifth neighbor
    N_neighbor = tree[N_level][N_box].neighbor[5];
    if(N_neighbor != -1) 
    {
        tree[N_lc][N_bc].inner[10] = tree[N_level][N_neighbor].children[0];
        tree[N_lc][N_bc].inner[9]  = tree[N_level][N_neighbor].children[1];
        tree[N_lc][N_bc].outer[14] = tree[N_level][N_neighbor].children[2];
        tree[N_lc][N_bc].outer[15] = tree[N_level][N_neighbor].children[3];             
    }

    //  Assign children of parent's sixth neighbor
    N_neighbor = tree[N_level][N_box].neighbor[6];
    if(N_neighbor != -1) 
    {
        tree[N_lc][N_bc].inner[12] = tree[N_level][N_neighbor].children[0];
        tree[N_lc][N_bc].inner[11] = tree[N_level][N_neighbor].children[1];
        tree[N_lc][N_bc].outer[16] = tree[N_level][N_neighbor].children[2];
        tree[N_lc][N_bc].outer[17] = tree[N_level][N_neighbor].children[3];             
    }

    //  Assign children of parent's seventh neighbor
    N_neighbor = tree[N_level][N_box].neighbor[7];
    if(N_neighbor != -1) 
    {
        tree[N_lc][N_bc].inner[14]   = tree[N_level][N_neighbor].children[0];
        tree[N_lc][N_bc].neighbor[7] = tree[N_level][N_neighbor].children[1];
        tree[N_lc][N_bc].neighbor[6] = tree[N_level][N_neighbor].children[2];
        tree[N_lc][N_bc].inner[13]   = tree[N_level][N_neighbor].children[3];             
    }
}

// Assigns the relations for child1 of a box
void FMM2DTree::assignChild1Relations(int N_level, int N_box) 
{
    // Level number for the child:
    int N_lc = N_level + 1;
    // Box number for the child:
    int N_bc = 4 * N_box + 1;
    // Neighbor number in consideration:
    int N_neighbor;

    // Assign siblings
    tree[N_lc][N_bc].neighbor[7] = N_bc - 1;
    tree[N_lc][N_bc].neighbor[5] = N_bc + 1;
    tree[N_lc][N_bc].neighbor[6] = N_bc + 2;

    // Assign children of parent's zeroth neighbor
    N_neighbor = tree[N_level][N_box].neighbor[0];
    if(N_neighbor != -1) 
    {
        tree[N_lc][N_bc].outer[23] = tree[N_level][N_neighbor].children[0];
        tree[N_lc][N_bc].inner[0]  = tree[N_level][N_neighbor].children[1];
        tree[N_lc][N_bc].inner[15] = tree[N_level][N_neighbor].children[2];
        tree[N_lc][N_bc].outer[22] = tree[N_level][N_neighbor].children[3];             
    }

    // Assign children of parent's first neighbor
    N_neighbor = tree[N_level][N_box].neighbor[1];
    if(N_neighbor != -1) 
    {
        tree[N_lc][N_bc].inner[1]    = tree[N_level][N_neighbor].children[0];
        tree[N_lc][N_bc].inner[2]    = tree[N_level][N_neighbor].children[1];
        tree[N_lc][N_bc].neighbor[1] = tree[N_level][N_neighbor].children[2];
        tree[N_lc][N_bc].neighbor[0] = tree[N_level][N_neighbor].children[3];             
    }

    // Assign children of parent's second neighbor
    N_neighbor = tree[N_level][N_box].neighbor[2];
    if(N_neighbor != -1) 
    {
        tree[N_lc][N_bc].inner[3]    = tree[N_level][N_neighbor].children[0];
        tree[N_lc][N_bc].inner[4]    = tree[N_level][N_neighbor].children[1];
        tree[N_lc][N_bc].inner[5]    = tree[N_level][N_neighbor].children[2];
        tree[N_lc][N_bc].neighbor[2] = tree[N_level][N_neighbor].children[3];             
    }

    // Assign children of parent's third neighbor
    N_neighbor = tree[N_level][N_box].neighbor[3];
    if(N_neighbor != -1)
    {
        tree[N_lc][N_bc].neighbor[3] = tree[N_level][N_neighbor].children[0];
        tree[N_lc][N_bc].inner[6]    = tree[N_level][N_neighbor].children[1];
        tree[N_lc][N_bc].inner[7]    = tree[N_level][N_neighbor].children[2];
        tree[N_lc][N_bc].neighbor[4] = tree[N_level][N_neighbor].children[3];             
    }

    // Assign children of parent's fourth neighbor
    N_neighbor = tree[N_level][N_box].neighbor[4];
    if(N_neighbor != -1) 
    {
        tree[N_lc][N_bc].inner[9]  = tree[N_level][N_neighbor].children[0];
        tree[N_lc][N_bc].inner[8]  = tree[N_level][N_neighbor].children[1];
        tree[N_lc][N_bc].outer[13] = tree[N_level][N_neighbor].children[2];
        tree[N_lc][N_bc].outer[14] = tree[N_level][N_neighbor].children[3];             
    }

    // Assign children of parent's fifth neighbor
    N_neighbor = tree[N_level][N_box].neighbor[5];
    if(N_neighbor != -1) 
    {
        tree[N_lc][N_bc].inner[11] = tree[N_level][N_neighbor].children[0];
        tree[N_lc][N_bc].inner[10] = tree[N_level][N_neighbor].children[1];
        tree[N_lc][N_bc].outer[15] = tree[N_level][N_neighbor].children[2];
        tree[N_lc][N_bc].outer[16] = tree[N_level][N_neighbor].children[3];             
    }

    // Assign children of parent's sixth neighbor
    N_neighbor = tree[N_level][N_box].neighbor[6];
    if(N_neighbor != -1) 
    {
        tree[N_lc][N_bc].outer[19] = tree[N_level][N_neighbor].children[0];
        tree[N_lc][N_bc].inner[12] = tree[N_level][N_neighbor].children[1];
        tree[N_lc][N_bc].outer[17] = tree[N_level][N_neighbor].children[2];
        tree[N_lc][N_bc].outer[18] = tree[N_level][N_neighbor].children[3];             
    }

    // Assign children of parent's seventh neighbor
    N_neighbor = tree[N_level][N_box].neighbor[7];
    if(N_neighbor != -1) 
    {
        tree[N_lc][N_bc].outer[21] = tree[N_level][N_neighbor].children[0];
        tree[N_lc][N_bc].inner[14] = tree[N_level][N_neighbor].children[1];
        tree[N_lc][N_bc].inner[13] = tree[N_level][N_neighbor].children[2];
        tree[N_lc][N_bc].outer[20] = tree[N_level][N_neighbor].children[3];             
    }
}

// Assigns the relations for child2 of a box
void FMM2DTree::assignChild2Relations(int N_level, int N_box)
{
    // Level number for the child:
    int N_lc = N_level + 1;
    // Box number for the child:
    int N_bc = 4 * N_box + 2;
    // Neighbor number in consideration:
    int N_neighbor;

    // Assign siblings
    tree[N_lc][N_bc].neighbor[0] = N_bc - 2;
    tree[N_lc][N_bc].neighbor[1] = N_bc - 1;
    tree[N_lc][N_bc].neighbor[7] = N_bc + 1;

    // Assign children of parent's zeroth neighbor
    N_neighbor = tree[N_level][N_box].neighbor[0];
    if(N_neighbor != -1) 
    {
        tree[N_lc][N_bc].outer[0]  = tree[N_level][N_neighbor].children[0];
        tree[N_lc][N_bc].outer[1]  = tree[N_level][N_neighbor].children[1];
        tree[N_lc][N_bc].inner[0]  = tree[N_level][N_neighbor].children[2];
        tree[N_lc][N_bc].outer[23] = tree[N_level][N_neighbor].children[3];             
    }

    //  Assign children of parent's first neighbor
    N_neighbor = tree[N_level][N_box].neighbor[1];
    if(N_neighbor != -1) 
    {
        tree[N_lc][N_bc].outer[2] = tree[N_level][N_neighbor].children[0];
        tree[N_lc][N_bc].outer[3] = tree[N_level][N_neighbor].children[1];
        tree[N_lc][N_bc].inner[2] = tree[N_level][N_neighbor].children[2];
        tree[N_lc][N_bc].inner[1] = tree[N_level][N_neighbor].children[3];             
    }

    // Assign children of parent's second neighbor
    N_neighbor = tree[N_level][N_box].neighbor[2];
    if(N_neighbor != -1) 
    {
        tree[N_lc][N_bc].outer[4] = tree[N_level][N_neighbor].children[0];
        tree[N_lc][N_bc].outer[5] = tree[N_level][N_neighbor].children[1];
        tree[N_lc][N_bc].inner[4] = tree[N_level][N_neighbor].children[2];
        tree[N_lc][N_bc].inner[3] = tree[N_level][N_neighbor].children[3];             
    }

    // Assign children of parent's third neighbor
    N_neighbor = tree[N_level][N_box].neighbor[3];
    if(N_neighbor != -1) 
    {
        tree[N_lc][N_bc].neighbor[2] = tree[N_level][N_neighbor].children[0];
        tree[N_lc][N_bc].inner[5]    = tree[N_level][N_neighbor].children[1];
        tree[N_lc][N_bc].inner[6]    = tree[N_level][N_neighbor].children[2];
        tree[N_lc][N_bc].neighbor[3] = tree[N_level][N_neighbor].children[3];             
    }

    // Assign children of parent's fourth neighbor
    N_neighbor = tree[N_level][N_box].neighbor[4];
    if(N_neighbor != -1) 
    {
        tree[N_lc][N_bc].neighbor[4] = tree[N_level][N_neighbor].children[0];
        tree[N_lc][N_bc].inner[7]    = tree[N_level][N_neighbor].children[1];
        tree[N_lc][N_bc].inner[8]    = tree[N_level][N_neighbor].children[2];
        tree[N_lc][N_bc].inner[9]    = tree[N_level][N_neighbor].children[3];             
    }

    // Assign children of parent's fifth neighbor
    N_neighbor = tree[N_level][N_box].neighbor[5];
    if(N_neighbor != -1) 
    {
        tree[N_lc][N_bc].neighbor[6] = tree[N_level][N_neighbor].children[0];
        tree[N_lc][N_bc].neighbor[5] = tree[N_level][N_neighbor].children[1];
        tree[N_lc][N_bc].inner[10]   = tree[N_level][N_neighbor].children[2];
        tree[N_lc][N_bc].inner[11]   = tree[N_level][N_neighbor].children[3];             
    }

    // Assign children of parent's sixth neighbor
    N_neighbor = tree[N_level][N_box].neighbor[6];
    if(N_neighbor != -1) 
    {
        tree[N_lc][N_bc].outer[20] = tree[N_level][N_neighbor].children[0];
        tree[N_lc][N_bc].inner[13] = tree[N_level][N_neighbor].children[1];
        tree[N_lc][N_bc].inner[12] = tree[N_level][N_neighbor].children[2];
        tree[N_lc][N_bc].outer[19] = tree[N_level][N_neighbor].children[3];             
    }

    // Assign children of parent's seventh neighbor
    N_neighbor = tree[N_level][N_box].neighbor[7];
    if(N_neighbor != -1) 
    {
        tree[N_lc][N_bc].outer[22] = tree[N_level][N_neighbor].children[0];
        tree[N_lc][N_bc].inner[15] = tree[N_level][N_neighbor].children[1];
        tree[N_lc][N_bc].inner[14] = tree[N_level][N_neighbor].children[2];
        tree[N_lc][N_bc].outer[21] = tree[N_level][N_neighbor].children[3];             
    }
}

// Assigns the relations for child2 of a box
void FMM2DTree::assignChild3Relations(int N_level, int N_box)
{
    // Level number for the child:
    int N_lc = N_level + 1;
    // Box number for the child:
    int N_bc = 4 * N_box + 3;
    // Neighbor number in consideration:
    int N_neighbor;

    //  Assign siblings
    tree[N_lc][N_bc].neighbor[1] = N_bc - 3;
    tree[N_lc][N_bc].neighbor[2] = N_bc - 2;
    tree[N_lc][N_bc].neighbor[3] = N_bc - 1;

    //  Assign children of parent's zeroth neighbor
    N_neighbor = tree[N_level][N_box].neighbor[0];
    if(N_neighbor != -1) 
    {
        tree[N_lc][N_bc].outer[1] = tree[N_level][N_neighbor].children[0];
        tree[N_lc][N_bc].outer[2] = tree[N_level][N_neighbor].children[1];
        tree[N_lc][N_bc].inner[1] = tree[N_level][N_neighbor].children[2];
        tree[N_lc][N_bc].inner[0] = tree[N_level][N_neighbor].children[3];             
    }

    //  Assign children of parent's first neighbor
    N_neighbor = tree[N_level][N_box].neighbor[1];
    if(N_neighbor != -1) 
    {
        tree[N_lc][N_bc].outer[3] = tree[N_level][N_neighbor].children[0];
        tree[N_lc][N_bc].outer[4] = tree[N_level][N_neighbor].children[1];
        tree[N_lc][N_bc].inner[3] = tree[N_level][N_neighbor].children[2];
        tree[N_lc][N_bc].inner[2] = tree[N_level][N_neighbor].children[3];             
    }

    //  Assign children of parent's second neighbor
    N_neighbor = tree[N_level][N_box].neighbor[2];
    if(N_neighbor != -1) 
    {
        tree[N_lc][N_bc].outer[5] = tree[N_level][N_neighbor].children[0];
        tree[N_lc][N_bc].outer[6] = tree[N_level][N_neighbor].children[1];
        tree[N_lc][N_bc].outer[7] = tree[N_level][N_neighbor].children[2];
        tree[N_lc][N_bc].inner[4] = tree[N_level][N_neighbor].children[3];             
    }

    //  Assign children of parent's third neighbor
    N_neighbor = tree[N_level][N_box].neighbor[3];
    if(N_neighbor != -1) 
    {
        tree[N_lc][N_bc].inner[5] = tree[N_level][N_neighbor].children[0];
        tree[N_lc][N_bc].outer[8] = tree[N_level][N_neighbor].children[1];
        tree[N_lc][N_bc].outer[9] = tree[N_level][N_neighbor].children[2];
        tree[N_lc][N_bc].inner[6] = tree[N_level][N_neighbor].children[3];             
    }

    //  Assign children of parent's fourth neighbor
    N_neighbor = tree[N_level][N_box].neighbor[4];
    if(N_neighbor != -1) 
    {
        tree[N_lc][N_bc].inner[7]  = tree[N_level][N_neighbor].children[0];
        tree[N_lc][N_bc].outer[10] = tree[N_level][N_neighbor].children[1];
        tree[N_lc][N_bc].outer[11] = tree[N_level][N_neighbor].children[2];
        tree[N_lc][N_bc].inner[8]  = tree[N_level][N_neighbor].children[3];             
    }

    //  Assign children of parent's fifth neighbor
    N_neighbor = tree[N_level][N_box].neighbor[5];
    if(N_neighbor != -1) 
    {
        tree[N_lc][N_bc].neighbor[5] = tree[N_level][N_neighbor].children[0];
        tree[N_lc][N_bc].neighbor[4] = tree[N_level][N_neighbor].children[1];
        tree[N_lc][N_bc].inner[9]    = tree[N_level][N_neighbor].children[2];
        tree[N_lc][N_bc].inner[10]   = tree[N_level][N_neighbor].children[3];             
    }

    //  Assign children of parent's sixth neighbor
    N_neighbor = tree[N_level][N_box].neighbor[6];
    if(N_neighbor != -1) 
    {
        tree[N_lc][N_bc].inner[13]   = tree[N_level][N_neighbor].children[0];
        tree[N_lc][N_bc].neighbor[6] = tree[N_level][N_neighbor].children[1];
        tree[N_lc][N_bc].inner[11]   = tree[N_level][N_neighbor].children[2];
        tree[N_lc][N_bc].inner[12]   = tree[N_level][N_neighbor].children[3];             
    }

    //  Assign children of parent's seventh neighbor
    N_neighbor = tree[N_level][N_box].neighbor[7];
    if(N_neighbor != -1) 
    {
        tree[N_lc][N_bc].inner[15]   = tree[N_level][N_neighbor].children[0];
        tree[N_lc][N_bc].neighbor[0] = tree[N_level][N_neighbor].children[1];
        tree[N_lc][N_bc].neighbor[7] = tree[N_level][N_neighbor].children[2];
        tree[N_lc][N_bc].inner[14]   = tree[N_level][N_neighbor].children[3];             
    }
}

void FMM2DTree::getM2LInteractions() 
{
    array M2L;

    // Getting outer interactions:
    for(unsigned i = 0; i < 6; i++)
    {
        nodes= af::join(1, this->standard_nodes(af::span, 0) + 2 * (i - 3),
                        this->standard_nodes(af::span, 1) - 6
                       )

        getM2L(this->standard_nodes, nodes, *(this->M_ptr), M2L);
        M2L.eval();
        this->M2L_outer[i] = M2L;
    }

    for(unsigned i = 0; i < 6; i++)
    {
        nodes= af::join(1, this->standard_nodes(af::span, 0) + 6,
                        this->standard_nodes(af::span, 1) + 2 * (i - 3)
                       )

        getM2L(this->standard_nodes, nodes, *(this->M_ptr), M2L);
        M2L.eval();
        this->M2L_outer[i + 6] = M2L;
    }

    for(unsigned i = 0; i < 6; i++)
    {
        nodes= af::join(1, this->standard_nodes(af::span, 0) + 2 * (3 - i),
                        this->standard_nodes(af::span, 1) + 6
                       )

        getM2L(this->standard_nodes, nodes, *(this->M_ptr), M2L);
        M2L.eval();
        this->M2L_outer[i + 12] = M2L;
    }

    for(unsigned i = 0; i < 6; i++)
    {
        nodes= af::join(1, this->standard_nodes(af::span, 0) - 6,
                        this->standard_nodes(af::span, 1) + 2 * (3 - i)
                       )

        getM2L(this->standard_nodes, nodes, *(this->M_ptr), M2L);
        M2L.eval();
        this->M2L_outer[i + 18] = M2L;
    }

    // Getting inner interactions:
    for(unsigned i = 0; i < 4; i++)
    {
        nodes= af::join(1, this->standard_nodes(af::span, 0) + 2 * (i - 2),
                        this->standard_nodes(af::span, 1) - 4
                       )

        getM2L(this->standard_nodes, nodes, *(this->M_ptr), M2L);
        M2L.eval();
        this->M2L_inner[i] = M2L;
    }

    for(unsigned i = 0; i < 4; i++)
    {
        nodes= af::join(1, this->standard_nodes(af::span, 0) + 4,
                        this->standard_nodes(af::span, 1) + 2 * (i - 2)
                       )

        getM2L(this->standard_nodes, nodes, *(this->M_ptr), M2L);
        M2L.eval();
        this->M2L_inner[i + 4] = M2L;
    }

    for(unsigned i = 0; i < 4; i++)
    {
        nodes= af::join(1, this->standard_nodes(af::span, 0) + 2 * (2 - i),
                        this->standard_nodes(af::span, 1) + 4
                       )

        getM2L(this->standard_nodes, nodes, *(this->M_ptr), M2L);
        M2L.eval();
        this->M2L_inner[i + 8] = M2L;
    }

    for(unsigned i = 0; i < 4; i++)
    {
        nodes= af::join(1, this->standard_nodes(af::span, 0) - 4,
                        this->standard_nodes(af::span, 1) + 2 * (2 - i)
                       )

        getM2L(this->standard_nodes, nodes, *(this->M_ptr), M2L);
        M2L.eval();
        this->M2L_inner[i + 12] = M2L;
    }
}

void FMM2DTree::getLeafLevelInteractions()
{
    array M2L;
    // Get neighbor interactions:
    r_x_leaf = this->tree[this->max_levels][0].r_x; // since its a uniform tree, 
                                                    // we don't bother about the box considered
    r_y_leaf = this->tree[this->max_levels][0].r_y;

    leaf_nodes = af::join(1, this->standard_nodes(af::span, 0) * r_x_leaf,
                          this->standard_nodes(af::span, 1) * r_y_leaf
                         );

    for(unsigned i = 0; i < 2; i++)
    {
        nodes= af::join(1, (this->standard_nodes(af::span, 0) + 2 * (i - 1)) * r_x_leaf,
                        (this->standard_nodes(af::span, 1) - 2) * r_y_leaf
                       );

        getM2L(leaf_nodes, nodes, *(this->M_ptr), M2L);
        M2L.eval();
        this->neighbor_interaction[i] = M2L;
    }

    for(unsigned i = 0; i < 2; i++)
    {
        nodes= af::join(1, (this->standard_nodes(af::span, 0) + 2) * r_x_leaf,
                        (this->standard_nodes(af::span, 1) + 2 * (i - 1)) * r_y_leaf
                       );
        getM2L(leaf_nodes, nodes, *(this->M_ptr), M2L);
        M2L.eval();
        this->neighbor_interaction[i + 2] = M2L;
    }


    for(unsigned i = 0; i < 2; i++)
    {
        nodes= af::join(1, (this->standard_nodes(af::span, 0) + 2 * (1 - i)) * r_x_leaf,
                        (this->standard_nodes(af::span, 1) + 2) * r_y_leaf
                       )

        getM2L(this->standard_nodes, nodes, *(this->M_ptr), M2L);
        M2L.eval();
        this->neighbor_interaction[i + 4] = M2L;
    }

    for(unsigned i = 0; i < 2; i++)
    {
        nodes= af::join(1, (this->standard_nodes(af::span, 0) - 2) * r_x_leaf,
                        (this->standard_nodes(af::span, 1) + 2 * (1 - i)) * r_y_leaf
                       )

        getM2L(this->standard_nodes, nodes, *(this->M_ptr), M2L);
        M2L.eval();
        this->neighbor_interaction[i + 6] = M2L;
    }

    // Getting self-interaction:
    getM2L(leaf_nodes, leaf_nodes, *(this->M_ptr), this->self_interaction);
    this->self_interaction.eval();
}

#endif
