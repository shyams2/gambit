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
    bool is_translation_invariant;

    std::array<af::array, 4> L2L;            // L2L operators from the parent to the four children
    std::array<af::array, 4> M2M;            // M2M operators from the children to the parent
    std::vector<std::vector<FMM2DBox>> tree; // The tree storing all the information.

    // ===== FUNCTIONS USED IN TREE BUILDING AND PRECOMPUTING OPERATIONS =====
    // Create Tree method:
    void buildTree();
    // Assigns relations amongst boxes:
    void assignChild0Relations(int N_level, int N_box);
    void assignChild1Relations(int N_level, int N_box);
    void assignChild2Relations(int N_level, int N_box);
    void assignChild3Relations(int N_level, int N_box);
    // Uses above functions to assign to all boxes in tree:
    void assignTreeRelations();
    // Transfer from parent to child:
    void getTransferMatrices();

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

// So when N_levels and radius is passed to the constructor, then we'll create
// the tree assuming that at the leaf level, the particles are at the leaf nodes
// This would mean that there is no need for the P2M and L2P operators:
FMM2DTree::FMM2DTree(MatrixData &M, unsigned N_nodes, std::string nodes_type,
                     unsigned N_levels = 0, double radius = -1
                    )
{   
    this->M_ptr      = &M;
    this->N_nodes    = N_nodes;
    this->nodes_type = nodes_type;
    this->rank       = N_nodes * N_nodes;
    this->potentials = af::constant(0, this->M_ptr->getNumCols(), f64);
    this->max_levels = N_levels;
    
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

    // Finding out the radius and center for the domain of consideration:
    determineCenterAndRadius((*(this->M_ptr->getSourceCoordsPtr()))(af::span, 0), this->c_x, this->r_x);
    determineCenterAndRadius((*(this->M_ptr->getSourceCoordsPtr()))(af::span, 1), this->c_y, this->r_y);

    // If we use the above obtained radii, we will get some particles on the edge of the domain
    // To account for these edge-cases, we will slightly increase the radius of the domain:
    this->r_x = (1 + 1e-14) * this->r_x;
    this->r_y = (1 + 1e-14) * this->r_y;

    // Ensuring that the problem is translation invariant:
    this->is_translation_invariant = this->M_ptr->isTranslationInvariant(2);
    if(this->is_translation_invariant == false)
    {
        cout << "Non translation-invariant Kernels not supported!!!" << endl;
        exit(1);
    }

    // Getting Transfer Matrices(that is M2M and L2L):
    // cout << "Getting Transfer Matrices..." << endl;
    FMM2DTree::getTransferMatrices();
    // cout << "Building Tree..." << endl;
    FMM2DTree::buildTree();

    // cout << "Assigning Relations Amongst Boxes in the tree..." << endl;
    FMM2DTree::assignTreeRelations();
}

void FMM2DTree::getTransferMatrices()
{
    // Getting the standard nodes for the children:
    array child_nodes;
    std::array<array> standard_nodes_child(4);

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
    standard_nodes_child[0] = child_nodes;

    // For child1:
    child_nodes = af::join(1,
                           0.5 * (this->standard_nodes(af::span, 0) + 1),
                           0.5 * (this->standard_nodes(af::span, 1) - 1)
                          );
    child_nodes.eval();
    standard_nodes_child[1] = child_nodes;

    // For child2:
    child_nodes = 0.5 * (this->standard_nodes + 1);
    child_nodes.eval();
    standard_nodes_child[2] = child_nodes;

    // For child3:
    child_nodes = af::join(1,
                           0.5 * (this->standard_nodes(af::span, 0) - 1),
                           0.5 * (this->standard_nodes(af::span, 1) + 1)
                          );
    child_nodes.eval();
    standard_nodes_child[3] = child_nodes;

    array nodes_x, nodes_y, L2L_array;
    // Initializing L2L array
    L2L_array = array(N_nodes * N_nodes, N_nodes * N_nodes, f64);

    // L2L is used to tranfer the information from the parent to the child
    for(unsigned short N_child = 0; N_child < 4; N_child++)
    {
        nodes_x = (standard_nodes_child[N_child])(af::span, 0);
        nodes_y = (standard_nodes_child[N_child])(af::span, 1);

        getL2L2D(nodes_x, nodes_y, this->standard_nodes_1d, L2L_array);
        L2L_array.eval();
        this->L2L[N_child] = L2L_array;
    }

    // M2M is used to transfer the information from the children to the parent:
    // Hence, M2M is going to be the transpose of the L2L operators:
    for(unsigned short N_child = 0; N_child < 4; N_child++)
    {
        this->M2M[N_child] = this->L2L[N_child].T();
        this->M2M[N_child].eval();
    }
}

void FMM2DTree::buildTree()
{
    // Box at the root level. This box will contain all the particles
    // under consideration. This is the box that will now be recursively
    // subdivided until the smallest box only contains 4 * N_nodes^2 particles
    FMM2DBox root;

    root.N_box       =  0; // only box on it's level
    root.N_level     =  0; // root is always on level 0
    root.parent      = -1; // since it doesn't have a parent
    root.is_assigned = true

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

    // Assigning the radius and the centers for root:
    root.c_x = this->c_x;
    root.c_y = this->c_y;
    root.r_x = this->r_x;
    root.r_y = this->r_y;

    // Pushing back the data onto the tree variable:
    std::vector<FMM2DBox> root_level;
    root_level.push_back(root);
    tree.push_back(root_level);

    // Pushing for the root level:
    this->number_of_boxes.push_back(1);
    
    // Using this flag to check whether we have reached leaf level:
    bool reached_leaf = false;
    unsigned N_level = 0;

    while(!reached_leaf)
    {
        N_level++;
        this->number_of_boxes.push_back(4 * this->number_of_boxes[N_level - 1]);
        std::vector<FMM2DBox> level;
        level.reserve(this->number_of_boxes[N_level]);

        for(unsigned i = 0; i < this->number_of_boxes[N_level]; i++) 
        {
            FMM2DBox box;
            box.N_level = N_level;
            box.N_box   = i;
            box.parent  = i / 4;
            // Initializing the value for the potentials:
            box.node_potentials = af::constant(0, this->rank, f64);

            FMM2DBox &parent_node = (tree.back())[box.parent];
            // Assigning the new centers and radii:
            box.c_x = parent_node.c_x + (((i - 4 * box.parent) % 2) - 0.5) * parent_node.r_x;
            box.c_y = parent_node.c_y + (((i - 4 * box.parent) / 2) - 0.5) * parent_node.r_y;

            if(i % 4 > 1)
                box.c_x = parent_node.c_x - (((i - 4 * box.parent) % 2) - 0.5) * parent_node.r_x;

            else
                box.c_x = parent_node.c_x + (((i - 4 * box.parent) % 2) - 0.5) * parent_node.r_x;

            box.r_x = 0.5 * parent_node.r_x;
            box.r_y = 0.5 * parent_node.r_y;

            if(parent_node.is_leaf == true)
                box.is_assigned = false;

            // When the charges at leaf level are at non-standard locations:
            if(this->max_levels == 0)
            {
                if(box.is_assigned = true)
                {
                    // Locations local to the parent box:
                    array locations_local_x = (*(this->M_ptr->getSourceCoordsPtr()))(parent_node.inds_in_box, 0); 
                    array locations_local_y = (*(this->M_ptr->getSourceCoordsPtr()))(parent_node.inds_in_box, 1); 

                    // Finding the indices in the children:
                    box.inds_in_box =
                    parent_node.inds_in_box(af::where(   locations_local_x <  box.c_x + box.r_x
                                                      && locations_local_x >= box.c_x - box.r_x
                                                      && locations_local_y <  box.c_y + box.r_y
                                                      && locations_local_y >= box.c_y - box.r_y
                                                     )
                                           );

                    box.inds_in_box.eval();

                    // If the box shows that it contains < 4 * rank particles, that's a leaf box:
                    if(box.inds_in_box.elements()<= 4 * this->rank)
                    {
                        box.is_leaf = true;
                    }
                }
            }

            // Assigning children for the box:
            #pragma omp parallel for
            for(unsigned j = 0; j < 4; j++) 
            {
                box.children[j] = 4 * i + j;
            }
            level.push_back(box);
        }

        // Terminating if all the boxes have reached the leaf criterion:
        reached_leaf = true
        for(unsigned i = 0; i < this->number_of_boxes[N_level]; i++) 
        {
            FMM2DBox &box = level[i];
            if(box.is_leaf == false)
                reached_leaf = false;
        }
        
        tree.push_back(level);
    }

    if(this->max_levels == 0)
        this->max_levels = N_level;
}

// Assigns the relations for the children all boxes in the tree
void FMM2DTree::assignTreeRelations() 
{
    // DO NOT USE OPENMP HERE!!!
    for(unsigned N_level = 0; N_level < this->max_levels; N_level++) 
    {
        #pragma omp parallel for
        for(int N_box = 0; N_box < this->number_of_boxes[N_level]; N_box++) 
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

void FMM2DTree::assignLeafCharges()
{
    for(unsigned N_level = 0; N_level < this->max_levels; N_level++)
    {
        for(unsigned i = 0; i < this->number_of_boxes[N_level]; i++)
        {   
            // Getting the box that is in consideration:
            FMM2DBox &box = this->tree[this->max_levels][i];

            if(box.is_leaf == true)
            {
                array std_locations_x, std_locations_y;
                // Mapping onto the standard interval of [-1, 1]:
                scalePoints(box.c_x, box.r_x, (*(this->M_ptr->getSourceCoordsPtr()))(box.inds_in_box, 0),
                            0, 1, std_locations_x
                           );

                scalePoints(box.c_y, box.r_y, (*(this->M_ptr->getSourceCoordsPtr()))(box.inds_in_box, 1),
                            0, 1, std_locations_y
                           );

                // Getting the P2M operator which interpolates information from the particles to nodes:
                array P2M_array = array(box.inds_in_box.elements(), N_nodes * N_nodes, f64);
                getL2L2D(std_locations_x, std_locations_y, this->standard_nodes_1d, P2M_array);
                box.node_charges = af::matmul(P2M_array.T(), 
                                              (*this->charges_ptr)(box.inds_in_box)
                                             );
                box.node_charges.eval();
            }
        }
    }
}

void FMM2DTree::upwardTraveral()
{
    // Starting at the level just above the leaf level:
    for (int N_level = this->max_levels - 1; N_level > 1; N_level--) 
    {
        // Level number for the child:
        int N_lc = N_level + 1;

        int N_bc0, N_bc1, N_bc2, N_bc3;
        for(unsigned N_box = 0; N_box < this->number_of_boxes[N_level]; N_box++) 
        {
            if(tree[N_level][N_box].is_assigned = true)
            {
                // Box number for child0:
                N_bc0 = 4 * N_box;
                // Box number for child1:
                N_bc1 = 4 * N_box + 1;
                // Box number for child2:
                N_bc2 = 4 * N_box + 2;
                // Box number for child3:
                N_bc3 = 4 * N_box + 3;

                tree[N_level][N_box].node_charges =  af::matmul(this->M2M[0], tree[N_lc][N_bc0].node_charges)
                                                   + af::matmul(this->M2M[1], tree[N_lc][N_bc1].node_charges)
                                                   + af::matmul(this->M2M[2], tree[N_lc][N_bc2].node_charges)
                                                   + af::matmul(this->M2M[3], tree[N_lc][N_bc3].node_charges);
            }
        }
    }
}
