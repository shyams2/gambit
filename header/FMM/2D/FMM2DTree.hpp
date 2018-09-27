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
    std::vector<double> box_radius;
    std::vector<double> box_log_homog_radius;
    std::vector<double> box_homog_radius;
    double alpha; // degree of homog

    std::array<af::array, 4> standard_nodes_child;    // Nodes of the child given that parent has standard nodes in [-1, 1]
    std::array<af::array, 4> L2L;                     // L2L operators from the parent to the four children
    std::array<af::array, 4> M2M;                     // M2M operators from the children to the parent
    std::vector<std::array<af::array, 16>> M2L_inner; // M2L operators from the inner well separated clusters
    std::vector<std::array<af::array, 24>> M2L_outer; // M2L operators from the outer well separated clusters
    std::vector<std::vector<FMM2DBox>> tree;          // The tree storing all the information.

    // The following operators are used when the charges at 
    // the leaf level are placed at the standard nodes:
    // NOTE: These operations only need to be performed at the leaf level:
    std::array<af::array, 8> neighbor_interaction; // Operators used for neighbor interaction
    array self_interaction;                        // Self-interaction operator 

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
    // Getting M2L operators:
    void getM2LInteractions();
    // Getting M2L operators when kernel is homogeneous:
    void getM2LInteractionsHomogeneous();
    // Getting the neighbour and self interactions:
    void getNeighborSelfInteractions();

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

array& FMM2DTree::getPotential(const array &charges)
{
    this->charges_ptr = &charges;
    // If there's only level, then only leaf-level evaluations are performed:
    if(this->max_levels > 1)
    {
        cout << "Getting the charges at the nodes of the leaf level" << endl;
        FMM2DTree::assignLeafCharges();
        cout << "Performing upward sweep to get charges" << endl;
        FMM2DTree::upwardTraveral();
        cout << "Performing M2L..." << endl;
        if(this->is_homogeneous || this->is_log_homogeneous)
            FMM2DTree::evaluateAllM2LHomogeneous();
        else
            FMM2DTree::evaluateAllM2L();
        cout << "Performing downward sweep" << endl;
        FMM2DTree::downwardTraversal();
        cout << "Getting potentials for particles using L2P" << endl;
        FMM2DTree::evaluateL2P();
    }

    cout << "Evaluation the direct interactions at the leaf levels" << endl;
    FMM2DTree::evaluateLeafLevelInteractions();

    this->potentials.eval();
    return this->potentials;
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
    this->potentials = af::array(this->M_ptr->getNumCols(), f64);
    this->max_levels = N_levels;

    // Depending upon whether the points are at the location
    // of the nodes on the leaf level, we'll set the following flag:
    if(N_levels == 0)
        this->user_def_locations = true;

    else
        this->user_def_locations = false;
    
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
    af_print(child_nodes);
    this->standard_nodes_child[2] = child_nodes;

    // For child3:
    child_nodes = af::join(1,
                           0.5 * (this->standard_nodes(af::span, 0) - 1),
                           0.5 * (this->standard_nodes(af::span, 1) + 1)
                          );
    child_nodes.eval();
    this->standard_nodes_child[3] = child_nodes;

    if(radius == -1)
    {
        // Finding out the radius and center for the domain of consideration:
        determineCenterAndRadius((*(this->M_ptr->getSourceCoordsPtr()))(af::span, 0), this->c_x, this->r_x);
        determineCenterAndRadius((*(this->M_ptr->getSourceCoordsPtr()))(af::span, 1), this->c_y, this->r_y);

        // If we use the above obtained radii, we will get some particles on the edge of the domain
        // To account for these edge-cases, we will slightly increase the radius of the domain:
        this->r_x = (1 + 1e-14) * this->r_x;
        this->r_y = (1 + 1e-14) * this->r_y;
    }

    // If the set of points isn't given then we are solving in a domain that uses the radius
    // that is given by the user. Additionally, it is assumed that the box is centered at (0, 0)
    else
    {
        this->c_x = this->c_y = 0;
        this->r_x = this->r_y = radius;
    }


    // Determining nature of the underlying problem:
    this->is_translation_invariant = this->M_ptr->isTranslationInvariant(2);
    if(this->is_translation_invariant == false)
    {
        cout << "Non translation-invariant Kernels not supported!!!" << endl;
        exit(1);
    }

    // Checking that the radius of r_x and r_y are approximately the same:
    // TODO: Need to see if this is a good enough treshold:
    if(2 * fabs(this->r_x - this->r_y) / (this->r_x + this->r_y) < 1e-3)
    {
        this->is_homogeneous     = this->M_ptr->isHomogeneous(2);
        this->is_log_homogeneous = this->M_ptr->isLogHomogeneous(2);
        
        if(this->is_homogeneous || this->is_log_homogeneous)
        {
            this->alpha = this->M_ptr->getDegreeOfHomog(2);
        } 
    }

    else
    {
        this->is_homogeneous     = false;
        this->is_log_homogeneous = false;
    }

    // Getting Transfer Matrices(that is M2M and L2L):
    cout << "Getting Transfer Matrices..." << endl;
    FMM2DTree::getTransferMatrices();
    cout << "Building Tree..." << endl;
    FMM2DTree::buildTree();
    cout << "Number of Levels in the Tree: " << this->max_levels << endl;
    cout << "Assigning Relations Amongst Boxes in the tree..." << endl;
    FMM2DTree::assignTreeRelations();
    cout << "Getting M2L interactions" << endl;
    
    if(this->is_homogeneous || this->is_log_homogeneous)
        FMM2DTree::getM2LInteractionsHomogeneous();
    else
        FMM2DTree::getM2LInteractions();
    
    if(this->user_def_locations == false)
        FMM2DTree::getNeighborSelfInteractions();
}

// ======================== END OF PUBLIC FUNCTIONS =======================
// ===== FUNCTIONS USED IN TREE BUILDING AND PRECOMPUTING OPERATIONS INVOLVED =====
void FMM2DTree::getTransferMatrices()
{
    array nodes_x, nodes_y, L2L_array;
    // Initializing L2L array
    L2L_array = array(N_nodes * N_nodes, N_nodes * N_nodes, f64);

    // L2L is used to tranfer the information from the parent to the child
    for(unsigned short N_child = 0; N_child < 4; N_child++)
    {
        nodes_x = (this->standard_nodes_child[N_child])(af::span, 0);
        nodes_y = (this->standard_nodes_child[N_child])(af::span, 1);

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
    if(this->is_homogeneous || this->is_log_homogeneous)
    {
        this->box_radius.push_back(root.r_x);
        this->box_homog_radius.push_back(pow(root.r_x, this->alpha));
        this->box_log_homog_radius.push_back(this->alpha * log(root.r_x));
    }

    // Using this flag to check whether we have reached leaf level:
    bool reached_leaf = false;
    unsigned N_level = 0;
    while(!reached_leaf)
    {
        N_level++;
        this->number_of_boxes.push_back(4 * this->number_of_boxes[N_level - 1]);
        
        if(this->is_homogeneous || this->is_log_homogeneous)
        {
            this->box_radius.push_back(0.5 * this->box_radius[N_level - 1]);
            this->box_homog_radius.push_back(pow(0.5, this->alpha) * this->box_homog_radius[N_level - 1]);
            this->box_log_homog_radius.push_back(  this->box_log_homog_radius[N_level - 1]
                                                 - this->alpha * log(2)
                                                );
        }

        std::vector<FMM2DBox> level;
        level.reserve(this->number_of_boxes[N_level]);
        for(unsigned i = 0; i < this->number_of_boxes[N_level]; i++) 
        {
            FMM2DBox box;
            box.N_level = N_level;
            box.N_box   = i;
            box.parent  = i / 4;

            FMM2DBox &parent_node = (tree.back())[box.parent];

            // Assigning the new centers and radii:
            box.c_x = parent_node.c_x + (((i - 4 * box.parent) % 2) - 0.5) * parent_node.r_x;
            box.c_y = parent_node.c_y + (((i - 4 * box.parent) / 2) - 0.5) * parent_node.r_y;
            box.r_x = 0.5 * parent_node.r_x;
            box.r_y = 0.5 * parent_node.r_y;

            // When the charges at leaf level are at non-standard locations:
            if(this->max_levels == 0)
            {
                // Locations local to the parent box:
                array locations_local_x = (*(this->M_ptr->getSourceCoordsPtr()))(parent_node.inds_in_box, 0); 
                array locations_local_y = (*(this->M_ptr->getSourceCoordsPtr()))(parent_node.inds_in_box, 1); 

                // Finding the indices in the children:
                box.inds_in_box =
                parent_node.inds_in_box(af::where(   locations_local_x < box.c_x + box.r_x
                                                  && locations_local_x >= box.c_x - box.r_x
                                                  && locations_local_y < box.c_y + box.r_y
                                                  && locations_local_y >= box.c_y - box.r_y
                                                 )
                                       );

                box.inds_in_box.eval();

                // If even a single box shows that it contains < 4 * rank particles, we terminate:
                if(box.inds_in_box.elements()<= 4 * this->rank)
                {
                    reached_leaf = true;
                }
            }

            else
            {
                if(N_level == this->max_levels)
                    reached_leaf = true;
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

void FMM2DTree::getM2LInteractionsHomogeneous() 
{
    std::array<array, 24> M2L_outer_homogeneous;
    std::array<array, 16> M2L_inner_homogeneous;

    array M2L, nodes;
    // The M2L array transfers information from the
    // surrounding box to the box of concern:
    // ============= Getting outer interations ============================
    // Getting interaction between boxes on the lower edge:
    for(unsigned i = 0; i < 6; i++)
    {
        nodes= af::join(1, this->standard_nodes(af::span, 0) + 2 * (i - 3),
                        this->standard_nodes(af::span, 1) - 6
                       );

        getM2L(this->standard_nodes, nodes, *(this->M_ptr), M2L);
        M2L.eval();
        M2L_outer_homogeneous[i] = M2L;
    }

    // Getting interaction between boxes on the right edge:
    for(unsigned i = 0; i < 6; i++)
    {
        nodes= af::join(1, this->standard_nodes(af::span, 0) + 6,
                        this->standard_nodes(af::span, 1) + 2 * (i - 3)
                       );

        getM2L(this->standard_nodes, nodes, *(this->M_ptr), M2L);
        M2L.eval();
        M2L_outer_homogeneous[i + 6] = M2L;
    }

    // Getting interaction between boxes on the top edge:
    for(unsigned i = 0; i < 6; i++)
    {
        nodes= af::join(1, this->standard_nodes(af::span, 0) + 2 * (3 - i),
                        this->standard_nodes(af::span, 1) + 6
                       );

        getM2L(this->standard_nodes, nodes, *(this->M_ptr), M2L);
        M2L.eval();
        M2L_outer_homogeneous[i + 12] = M2L;
    }

    // Getting interaction between boxes on the left edge:
    for(unsigned i = 0; i < 6; i++)
    {
        nodes= af::join(1, this->standard_nodes(af::span, 0) - 6,
                        this->standard_nodes(af::span, 1) + 2 * (3 - i)
                       );

        getM2L(this->standard_nodes, nodes, *(this->M_ptr), M2L);
        M2L.eval();
        M2L_outer_homogeneous[i + 18] = M2L;
    }

    // ============= Getting inner interations ============================
    // Getting interaction between boxes on the lower edge:
    for(unsigned i = 0; i < 4; i++)
    {
        nodes= af::join(1, this->standard_nodes(af::span, 0) + 2 * (i - 2),
                        this->standard_nodes(af::span, 1) - 4
                       );

        getM2L(this->standard_nodes, nodes, *(this->M_ptr), M2L);
        M2L.eval();
        M2L_inner_homogeneous[i] = M2L;
    }

    // Getting interaction between boxes on the right edge:
    for(unsigned i = 0; i < 4; i++)
    {
        nodes= af::join(1, this->standard_nodes(af::span, 0) + 4,
                        this->standard_nodes(af::span, 1) + 2 * (i - 2)
                       );

        getM2L(this->standard_nodes, nodes, *(this->M_ptr), M2L);
        M2L.eval();
        M2L_inner_homogeneous[i + 4] = M2L;
    }

    // Getting interaction between boxes on the top edge:
    for(unsigned i = 0; i < 4; i++)
    {
        nodes= af::join(1, this->standard_nodes(af::span, 0) + 2 * (2 - i),
                        this->standard_nodes(af::span, 1) + 4
                       );

        getM2L(this->standard_nodes, nodes, *(this->M_ptr), M2L);
        M2L.eval();
        M2L_inner_homogeneous[i + 8] = M2L;
    }

    // Getting interaction between boxes on the left edge:
    for(unsigned i = 0; i < 4; i++)
    {
        nodes= af::join(1, this->standard_nodes(af::span, 0) - 4,
                        this->standard_nodes(af::span, 1) + 2 * (2 - i)
                       );

        getM2L(this->standard_nodes, nodes, *(this->M_ptr), M2L);
        M2L.eval();
        M2L_inner_homogeneous[i + 12] = M2L;
    }

    this->M2L_outer.push_back(M2L_outer_homogeneous);
    this->M2L_inner.push_back(M2L_inner_homogeneous);
}

void FMM2DTree::getM2LInteractions() 
{
    // M2L inner and outer operators for a single level:
    std::array<array, 24> M2L_outer_level;
    std::array<array, 16> M2L_inner_level;

    array M2L, nodes, scaled_standard_nodes_x, scaled_standard_nodes_y, scaled_standard_nodes;
    // Starting at L2 since only from then onwards do we have interaction lists:
    for(unsigned i = 2; i <= this->max_levels; i++)
    {
        // Finding the radii of the boxes at this level:
        // Note that although we've considered the box number 0, any box number maybe considered
        // (Box number doesn't matter since we're using a uniform tree)
        double r_x = this->tree[i][0].r_x;
        double r_y = this->tree[i][0].r_y;

        // Scaling to the current level's box radius
        scalePoints(0, 1, this->standard_nodes(af::span, 0),
                    0, r_x, scaled_standard_nodes_x
                   );

        scalePoints(0, 1, this->standard_nodes(af::span, 1),
                    0, r_y, scaled_standard_nodes_y
                   );

        scaled_standard_nodes = af::join(1, scaled_standard_nodes_x, scaled_standard_nodes_y);

        // The M2L array transfers information from the
        // surrounding box to the box of concern:
        // ============= Getting outer interations ============================
        // Getting interaction between boxes on the lower edge:
        for(unsigned i = 0; i < 6; i++)
        {
            nodes= af::join(1, scaled_standard_nodes_x + 2 * r_x * (i - 3),
                            scaled_standard_nodes_y - 6 * r_y
                           );

            getM2L(scaled_standard_nodes, nodes, *(this->M_ptr), M2L);
            M2L.eval();
            M2L_outer_level[i] = M2L;
        }

        // Getting interaction between boxes on the right edge:
        for(unsigned i = 0; i < 6; i++)
        {
            nodes= af::join(1, scaled_standard_nodes_x + 6 * r_x,
                            scaled_standard_nodes_y + 2 * r_y * (i - 3)
                           );

            getM2L(scaled_standard_nodes, nodes, *(this->M_ptr), M2L);
            M2L.eval();
            M2L_outer_level[i + 6] = M2L;
        }

        // Getting interaction between boxes on the top edge:
        for(unsigned i = 0; i < 6; i++)
        {
            nodes= af::join(1, scaled_standard_nodes_x + 2 * r_x * (3 - i),
                            scaled_standard_nodes_y + 6 * r_y
                           );

            getM2L(scaled_standard_nodes, nodes, *(this->M_ptr), M2L);
            M2L.eval();
            M2L_outer_level[i + 12] = M2L;
        }

        // Getting interaction between boxes on the left edge:
        for(unsigned i = 0; i < 6; i++)
        {
            nodes= af::join(1, scaled_standard_nodes_x - 6 * r_x,
                            scaled_standard_nodes_y + 2 * r_y * (3 - i)
                           );

            getM2L(scaled_standard_nodes, nodes, *(this->M_ptr), M2L);
            M2L.eval();
            M2L_outer_level[i + 18] = M2L;
        }

        // ====================================================================

        // ============= Getting inner interations ============================
        // Getting interaction between boxes on the lower edge:
        for(unsigned i = 0; i < 4; i++)
        {
            nodes= af::join(1, scaled_standard_nodes_x + 2 * r_x * (i - 2),
                            scaled_standard_nodes_y - 4 * r_y
                           );

            getM2L(scaled_standard_nodes, nodes, *(this->M_ptr), M2L);
            M2L.eval();
            M2L_inner_level[i] = M2L;
        }

        // Getting interaction between boxes on the right edge:
        for(unsigned i = 0; i < 4; i++)
        {
            nodes= af::join(1, scaled_standard_nodes_x + 4 * r_x,
                            scaled_standard_nodes_y + 2 * r_y * (i - 2)
                           );

            getM2L(scaled_standard_nodes, nodes, *(this->M_ptr), M2L);
            M2L.eval();
            M2L_inner_level[i + 4] = M2L;
        }

        // Getting interaction between boxes on the top edge:
        for(unsigned i = 0; i < 4; i++)
        {
            nodes= af::join(1, scaled_standard_nodes_x + 2 * r_x * (2 - i),
                            scaled_standard_nodes_y + 4 * r_y
                           );

            getM2L(scaled_standard_nodes, nodes, *(this->M_ptr), M2L);
            M2L.eval();
            M2L_inner_level[i + 8] = M2L;
        }

        // Getting interaction between boxes on the left edge:
        for(unsigned i = 0; i < 4; i++)
        {
            nodes= af::join(1, scaled_standard_nodes_x - 4 * r_x,
                            scaled_standard_nodes_y + 2 * r_y * (2 - i)
                           );

            getM2L(scaled_standard_nodes, nodes, *(this->M_ptr), M2L);
            M2L.eval();
            M2L_inner_level[i + 12] = M2L;
        }
    
        this->M2L_outer.push_back(M2L_outer_level);
        this->M2L_inner.push_back(M2L_inner_level);
    }
}


void FMM2DTree::getNeighborSelfInteractions() 
{
    array M2L, nodes, leaf_nodes;
    
    // Get neighbor interactions:
    double r_x_leaf = this->tree[this->max_levels][0].r_x; // since its a uniform tree, 
                                                           // we don't bother about the box considered
    double r_y_leaf = this->tree[this->max_levels][0].r_y;

    leaf_nodes = af::join(1, this->standard_nodes(af::span, 0) * r_x_leaf,
                          this->standard_nodes(af::span, 1) * r_y_leaf
                         );

    // Finding interaction with boxes on the 
    for(unsigned short i = 0; i < 2; i++)
    {
        nodes= af::join(1, (this->standard_nodes(af::span, 0) + 2 * (i - 1)) * r_x_leaf,
                        (this->standard_nodes(af::span, 1) - 2) * r_y_leaf
                       );

        getM2L(leaf_nodes, nodes, *(this->M_ptr), M2L);
        M2L.eval();
        this->neighbor_interaction[i] = M2L;
    }

    for(unsigned short i = 0; i < 2; i++)
    {
        nodes= af::join(1, (this->standard_nodes(af::span, 0) + 2) * r_x_leaf,
                        (this->standard_nodes(af::span, 1) + 2 * (i - 1)) * r_y_leaf
                       );
        getM2L(leaf_nodes, nodes, *(this->M_ptr), M2L);
        M2L.eval();
        this->neighbor_interaction[i + 2] = M2L;
    }


    for(unsigned short i = 0; i < 2; i++)
    {
        nodes= af::join(1, (this->standard_nodes(af::span, 0) + 2 * (1 - i)) * r_x_leaf,
                        (this->standard_nodes(af::span, 1) + 2) * r_y_leaf
                       );

        getM2L(this->standard_nodes, nodes, *(this->M_ptr), M2L);
        M2L.eval();
        this->neighbor_interaction[i + 4] = M2L;
    }

    for(unsigned short i = 0; i < 2; i++)
    {
        nodes= af::join(1, (this->standard_nodes(af::span, 0) - 2) * r_x_leaf,
                        (this->standard_nodes(af::span, 1) + 2 * (1 - i)) * r_y_leaf
                       );

        getM2L(this->standard_nodes, nodes, *(this->M_ptr), M2L);
        M2L.eval();
        this->neighbor_interaction[i + 6] = M2L;
    }

    // Getting self-interaction:
    getM2L(leaf_nodes, leaf_nodes, *(this->M_ptr), this->self_interaction);
    this->self_interaction.eval();
}

// ===== END OF FUNCTIONS USED IN TREE BUILDING AND PRECOMPUTING OPERATIONS INVOLVED =====

// ========== FUNCTIONS USED IN EVALUATION OF POTENTIAL =============

void FMM2DTree::assignLeafCharges()
{
    // Looping over the boxes at the leaf level:
    for(unsigned i = 0; i < this->number_of_boxes[this->max_levels]; i++)
    {   
        // Getting the box that is in consideration:
        FMM2DBox &node = this->tree[this->max_levels][i];

        array std_locations_x, std_locations_y;
        // Mapping onto the standard interval of [-1, 1]:
        scalePoints(node.c_x, node.r_x, (*(this->M_ptr->getSourceCoordsPtr()))(node.inds_in_box, 0),
                    0, 1, std_locations_x
                   );

        scalePoints(node.c_x, node.r_x, (*(this->M_ptr->getSourceCoordsPtr()))(node.inds_in_box, 1),
                    0, 1, std_locations_y
                   );

        // Getting the P2M operator which interpolates information from the particles to nodes:
        array P2M_array = array(node.inds_in_box.elements(), N_nodes * N_nodes, f64);
        getL2L2D(std_locations_x, std_locations_y, this->standard_nodes_1d, P2M_array);
        tree[this->max_levels][i].node_charges = af::matmul(P2M_array.T(), (*this->charges_ptr)(node.inds_in_box));
        tree[this->max_levels][i].node_charges.eval();
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

void FMM2DTree::evaluateAllM2LHomogeneous()
{
    // We start from level 2 since L1 doesn't have any boxes in its interaction list
    for(unsigned N_level = 2; N_level <= this->max_levels; N_level++) 
    {
        for(unsigned N_box = 0; N_box < this->number_of_boxes[N_level]; N_box++)
        {
            // Getting the box in consideration:
            FMM2DBox &box = this->tree[N_level][N_box];
            // Initializing the value for the potentials:
            box.node_potentials = af::constant(0, this->rank, f64);

            if(this->is_homogeneous)
            {
                // Inner well-separated clusters
                for(unsigned short i = 0; i < 16; i++) 
                {
                    int N_inner = box.inner[i];
                    if(N_inner > -1) 
                    {
                        box.node_potentials += af::matmul(this->M2L_inner[0][i],
                                                          this->tree[N_level][N_inner].node_charges
                                                         );
                    }
                }

                // Outer well-separated clusters
                for(unsigned short i = 0; i < 24; i++) 
                {
                    int N_outer = box.outer[i];
                    if(N_outer > -1) 
                    {
                        box.node_potentials += af::matmul(this->M2L_outer[0][i],
                                                          this->tree[N_level][N_outer].node_charges
                                                         );
                    }
                }

                // Applying the scaling factor to account for level:
                box.node_potentials *= this->box_homog_radius[N_level];                   
            }

            else if(this->is_log_homogeneous)
            {
                // Inner well-separated clusters
                for(unsigned short i = 0; i < 16; i++) 
                {
                    int N_inner = box.inner[i];
                    if(N_inner > -1) 
                    {
                        box.node_potentials += af::matmul(this->M2L_inner[0][i],
                                                          this->tree[N_level][N_inner].node_charges
                                                         );

                        // Here the scaling factor to account for level shows up as an addition operator:
                        box.node_potentials +=   this->box_log_homog_radius[N_level]
                                               * af::tile(af::sum(this->tree[N_level][N_inner].node_charges),
                                                          this->rank
                                                         );
                    }
                }

                //  Outer well-separated clusters
                for(unsigned short i = 0; i < 24; i++) 
                {
                    int N_outer = box.outer[i];
                    if(N_outer > -1) 
                    {
                        box.node_potentials += af::matmul(this->M2L_outer[0][i],
                                                          this->tree[N_level][N_outer].node_charges
                                                         );

                        // Here the scaling factor to account for level shows up as an addition operator:
                        box.node_potentials +=   this->box_log_homog_radius[N_level]
                                               * af::tile(af::sum(this->tree[N_level][N_outer].node_charges),
                                                          this->rank
                                                         );
                    }
                }
            }

            // Performing eval before proceeding to next box
            box.node_potentials.eval();
        }
    }
}

void FMM2DTree::evaluateAllM2L()
{
    // We start from level 2 since L1 doesn't have any boxes in its interaction list
    for(unsigned N_level = 2; N_level <= this->max_levels; N_level++) 
    {
        for(unsigned N_box = 0; N_box < this->number_of_boxes[N_level]; N_box++)
        {
            // Getting the box in consideration:
            FMM2DBox &box = this->tree[N_level][N_box];
            // Initializing the value for the potentials:
            box.node_potentials = af::constant(0, this->rank, f64);
            // Inner well-separated clusters
            for(unsigned short i = 0; i < 16; i++) 
            {
                int N_inner = box.inner[i];
                if(N_inner > -1) 
                {   
                    box.node_potentials += matmul(this->M2L_inner[N_level-2][i],
                                                  this->tree[N_level][N_inner].node_charges
                                                 );
                }
            }

            // Outer well-separated clusters
            for(unsigned short i = 0; i < 24; i++) 
            {
                int N_outer = box.outer[i];
                if(N_outer > -1) 
                {
                    box.node_potentials += matmul(this->M2L_outer[N_level-2][i],
                                                  this->tree[N_level][N_outer].node_charges
                                                 );
                }
            }

            // Performing eval before proceeding to next box
            box.node_potentials.eval();
        }
    }
}

void FMM2DTree::downwardTraversal()
{
    for(unsigned N_level = 2; N_level < this->max_levels; N_level++) 
    {
        // Level number for the child:
        int N_lc = N_level + 1;
        int N_bc0;

        for(unsigned N_box = 0; N_box < this->number_of_boxes[N_level]; N_box++) 
        {
            // Box number for child0. Box numbers for the remaining children are defined w.r.t to this:
            N_bc0 = 4 * N_box;

            for(unsigned short N_child = 0; N_child < 4; N_child++)
            {
                tree[N_lc][N_bc0 + N_child].node_potentials += af::matmul(L2L[N_child],
                                                                          tree[N_level][N_box].node_potentials
                                                                         );
                tree[N_lc][N_bc0 + N_child].node_potentials.eval();
            }
        }
    }
}

void FMM2DTree::evaluateL2P()
{
    array std_locations_x, std_locations_y;
    // Looping over the boxes at the leaf level:
    for(unsigned i = 0; i < this->number_of_boxes[this->max_levels]; i++)
    {   
        // Getting the box that is in consideration:
        FMM2DBox &node = this->tree[this->max_levels][i];

        // Mapping onto the standard interval of [-1, 1]:
        scalePoints(node.c_x, node.r_x, (*(this->M_ptr->getSourceCoordsPtr()))(node.inds_in_box, 0),
                    0, 1, std_locations_x
                   );

        scalePoints(node.c_x, node.r_x, (*(this->M_ptr->getSourceCoordsPtr()))(node.inds_in_box, 1),
                    0, 1, std_locations_y
                   );

        // Getting the P2M operator:
        array L2P = array(node.inds_in_box.elements(), N_nodes * N_nodes, f64);
        getL2L2D(std_locations_x, std_locations_y, this->standard_nodes_1d, L2P);
        this->potentials(node.inds_in_box) = af::matmul(L2P, node.node_potentials);
    }

    this->potentials.eval();
}

void FMM2DTree::evaluateLeafLevelInteractions()
{
    // Looping over the boxes at the leaf level:
    for(unsigned i = 0; i < this->number_of_boxes[this->max_levels]; i++)
    {   
        // Getting the box that is in consideration:
        FMM2DBox &box = this->tree[this->max_levels][i];
        int N_neighbor;
        array targets = (*(this->M_ptr->getSourceCoordsPtr()))(box.inds_in_box, af::span);

        for(unsigned j = 0; j < 8; j++) 
        {
            N_neighbor = this->tree[this->max_levels][i].neighbor[j];
    
            if(N_neighbor > -1) 
            {
                FMM2DBox &neighbor_box = tree[this->max_levels][N_neighbor];
                array sources = (*(this->M_ptr->getSourceCoordsPtr()))(neighbor_box.inds_in_box, af::span);
                array charges = (*(this->charges_ptr))(neighbor_box.inds_in_box);
                // Updating the potentials for the particles in the box of consideration
                // by evaluating the direct interaction with the particles in the neighbor box
                this->potentials(box.inds_in_box) +=
                af::matmul(this->M_ptr->buildArray(targets, sources), charges);
            }
        }

        // Self interaction:
        array charges = (*(this->charges_ptr))(box.inds_in_box);
        this->potentials(box.inds_in_box) +=
        af::matmul(this->M_ptr->buildArray(targets, targets),charges);
    }
}

// ========== END OF FUNCTIONS USED IN EVALUATION OF POTENTIAL =============

#endif
