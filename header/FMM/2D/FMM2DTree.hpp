#ifndef __2DTree_hpp__
#define __2DTree_hpp__

#include <cstdlib>
#include "MatrixData.hpp"
#include "FMM2DBox.hpp"
#include "utils/getStandardNodes.hpp"
#include "utils/determineCenterAndRadius.hpp"
#include "utils/scalePoints.hpp"
#include "MatrixFactorizations/getInterpolation.hpp"

class FMM2DTree 
{
private:
    MatrixData *M_ptr;       // Pointer to the structure describing the underlying problem
    unsigned N;              // Total number of source particles
    unsigned N_nodes;        // Number of nodes used for interpolation
    unsigned rank;           // Rank of the low-rank interaction
    unsigned max_levels;     // Number of levels in the tree.
                       
    const array *charges_ptr;  // Holds the information for the charges of the points
    array standard_nodes_1d;   // Standard nodes alloted as given by user choice in [-1, 1]
    array standard_nodes;      // Gets the points in 2D
    std::string nodes_type;    // Whether Chebyshev, or Legendre or Equispaced

    double c_x, c_y, r_x, r_y; // Radii and the center coordinates for the domain
    // Flags which determine structure for the underlying kernel:
    bool is_homogeneous, is_log_homogeneous, is_translation_invariant;

    std::vector<size_t> number_of_boxes; // Number of boxes at each level
    // The following functions will only be used if the function is homog / loghomog
    std::vector<double> box_radius;
    std::vector<double> box_log_homog_radius;
    std::vector<double> box_homog_radius;
    double alpha; // degree of homog

    std::array<af::array, 4> standard_nodes_child; // Nodes of the child given that parent has standard nodes in [-1, 1]
    std::vector<af::array> L2L;                    // L2L operators from the parent to the four children
    std::vector<af::array> M2M;                    // M2M operators from the children to the parent
    std::array<af::array, 16> M2L_inner;           // M2L operators from the inner well separated clusters
    std::array<af::array, 24> M2L_outer;           // M2L operators from the outer well separated clusters
    std::vector<std::vector<FMM2DBox>> tree;       // The tree storing all the information.

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
    void getM2LInteractions();
    // Gets the charges at the leaf level:
    void assignLeafCharges();
    // Computes the multipoles for all the boxes using upward sweep:
    void upwardTraveral();
    // Computing M2Ls for all the boxes:
    void evaluateAllM2L();
    // Performs the L2L using a downwardsweep:
    void downwardTraversal();
    // Performs the leaf level interactions:
    void evaluateLeafLevelInteractions();
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

    // Finding out the radius and center for the domain of consideration:
    determineCenterAndRadius((*(this->M_ptr->getSourceCoordsPtr()))(af::span, 0), this->c_x, this->r_x);
    determineCenterAndRadius((*(this->M_ptr->getSourceCoordsPtr()))(af::span, 1), this->c_y, this->r_y);

    // Determining nature of the underlying problem:
    this->is_translation_invariant = this->M_ptr->isTranslationInvariant();
    // Checking that the radius of r_x and r_y are approximately the same:
    // TODO: Need to see if this is a good enough treshold:
    if(2 * fabs(this->r_x - this->r_y) / (this->r_x + this->r_y) < 1e-3)
    {
        this->is_homogeneous     = this->M_ptr->isHomogeneous();
        this->is_log_homogeneous = this->M_ptr->isLogHomogeneous();
        
        if(this->is_homogeneous || this->is_log_homogeneous)
        {
            this->alpha = this->M_ptr->getDegreeOfHomog();
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
    cout << "Getting M2L and neighbor interactions" << endl;
    FMM2DTree::getM2LInteractions();
    cout << "Getting the charges at the nodes of the leaf level" << endl;
    FMM2DTree::assignLeafCharges();
    cout << "Performing upward sweep to get charges" << endl;
    FMM2DTree::upwardTraveral();
}

void FMM2DTree::getTransferMatrices()
{
    array nodes_x, nodes_y, L2L_array;
    // Initializing L2L array
    L2L_array = array(N_nodes * N_nodes, N_nodes * N_nodes, f64);

    // Reserving space on the vector:
    // this->L2L.reserve(4);
    // this->M2M.reserve(4);

    for(unsigned short N_child = 0; N < 4; N_child++)
    {
        nodes_x = (this->standard_nodes_child[N_child])(af::span, 0);
        nodes_y = (this->standard_nodes_child[N_child])(af::span, 1);

        getL2L2D(nodes_x, nodes_y, this->standard_nodes_1d, L2L_array);
        L2L_array.eval();
        this->L2L.push_back(L2L_array);
    }

    cout << "Type for nodes" << endl;
    cout << this->L2L[0].type() << endl;

    // M2M is going to be the transpose of the L2L operators:
    for(unsigned short N_child = 0; N < 4; N_child++)
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
    while(!reached_leaf)
    {
        this->max_levels++;
        this->number_of_boxes.push_back(4 * this->number_of_boxes[this->max_levels - 1]);
        
        if(this->is_homogeneous || this->is_log_homogeneous)
        {
            this->box_radius.push_back(0.5 * this->box_radius[this->max_levels - 1]);
            this->box_homog_radius.push_back(pow(0.5, this->alpha) * this->box_homog_radius[this->max_levels - 1]);
            this->box_log_homog_radius.push_back(  this->box_log_homog_radius[this->max_levels - 1]
                                                 - this->alpha * log(2)
                                                );
        }

        std::vector<FMM2DBox> level;
        
        level.reserve(this->number_of_boxes[this->max_levels]);
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
        if(N_leaf_criterion == this->number_of_boxes[this->max_levels])
            reached_leaf = true;
    }
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

void FMM2DTree::getM2LInteractions() 
{
    array M2L, nodes;

    // ============= Getting outer interations ============================
    // Getting interaction between boxes on the lower edge:
    for(unsigned i = 0; i < 6; i++)
    {
        nodes= af::join(1, this->standard_nodes(af::span, 0) + 2 * (i - 3),
                        this->standard_nodes(af::span, 1) - 6
                       );

        getM2L(this->standard_nodes, nodes, *(this->M_ptr), M2L);
        M2L.eval();
        this->M2L_outer[i] = M2L;
    }

    // Getting interaction between boxes on the right edge:
    for(unsigned i = 0; i < 6; i++)
    {
        nodes= af::join(1, this->standard_nodes(af::span, 0) + 6,
                        this->standard_nodes(af::span, 1) + 2 * (i - 3)
                       );

        getM2L(this->standard_nodes, nodes, *(this->M_ptr), M2L);
        M2L.eval();
        this->M2L_outer[i + 6] = M2L;
    }

    // Getting interaction between boxes on the top edge:
    for(unsigned i = 0; i < 6; i++)
    {
        nodes= af::join(1, this->standard_nodes(af::span, 0) + 2 * (3 - i),
                        this->standard_nodes(af::span, 1) + 6
                       );

        getM2L(this->standard_nodes, nodes, *(this->M_ptr), M2L);
        M2L.eval();
        this->M2L_outer[i + 12] = M2L;
    }

    // Getting interaction between boxes on the left edge:
    for(unsigned i = 0; i < 6; i++)
    {
        nodes= af::join(1, this->standard_nodes(af::span, 0) - 6,
                        this->standard_nodes(af::span, 1) + 2 * (3 - i)
                       );

        getM2L(this->standard_nodes, nodes, *(this->M_ptr), M2L);
        M2L.eval();
        this->M2L_outer[i + 18] = M2L;
    }

    // ====================================================================

    // ============= Getting inner interations ============================
    // Getting interaction between boxes on the lower edge:
    for(unsigned i = 0; i < 4; i++)
    {
        nodes= af::join(1, this->standard_nodes(af::span, 0) + 2 * (i - 2),
                        this->standard_nodes(af::span, 1) - 4
                       );

        getM2L(this->standard_nodes, nodes, *(this->M_ptr), M2L);
        M2L.eval();
        this->M2L_inner[i] = M2L;
    }

    // Getting interaction between boxes on the right edge:
    for(unsigned i = 0; i < 4; i++)
    {
        nodes= af::join(1, this->standard_nodes(af::span, 0) + 4,
                        this->standard_nodes(af::span, 1) + 2 * (i - 2)
                       );

        getM2L(this->standard_nodes, nodes, *(this->M_ptr), M2L);
        M2L.eval();
        this->M2L_inner[i + 4] = M2L;
    }

    // Getting interaction between boxes on the top edge:
    for(unsigned i = 0; i < 4; i++)
    {
        nodes= af::join(1, this->standard_nodes(af::span, 0) + 2 * (2 - i),
                        this->standard_nodes(af::span, 1) + 4
                       );

        getM2L(this->standard_nodes, nodes, *(this->M_ptr), M2L);
        M2L.eval();
        this->M2L_inner[i + 8] = M2L;
    }

    // Getting interaction between boxes on the left edge:
    for(unsigned i = 0; i < 4; i++)
    {
        nodes= af::join(1, this->standard_nodes(af::span, 0) - 4,
                        this->standard_nodes(af::span, 1) + 2 * (2 - i)
                       );

        getM2L(this->standard_nodes, nodes, *(this->M_ptr), M2L);
        M2L.eval();
        this->M2L_inner[i + 12] = M2L;
    }
}

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

        // Getting the P2M operator:
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

            cout << this->M2M[0].type() << endl;
            cout << tree[N_lc][N_bc0].node_charges.type() << endl;

            tree[N_level][N_box].node_charges =  af::matmul(this->M2M[0], tree[N_lc][N_bc0].node_charges)
                                               + af::matmul(this->M2M[1], tree[N_lc][N_bc1].node_charges)
                                               + af::matmul(this->M2M[2], tree[N_lc][N_bc2].node_charges)
                                               + af::matmul(this->M2M[3], tree[N_lc][N_bc3].node_charges);
        }
    }
}

// void FMM2DTree::evaluateAllM2L()
// {
//     for(int N_level = 2; N_level <= this->max_levels; N_level++) 
//     {
//         for(unsigned N_box = 0; N_box < pow(4, N_level); N_box++)
//         {
//             // Getting the box in consideration:
//             FMM2DBox &box = this->tree[N_level][N_box];
//             // Initializing the value for the locals:
//             box.locals = af::constant(0, this->rank, f64);

//             if(this->is_homogeneous = true)
//             {
//                 // Inner well-separated clusters
//                 for(unsigned short i = 0; i < 16; i++) 
//                 {
//                     int N_inner = box.inner[i];
//                     if(N_inner > -1) 
//                     {
//                         box.locals += this->M2L_inner[i] * this->tree[N_level][N_inner].multipoles;
//                     }
//                 }

//                 // Outer well-separated clusters
//                 for(unsigned short i = 0; i < 24; i++) 
//                 {
//                     int N_outer = box.outer[i];
//                     if(N_outer > -1) 
//                     {
//                         box.locals += this->M2L_outer[i] * this->tree[N_level][N_outer].multipoles;
//                     }
//                 }

//             // Applying the scaling factor to account for level:
//             box.locals *= this->box_homog_radius[N_level];                   
//             }

//             else if(this->is_log_homogeneous = true)
//             {
//                 // Inner well-separated clusters
//                 for(unsigned short i = 0; i < 16; i++) 
//                 {
//                     int N_inner = box.inner[i];
//                     if(N_inner > -1) 
//                     {
//                         box.locals += this->M2L_inner[i] * this->tree[N_level][N_inner].multipoles;
//                         // Here the scaling factor to account for level shows up as an addition operator:
//                         box.locals +=   this->box_log_homog_radius[N_level]
//                                       * af::tile(af::sum(this->tree[N_level][N_inner].multipoles), this->rank);
//                     }
//                 }

//                 //  Outer well-separated clusters
//                 for(unsigned short i = 0; i < 24; i++) 
//                 {
//                     int N_outer = box.outer[i];
//                     if(N_outer > -1) 
//                     {
//                         box.locals += this->M2L_outer[i] * this->tree[N_level][N_outer].multipoles;
//                         // Here the scaling factor to account for level shows up as an addition operator:
//                         box.locals +=   this->box_log_homog_radius[N_level]
//                                       * af::tile(af::sum(this->tree[N_level][N_inner].multipoles), this->rank);

//                     }
//                 }
//             }

//             else
//             {
//                 cout << "Feature still in progress!!" << endl;
//                 exit(1);
//             }
//         }
//     }
// }

// void FMM2DTree::downwardTraversal()
// {
//     for(int N_level = 2; N_level < pow(4, this->max_levels); N_level++) 
//     {
//         // Level number for the child:
//         int N_lc = N_level + 1;
//         int N_bc0, N_bc1, N_bc2, N_bc3;
        
//         for(unsigned N_box = 0; N_box < pow(4, N_level); N_box++) 
//         {
//             // Box number for child0:
//             N_bc0 = 4 * N_box;
//             // Box number for child1:
//             N_bc1 = 4 * N_box + 1;
//             // Box number for child2:
//             N_bc2 = 4 * N_box + 2;
//             // Box number for child3:
//             N_bc3 = 4 * N_box + 3;

//             tree[N_lc][N_bc3].locals += L2L[0] * tree[N_level][N_box].locals;
//             tree[N_lc][N_bc2].locals += L2L[1] * tree[N_level][N_box].locals;
//             tree[N_lc][N_bc2].locals += L2L[2] * tree[N_level][N_box].locals;
//             tree[N_lc][N_bc3].locals += L2L[3] * tree[N_level][N_box].locals;
//         }
//     }
// }

// void FMM2DTree::evaluateLeafLevelInteractions() 
// {
//     if(this->max_levels < 2) 
//     {
//         for(unsigned N_box = 0; N_box < pow(4, this->max_levels); N_box++) 
//         {
//             tree[this->max_levels][N_box].locals = af::constant(0, rank, f64);
//         }
//     }
    
//     for(unsigned N_box = 0; N_box < pow(4, this->max_levels); N_box++) 
//     {
//         for(int i = 0; i < 8; i++) 
//         {
//             int N_neighbor = tree[this->max_levels][N_box].neighbor[i];
//             if(N_neighbor > -1) 
//             {
//                 tree[this->max_levels][i].locals +=   this->neighbor_interaction[i]
//                                                     * tree[this->max_levels][N_neighbor].multipoles;
//             }
//         }

//         tree[this->max_levels][N_box].locals += self_interaction * tree[this->max_levels][N_box].multipoles;
//     }
// }

#endif
