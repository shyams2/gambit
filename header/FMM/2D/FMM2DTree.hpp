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

// The following nomenclature is used to describe the child cell number:
// ============================
// ||            |           ||
// ||            |           ||
// ||     3      |      4    ||
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

void FMM2DTree::assignChildrenInTree(FMM2DBox *node)
{
    // When number of particles in the box is zero:
    if(node->N == 0)
    {
        node->is_leaf  = true;
        node->is_empty = true;
    }
    
    else
    {
        node.is_leaf  = false;
        node.is_empty = false;
        // Get the scaled nodes for the radius and center of this node:
        scalePoints(0, 1, this->standard_nodes, node->c_x, node->r_x, node->scaled_nodes_x);
        scalePoints(0, 1, this->standard_nodes, node->c_y, node->r_y, node->scaled_nodes_y);

        // The following are operations on leaf cells:
        // We are saying that if N < 4 * N_nodes^2, then terminate
        if(node->N <= 4 * this->rank)
        {
            node.is_leaf = true;

            // Updating max_levels if needed:
            if(this->max_levels < node->N_level)
            {
                this->max_levels = node->N_level;
            }
        }
    
        // If it's not a leaf, then create children:
        else
        {
            for(unsigned i = 0, i < 4, i++)
            {
                FMM2DBox *child = new FMM2DBox;
                child->N_level  = node->N_level + 1;
                child->parent   = node

                // Assigning the new centers and radii:
                child->c_x = node->c_x + ((k % 2) - 0.5) * node->r_x;
                child->c_y = node->c_y + ((k / 2) - 0.5) * node->r_y;
                child->r_x = 0.5 * node->r_x;
                child->r_y = 0.5 * node->r_y;

                // Locations local to this current box:
                array locations_local_x = this->locationGlobal(node->inds_in_box, 0); 
                array locations_local_y = this->locationGlobal(node->inds_in_box, 1); 
                
                // Finding the indices in the children:
                child->inds_in_box = node->inds_in_box(af::where(   location_local < node->c_x + node->r_x
                                                                 && location_local < node->c_x + node->r_x
                                                                 && location_local < node->c_x + node->r_x
                                                                 && location_local < node->c_x + node->r_x
                                                                )
                                                      )
            }
        
        // Now recursively perform the routine for the
        // individual children as well until we hit a leaf node:
        for(unsigned i = 0, i < 4, i++)
            FMM2DBox::assignChildrenInTree(node->child.at(i));
        }
}

void buildTree(FMM2DBox *node)
{
    if(!node.is_empty)
        if(!node.is_eaf)
            assignSiblings(node)
}


void FMM2DTree::assignSiblings(std::shared_ptr<FMM2DBox> node)
{
    // Consider the following numbering scheme for neighbors:
    // ========================
    // ||     ||      ||     ||
    // ||  6  ||   5  ||  4  ||
    // ||     ||      ||     ||
    // ========================
    // ||     ||      ||     ||
    // ||  7  ||   x  ||  3  ||
    // ||     ||      ||     ||
    // ========================
    // ||     ||      ||     ||
    // ||  0  ||   1  ||  2  ||
    // ||     ||      ||     ||
    // ========================

    // Assigning the siblings to child[0]:
    ((node->child.at(0))->neighbor).at(3)      = node->child.at(1);
    ((node->child.at(0))->neighbor).at(5)      = node->child.at(2);
    ((node->child.at(0))->neighbor).at(4)      = node->child.at(3);
    (node->child.at(0))->N_neighbor_allocated += 3;

    // Assigning the siblings to child[1]:
    ((node->child.at(1))->neighbor).at(7)      = node->child.at(0);
    ((node->child.at(1))->neighbor).at(6)      = node->child.at(2);
    ((node->child.at(1))->neighbor).at(5)      = node->child.at(3);
    (node->child.at(0))->N_neighbor_allocated += 3;

    // Assigning the siblings to child[2]:
    ((node->child.at(2))->neighbor).at(1)      = node->child.at(0);
    ((node->child.at(2))->neighbor).at(2)      = node->child.at(1);
    ((node->child.at(2))->neighbor).at(3)      = node->child.at(3);
    (node->child.at(0))->N_neighbor_allocated += 3;

    // Assigning the siblings to child[3]:
    ((node->child.at(3))->neighbor).at(0)      = node->child.at(0);
    ((node->child.at(3))->neighbor).at(1)      = node->child.at(1);
    ((node->child.at(3))->neighbor).at(7)      = node->child.at(2);
    (node->child.at(0))->N_neighbor_allocated += 3;
}

void FMM2DTree::assignFourWellSeparatedCousins(std::shared_ptr<FMM2DBox> node, 
                                               unsigned neighbor_number,
                                               unsigned child_number
                                              )
{
      (node->child.at(child_number)->interaction).at(node->child.at(child_number)->N_interaction_allocated + 0) 
    = (node->neighbor.at(neighbor_number))->child.at(0);

      (node->child.at(child_number)->interaction).at(node->child.at(child_number)->N_interaction_allocated + 1) 
    = (node->neighbor.at(neighbor_number))->child.at(1);

      (node->child.at(child_number)->interaction).at(node->child.at(child_number)->N_interaction_allocated + 2) 
    = (node->neighbor.at(neighbor_number))->child.at(2);

      (node->child.at(child_number)->interaction).at(node->child.at(child_number)->N_interaction_allocated + 3) 
    = (node->neighbor.at(neighbor_number))->child.at(3);

    node->child.at(child_number)->N_interaction_allocated += 4;
}

void FMM2DTree::assignCousins(std::shared_ptr<FMM2DBox> node)
{
    // Assigning children of neighbour 0
    {
        // Assigning cousins to child[0]:
        // Three well-separated cousins:
          (node->child.at(0)->interaction).at(node->child.at(0)->N_interaction_allocated + 0) 
        = (node->neighbor.at(0))->child.at(0);

          (node->child.at(0)->interaction).at(node->child.at(0)->N_interaction_allocated + 1) 
        = (node->neighbor.at(0))->child.at(1);

          (node->child.at(0)->interaction).at(node->child.at(0)->N_interaction_allocated + 2) 
        = (node->neighbor.at(0))->child.at(2);
        
        // One neighbor:
        (node->child.at(0)->neighbor).at(0) = (node->neighbor.at(0))->child.at(3);
        node->child.at(0)->N_interaction_allocated += 3;
        node->child.at(0)->N_neighbor_allocated    += 1;

        // Assigning cousins to child[1]:
        // Four well-separated cousins:
        FMM2DTree::assignFourWellSeparatedCousins(node, 0, 1);

        // Assigning cousins to child[2]:
        // Four well-separated cousins:
        FMM2DTree::assignFourWellSeparatedCousins(node, 0, 2);

        // Assigning cousins to child[3]:
        // Four well-separated cousins:
        FMM2DTree::assignFourWellSeparatedCousins(node, 0, 3);
    }

    // Assigning children of neighbour 1
    {
        // Assigning cousins to child[0]:
        // Two well-separated cousins:
          (node->child.at(0)->interaction).at(node->child.at(0)->N_interaction_allocated + 0) 
        = (node->neighbor.at(1))->child.at(0);

          (node->child.at(0)->interaction).at(node->child.at(0)->N_interaction_allocated + 1) 
        = (node->neighbor.at(1))->child.at(1);

        // Two neighbors:
        (node->child.at(0)->neighbor).at(1) = (node->neighbor.at(1))->child.at(2);
        (node->child.at(0)->neighbor).at(2) = (node->neighbor.at(1))->child.at(3);

        node->child.at(0)->N_interaction_allocated += 2;
        node->child.at(0)->N_neighbor_allocated    += 2;

        // Assigning cousins to child[1]:
        // Two well-separated cousins:
          (node->child.at(1)->interaction).at(node->child.at(1)->N_interaction_allocated + 0) 
        = (node->neighbor.at(1))->child.at(0);

          (node->child.at(1)->interaction).at(node->child.at(1)->N_interaction_allocated + 1) 
        = (node->neighbor.at(1))->child.at(1);

        // Two neighbors:
        (node->child.at(1)->neighbor).at(0) = (node->neighbor.at(1))->child.at(2);
        (node->child.at(1)->neighbor).at(1) = (node->neighbor.at(1))->child.at(3);

        node->child.at(1)->N_interaction_allocated += 2;
        node->child.at(1)->N_neighbor_allocated    += 2;

        // Assigning cousins to child[2]:
        // Four well-separated cousins:
        FMM2DTree::assignFourWellSeparatedCousins(node, 1, 2);

        // Assigning cousins to child[3]:
        // Four well-separated cousins:
        FMM2DTree::assignFourWellSeparatedCousins(node, 1, 3);
    }

    // Assigning children of neighbour 2
    {
        // Assigning cousins to child[0]:
        // Four well-separated cousins:
        FMM2DTree::assignFourWellSeparatedCousins(node, 2, 0);
        // Assigning cousins to child[1]:
        // Three well-separated cousins:
          (node->child.at(1)->interaction).at(node->child.at(1)->N_interaction_allocated + 0) 
        = (node->neighbor.at(2))->child.at(0);

          (node->child.at(1)->interaction).at(node->child.at(1)->N_interaction_allocated + 1) 
        = (node->neighbor.at(2))->child.at(1);

          (node->child.at(1)->interaction).at(node->child.at(1)->N_interaction_allocated + 2) 
        = (node->neighbor.at(2))->child.at(3);

        // One neighbor:
        (node->child.at(1)->neighbor).at(2) = (node->neighbor.at(2))->child.at(2);
  
        node->child.at(1)->N_interaction_allocated += 3;
        node->child.at(1)->N_neighbor_allocated    += 1;

        // Assigning cousins to child[2]:
        // Four well-separated cousins:
        FMM2DTree::assignFourWellSeparatedCousins(node, 2, 2);

        // Assigning cousins to child[3]:
        // Four well-separated cousins:
        FMM2DTree::assignFourWellSeparatedCousins(node, 2, 3);
    }
        node.child[3].nInteraction += 2

        # Update neighbor count.
        node.child[1].nNeighbor += 2
        node.child[3].nNeighbor += 2

    // Assigning children of neighbour 3
    {
        // Assigning cousins to child[0]:
        // Four well-separated cousins:
        FMM2DTree::assignFourWellSeparatedCousins(node, 3, 0);

        // Assigning cousins to child[1]:
        // Two well-separated cousins:
          (node->child.at(1)->interaction).at(node->child.at(1)->N_interaction_allocated + 0) 
        = (node->neighbor.at(3))->child.at(1);

          (node->child.at(1)->interaction).at(node->child.at(1)->N_interaction_allocated + 1) 
        = (node->neighbor.at(3))->child.at(3);

        // Two neighbors:
        (node->child.at(1)->neighbor).at(3) = (node->neighbor.at(3))->child.at(0);
        (node->child.at(1)->neighbor).at(4) = (node->neighbor.at(3))->child.at(2);
  
        node->child.at(1)->N_interaction_allocated += 2;
        node->child.at(1)->N_neighbor_allocated    += 2;

        // Assigning cousins to child[2]:
        // Four well-separated cousins:
        FMM2DTree::assignFourWellSeparatedCousins(node, 3, 2);

        // Assigning cousins to child[3]:
        // Two well-separated cousins:
          (node->child.at(3)->interaction).at(node->child.at(3)->N_interaction_allocated + 0) 
        = (node->neighbor.at(3))->child.at(1);

          (node->child.at(3)->interaction).at(node->child.at(3)->N_interaction_allocated + 1) 
        = (node->neighbor.at(3))->child.at(3);

        // Two neighbors:
        (node->child.at(3)->neighbor).at(2) = (node->neighbor.at(3))->child.at(0);
        (node->child.at(3)->neighbor).at(3) = (node->neighbor.at(3))->child.at(2);

        node->child.at(3)->N_interaction_allocated += 2;
        node->child.at(3)->N_neighbor_allocated    += 2;
    }

    // Assigning children of neighbour 4
    {
        // Assigning cousins to child[0]:
        // Four well-separated cousins:
        FMM2DTree::assignFourWellSeparatedCousins(node, 4, 0);

        // Assigning cousins to child[1]:
        // Four well-separated cousins:
        FMM2DTree::assignFourWellSeparatedCousins(node, 4, 1);

        // Assigning cousins to child[2]:
        // Four well-separated cousins:
        FMM2DTree::assignFourWellSeparatedCousins(node, 4, 2);

        // Assigning cousins to child[3]:
        // Three well-separated cousins:
          (node->child.at(3)->interaction).at(node->child.at(3)->N_interaction_allocated + 0) 
        = (node->neighbor.at(4))->child.at(1);

          (node->child.at(3)->interaction).at(node->child.at(3)->N_interaction_allocated + 1) 
        = (node->neighbor.at(4))->child.at(2);

          (node->child.at(3)->interaction).at(node->child.at(3)->N_interaction_allocated + 2) 
        = (node->neighbor.at(4))->child.at(3);

        // One neighbor:
        (node->child.at(3)->neighbor).at(5) = (node->neighbor.at(4))->child.at(0);

        node->child.at(3)->N_interaction_allocated += 3;
        node->child.at(3)->N_neighbor_allocated    += 1;
    }


    // Assigning children of neighbour 5
    {
        // Assigning cousins to child[0]:
        // Four well-separated cousins:
        FMM2DTree::assignFourWellSeparatedCousins(node, 5, 0);

        // Assigning cousins to child[1]:
        // Four well-separated cousins:
        FMM2DTree::assignFourWellSeparatedCousins(node, 5, 1);

        // Assigning cousins to child[2]:
        // Two well-separated cousins:
          (node->child.at(2)->interaction).at(node->child.at(2)->N_interaction_allocated + 0) 
        = (node->neighbor.at(5))->child.at(2);

          (node->child.at(2)->interaction).at(node->child.at(2)->N_interaction_allocated + 1) 
        = (node->neighbor.at(5))->child.at(3);

        // Two neighbors:
        (node->child.at(2)->neighbor).at(5) = (node->neighbor.at(5))->child.at(0);
        (node->child.at(2)->neighbor).at(4) = (node->neighbor.at(5))->child.at(1);

        // Assigning cousins to child[3]:
        // Two well-separated cousins:
          (node->child.at(3)->interaction).at(node->child.at(3)->N_interaction_allocated + 0) 
        = (node->neighbor.at(5))->child.at(2);

          (node->child.at(3)->interaction).at(node->child.at(3)->N_interaction_allocated + 1) 
        = (node->neighbor.at(5))->child.at(3);

        // Two neighbors:
        (node->child.at(3)->neighbor).at(5) = (node->neighbor.at(5))->child.at(1);
        (node->child.at(3)->neighbor).at(6) = (node->neighbor.at(5))->child.at(0);

        node->child.at(3)->N_interaction_allocated += 2;
        node->child.at(3)->N_neighbor_allocated    += 2;
    }
    
    // Assigning children of neighbour 6
    {
        // Assigning cousins to child[0]:
        // Four well-separated cousins:
        FMM2DTree::assignFourWellSeparatedCousins(node, 6, 0);

        // Assigning cousins to child[1]:
        // Four well-separated cousins:
        FMM2DTree::assignFourWellSeparatedCousins(node, 6, 1);

        // Assigning cousins to child[2]:
        // Three well-separated cousins:
          (node->child.at(2)->interaction).at(node->child.at(2)->N_interaction_allocated + 0) 
        = (node->neighbor.at(6))->child.at(0);

          (node->child.at(2)->interaction).at(node->child.at(2)->N_interaction_allocated + 1) 
        = (node->neighbor.at(6))->child.at(2);

          (node->child.at(2)->interaction).at(node->child.at(2)->N_interaction_allocated + 2) 
        = (node->neighbor.at(6))->child.at(3);

        // One neighbor:
        (node->child.at(2)->neighbor).at(6) = (node->neighbor.at(6))->child.at(1);

        node->child.at(2)->N_interaction_allocated += 1;
        node->child.at(2)->N_neighbor_allocated    += 3;

        // Assigning cousins to child[3]:
        // Four well-separated cousins:
        FMM2DTree::assignFourWellSeparatedCousins(node, 6, 3);
    }

    // Assigning children of neighbour 7
    {
        // Assigning cousins to child[0]:
        // Two well-separated cousins:
          (node->child.at(0)->interaction).at(node->child.at(0)->N_interaction_allocated + 0) 
        = (node->neighbor.at(7))->child.at(0);

          (node->child.at(0)->interaction).at(node->child.at(0)->N_interaction_allocated + 1) 
        = (node->neighbor.at(7))->child.at(2);

        // One neighbor:
        (node->child.at(0)->neighbor).at(6) = (node->neighbor.at(7))->child.at(3);
        (node->child.at(0)->neighbor).at(7) = (node->neighbor.at(7))->child.at(1);

        node->child.at(0)->N_interaction_allocated += 2;
        node->child.at(0)->N_neighbor_allocated    += 2;

        // Assigning cousins to child[1]:
        // Four well-separated cousins:
        FMM2DTree::assignFourWellSeparatedCousins(node, 7, 1);

        // Assigning cousins to child[2]:
        // Two well-separated cousins:
          (node->child.at(2)->interaction).at(node->child.at(2)->N_interaction_allocated + 0) 
        = (node->neighbor.at(7))->child.at(0);

          (node->child.at(2)->interaction).at(node->child.at(2)->N_interaction_allocated + 1) 
        = (node->neighbor.at(7))->child.at(2);

        // Two neighbors:
        (node->child.at(2)->neighbor).at(0) = (node->neighbor.at(7))->child.at(1);
        (node->child.at(2)->neighbor).at(7) = (node->neighbor.at(7))->child.at(3);

        node->child.at(0)->N_interaction_allocated += 2;
        node->child.at(0)->N_neighbor_allocated    += 2;

        // Assigning cousins to child[3]:
        // Four well-separated cousins:
        FMM2DTree::assignFourWellSeparatedCousins(node, 7, 3);
    }
}
