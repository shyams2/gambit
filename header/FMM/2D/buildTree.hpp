#ifndef __buildTree_hpp__
#define __buildTree_hpp__

#include "FMM2DBox.hpp"
#include "FMM2DTree.hpp"
#include "utils/determineCenterAndRadius.hpp"
#include <math.h>

// Assigns the relations for the children all boxes in the tree
void FMM2DTree::assignTreeRelations() 
{
    for(int N_level = 0; N_level < this->max_levels; N_level++) 
    {
        #pragma omp parallel for
        for(int N_box = 0; k < pow(4, N_level); N_box++) 
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
void assignChild0Relations(int N_level, int N_box) 
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
void assignChild1Relations(int N_level, int N_box) 
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
void assignChild2Relations(int N_level, int N_box)
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
void assignChild3Relations(int N_level, int N_box)
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

#endif
