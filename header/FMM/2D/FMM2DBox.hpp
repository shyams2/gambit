#ifndef __2DBox_hpp__
#define __2DBox_hpp__

#include <iostream>
#include <arrayfire.h>
#include <vector>

using af::array;

class FMM2DBox 
{
public:
    double c_x, c_y, r_x, r_y; // centers and radii for this box
    array inds_in_box;         // indices of particles in this box

    // NOTE:
    // Box numbers are provided w.r.t a single level
    // That is, the box numbers start again from zero at each level
    // So, at L0:
    // ================================
    // ||                            ||
    // ||                            ||
    // ||                            ||
    // ||                            ||
    // ||                            ||
    // ||             0              ||
    // ||                            ||
    // ||                            ||
    // ||                            ||
    // ||                            ||
    // ||                            ||
    // ||                            ||
    // ================================
    // At L1:
    // ================================
    // ||             |              ||
    // ||             |              ||
    // ||      2      |      3       ||
    // ||             |              ||
    // ||             |              ||
    // ||-------------|--------------||
    // ||             |              ||
    // ||             |              ||
    // ||      0      |      1       ||
    // ||             |              ||
    // ||             |              ||
    // ||             |              ||
    // ================================
    // At L2:
    // ================================
    // ||      |      |       |      ||
    // ||  10  |  11  |   14  |  15  ||
    // ||------|------|-------|------||
    // ||      |      |       |      ||
    // ||   8  |  9   |   12  |  13  ||
    // ||-------------|--------------||
    // ||      |      |       |      ||
    // ||   2  |   3  |    5  |   7  ||
    // ||------|------|-------|------||
    // ||      |      |       |      ||
    // ||   0  |  1   |    4  |  6   ||
    // ||      |      |       |      ||
    // ================================


    int N_level;      // level in the tree
    int N_box;        // box number assigned to the current instance
    int parent;       // box number of the parent
    int children[4];  // box numbers of the children
    int neighbor[8];  // box numbers of the neighbors

    // Used in adaptive, level restricted tree:
    bool is_assigned;
    bool is_leaf;
    
    // Members used in a level restricted tree:
    int neighbor_fine[12];
    int sep_neighbor_fine[20];
    int neighbor_coarse[12];
    
    // Only currently being used in uniform tree:
    int inner[16];    // box numbers of the inner boxes
    int outer[24];    // box numbers of the outer boxes

    array nodes;
    array node_charges, node_potentials;
    // This is only allocated when performing checks:
    array exact_potentials;
    // Constructor for the class
    FMM2DBox();
    // Destructor for the class:
    ~FMM2DBox(){};
    // Prints details about the box:
    void printBoxDetails();
};

void FMM2DBox::printBoxDetails()
{
    cout << "Level Number       :" << this->N_level << endl;
    cout << "Box Number         :" << this->N_box << endl;
    cout << "Parent             :" << this->parent << endl;
    cout << "Radius             :" << this->r_x << ", " << this->r_y << endl;
    cout << "Center             :" << this->c_x << ", " << this->c_y << endl;
    cout << "Number of Particles:" << this->inds_in_box.elements() << endl;
    
    cout << "Children           :";
    for(int i = 0; i < 4; i++)
        cout << this->children[i] << ", ";
    cout << endl;
    
    cout << "Neighbor           :";
    for(int i = 0; i < 8; i++)
        cout << this->neighbor[i] << ", ";
    cout << endl;
    
    cout << "Inner              :";
    for(int i = 0; i < 16; i++)
        cout << this->inner[i] << ", ";
    cout << endl;
    
    cout << "Outer              :";
    for(int i = 0; i < 24; i++)
        cout << this->outer[i] << ", ";
    cout << endl;
}

FMM2DBox::FMM2DBox()
{
    // Setting all values to -1 at initialization:
    this->N_level     = -1;
    this->N_box       = -1;
    this->parent      = -1;
    this->is_assigned = true;
    this->is_leaf     = false;

    #pragma omp parallel for
    for(int i = 0; i < 4; i++)
        this->children[i] = -1;

    #pragma omp parallel for
    for(int i = 0; i < 8; i++)
        this->neighbor[i] = -1;

    #pragma omp parallel for
    for(int i = 0; i < 16; i++)
        this->inner[i] = -1;

    #pragma omp parallel for
    for(int i = 0; i < 24; i++)
        this->outer[i] = -1;
}

#endif
