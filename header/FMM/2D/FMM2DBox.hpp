#ifndef __2DBox_hpp__
#define __2DBox_hpp__

#include <iostream>
// #include <arrayfire.h>
#include <vector>

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
    int inner[16];    // box numbers of the inner boxes
    int outer[24];    // box numbers of the outer boxes

    array nodes;      // nodes at the leaf level
    array multipoles, locals;

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
    for(unsigned k = 0; k < 4; k++)
        cout << this->children[k] << ", ";
    cout << endl;
    cout << "Neighbor           :";
    for(unsigned k = 0; k < 8; k++)
        cout << this->neighbor[k] << ", ";
    cout << endl;
    cout << "Inner              :";
    for(unsigned k = 0; k < 16; k++)
        cout << this->inner[k] << ", ";
    cout << endl;
    cout << "Outer              :";
    for(unsigned k = 0; k < 24; k++)
        cout << this->outer[k] << ", ";
    cout << endl;
}

FMM2DBox::FMM2DBox()
{
    // Setting all values to -1 at initialization:
    this->N_level = -1;
    this->N_box   = -1;
    this->parent  = -1;

    #pragma omp parallel for
    for(unsigned i = 0; i < 24; i++)
    {
        this->children[(int)(i / 6)]  = -1;
        this->neighbor[(int)(i / 3)]  = -1;
        this->inner[(int)(2 * i / 3)] = -1;
        this->outer[i]                = -1;
    }
}

#endif
