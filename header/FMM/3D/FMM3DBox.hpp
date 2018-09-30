#ifndef __3DBox_hpp__
#define __3DBox_hpp__

#include <iostream>
#include <arrayfire.h>
#include <vector>

using af::array;

class FMM3DBox 
{
public:
    double c_x, c_y, c_z, r_x, r_y, r_z; // centers and radii for this box
    array inds_in_box;                   // indices of particles in this box

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
    // back face:
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
    // front face:
    // ================================
    // ||             |              ||
    // ||             |              ||
    // ||      7      |      6       ||
    // ||             |              ||
    // ||             |              ||
    // ||-------------|--------------||
    // ||             |              ||
    // ||             |              ||
    // ||      4      |      5       ||
    // ||             |              ||
    // ||             |              ||
    // ||             |              ||
    // ================================
    // At L2:
    // back face:
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
    // second from back:
    // ================================
    // ||      |      |       |      ||
    // ||  27  |  26  |   31  |  30  ||
    // ||------|------|-------|------||
    // ||      |      |       |      ||
    // ||  24  |  25  |   28  |  29  ||
    // ||-------------|--------------||
    // ||      |      |       |      ||
    // ||   19 |  18  |   23  |  22  ||
    // ||------|------|-------|------||
    // ||      |      |       |      ||
    // ||   16 |  17  |   20  |   21 ||
    // ||      |      |       |      ||
    // ================================
    // ...
    // ...
    // and so on ...

    int N_level;      // level in the tree
    int N_box;        // box number assigned to the current instance
    int parent;       // box number of the parent
    int children[8];  // box numbers of the children
    int neighbor[26]; // box numbers of the neighbors
    int inner[98];    // box numbers of the inner boxes
    int outer[218];   // box numbers of the outer boxes

    array nodes;      // nodes at the leaf level
    array node_charges, node_potentials;

    // Constructor for the class
    FMM3DBox();
    // Destructor for the class:
    ~FMM3DBox(){};
    // Prints details about the box:
    void printBoxDetails();
};

void FMM3DBox::printBoxDetails()
{
    cout << "Level Number       :" << this->N_level << endl;
    cout << "Box Number         :" << this->N_box << endl;
    cout << "Parent             :" << this->parent << endl;
    cout << "Radius             :" << this->r_x << ", " << this->r_y << endl;
    cout << "Center             :" << this->c_x << ", " << this->c_y << endl;
    cout << "Number of Particles:" << this->inds_in_box.elements() << endl;
    cout << "Children           :";
    for(unsigned k = 0; k < 8; k++)
        cout << this->children[k] << ", ";
    cout << endl;
    cout << "Neighbor           :";
    for(unsigned k = 0; k < 26; k++)
        cout << this->neighbor[k] << ", ";
    cout << endl;
    cout << "Inner              :";
    for(unsigned k = 0; k < 98; k++)
        cout << this->inner[k] << ", ";
    cout << endl;
    cout << "Outer              :";
    for(unsigned k = 0; k < 218; k++)
        cout << this->outer[k] << ", ";
    cout << endl;
}

FMM3DBox::FMM3DBox()
{
    // Setting all values to -1 at initialization:
    this->N_level = -1;
    this->N_box   = -1;
    this->parent  = -1;

    #pragma omp parallel for
    for(unsigned i = 0; i < 8; i++)
    {
        this->children[i] = -1;
    }

    #pragma omp parallel for
    for(unsigned i = 0; i < 26; i++)
    {
        this->neighbor[i] = -1;
    }
 
    #pragma omp parallel for
    for(unsigned i = 0; i < 98; i++)
    {
        this->inner[i] = -1;
    }

    #pragma omp parallel for
    for(unsigned i = 0; i < 218; i++)
    {
        this->outer[i] = -1;
    }
}

#endif
