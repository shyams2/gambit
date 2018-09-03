#include <iostream>
#include <arrayfire.h>
#include <vector>;
#include "MatrixData.hpp"

class FMM2DBox 
{
private:
    double c_x, c_y, r_x, r_y; // centers and radii for this box
    bool is_root, is_leaf;    
    array inds_in_box;

    int box_number;
    int parent_number;
    int children_numbers[4];
    int neighbor_numbers[8];
    int inner_numbers[16];
    int outer_numbers[24];

public:
    // Constructor for the class
    FMM2DBox(MatrixData M);
    // Alternate constructor:
    FMM2DBox(double c_x, double c_y, double r_x, double r_y,
             array inds_in_box, n_level
            );
    // Destructor for the class:
    ~FMM2DBox(){};

    void createChildren();
    void getRadius(double &r_x, double &r_y);
    void getCenter(double &c_x, double &c_y);
};

// When declared with MatrixData it means root box:
FMM2DBox::FMM2DBox(MatrixData M)
{
    this->box_number = 0;
    this->n_level = 0;
    this->is_root = true;
    this->is_leaf = false;
    // Determining the centers and radii of sources:
    determineCenterAndRadius(M.getSourceCoords(af::span, 0), this->c_x, this->r_x);
    determineCenterAndRadius(M.getSourceCoords(af::span, 1), this->c_y, this->r_y);
}

void FMM2DBox::createChildren()
{
    
}