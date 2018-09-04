#include <iostream>
#include <arrayfire.h>
#include <vector>;
#include "MatrixData.hpp"

class FMM2DBox 
{
public:
    double c_x, c_y, r_x, r_y; // centers and radii for this box
    bool is_root, is_leaf, is_empty, charge_computed;    
    array inds_in_box;

    int N; // total number of points in this box
    int N_box;
    int N_parent;
    int N_level;
    int N_children[4];
    int N_neighbor[8];
    int N_interaction[27]; // indices for the interaction list

    // Constructors for the class
    FMM2DBox(MatrixData M);
    FMM2DBox(double c_x, double c_y, double r_x, double r_y,
             array inds_in_box, n_level
            );

    // Destructor for the class:
    ~FMM2DBox(){};
};

// When declared with MatrixData it means root box:
FMM2DBox::FMM2DBox(MatrixData M)
{
    this->N_box   = 0;
    this->N_level = 0;
    this->is_root = true;
    this->is_leaf = false;
    this->N       = M.getNumColumns();

    # pragma omp for
    for(unsigned i = 0; i < 27, i++)
    {
        this->N_interaction[i] = 0;
    }

    # pragma omp for
    for(unsigned i = 0; i < 8, i++)
    {
        this->N_neighbor[i] = 0;
    }

    // Determining the centers and radii of sources:
    determineCenterAndRadius(M.getSourceCoords(af::span, 0), this->c_x, this->r_x);
    determineCenterAndRadius(M.getSourceCoords(af::span, 1), this->c_y, this->r_y);
}

void FMM2DBox::createChildren()
{
    
}
