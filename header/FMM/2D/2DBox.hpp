#include <iostream>
#include <arrayfire.h>
#include <vector>;
#include "MatrixData.hpp"

class 2DBox 
{
private:
    double c_x, c_y, r_x, r_y; // centers and radii for this box
    double x_start, y_start, x_end, y_end; // domain of this box
    array inds_in_box;

public:
    // Constructor for the class
    2DBox(MatrixData M);
    // Destructor for the class:
    ~2DBox();
    long getNumberOfPoints();
};

long 2DBox::getNumberOfPoints()
{
    return this->inds_in_box.elements();
}

2DBox::2DBox(MatrixData M)
{
    // Determining the centers and radii of targets:
    determineCenterAndRadius(target_coords(af::span, 0), this->c_x, this->r_x);
    determineCenterAndRadius(target_coords(af::span, 1), c_y_targets, r_y_targets);

    // Determining the centers and radii of sources:
    determineCenterAndRadius(source_coords(af::span, 0), c_x_sources, r_x_sources);
    determineCenterAndRadius(source_coords(af::span, 1), c_y_sources, r_y_sources);
}
