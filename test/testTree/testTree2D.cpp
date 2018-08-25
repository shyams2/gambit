// This file is used to check different aspects of the Quadtree implementation 
// that has been carried out to implement to the Fast Multipole Method.

#include "MatrixData.hpp"

// Log(R) kernel:
// K(r) = log(R):
array interaction_kernel(const array &i, const array &j, const array &targets, const array &sources)
{   
    array x_targets = targets(af::span, 0);
    array x_sources = sources(af::span, 0);

    array y_targets = targets(af::span, 1);
    array y_sources = sources(af::span, 1);

    array x_diff = x_targets(i) - (x_sources.T())(j);
    array y_diff = y_targets(i) - (y_sources.T())(j);
    
    array r_squared = x_diff * x_diff + y_diff * y_diff;
    return(0.5 * af::log(r_squared));
}

int main(int argc, char** argv)
{
    // Printing backend information:
    af::info();
    cout << endl;

    int size = atoi(argv[1]);
    int rank = atoi(argv[2]);

    // Initializing the array which we need to approximate:
    // Location of points in 2D:
    // (p1(:, 0) is x-coords);(p1(:, 1) is y-coords) for p1
    array p1 = -0.5 - af::randu(size, 2, f64); // r = 0.5 c = -1
    array p2 =  0.5 + af::randu(size, 2, f64); // r = 0.5 c = 1
    // Creating an instance of MatrixData:
    MatrixData M(interaction_kernel, p1, p2);
}
