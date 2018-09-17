// This file is used to check different aspects of the Quadtree implementation 
// that has been carried out to implement to the Fast Multipole Method.

#include "MatrixData.hpp"
#include "FMM/2D/FMM2DTree.hpp"

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

    // If r ~ 0:
    //    return 0
    // Else:
    //    return 0.5 * log(R^2)
    return(af::select(r_squared < 1e-14, 0, 0.5 * af::log(r_squared)));
}

int main(int argc, char** argv)
{
    unsigned size = atoi(argv[1]);
    // Printing backend information:
    af::info();
    cout << endl;

    // Initializing the set of points:
    // (p1(:, 0) is x-coords);(p1(:, 1) is y-coords) for p1
    array p1 = -0.5 - af::randu(size, 2, f64); // r = 0.5 c = -1
    array p2 =  0.5 + af::randu(size, 2, f64); // r = 0.5 c = 1
    // Creating an instance of MatrixData:
    MatrixData M(interaction_kernel, p1, p2);

    // Array for the charges:
    array charges = 2 * (af::randu(size, f64) - 1);

    // We then will pass these set of points to the FMM2D tree class:
    FMM2DTree T(M, 3, "CHEBYSHEV", charges);
}
