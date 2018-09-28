// This file is used to check the implementation of the FMM method implemented
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
    // Printing backend information:
    af::info();
    cout << endl;

    // Creating an instance of MatrixData:
    MatrixData M(interaction_kernel);

    // We then will pass this to the FMM2D tree class:
    FMM2DTree T(M, 8, "CHEBYSHEV", 5, 1);
    T.checkPotentialInBox(0);

    // // Array for the charges:
    // array charges = 2 * (af::randn(size, f64) - 1);
    // array &potential = T.getPotential(charges);

    return 0;
}
