// This file is used to check the implementation of the FMM method implemented
#include "MatrixData.hpp"
#include "FMM/2D/FMM2DTree.hpp"
#include "../interactionKernels.hpp"

int main(int argc, char** argv)
{
    unsigned size = atoi(argv[1]);
    // Printing backend information:
    af::info();
    cout << endl;

    // Initializing the set of points:
    array p = 2 * af::randu(size, 2, f64) - 1; // r = 1, c = 0
    // For the moment, we shall take the source and target
    // particles to be the same:
    // Creating an instance of MatrixData:
    MatrixData M(laplaceSingleLayer, p, p);

    // We then will pass these set of points to the FMM2D tree class:
    FMM2DTree T(M, 12, "LEGENDRE");

    // Array for the charges:
    array charges = af::constant(1, size, f64);
    array &potential = T.getPotential(charges);
    
    // Direct evaluation:
    cout << "Performing a direct evaluation using MatVec multiplication:" << endl;
    array potential_direct = af::matmul(M.getArray(), charges);
    potential_direct.eval();

    cout << "=============ERROR IN THE CALCULATED POTENTIAL=============" << endl;
    cout << "Absolute Error: " << af::norm(potential_direct - potential) << endl;
    cout << "Relative Error: " << af::norm(potential_direct - potential) / af::norm(potential_direct) << endl;
    return 0;
}
