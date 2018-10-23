// This file is used to check the implementation of the FMM method implemented
#include "MatrixData.hpp"
#include "FMM/2D/FMM2DTree_adaptive.hpp"
#include "../interactionKernels.hpp"

int main(int argc, char** argv)
{
    unsigned N_nodes  = atoi(argv[1]);
    unsigned N_levels = atoi(argv[2]);

    // Printing backend information:
    af::info();
    cout << endl;

    // Initializing the set of points:
    array p = 2 * af::randn(size, 2, f64);
    // Creating an instance of MatrixData:
    MatrixData M(gaussianKernel, p, p);
    // We then will pass this to the FMM2D tree class:
    FMM2DTree T(M, N_nodes, "CHEBYSHEV", N_levels, 1);

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
