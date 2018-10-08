// In this file, we check that the result given by the implemented interpolation algorithm 
// is correct. We pass the matrix to the approximate to the getInterpolation function which 
// then provides an approximation A = U * S * V. We check that if the approximation given 
// reduces with increasing the rank of the approximation.

// This file tests that the algorithm is able to handle when passed with points in 3D.
// That is x, y and z coordinates are passed for both the source and the targets

#include "MatrixData.hpp"
#include "MatrixFactorizations/getInterpolation.hpp"
#include "utils/computeError.hpp"
#include "../interactionKernels.hpp"

int main(int argc, char** argv)
{
    // Printing backend information:
    af::info();
    cout << endl;

    int size = atoi(argv[1]);
    int rank = atoi(argv[2]);

    // Initializing the array which we need to approximate:
    // Location of points in 2D:
    // (p1(:, 0) is x-coords);(p1(:, 1) is y-coords);(p1(:, 2) is z-coords) for p1
    array p1 = -0.5 - af::randu(size, 3, f64); // r = 0.5 c = -1
    array p2 =  0.5 + af::randu(size, 3, f64); // r = 0.5 c = 1
    // Creating an instance of MatrixData:
    MatrixData M(stokesSingleLayer, p1, p2);

    // Initializing the arrays U, S, V:
    array U, S, V;
    
    cout << "Using Chebyshev Nodes" << endl;
    MatrixFactorizer::getInterpolation(U, S, V, rank, "CHEBYSHEV", M);
    // Finding Z_approx:
    array Z_approx = af::matmul(af::tile(U, 1, 1, S.dims(2)), S, af::tile(V.T(), 1, 1, S.dims(2)));        
    printAllErrorNorms(Z_approx, M.getArray());

    cout << "Using Legendre Nodes" << endl;
    MatrixFactorizer::getInterpolation(U, S, V, rank, "LEGENDRE", M);
    // Finding Z_approx:
    Z_approx = af::matmul(af::tile(U, 1, 1, S.dims(2)), S, af::tile(V.T(), 1, 1, S.dims(2)));        
    printAllErrorNorms(Z_approx, M.getArray());

    cout << "Using Equispaced Nodes" << endl;
    MatrixFactorizer::getInterpolation(U, S, V, rank, "EQUISPACED", M);
    // Finding Z_approx:
    Z_approx = af::matmul(af::tile(U, 1, 1, S.dims(2)), S, af::tile(V.T(), 1, 1, S.dims(2)));        
    printAllErrorNorms(Z_approx, M.getArray());
}
