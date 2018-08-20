// In this file, we check that the result given by the implemented interpolation algorithm 
// is correct. We pass the matrix to the approximate to the getInterpolation function which 
// then provides an approximation A = U * S * V. We check that if the approximation given 
// reduces with increasing the rank of the approximation.

// This file tests that the algorithm is able to handle when passed with points in 2D.
// That is both x, and y coordinates are passed for both the source and the targets

#include "MatrixData.hpp"
#include "MatrixFactorizations/getInterpolation.hpp"
#include "utils/computeError.hpp"

// Inverse Quadric kernel:
// K(r) = 1 / (1 + r^2)
array interaction_kernel(const array &i, const array &j, const array &targets, const array &sources)
{   
    array x_targets = targets(af::span, 0);
    array x_sources = sources(af::span, 0);

    array y_targets = targets(af::span, 1);
    array y_sources = sources(af::span, 1);

    array x_diff = x_targets(i) - (x_sources.T())(j);
    array y_diff = y_targets(i) - (y_sources.T())(j);
    
    array r_squared = x_diff * x_diff + y_diff * y_diff;
    return(1 / (1 + r_squared));
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

    // array p1 = -0.5 - (0.5 + af::range(size).as(f64)) * 1 / size; //-0.5  - af::randu(size, 2, f64); // r = 0.5 c = -1
    // array p2 =  0.5 + (0.5 + af::range(size).as(f64)) * 1 / size; //, 0.5  + af::randu(size, 2, f64); // r = 0.5 c = 1

    // p1 = af::join(1, p1, p1);
    // p2 = af::join(1, p2, p2);
    
    p1 = af::sort(p1);
    p2 = af::sort(p2);

    af_print(p1);
    af_print(p2);

    // Creating an instance of MatrixData:
    MatrixData M(interaction_kernel, p1, p2);
    af_print(M.getArray());

    // Initializing the arrays U, S, V:
    array U, S, V;
    cout << "Using Chebyshev Nodes" << endl;
    MatrixFactorizer::getInterpolation(U, S, V, rank, "CHEBYSHEV", M);
    cout << "Printing shape of U, S and V" << endl;
    cout << "For U = " << U.dims(0) << "," << U.dims(1) << endl;
    cout << "For S = " << S.dims(0) << "," << S.dims(1) << endl;
    cout << "For V = " << V.dims(0) << "," << V.dims(1) << endl;
    // Finding Z_approx:
    array Z_approx = af::matmul(U, S, V);
    af_print(Z_approx);
    printAllErrorNorms(Z_approx, M.getArray());

    // cout << "Using Legendre Nodes" << endl;
    // MatrixFactorizer::getInterpolation(U, S, V, rank, "LEGENDRE", M);
    // // Finding Z_approx:
    // Z_approx = af::matmul(U, S, V);
    // printAllErrorNorms(Z_approx, M.getArray());

    // cout << "Using Equispaced Nodes" << endl;
    // MatrixFactorizer::getInterpolation(U, S, V, rank, "EQUISPACED", M);
    // // Finding Z_approx:
    // Z_approx = af::matmul(U, S, V);
    // printAllErrorNorms(Z_approx, M.getArray());
}
