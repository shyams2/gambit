// In this file, we check that the result given by the implemented interpolation algorithm 
// is correct. We pass the matrix to the approximate to the getInterpolation function which 
// then provides an approximation A = U * S * V. We check that if the approximation given 
// reduces with increasing the rank of the approximation.

// This file tests that the algorithm is able to handle when passed with points in 2D.
// That is both x, and y coordinates are passed for both the source and the targets

#include "MatrixData.hpp"
#include "MatrixFactorizations/getInterpolation.hpp"
#include "utils/computeError.hpp"

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

    // Initializing the arrays U, S, V:
    array U, S, V;
    cout << "Using Chebyshev Nodes" << endl;
    MatrixFactorizer::getInterpolation(U, S, V, rank, "CHEBYSHEV", M);
    // Finding Z_approx:
    array Z_approx = af::matmul(U, S, V);
    printAllErrorNorms(Z_approx, M.getArray());

    cout << "Using Legendre Nodes" << endl;
    MatrixFactorizer::getInterpolation(U, S, V, rank, "LEGENDRE", M);
    // Finding Z_approx:
    Z_approx = af::matmul(U, S, V);
    printAllErrorNorms(Z_approx, M.getArray());

    cout << "Using Equispaced Nodes" << endl;
    MatrixFactorizer::getInterpolation(U, S, V, rank, "EQUISPACED", M);
    // Finding Z_approx:
    Z_approx = af::matmul(U, S, V);
    printAllErrorNorms(Z_approx, M.getArray());
}
