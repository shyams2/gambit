// In this file, we check that the result given by the implemented interpolation algorithm 
// is correct. We pass the matrix to the approximate to the getInterpolation function which 
// then provides an approximation A = U * S * V. We check that if the approximation given 
// reduces with increasing the rank of the approximation.

// In this file, we test the implementation when the points given are in 1D

#include "MatrixData.hpp"
#include "MatrixFactorizations/getInterpolation.hpp"
#include "utils/computeError.hpp"
#include "utils/printElementsOfArray.hpp"

// Inverse Quadric kernel:
// K(r) = 1 / (1 + r^2)
array interaction_kernel(array i, array j, array targets, array sources)
{   
    array r = targets(i) - (sources.T())(j);
    return(1 / (1 + r * r));
}

int main(int argc, char** argv)
{
    // Printing backend information:
    af::info();
    cout << endl;

    int size = atoi(argv[1]);
    int rank = atoi(argv[2]);

    // Initializing the array which we need to approximate:
    // Location of points:
    af::array x1 =  1.5  - af::randu(size, f64); // r = 0.5 c = 1
    af::array x2 = -0.5  + af::randu(size, f64); // r = 0.5 c = 0

    // Creating an instance of MatrixData:
    MatrixData M(interaction_kernel, x1, x2);

    // Initializing the arrays U, S, V:
    array U, S, V;
    cout << "Using Chebyshev Nodes" << endl;
    MatrixFactorizer::getInterpolation(U, S, V, rank, "CHEBYSHEV", M);
    // Finding Z_approx:
    array Z_approx = af::matmul(U, S, V.T());
    printAllErrorNorms(Z_approx, M.getArray());

    cout << "Using Legendre Nodes" << endl;
    MatrixFactorizer::getInterpolation(U, S, V, rank, "LEGENDRE", M);
    // Finding Z_approx:
    Z_approx = af::matmul(U, S, V.T());
    printAllErrorNorms(Z_approx, M.getArray());

    cout << "Using Equispaced Nodes" << endl;
    MatrixFactorizer::getInterpolation(U, S, V, rank, "EQUISPACED", M);
    // Finding Z_approx:
    Z_approx = af::matmul(U, S, V.T());
    printAllErrorNorms(Z_approx, M.getArray());
}

// Alternate form to be worked upon:
// Inverse Quadric kernel:
// K(r) = 1 / (1 + r^2)
// double interaction_kernel(unsigned i, unsigned j, 
//                           const double *targets, const size_t n_targets, 
//                           const double *sources, const size_t n_sources,
//                           const size_t n_dim
//                          )
// {
//     // In 1D: ||r|| = |x2 - x1|
//     double r = targets[i] - sources[j];
//     return (1 / (1 + r * r));
// }

// int main(int argc, char** argv)
// {
//     // Printing backend information:
//     af::info();
//     cout << endl;

//     int size = atoi(argv[1]);
//     int rank = atoi(argv[2]);

//     // Initializing the array which we need to approximate:
//     // Location of points:
//     af::array x1 =  1.5  - af::randu(size, f64); // r = 0.5 c = 1
//     af::array x2 = -0.5  + af::randu(size, f64); // r = 0.5 c = 0

//     // Getting it in double:
//     double *x_targets = new double[size];
//     double *x_sources = new double[size];

//     x1.host(x_targets);
//     x2.host(x_sources);

//     // Creating an instance of MatrixData:
//     MatrixData M(interaction_kernel, x_targets, size, x_sources, size, 1);
// }
