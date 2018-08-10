// In this file, we check that the result given by the implemented interpolation algorithm 
// is correct. We pass the matrix to the approximate to the getInterpolation function which 
// then provides an approximation A = U * S * V. We check that if the approximation given 
// reduces with increasing the rank of the approximation.

// In this file, we test the implementation when the points given are in 1D

#include "MatrixData.hpp"
#include "MatrixFactorizer.hpp"

struct Point
{
    array x;
};

// K(r) = 1 / (1 + r^2)
array interaction_kernel(array i, array j, Point targets, Point sources)
{   
    array r = targets.x(i) - sources.x(j);
    return(1 / (1 + r * r));
}

void compute_error(array Z_approx, array Z)
{
    // Calculating the error using L1-norm:
    double abs_error = af::sum<double>(af::abs(Z_approx - Z));
    double rel_error = af::sum<double>(af::abs(Z_approx - Z)) / af::sum<double>(af::abs(Z));
    cout << "Using L1-norm:" << endl;
    cout << "Absolute Error:" << abs_error << endl;
    cout << "Relative Error:" << rel_error << endl << endl;

    // Calculating the error using L2-norm:
    abs_error = af::norm(Z_approx - Z);
    rel_error = af::norm(Z_approx - Z) / af::norm(Z);
    cout << "Using L2-norm:" << endl;
    cout << "Absolute Error:" << abs_error << endl;
    cout << "Relative Error:" << rel_error << endl << endl;

    // Calculating the error using L∞-norm:
    abs_error = af::max<double>(Z_approx - Z);
    rel_error = af::max<double>(Z_approx - Z) / af::max<double>(Z);
    cout << "Using max-norm:" << endl;
    cout << "Absolute Error:" << abs_error << endl;
    cout << "Relative Error:" << rel_error << endl << endl;
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
    array x1 =  1.5  - af::randu(size, f64); // r = 0.5 c = 1
    array x2 = -1.5  + af::randu(size, f64); // r = 0.5 c = -1

    // Creating an instance of MatrixData:
    MatrixData M(interaction_kernel, x1, x2);

    // Initializing the arrays U, S, V:
    array U, S, V;
    cout << "Using Chebyshev Nodes" << endl;
    MatrixFactorizer::getInterpolation(U, S, V, rank, "CHEBYSHEV", M);
    // Finding Z_approx:
    array Z_approx = af::matmul(U, S, V);
    compute_error(Z_approx, M.getArray());

    cout << "Using Legendre Nodes" << endl;
    MatrixFactorizer::getInterpolation(U, S, V, rank, "LEGENDRE", M);
    // Finding Z_approx:
    Z_approx = af::matmul(U, S, V);
    compute_error(Z_approx, M.getArray());

    cout << "Using Equispaced Nodes" << endl;
    MatrixFactorizer::getInterpolation(U, S, V, rank, "EQUISPACED", M);
    // Finding Z_approx:
    Z_approx = af::matmul(U, S, V);
    compute_error(Z_approx, M.getArray());
}
