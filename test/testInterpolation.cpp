// In this file, we check that the result given by the implemented Chebyshev interpolation algorithm 
// is correct. We pass the matrix to the approximate to the getChebyshev function which then provides 
// an approximation A = U * S * V. We check that if the approximation given reduces with increasing the
// rank of the approximation.

#include "MatrixData.hpp"
#include "MatrixFactorizer.hpp"

// K(r) = 1 / (1 + r^2)
array interaction_kernel(array i, array j, array targets, array sources)
{   
    array r = targets(i) - sources(j);
    return(1 / (1 + r * r));
}

int main(int argc, char** argv)
{
    // Printing backend information:
    af::info();
    cout << endl;

    int size = atoi(argv[1]);
    int rank = atoi(argv[2]);

    //Initializing the array  which we need to approximate:
    // Location of points:
    array x1 =  1.5  - af::randu(size, f64); // r = 0.5 c = 1
    array x2 = -1.5  + af::randu(size, f64); // r = 0.5 c = -1

    // Creating an instance of MatrixData:
    MatrixData M(interaction_kernel, x1, x2);

    // Initializing the arrays U, S, V:
    array U, S, V;
    MatrixFactorizer::getInterpolation(U, S, V, rank, M);
    // Finding Z_approx:
    array approx = af::matmul(U, S, V);

    // Calculating the error:
    double abs_error = af::norm(approx - M.getArray());
    double rel_error = af::norm(approx - M.getArray()) / af::norm(M.getArray());

    cout << "Absolute Error:" << abs_error << endl;
    cout << "Relative Error:" << rel_error << endl;
}
