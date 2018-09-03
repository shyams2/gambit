// In this file, we check that the result given by the implemented ACA algorithm is correct
// We pass the matrix to the approximate to the ACA function which then provides an approximation
// Z' = U * V. We check that if the approximation given is accurate to the requested tolerance.

#include "MatrixData.hpp"
#include "MatrixFactorizer.hpp"

int main(int argc, char** argv)
{
    // Printing backend information:
    af::info();
    cout << endl;

    int size   = atoi(argv[1]);
    double tol = strtod(argv[2], NULL); // tolerance
    
    //Initializing the array Z which we need to approximate:
    // Location of points:
    array r = af::randu(size, f64);
    // Distance between every pair of points:
    r       = af::tile(r, 1, size) - af::tile(r.T(), size);
    array Z = 1 / (1 + af::pow(r, 2));

    // Estimated Rank of the matrix:
    uint k;

    // Initializing the arrays U, V:
    array U, V;
    MatrixFactorizer::getACA(U, V, k, tol, MatrixData(Z));

    // Finding Z_approx:
    array Z_approx = af::matmul(U, V);

    // Printing rank of the approximation:
    cout << "Rank of Approximation:" << k << endl;
    // Calculating the error:
    double abs_error = af::norm(Z_approx - Z);
    double rel_error = af::norm(Z_approx - Z) / af::norm(Z);

    cout << "Absolute Error:" << abs_error << endl;
    cout << "Relative Error:" << rel_error << endl;
}
