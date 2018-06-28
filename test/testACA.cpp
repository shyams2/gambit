// In this file, we check that the result given by the implemented ACA algorithm is correct
// We pass the matrix to the approximate to the ACA function which then provides an approximation
// Z' = U * V. We check that if the approximation given is accurate to the requested tolerance.

#include <iostream>
#include <arrayfire.h>
#include "../header/aca.hpp"

int main(int argc, char** argv)
{
    // Printing backend information:
    af::info();
    std::cout << std::endl;

    int size   = atoi(argv[1]);
    double tol = strtod(argv[2], NULL); // tolerance
    
    //Initializing the array Z which we need to approximate:
    // Location of points:
    af::array r = af::randu(size, f64);
    // Distance between every pair of points:
    r           = af::tile(r, 1, size) - af::tile(r.T(), size);
    af::array Z = 1 / (1 + af::pow(r, 2));

    // Estimated Rank of the matrix:
    long int k;

    // Initializing the arrays U, V:
    af::array U, V;
    aca(U, V, k, tol, Z);

    // Finding Z_approx:
    af::array Z_approx = af::matmul(U, V);

    // Printing rank of the approximation:
    std::cout << "Rank:" << k << std::endl;
    // Calculating the error:
    double abs_error = af::norm(Z_approx - Z);
    double rel_error = af::norm(Z_approx - Z) / af::norm(Z);

    std::cout << "Absolute Error:" << abs_error << std::endl;
    std::cout << "Relative Error:" << rel_error << std::endl;
}
