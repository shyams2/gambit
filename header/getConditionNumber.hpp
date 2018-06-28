// Contains the definition of function that returns the condition
// number κ for a given matrix A.

#ifndef __getConditionNumber_hpp__
#define __getConditionNumber_hpp__

#include <arrayfire.h>
#include "pinv.hpp"

double getConditionNumber(af::array A) 
{
    double kappa; // Condition number

    // Approach 2 given below seems to give better results
    // TODO: Look into further
    // Approach 1: computing κ = ||A|| ||A^{-1}||

    // Getting the dimensions of the matrix:
    // int m = A.dims(0);
    // int n = A.dims(1);

    // Declaring the array for the inverse of A:
    // af::array Ainv;

    // If square then the regular inverse is used:
    // if(m == n)
    // {
    //     Ainv = af::inverse(A);
    // }

    // When non-square, we will use the Moore-Penrose inverse:
    // else
    // {
    //     Ainv = pinv(A);
    // }

    // kappa = af::norm(A) * af::norm(Ainv);

    // Approach 2: Alternatively using SVD for taking the ratio of max / min singular values:
    af::array U, S, Vh;
    af::svd(U, S, Vh, A);

    kappa = (af::max(S) / af::min(S)).scalar<double>();
    return kappa;
}

#endif
