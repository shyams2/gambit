// Returns a Random Matrix with the condition number κ that is passed to the function
// Original Source:
// https://ieee754wonderland.quora.com/Numerical-Stability-and-Orthogonalization

#ifndef __getRandomMatrix_hpp__
#define __getRandomMatrix_hpp__

#include <iostream>
#include <arrayfire.h>
#include <cmath>

void getRandomMatrix(const int m, const int n, double kappa, af::array& A) 
{
    af::array Q1, Q2, R1, R2, T1, T2;
    int k = std::min(m, n);

    af::qr(Q1, R1, T1, af::randn(m, n, f64));
    af::qr(Q2, R2, T2, af::randn(m, n, f64));
    
    af::array U = Q1(af::seq(m), af::seq(k));
    af::array V = Q2(af::seq(k), af::seq(n));

    double l = std::pow(kappa, ((double) 1 / (double) (k-1)));
    
    af::array S = af::diag(af::pow(l, af::seq(0, -(k-1), -1)), 0, false);

    A = af::matmul(U, S.as(U.type()), V);
}

#endif
