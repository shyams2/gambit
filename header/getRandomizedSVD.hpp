#ifndef __getRandomizedSVD_hpp__
#define __getRandomizedSVD_hpp__

#include <iostream>
#include <arrayfire.h>
#include <cmath>

// Computes an orthonormal matrix whose range approximates the range of A
void randomizedRangeFinder(af::array& Q, af::array A, int size, int nIter) 
{
    Q = af::randn(A.dims(1), size, f64);

    // Arrays which are just place holders for values which aren't going to be used:
    af::array B, C;
    // Power iteration normalization(adds stability):
    for(int i = 0; i < nIter; i++)
    {   
        // af::lu(Q, B, C, af::matmul(A, Q));
        // af::lu(Q, B, C, af::matmul(A.T(), Q));

        af::qr(Q, B, C, af::matmul(A, Q));
        af::qr(Q, B, C, af::matmul(A.T(), Q));

        // Q = af::matmul(A, Q);
        // Q = af::matmul(A.T(), Q);
    }

    af::array temp = af::matmul(A, Q);
    af::qr(Q, B, C, temp);
}

void randomizedSVD(af::array& U, af::array& S, af::array& V, af::array A, int nComponents, int overSamples, int nIter) 
{
    af::array Q;
    int nRandom = nComponents + overSamples;
    randomizedRangeFinder(Q, A, nRandom, nIter);
    // Project to the (nComponents + overSamples) dimensional space using the basis vectors
    af::array B = af::matmul(Q.T(), A);
    
    af::array Uhat;
    // Computing SVD of the thin matrix:
    af::svdInPlace(Uhat, S, V, B);
    U = af::matmul(Q, Uhat);
    
    // Excluding oversamples:
    U = U(af::span, af::seq(nComponents));
    S = S(af::seq(nComponents));
    V = V(af::seq(nComponents), af::span);
    af::eval(U);
    af::eval(S);
    af::eval(V);
}

#endif
