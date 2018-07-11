// Getting the SVD decomposition using a randomized algorithm
// Follows Algorithm 4.3 of http://arxiv.org/pdf/0909.4061

// Other References: https://research.fb.com/fast-randomized-svd/
//                   https://relate.cs.illinois.edu/course/cs598apk-f17/f/demos/upload/02-tools-for-low-rank/Randomized%20SVD.html

#ifndef __getRandomizedSVD_hpp__
#define __getRandomizedSVD_hpp__

#include "MatrixData.hpp"

// Computes an orthonormal matrix whose range approximates the range of A
void randomizedRangeFinder(array& Q, array A, int size, int n_iter) 
{
    Q = af::randn(A.dims(1), size, f64);

    // Arrays which are just place holders for values which aren't going to be used:
    array B, C;
    // Power iteration normalization(adds stability):
    for(int i = 0; i < n_iter; i++)
    {   
        // Applying iterations of QR and LU is more stable:
        // af::lu(Q, B, C, af::matmul(A, Q));
        // af::lu(Q, B, C, af::matmul(A.T(), Q));

        // af::qr(Q, B, C, af::matmul(A, Q));
        // af::qr(Q, B, C, af::matmul(A.T(), Q));

        Q = af::matmul(A, Q);
        Q = af::matmul(A.T(), Q);
    }

    array temp = af::matmul(A, Q);
    af::qr(Q, B, C, temp);
}

// TODO: Find why this wasn't faster than regular SVD. Rerun benchmarks?
namespace MatrixFactorizer
{
    void getRandomizedSVD(array& U, array& S, array& V, 
                          int n_components, MatrixData M, 
                          int n_oversamples = 5, int n_iter = 4
                         ) 
    {
        array Q;
        int n_random = n_components + n_oversamples;
        array A = M.getArray();
        randomizedRangeFinder(Q, A, n_random, n_iter);
        // Project to the (n_components + n_oversamples) dimensional space 
        // using the basis vectors
        array B = af::matmul(Q.T(), A);
        
        array U_hat;
        // Computing SVD of the thin matrix:
        af::svdInPlace(U_hat, S, V, B);
        U = af::matmul(Q, U_hat);
        
        // Excluding oversamples:
        U = U(af::span, af::seq(n_components));
        S = S(af::seq(n_components));
        V = V(af::seq(n_components), af::span);
        U.eval(); S.eval(); V.eval();
    }
}

#endif
