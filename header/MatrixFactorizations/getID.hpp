#ifndef __getID_hpp__
#define __getID_hpp__

#include "MatrixData.hpp"
#include "getPseudoInverse.hpp"
#include "getRRQR.hpp"

using af::matmul;

// Interpolative decomposition of a matrix decomposes a matrix A as a product of an skeleton
// and an interpolation matrix:
// A = BP + E
// where ||E|| < σ_{k + 1} when performing a rank k decomposition. This error is comparative
// to the error that is obtained from SVD in addition to being cheaper in terms of computational
// cost.
// Reference: Slide No. 44-48 of https://andreask.cs.illinois.edu/cs598apk-f17/notes.pdf

// Algorithm:
// A Π = Q R (rank revealing QR)
// A Π = Q (R_{11} R_{12})
// To get A = B P
// Set  B = Q R_{11}
// Then P = (Id inv(R_{11}) R_{12}) Π.T

// NOTE: A very inefficient ID has been implemented
// TODO: change this
namespace MatrixFactorizer
{
    void getID(array& B, array& P, int rank, MatrixData M) 
    {
        array Q, R, Perm;
        MatrixFactorizer::getRRQR(Q, R, Perm, M);
        array R11 = R(af::span, af::seq(rank));
        array R12 = R(af::span, af::seq(rank, af::end));
        array Id  = af::identity(rank, rank, f64);
        // B = Q R_{11}
        B = matmul(Q, R11);
        // P = (Id inv(R_{11}) R_{12}) Π.T
        P = matmul(join(1, Id, matmul(MatrixFactorizer::getPseudoInverse(MatrixData(R11)),
                                      R12
                                     )
                       ),
                   Perm.T()
                  );

        B.eval(); P.eval();
    }
}

#endif
