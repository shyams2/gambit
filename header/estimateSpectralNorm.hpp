#ifndef __estimateSpectralNorm_hpp__
#define __estimateSpectralNorm_hpp__

#include <arrayfire.h>




void getCUR(af::array& C, af::array& U, af::array& R, int rowRank, int colRank, af::array A) 
{
    af::array pCol, pRow, rowInd, colInd;
    sampleProbability(pCol, pRow, A);

    sample(rowInd, rowRank, pRow);
    sample(colInd, colRank, pCol);

    C = af::matmul(A(af::span, colInd), af::diag(af::constant(1, colInd.elements(), f64), 0, false));        
    R = af::matmul(af::diag(af::constant(1, colInd.elements(), f64), 0, false), A(rowInd, af::span));

    U = af::matmul(pinv(C), A, pinv(R));
    af::eval(C);
    af::eval(R);
    af::eval(U);
}

#endif
