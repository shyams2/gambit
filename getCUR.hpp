#ifndef __getCUR_hpp__
#define __getCUR_hpp__

#include <arrayfire.h>

void sampleProbability(af::array& pCol, af::array& pRow, af::array A)
{
    af::array A2 = af::pow(A, 2);

    pCol = af::sum(A2, 1);
    pRow = af::flat(af::sum(A2, 0));

    // Normalizing:
    af::gforSet(true);
    pCol /= af::sum(pCol);
    pRow /= af::sum(pRow);
    af::gforSet(false);
}

void sample(af::array& indices, int size, af::array prob)
{
    af::array probRows = af::accum(prob);           
    af::array tempInd  = af::constant(0, size, s32);

    for(int i = 0; i < size; i++)
    {
        double v        = (double) rand() / (double) RAND_MAX;
        af::array tempI = af::where(probRows >= v);
        if(tempI.elements() == 0)
        {
            tempInd(i) = probRows.elements();
        }
        else
        {
            tempInd(i) = tempI(0);
        }
    }        
    indices = af::sort(tempInd);
}

/*
The following function returns the Moore Penrose Psuedo inverse
*/
af::array pinv(af::array A)
{
    af::array U, S, Vh;
    af::svd(U, S, Vh, A);
    af::array S_pinv = af::constant(0, A.dims(0), A.dims(1), f64);
    
    for (int i = 0; i < S.dims(0); i++) 
    {
        S_pinv(i, i) = 1. / S(i);
    }
    
    S_pinv(af::isInf(S_pinv)) = 0;
    af::array pinv = af::matmul(af::matmul(Vh.T(), S_pinv.T()), U.T());
    af::eval(pinv);
    return pinv;
}

/*!
 Obtains the CUR decomposition as A = C * U * R, where all the matrices are in af::array format.
 */
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
