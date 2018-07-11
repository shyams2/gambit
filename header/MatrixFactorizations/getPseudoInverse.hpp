#ifndef __getPseudoInverse_hpp__
#define __getPseudoInverse_hpp__

#include "MatrixData.hpp"

namespace MatrixFactorizer
{
    array getPseudoInverse(MatrixData M)
    {
        array U, S, Vh;
        af::svd(U, S, Vh, M.getArray());
        // Initializing S_pinv:
        array S_pinv = af::constant(0, M.getNumRows(), M.getNumColumns(), f64);

        for (int i = 0; i < S.dims(0); i++) 
        {
            S_pinv(i, i) = 1. / S(i);
        }
        
        // Changing entries which were divided by zero = inf to zero:
        S_pinv(af::isInf(S_pinv)) = 0;
        array pinv = af::matmul(af::matmul(Vh.T(), S_pinv.T()), U.T());
        af::eval(pinv);
        return pinv;
    }
}

#endif
