// The function contained in this file returns the Moore Penrose Psuedo inverse

#ifndef __pinv_hpp__
#define __pinv_hpp__

#include <arrayfire.h>

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

#endif
