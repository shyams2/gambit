#ifndef __kroneckerProduct_hpp__
#define __kroneckerProduct_hpp__

#include <iostream>
#include <arrayfire.h>
using af::array;

// Returns the Kronecker Product:
array kroneckerProduct(const array &A, const array &B)
{   
    array A_tiled = af::tile(af::flat(A), 1, B.elements());
    array B_tiled = af::tile(af::flat(B).T(), A.elements(), 1);
    array C       = af::moddims(A_tiled * B_tiled, A.dims(0), A.dims(1), B.dims(0), B.dims(1));
    C             = af::reorder(C, 0, 2, 1, 3);
    C             = af::moddims(C, A.dims(0) * B.dims(0), A.dims(1) * B.dims(1));

    C.eval();
    return C
}

#endif
