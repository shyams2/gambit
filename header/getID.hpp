#ifndef __getID_hpp__
#define __getID_hpp__

#include <arrayfire.h>
#include <assert.h>

/*
Interpolative decomposition of a matrix decomposes a matrix A as a product of an skeleton
and an interpolation matrix:
A = BP + E
where ||E|| < σ_{k + 1} when performing a rank k decomposition. This error is comparative
to the error that is obtained from SVD in addition to being cheaper in terms of computational
cost. The skeleton and interpolation matrices can further be expressed in terms of:

B = A * Π1
P = [I, T] * Π.T, wher Π = [Π1, Π2]
*/

/*
References for this implementation:
*/


/*!
 Obtains the ID as A = B * P, where all the matrices are in af::array format.
 */
void getID(af::array A, double rankOrTolerance, bool randomSampling = true) 
{
    // Checking that the input is always greater than zero:
    assert(rankOrTolerance >= 0);

    int m, n;
    m = A.dims(0);
    n = A.dims(1);


}

#endif
