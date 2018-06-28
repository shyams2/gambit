#ifndef __estimateRank_hpp__
#define __estimateRank_hpp__

#include <arrayfire.h>

/*
This header file contains the definition for the function that estimate matrix 
rank to a specified relative precision.
*/

/*
Returns the rank of the Matrix A
*/

/*
For now rank estimation has been carried out using SVD which is a slow algo. 
Will include faster methods in the future
*/

int estimateRank(af::array A, double tolerance) 
{
    int rank = 0;
    af::array U, S, Vh;
    af::svd(U, S, Vh, A);
    for(int i = 0; i < S.elements(); S++)
    {
        rank++;
        if(S(i).scalar<double>() < tolerance)
            break
    }
    return rank
}

#endif
