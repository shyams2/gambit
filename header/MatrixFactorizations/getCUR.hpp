// Getting the CUR decomposition of a matrix A as
// A = C * U * R
// Reference:http://infolab.stanford.edu/~ullman/mmds/ch11.pdf

#ifndef __getCUR_hpp__
#define __getCUR_hpp__

#include "MatrixData.hpp"
#include "getPseudoInverse.hpp"

void sampleProbability(array& p_col, array& p_row, array A)
{
    // Taking square of the input matrix
    array A_squared = af::pow(A, 2);

    // Summing along the axes to get the probabilities of selection:
    p_col = af::sum(A_squared, 1);
    p_row = af::flat(af::sum(A_squared, 0));

    // Normalizing:
    // gforSet has been used for broadcasting:
    af::gforSet(true);
    p_col /= af::sum(p_col);
    p_row /= af::sum(p_row);
    af::gforSet(false);
}

void sample(std::vector<uint>& indices, int size, array prob)
{
    // Taking cumulative sum to find allow us to pick with the prob
    // distribution that has been given
    array p_cumsum = af::accum(prob);
    uint ind_picked;

    int i = 0; // loop counter

    while(i < size)
    {
        // Flag to check if ind_picked already exists in indices:
        bool flag = false;
        // Generating random numbers:
        double r    = (double) rand() / (double) RAND_MAX;
        // Picking the column within which the random number falls
        // For instance if the probabilites are given by:
        //   p        = {0.1, 0.2, 0.4, 0.05, 0.15, 0.05, 0.05}
        //=> p_cumsum = {0.1, 0.3, 0.7, 0.75, 0.9 , 0.95, 1.  }
        // If random number generated is r = 0.06, then chosen index is 0
        // If random number generated is r = 0.16, then chosen index is 1
        // If random number generated is r = 0.26, then chosen index is 2
        // If random number generated is r = 0.36, then chosen index is 3
        // ...
        ind_picked = af::where(p_cumsum >= r)(0).scalar<uint>();

        for(uint j = 0; j < indices.size(); j++)
        {
            if(ind_picked == indices.at(j))
            {
                flag = true;
            }
        }

        if(flag == false)
        {
            indices.push_back(ind_picked);
            i++;
        }
    }
}

namespace MatrixFactorizer
{
    // Obtains the CUR decomposition as A = C * U * R
    void getCUR(array& C, array& U, array& R, int rank, MatrixData M) 
    {
        array p_col, p_row;
        std::vector<uint> row_ind, col_ind;
        array A = M.getArray();

        sampleProbability(p_col, p_row, A);

        sample(row_ind, rank, p_row);
        sample(col_ind, rank, p_col);

        // Converting the STL vector to array:
        array row_ind_af(rank, row_ind.data());
        array col_ind_af(rank, col_ind.data());
        
        // Getting matrices C and R:
        C = A(af::span, col_ind_af); 
        R = A(row_ind_af, af::span);

        // Constructing the W matrix which is used to get the U matrix:
        array W = A(row_ind_af, col_ind_af);
        U       = getPseudoInverse(MatrixData(W));

        C.eval(); R.eval(); U.eval();
    }
}

#endif
