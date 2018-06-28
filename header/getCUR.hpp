// Getting the CUR decomposition of a matrix A as
// A = C * U * R
// Reference:http://infolab.stanford.edu/~ullman/mmds/ch11.pdf

#ifndef __getCUR_hpp__
#define __getCUR_hpp__

#include <arrayfire.h>

void sampleProbability(af::array& p_col, af::array& p_row, af::array A)
{
    // Taking square of the input matrix
    af::array A_squared = af::pow(A, 2);

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

void sample(af::array& indices, int size, af::array prob)
{
    // Taking cumulative sum to find allow us to pick with the prob
    // distribution that has been given
    af::array p_cumsum = af::accum(prob);           
    af::array temp_ind = af::constant(0, size, s32);

    for(int i = 0; i < size; i++)
    {
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
        temp_ind(i) = af::where(p_cumsum >= r)(0);
    }

    indices = af::sort(temp_ind);
}

// Obtains the CUR decomposition as A = C * U * R, where all the matrices are in af::array format.
void getCUR(af::array& C, af::array& U, af::array& R, int rank, af::array A) 
{
    af::array p_col, p_row, row_ind, col_ind;
    sampleProbability(p_col, p_row, A);

    sample(row_ind, rank, p_row);
    sample(col_ind, rank, p_col);

    // Dividing the chosen columns and row by  sqrt(r * prob)
    // Using gforSet for broadcasting:
    af::gforSet(true);
    C = A(af::span, col_ind) / af::sqrt(rank * af::moddims(p_col(col_ind), af::dim4(1, col_ind.elements())));        
    R = A(row_ind, af::span) / af::sqrt(rank * p_row(row_ind));
    af::gforSet(false);

    // Constructing the W matrix which is used to get the U matrix:
    af::array W = A(row_ind, col_ind);

    // Taking SVD of the W matrix:
    // W = W1 * Σ *  W2
    af::array W_singular1, sigma, W_singular2; // Declaring arrays for the singular vectors and values
    af::svdInPlace(W_singular1, sigma, W_singular2, W); // Performing inPlaceSVD since we don't need W after

    // Taking Moore-Penrose inverse of the diagonal matrix:
    af::array sigma_prime = (1. / sigma);
    // Converting the inf caused by divide by zero to zero:
    sigma_prime(af::isInf(sigma_prime)) = 0;

    // Getting the middle matrix U of A = C * U * R by:
    // U = W2.T * (Σ')^2 * W1.T
    U = af::matmul(W_singular2.T(), af::diag(af::pow(sigma_prime, 2), 0, false), W_singular1.T());
    af::eval(C);
    af::eval(R);
    af::eval(U);
}

#endif
