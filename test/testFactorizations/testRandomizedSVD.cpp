#include "MatrixData.hpp"
#include "MatrixFactorizer.hpp"

int main(int argc, char **argv)
{
    int N    = atoi(argv[1]);
    int rank = atoi(argv[2]);

    // Creating an instance of MatrixData having size N X N and rank provided
    MatrixData M(N, N, rank);

    // We perform Randomized SVD uptil rank and check that the matrix
    // reconstructed from product of U, S and Vh is accurate enough:
    array U, S, Vh;
    MatrixFactorizer::getRandomizedSVD(U, S, Vh, rank, M);
    array reconstructed_array = af::matmul(U, af::diag(S, 0, false), Vh);

    cout << "Error in Reconstruction = " << af::norm(M.getArray() - reconstructed_array) << endl;
}
