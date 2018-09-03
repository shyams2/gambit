#include "MatrixData.hpp"
#include "MatrixFactorizer.hpp"

int main(int argc, char **argv)
{
    int N    = atoi(argv[1]);
    int rank = atoi(argv[2]);

    // Creating an instance of MatrixData having size N X N and rank provided
    MatrixData M(N, N, rank);

    // We perform interpolative decomposition for given rank and check the error of reconstructed array:
    array B, P;
    MatrixFactorizer::getID(B, P, rank, M);
    array reconstructed_array = af::matmul(B, P);

    cout << "Error in Reconstruction = " << af::norm(M.getArray() - reconstructed_array) << endl;
}
