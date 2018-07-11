#include "MatrixFactorizer.hpp"
#include "MatrixData.hpp"

int main(int argc, char **argv)
{
    int N    = atoi(argv[1]);
    int rank = atoi(argv[2]); // rank of approximation

    array A = af::randu(N, N, f64);
    array C, U, R;

    MatrixFactorizer::getCUR(C, U, R, rank, MatrixData(A));
    // Reconstucting using C, U and R:
    array A_rec = af::matmul(C, U, R);
    cout << "Error in Approximation:" << af::norm(A - A_rec) << endl;
    return 0;
}
