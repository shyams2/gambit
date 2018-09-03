#include "MatrixData.hpp"
#include "MatrixFactorizer.hpp"
#include <assert.h>

int main(int argc, char **argv)
{
    int N = atoi(argv[1]);
    // Array to factorize:
    array A = af::randu(N, N, f64);
    MatrixData M(A);

    array Q, R, P;
    MatrixFactorizer::getRRQR(Q, R, P, M);

    cout << "Error in Factorization = " << 
            af::norm(af::matmul(M.getArray(), P) - af::matmul(Q, R)) << endl;

    // Ensuring that the magnitude of values on diagonal is dropping:
    array R_diag = af::abs(af::diag(R));
    assert(af::norm(R_diag - af::sort(R_diag, true))<1e-14);            
}
