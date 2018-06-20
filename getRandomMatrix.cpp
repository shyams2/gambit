#include "getRandomMatrix.hpp"

int main()
{
    double kappa = 1e12;
    af::array A;
    getRandomMatrix(32, 32, kappa, A);
    af::array U, S, Vh;
    af::svd(U, S, Vh, A);
    // Checking that the singular values fall off rapidly
    af_print(S);
    return 0;
}
