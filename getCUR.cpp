#include "getRandomMatrix.hpp"
#include "getCUR.hpp"

int main()
{
    double kappa = 1e14;
    af::array A = af::randu(32, 32, f64);
    // getRandomMatrix(32, 32, kappa, A);
    af::array C, U, R;
    getCUR(C, U, R, 32, 32, A);
    af::array A_rec = af::matmul(C, U, R);
    std::cout << af::norm(A - A_rec);
    return 0;
}
