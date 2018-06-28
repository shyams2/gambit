#include "../header/getCUR.hpp"
#include <arrayfire.h>

int main()
{
    af::array A = af::randu(10000, 10000, f64);
    af::array C, U, R;
    // Performing a rank 32 decomposition:
    getCUR(C, U, R, 32, A);
    // Reconstucting using C, U and R:
    af::array A_rec = af::matmul(C, U, R);
    std::cout << "Error in Approximation:" << af::norm(A - A_rec) << std::endl;
    return 0;
}
