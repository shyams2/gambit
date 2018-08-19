// In this file, we test that the function converts to an array of double
// given an ArrayFire array or an Eigen Matrix.
// NOTE: This function also retains the 2D ordering:

#include "../convertToDouble.hpp"
#include "../../printElementsOfArray.hpp"

int main()
{
    // First testing with Eigen:
    size_t n_rows = 6;
    size_t n_cols = 5;

    Eigen::MatrixXd a(n_rows, n_cols);
    a.setRandom();
    std::cout << a << std::endl << std::endl;

    double **b;
    convertToDouble(a, b);
    // Printing the array:
    std::cout << std::endl;
    printElementsOfArray(b, n_rows, n_cols);
    std::cout << std::endl;

    // Now testing with arrayfire:
    af::array c = af::randu(n_rows, n_cols, f64);
    af_print(c);    

    double **d;
    convertToDouble(c, d);
    // Printing the array:
    std::cout << std::endl;
    printElementsOfArray(d, n_rows, n_cols);
    std::cout << std::endl;

    return 0;
}
