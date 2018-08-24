// In this file, we test that the function converts to an array of double
// given an ArrayFire array or an Eigen Matrix.

#include "../convertToDouble.hpp"
#include "../../printElementsOfArray.hpp"

int main()
{
    // First testing with Eigen:
    size_t n_rows = 6;
    size_t n_cols = 5;

    Eigen::MatrixXd a(n_rows, n_cols);
    a.setRandom();
    std::cout << "Printing Matrix stored under Eigen::MatrixXd:" << std::endl;
    std::cout << a << std::endl << std::endl;

    // Converting to 2D double array:
    double **a_double_2d;
    convertToDouble(a, a_double_2d);
    // Printing the array:
    std::cout << std::endl;
    printElementsOfArray(a_double_2d, n_rows, n_cols);
    std::cout << std::endl;

    // Converting to 1D double array:
    double *a_double_1d;
    convertToDouble(a, a_double_1d);
    // Printing the array:
    std::cout << std::endl;
    printElementsOfArray(a_double_1d, n_rows * n_cols);
    std::cout << std::endl;

    // Now testing with arrayfire:
    af::array a_af = af::randu(n_rows, n_cols, f64);
    std::cout << "Printing Matrix stored under af::Array:" << std::endl;
    af_print(a_af);    

    // Converting to 2D double array:
    convertToDouble(a_af, a_double_2d);
    // Printing the array:
    std::cout << std::endl;
    printElementsOfArray(a_double_2d, n_rows, n_cols);
    std::cout << std::endl;

    // Converting to 1D double array:
    convertToDouble(a_af, a_double_1d);
    // Printing the array:
    std::cout << std::endl;
    printElementsOfArray(a_double_1d, n_rows * n_cols);
    std::cout << std::endl;

    return 0;
}
