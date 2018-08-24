// In this file, we test that the function converts to the Eigen Matrix in
// the format the we actually desire. For this we print the array and 
// correlate to see if it's in the form that is needed.

#include "../convertToEigen.hpp"

int main()
{
    // First testing with a 2D double:
    double a[6][5] = {{ 1,  2,  3,  4,  5}, 
                      { 6,  7,  8,  9, 10},
                      {11, 12, 13, 14, 15},
                      {16, 17, 18, 19, 20},
                      {21, 22, 23, 24, 25},
                      {26, 27, 28, 29, 30},
                     };

    Eigen::MatrixXd a_E;
    convertToEigen(a, a_E); 
    std::cout << a_E << std::endl << std::endl << std::endl;

    // Testing with single dimensional array:
    double b[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
    Eigen::MatrixXd b_E;
    convertToEigen(b, 3, 2, b_E); 
    std::cout << b_E << std::endl << std::endl << std::endl;

    // Now with AF:
    af::array a_af = af::randu(6, 5, f64);

    std::cout << "Original Matrix Stored in af::array" << std::endl;
    af_print(a_af);
    convertToEigen(a_af, a_E); 
    std::cout << "After converting to Eigen:" << std::endl;
    std::cout << a_E << std::endl << std::endl << std::endl;
    return 0;
}
