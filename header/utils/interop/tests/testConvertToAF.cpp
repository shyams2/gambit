// In this file, we test that the function converts to the AF array in
// the format the we actually desire. For this we print the array and 
// correlate to see if it's in the form that is needed.

#include "../convertToAF.hpp"

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

    af::array a_af;
    convertToAF(a, a_af); 
    af_print(a_af);

    // Testing with single dimensional array:
    double b[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
    af::array b_af;
    convertToAF(b, 3, 2, b_af); 
    af_print(b_af);

    // Now with Eigen:
    Eigen::MatrixXd c(6 , 5);
    c.setRandom();
    af::array c_af;

    std::cout << "Original Matrix Stored in Eigen::MatrixXd:" << std::endl;
    std::cout << c << std::endl << std::endl;
    convertToAF(c, c_af);
    std::cout << "After converting to AF:" << std::endl;
    af_print(c_af);
    return 0;
}
