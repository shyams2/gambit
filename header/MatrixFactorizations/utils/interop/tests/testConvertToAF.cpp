// In this file, we test that the function converts to the AF array in
// the format the we actually desire. For this we print the array and 
// correlate to see if it's in the form that is needed.

#include "../convertToAF.hpp"

int main()
{
    // First testing with double:
    double a[6][5] = {{ 1,  2,  3,  4,  5}, 
                      { 6,  7,  8,  9, 10},
                      {11, 12, 13, 14, 15},
                      {16, 17, 18, 19, 20},
                      {21, 22, 23, 24, 25},
                      {26, 27, 28, 29, 30},
                     };

    af_print(convertToAF(a));

    // Now with Eigen:
    Eigen::MatrixXd b(6 , 5);
    b.setRandom();

    std::cout << b << std::endl;
    std::cout << "After converting to AF:" << std::endl;
    af_print(convertToAF(b));
    return 0;
}
