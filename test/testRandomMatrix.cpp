// This file checks whether the getRandomMatrix and getConditionNumber 
// functions are functioning properly

#include "../header/getRandomMatrix.hpp"
#include "../header/getConditionNumber.hpp"
#include <arrayfire.h>

int main(int argc, char** argv)
{
    af::info();
    int N              = atoi(argv[1]); //size of matrix
    double input_kappa = strtod(argv[2], NULL); //condition number
    
    af::array A;
    getRandomMatrix(N, N, input_kappa, A);
    double estimated_kappa = getConditionNumber(A);

    // Printing the results:
    std::cout << "Input Condition Number:"     << input_kappa << std::endl;
    std::cout << "Estimated Condition Number:" << estimated_kappa << std::endl;
    return 0;
}
