// In this file, we test that the function converts to an array of double
// given an ArrayFire array or an Eigen Matrix.
// NOTE: This function also retains the 2D ordering:

#include "../convertToDouble.hpp"
// #include <iostream>
// #include "Eigen/Dense"

void printElementsOfArray(double *x, int n)
{
    for(int i = 0; i < n; i++)
    {
        std::cout << x[i] << std::endl;
    }
}

int main()
{
    // First testing with Eigen:
    Eigen::MatrixXd a(6 , 5);
    a.setRandom();

    std::cout << a << std::endl << std::endl;

    // // printElementsOfArray(a.data(), 30);
    // double **b = new double*[6];

    // for(int i = 0; i < 6; i++)
    // {   
    //     b[i] = new double[5];
    // }

    // for(int j = 0;j < 5; j++)
    // {
    //     for(int i = 0;i < 6; i++)
    //     {   
    //          b[i][j] = a.data()[i + j * 6];
    //     }
    // }

    // for(int i = 0; i < 6; i++)
    // {
    //     for(int j = 0; j < 5; j++)
    //     {
    //         std::cout << b[i][j] << "  ";
    //     }
    //     std::cout << std::endl;
    // }

    //     b[i][j] = a(i, j);
    // }
    double **b = convertToDouble(a);
    return 0;
}
