#include <iostream>
#include <iomanip>
#include <vector>

void printElementsOfArray(double *x, int n)
{
    std::cout << std::setprecision(4) << std::fixed;
    for(int i = 0; i < n; i++)
    {
        std::cout << x[i] << std::endl;
    }
}

void printElementsOfArray(double **x, size_t n_rows, size_t n_cols)
{
    std::cout << std::setprecision(4) << std::fixed;
    for(int i = 0; i < n_rows; i++)
    {
        for(int j = 0; j < n_cols; j++)
        {
            std::cout << x[i][j] << "  ";
        }
        std::cout << std::endl;
    }
}
