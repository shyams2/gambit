#include <iostream>
#include <iomanip>
#include <vector>

void printElementsOfArray(double *x, size_t n, unsigned precision = 4)
{
    std::cout << std::setprecision(precision) << std::fixed;
    for(unsigned i = 0; i < n; i++)
    {
        std::cout << x[i] << std::endl;
    }
}

void printElementsOfArray(double **x, size_t n_rows, size_t n_cols, unsigned precision = 4)
{
    std::cout << std::setprecision(precision) << std::fixed;
    for(unsigned i = 0; i < n_rows; i++)
    {
        for(unsigned j = 0; j < n_cols; j++)
        {
            std::cout << x[i][j] << "  ";
        }
        std::cout << std::endl;
    }
}
