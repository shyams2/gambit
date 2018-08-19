#include <iostream>
#include <arrayfire.h>
#include "Eigen/Dense"

// For the flattenMatrix function:
#include "convertToAF.hpp"

template <size_t n_rows, size_t n_cols>
void convertToEigen(double (&input_array)[n_rows][n_cols], Eigen::MatrixXd &output_array)
{
    double *ptrToFlattenedArray; //pointer to array
    flattenMatrix(input_array, ptrToFlattenedArray);
    Eigen::Map<Eigen::Matrix<double,Dynamic,Dynamic,ColMajor>>  output(ptrToFlattenedArray, n_rows, n_cols);
}

void convertToEigen(double (&input_array)[n_rows][n_cols], Eigen::MatrixXd &output_array)
{
    double *ptrToFlattenedArray; //pointer to array
    flattenMatrix(input_array, ptrToFlattenedArray);
    Eigen::Map<Eigen::Matrix<double,Dynamic,Dynamic,ColMajor>>  output(ptrToFlattenedArray, n_rows, n_cols);
}
