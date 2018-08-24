#include <iostream>
#include <arrayfire.h>
#include "Eigen/Dense"
// For the flattenMatrix function:
#include "convertToAF.hpp"

using Eigen::Dynamic;

template <size_t n_rows, size_t n_cols>
void convertToEigen(const double (&input_array)[n_rows][n_cols], Eigen::MatrixXd &output_array)
{
    double *ptr_to_flattened_array; //pointer to array
    flattenMatrix(input_array, ptr_to_flattened_array);
    output_array = Eigen::Map<Eigen::Matrix<double,Dynamic,Dynamic>>(ptr_to_flattened_array, n_rows, n_cols);
}

// When the data passed is a 1D array:
void convertToEigen(double *ptr_to_array, const size_t n_rows, const size_t n_cols, Eigen::MatrixXd &output_array)
{
    output_array = Eigen::Map<Eigen::Matrix<double,Dynamic,Dynamic>>(ptr_to_array, n_rows, n_cols);
}

// When converting from af::Array():
void convertToEigen(const af::array &input_array, Eigen::MatrixXd &output_array)
{
    double *ptr_to_data = input_array.host<double>();
    size_t n_rows       = input_array.dims(0);
    size_t n_cols       = input_array.dims(1);

    output_array = Eigen::Map<Eigen::Matrix<double,Dynamic,Dynamic>>(ptr_to_data, n_rows, n_cols);
    af::freeHost(ptr_to_data);
}
