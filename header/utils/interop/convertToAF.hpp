#include <iostream>
#include <arrayfire.h>
#include "Eigen/Dense"

// Flattened using column major format:
template <const size_t n_rows, const size_t n_cols>
void flattenMatrix(const double (&input_array)[n_rows][n_cols], double *&flattened_array)
{
    flattened_array = new double[n_rows * n_cols];

    for(int j = 0;j < n_cols; j++)
    {
        for(int i = 0;i < n_rows; i++)
        {   
            flattened_array[i + j * n_rows] = input_array[i][j];
        }
    }
}

// For 2 dimensional arrays which are given as an input:
template <size_t n_rows, size_t n_cols>
void convertToAF(const double (&input_array)[n_rows][n_cols], af::array &output_array)
{
    double *ptr_to_flattened_array; //pointer to array
    flattenMatrix(input_array, ptr_to_flattened_array);
    output_array = af::array(n_rows, n_cols, ptr_to_flattened_array);
}

// When the data passed is a 1D array:
void convertToAF(const double *ptr_to_array, const size_t n_rows, const size_t n_cols, af::array &output_array)
{
    output_array = af::array(n_rows, n_cols, ptr_to_array);
}

// When converting from Eigen::MatrixXd
void convertToAF(const Eigen::MatrixXd &input_array, af::array &output_array)
{
    output_array = af::array(input_array.rows(), input_array.cols(), input_array.data());
}
