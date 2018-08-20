#include <iostream>
#include <arrayfire.h>
#include "Eigen/Dense"

// Flattened using column major format:
template <size_t n_rows, size_t n_cols>
void flattenMatrix(double (&input_array)[n_rows][n_cols], double *&flattenedArray)
{
    flattenedArray = new double[n_rows * n_cols];

    for(int j = 0;j < n_cols; j++)
    {
        for(int i = 0;i < n_rows; i++)
        {   
            flattenedArray[i + j * n_rows] = input_array[i][j];
        }
    }
}

template <size_t n_rows, size_t n_cols>
void convertToAF(const double (&input_array)[n_rows][n_cols], af::array &output_array)
{
    double *ptrToFlattenedArray; //pointer to array
    flattenMatrix(input_array, ptrToFlattenedArray);
    output_array = af::array(n_rows, n_cols, ptrToFlattenedArray);
}

void convertToAF(Eigen::MatrixXd &input_array, af::array &output_array)
{
    output_array = af::array(input_array.rows(), input_array.cols(), input_array.data());
}
