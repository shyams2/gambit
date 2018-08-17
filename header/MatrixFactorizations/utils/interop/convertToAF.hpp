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
af::array convertToAF(double (&input_array)[n_rows][n_cols])
{
    double *ptrToFlattenedArray; //pointer to array
    flattenMatrix(input_array, ptrToFlattenedArray);
    return af::array(n_rows, n_cols, ptrToFlattenedArray);
}

af::array convertToAF(Eigen::MatrixXd &input_array)
{
    return af::array(input_array.rows(), input_array.cols(), input_array.data());
}
