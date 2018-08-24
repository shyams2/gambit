#include <iostream>
#include <arrayfire.h>
#include "Eigen/Dense"

// Reshaping using column major format:
void resizeArray(const double *input_array, double **&reshaped_array, size_t n_rows, size_t n_cols)
{
    for(int j = 0;j < n_cols; j++)
    {
        for(int i = 0;i < n_rows; i++)
        {   
            reshaped_array[i][j] = input_array[i + j * n_rows];
        }
    }
}

// Conversion from Eigen::MatrixXd:
// Converting to 2D array:
void convertToDouble(const Eigen::MatrixXd &input_array, double **&ptr_to_output_array)
{
    const size_t n_rows = input_array.rows();
    const size_t n_cols = input_array.cols();

    ptr_to_output_array = new double*[input_array.rows()];
    for(int i; i < input_array.rows(); i++)
    {
        ptr_to_output_array[i] = new double[input_array.cols()];
    }

    resizeArray(input_array.data(), ptr_to_output_array, n_rows, n_cols);
}

// Converting to 1D array:
void convertToDouble(Eigen::MatrixXd &input_array, double *&ptr_to_output_array)
{
    ptr_to_output_array = input_array.data();
}

// Conversion from af:
// Converting to 2D array:
void convertToDouble(const af::array &input_array, double **&ptr_to_output_array)
{
    const size_t n_rows = input_array.dims(0);
    const size_t n_cols = input_array.dims(1);

    ptr_to_output_array = new double*[n_rows];
    for(int i; i < n_rows; i++)
    {
        ptr_to_output_array[i] = new double[n_cols];
    }

    // Getting raw_ptr of 1D data from af:
    double *data_1d = new double[input_array.elements()];
    input_array.host(data_1d);
    resizeArray(data_1d, ptr_to_output_array, n_rows, n_cols);
    delete[] data_1d;
}

// Converting to 1D array:
void convertToDouble(const af::array &input_array, double *&ptr_to_output_array)
{
    ptr_to_output_array = input_array.host<double>();
}
