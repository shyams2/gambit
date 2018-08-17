#include <iostream>
#include <arrayfire.h>
#include "Eigen/Dense"

// Reshaping using column major format:
void resizeArray(double *input_array, double **reshaped_array, size_t n_rows, size_t n_cols)
{
    // reshaped_array = new double*[n_rows];
    // for(int i; i < n_rows; i++)
    // {
    //     reshaped_array[i] = new double[n_cols];
    // }

    for(int i = 0; i < 6; i++)
    {
        for(int j = 0; j < 5; j++)
        {
            std::cout << reshaped_array[i][j] << "  ";
        }
        std::cout << std::endl;
    }

    for(int j = 0;j < n_cols; j++)
    {
        for(int i = 0;i < n_rows; i++)
        {   
            std::cout << input_array[i + j * n_rows] << std::endl;
            reshaped_array[i][j] = 1.0; //input_array[i + j * n_rows];
        }
    }
}

double** convertToDouble(Eigen::MatrixXd &input_array)
{
    // std::cout << input_array << std::endl;
    std::cout << "DONE!" << std::endl;
    const size_t n_rows = input_array.rows();
    const size_t n_cols = input_array.cols();
    std::cout << "DONE!" << std::endl;

    // printElementsOfArray(a.data(), 30);
    double **b = new double*[6];
    for(int i = 0; i < 6; i++)
    {   
        b[i] = new double[5];
    }
    std::cout << "DONE!" << std::endl;


    // for(int j = 0;j < 5; j++)
    // {
    //     for(int i = 0;i < 6; i++)
    //     {   
    //          b[i][j] = input_array.data()[i + j * 6];
    //     }
    // }


    double **reshaped_array = new double*[input_array.rows()];
    for(int i; i < input_array.rows(); i++)
    {
        reshaped_array[i] = new double[input_array.cols()];
    }

    resizeArray(input_array.data(), reshaped_array, input_array.rows(), input_array.cols());
    return reshaped_array;
}
