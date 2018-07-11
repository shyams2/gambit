#ifndef __getRRQR_hpp__
#define __getRRQR_hpp__

#include "MatrixData.hpp"
#include <Eigen/Dense>

// TODO: Need to implement this from scratch in AF itself for speed
namespace MatrixFactorizer
{
    void getRRQR(array& Q, array& R, array &P, MatrixData M) 
    {
        // Getting sizes:
        size_t row_size = M.getNumRows();
        size_t col_size = M.getNumColumns();

        double *temp = new double[row_size * col_size];
        M.getArray().host(temp);
        // Converting to Eigen to be able to do column-pivoted QR:
        Eigen::Map<Eigen::MatrixXd> A_eig(temp, row_size, col_size);

        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(A_eig);
        Eigen::MatrixXd Q_eig = qr.householderQ();
        Eigen::MatrixXd R_eig = qr.matrixQR().triangularView<Eigen::Upper>();
        Eigen::MatrixXd P_eig = qr.colsPermutation();

        // Deleting the temporary variable used for transfer:
        delete[] temp;

        Q = array(row_size, col_size, Q_eig.data());
        R = array(row_size, col_size, R_eig.data());
        P = array(row_size, col_size, P_eig.data());
    }
}

#endif
