// The function contained in this file returns the ACA 
// decomposition of the input matrix given
// Reference: https://ieeexplore.ieee.org/document/1580747/

#ifndef __getACA_hpp__
#define __getACA_hpp__

#include "MatrixData.hpp"

// ||Z^{i + 1}|| = ||Z^{i}|| + ||U||**2 * ||V||**2 + term
// Where term = 2 * Σ(from j = 1 to k - 1)((u_j.T * u_k) * (v_k * v_j.T))(NOTE the correction)
double term(array& U, array& V, uint rank)
{
    // For rank = 1:
    // ||Z^{2}|| = ||Z^{1}|| + ||U||**2 * ||V||**2
    if(rank == 1)
    {
        return 0;
    }

    else
    {
        double sum = 0;
        for (uint i = 0; i < (rank-1);  i++)
        {
            sum = sum + af::abs(U(0, i) * U(0, rank - 1) * V(rank - 1, 0) * V(i, 0)).scalar<double>();
        }
        return(2 * sum);
    }
}

// We write a separate function for max_index instead of using the functions
// provided by the library since we'd also need to ensure that the indices selected
// are unique(non-repetitive)

// For the following function:
// A   - This would either be a row / column of R the approximate error matrix.
// idx - Vector containing the indices
int max_index(const array& A, const std::vector<int>& idx)
{
    uint i, j; // used in the loops
    uint max_idx = 0;
    
    // If no suitable index if found, j = -1 is returned which breaks the iterations loop:
    j = -1;
    
    double max = 0;

    for(i = 0; i < A.elements(); i++)
    {
        // Skip checking of max if i exists in idx:
        bool flag = true;
        for(j = 0; j < idx.size(); j++)
        {
            if(idx.at(j) == i)
            {
                flag = false;
            }
        }

        if(flag && abs(A(i).scalar<double>()) >= max)
        {
            max     = abs(A(i).scalar<double>());
            max_idx = i;
        }
    }
    return max_idx;
}

namespace MatrixFactorizer
{
    void getACA(array& U, array& V, uint &k, double tol, MatrixData M)
    {
        // Initializing k = 0:
        k = 0;

        // Getting row_size and column_size:
        size_t row_size    = M.getNumRows();
        size_t column_size = M.getNumColumns();
        array Z            = M.getArray();

        // Initialization of approximate error matrix:
        array R = af::constant(0, row_size, column_size, f64);

        // Norm of the approximation for Z:
        double z = 0;

        // I vector stores the row indices of the pivots:
        // J vector stores the column indices of the pivots:
        std::vector<int> I, J;

        // Z will be decomposed into: Z = U * V:
        // If A is of shape m X n. Then, for a rank p approximation
        // U will be of shape m X p
        // V will be of shape p X n

        // Initializing U to be a column array:
        // Columns will be appended to this array:
        U = af::constant(0, column_size, f64);
        // Initializing V to be a row array:
        // Rows will be appended to this array:
        V = af::constant(0, 1, row_size, f64);

        // Variable used to check if we've reached expected tolerance
        // Initializing with a high value:
        double epsilon = 100;

        while(abs(epsilon) > tol && k < std::min(row_size, column_size))
        {
            cout << "Iteration Number:" << (k+1) << endl;
        
            // When k = 0, we perform the intialization:
            // When k = 0, then i is chosen arbitrarily to be 0
            if(k == 0)
            {
                I.push_back(k);
            }
            
            // Finding row index:i S.T R(I, J) = MAX(R(:, J))
            else
            {
                I.push_back(max_index(R.col(J.at(k-1)), I));
            }

            // Terminating due to lack of finding an index
            if(I.at(k) == -1)
            { 
                break;  
            }

            else
            {
                // R(I, :) = Z(I, :) - Σ(from l = 0 to k-1)(U(I, l) * V(l, :))
                R.row(I.at(k)) = Z.row(I.at(k));
                for(uint i = 0; i < k; i++)
                { 
                    R.row(I.at(k)) = R.row(I.at(k)) - af::matmul(U(I.at(k),i), V.row(i));
                }
            }
            
            // Finding column index:j S.T R(I, J) = MAX(R(I, :))
            J.push_back(max_index(R.row(I.at(k)), J));
            
            // Terminating due to lack of finding an index
            if(J.at(k) == -1)
            {
                break;
            }

            else
            {
                // When k = 0, we assign the values to the initialized row matrix:
                if(k == 0)
                {
                    V = R.row(I.at(k))/ R(I.at(k), J.at(k)).scalar<double>();
                }
                
                // Appending to the existing row matrix:
                else
                {
                    V = af::join(0, V, R.row(I.at(k))/ R(I.at(k), J.at(k)).scalar<double>());
                }

                // R(:, J) = Z(:, J) - Σ(from l = 0 to k-1)(V(l, J) * U(:, l))
                R.col(J.at(k)) = Z.col(J.at(k));
                for(uint i = 0; i < k; i++)
                {
                    R.col(J.at(k)) = R.col(J.at(k)) - V(i, J.at(k)).scalar<double>() * U.col(i);
                }
            }
            
            // When k = 0, we assign the values to the initialized column matrix:
            if(k == 0)
            {
                U = R.col(J.at(k));
            }

            // Appending to the existing column matrix:
            else
            {
                U = af::join(1, U, R.col(J.at(k)));
            }
            
            // Updating ||Z||^2:
            z = z + std::pow(af::norm(V.row(k)), 2) * std::pow(af::norm(U.col(k)), 2) + term(U, V, k+1);
            // ε = ||u_k|| * ||v_k|| / ||Z^k||
            epsilon = (af::norm(U.col(k)) * af::norm(V.row(k))) / std::sqrt(z);
            
            cout << "Epsilon = " << epsilon << endl << endl;
            k++;
        }
        U.eval(); V.eval();
    }
}

#endif
