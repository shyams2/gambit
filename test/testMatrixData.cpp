#include "MatrixData.hpp"

// K(r) = 1 / (1 + r^2)
array interaction_kernel(array i, array j, array targets, array sources)
{   
    array r = targets(i) - sources(j);
    return(1 / (1 + r * r));
}

int main(int argc, char** argv)
{
    // Printing backend info:
    af::info();

    cout << endl << "Method 1 of creating an instance of MatrixData" << endl;
    cout << "Array created by user, and then passed to the MatrixData" << endl;
    // Creating the array to pass to the MatrixData class:
    array A = af::randu(100, 100, f64);
    MatrixData M1(A);
    cout << "Checking the data stored in the object is the same" << endl;
    cout << "||A - M1->A|| = " << af::norm(A - M1.getArray()) << endl;

    cout << endl << "Testing File-Writing" << endl;
    M1.dumpArray("data.h5");
    cout << "DONE!" << endl;

    cout << endl << "Testing File-Loading" << endl;
    M1.loadArray("data.h5");
    cout << "DONE!" << endl << endl;

    cout << "Checking that the data is the same even after writing and loading" << endl;
    cout << "||A - M1->A|| = " << af::norm(A - M1.getArray()) << endl;


    cout << endl << "Method 2 of creating an instance of MatrixData" << endl;
    cout << "Providing the dimensions and rank of the matrix data" << endl;
    MatrixData M2(50, 100, 20);
    cout << "Checking rank of the created matrix..." << endl;
    cout << "Rank of the generated matrix is = " << M2.estimateRank() << endl;

    cout << endl << "Method 3 of creating an instance of MatrixData" << endl;
    cout << "Providing the dimensions and range for singular values" << endl;
    MatrixData M3(50, 100, 1e7, 1e-7);
    cout << "Checking singular values of the created matrix..." << endl;
    af::array U, S, Vh;
    af::svd(U, S, Vh, M3.getArray());
    cout << "Maximum Singular Value = " << af::max<double>(S) << endl;
    cout << "We'll check the same using the power iteration based estimateSpectralNorm()" << endl;
    cout << "Maximum Singular Value = " << M3.estimateSpectralNorm() << endl;
    cout << "Minimum Singular Value = " << af::min<double>(S) << endl;
    cout << "Condition Number       = " << M3.getConditionNumber() << endl;

    cout << endl << "Method 4 of creating an instance of MatrixData" << endl;
    cout << "Providing the blueprint for the way the matrix is generated via a function" << endl;
    array target_coords = 1 + af::randu(100, f64);
    array source_coords = 3 + af::randu(200, f64);
    MatrixData M4(interaction_kernel, target_coords, source_coords);
    cout << "Getting the Matrix by using the getArray() method..." << endl;
    A = M4.getArray();
    cout << "Verifying that the expected matrix is returned by the method" << endl;
    // Allowing broadcasting:
    af::gforSet(true);
    af::array A_expected = interaction_kernel(af::seq(100), af::seq(200), 
                                              target_coords, source_coords.T()
                                             );
    af::gforSet(false);
    cout << "||A - A_expected|| = " << af::mean<double>(A - A_expected) << endl;

    // Checking the getRow and getColumn methods:
    // NOTE: 27 and 47 are arbitrary choices
    cout << "||M.getRow(27) - A_expected.row(27)|| = " << 
            af::mean<double>(M4.getRow(27) - A_expected.row(27)) << endl;
    cout << "||M.getColumn(47) - A_expected.row(47)|| = " << 
            af::mean<double>(M4.getColumn(47) - A_expected.col(47)) << endl;
}
