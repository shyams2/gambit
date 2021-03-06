// This file is used to check the implementation of the FMM method implemented
#include "MatrixData.hpp"
#include "FMM/2D/FMM2DTreeUniform.hpp"
#include "../interactionKernels.hpp"

int main(int argc, char** argv)
{
    unsigned N_nodes  = atoi(argv[1]);
    unsigned N_levels = atoi(argv[2]);

    // Printing backend information:
    af::info();
    cout << endl;
    // Creating an instance of MatrixData:
    MatrixData M(gaussianKernel);
    // We then will pass this to the FMM2D tree class:
    FMM2DTree T(M, N_nodes, "CHEBYSHEV", N_levels, 1);

    af::timer tic = af::timer::start();
    T.getPotential(array());
    af::sync();
    double total_time = af::timer::stop(tic);
    cout << "Time taken by FMM:" << total_time << endl;

    tic = af::timer::start();
    // T.evaluateExactPotentials();
    af::sync();
    total_time = af::timer::stop(tic);
    cout << "Time taken by exact:" << total_time << endl << endl;

    // Checking for random boxes:
    int box_1 = pow(4, N_levels) * (double) rand() / RAND_MAX;
    int box_2 = pow(4, N_levels) * (double) rand() / RAND_MAX;
    int box_3 = pow(4, N_levels) * (double) rand() / RAND_MAX;
    int box_4 = pow(4, N_levels) * (double) rand() / RAND_MAX;
    int box_5 = pow(4, N_levels) * (double) rand() / RAND_MAX;
    
    T.checkPotentialsInBox(box_1);
    T.checkPotentialsInBox(box_2);
    T.checkPotentialsInBox(box_3);
    T.checkPotentialsInBox(box_4);
    T.checkPotentialsInBox(box_5);

    return 0;
}

