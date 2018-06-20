#include "getRandomMatrix.hpp"
#include "getRandomizedSVD.hpp"

int main()
{
    // Ratio of largest to smallest singular value:
    double kappa = 1e14;
    const int N[] = {32, 48, 64, 96, 128, 192, 256};
    double timing_data_regular_svd[sizeof(N) / sizeof(N[0])];
    double timing_data_randomized_svd[sizeof(N) / sizeof(N[0])];
    int i, j; // loop counters
    // Input Array used to test:
    af::array A;
    // Arrays returned by SVD:
    af::array U, S, Vh;
    
    for(i = 0; i < sizeof(N) / sizeof(N[0]); i++)
    {
        getRandomMatrix(N[i], N[i], kappa, A);
        
        // Kernel Compilation:
        af::svd(U, S, Vh, A);
        af::eval(U);
        af::eval(S);
        af::eval(Vh);
        af::sync();

        af::timer::start();
        for(j = 0; j < 100; j++)
        {
            af::svd(U, S, Vh, A);
            af::eval(U);
            af::eval(S);
            af::eval(Vh);
        }
        af::sync();
        timing_data_regular_svd[i] = af::timer::stop();
    }

    std::cout << "Regular SVD:" << std::endl;
    for(i = 0; i < sizeof(N) / sizeof(N[0]); i++)
    {
        std::cout << (int) N[i] << " X " << (int) N[i] << " took " << timing_data_regular_svd[i] << " seconds" << std::endl;
    }

    // Arrays returned by the randomized SVD:
    af::array Ua, Sa, Vha;
    for(i = 0; i < sizeof(N) / sizeof(N[0]); i++)
    {
        getRandomMatrix(N[i], N[i], kappa, A);

        // Kernel Compilation:
        randomizedSVD(Ua, Sa, Vha, A, 5, 5, 0);
        af::eval(Ua);
        af::eval(Sa);
        af::eval(Vha);
        af::sync();

        af::timer::start();
        for(j = 0; j < 100; j++)
        {
            randomizedSVD(Ua, Sa, Vha, A, 5, 5, 0);
            af::eval(Ua);
            af::eval(Sa);
            af::eval(Vha);
        }
        af::sync();
        timing_data_randomized_svd[i] = af::timer::stop();
    }

    std::cout << std::endl << "Randomized SVD:" << std::endl;
    for(i = 0; i < sizeof(N) / sizeof(N[0]); i++)
    {
        std::cout << (int) N[i] << " X " << (int) N[i] << " took " << timing_data_randomized_svd[i] << " seconds" << std::endl;
    }

    return 0;
}
