#include <iostream>
#include <vector>
#include <arrayfire.h>

/*!
 Obtains the SVD decomposition as A = U * S * Vh, where all the matrices are in af::array format.
 */
void getSVD(const af::array A, const double tolerance, af::array& U, af::array& S, af::array& Vh, int& rank) 
{
    af::svd(U, S, Vh, A);
    rank = 0;
    while ((rank < S.elements()) && (S(rank).scalar<double>() > tolerance))
    {
        rank++;
    }
    U  = U(af::span, af::seq(0, rank));
    S  = S(af::seq(0,rank));
    Vh = Vh(af::seq(0, rank), af::span);
}

// For now, we will add the segment of code used to test this function here itself:
int main()
{
    const int N[] = {32, 48, 64, 96, 128, 192, 256, 384, 512};
    double timing_data[sizeof(N) / sizeof(N[0])];
    int i, j; // loop counters
    af::array U, S, Vh;
    
    for(i = 0; i < sizeof(N) / sizeof(N[0]); i++)
    {
        af::array A = af::randu(N[i], N[i], f64);
        
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
        timing_data[i] = af::timer::stop();
    }

    for(i = 0; i < sizeof(N) / sizeof(N[0]); i++)
    {
        std::cout << (int) N[i] << " X " << (int) N[i] << " took " << timing_data[i] << " seconds" << std::endl;
    }

    return 0;
}
