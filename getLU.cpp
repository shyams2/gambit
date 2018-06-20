#include <iostream>
#include <vector>
#include <arrayfire.h>

int main()
{

    const int N[] = {32, 48, 64, 96, 128, 192, 256, 384, 512};
    double timing_data[sizeof(N) / sizeof(N[0])];
    int i, j; // loop counters
    af::array L, U, P;
    
    for(i = 0; i < sizeof(N) / sizeof(N[0]); i++)
    {
        af::array A = af::randu(N[i], N[i], f64);
        
        // Kernel Compilation:
        af::lu(L, U, P, A);
        af::eval(L);
        af::eval(U);
        af::eval(P);
        af::sync();

        af::timer::start();
        for(j = 0; j < 100; j++)
        {
            af::qr(L, U, P, A);
            af::eval(L);
            af::eval(U);
            af::eval(P);
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
