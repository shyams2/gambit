// This function gets the matrix matrix products for different sizes 
// of arrays and outputs the timing data to file for post-processing.

#include <iostream>
#include <fstream>
#include <arrayfire.h>

int main()
{
    std::cout << "Device Details" << std::endl;
    af::info();
    std::cout << std::endl;

    const int N[] = {32, 64, 128, 256, 512, 1024, 2048, 4096, 9192};
    double timing_data[sizeof(N) / sizeof(N[0])];
    int i, j; // loop counters

    for(i = 0; i < sizeof(N) / sizeof(N[0]); i++)
    {
        std::cout << "Benchmarking for N = " << N[i] << std::endl;
        af::array A = af::randu(N[i], N[i], f64);
        af::array B = af::randu(N[i], N[i], f64);
        
        // Declaration of the C which will store A * B:
        af::array C;

        // Kernel Compilation:
        C = af::matmul(A, B);
        af::eval(C)
        af::sync();

        af::timer::start();
        // Performing over 10 runs and then taking average
        for(j = 0; j < 10; j++)
        {
            C = af::matmul(A, B);
            af::eval(C)
        }
        af::sync();
        timing_data[i] = af::timer::stop() / 10; // dividing by 10 to average over iterations
    }

    // Writing the data to file:
    std::ofstream file;
    file.open("data.txt");
    for(i = 0; i < sizeof(N) / sizeof(N[0]); i++)
    {
        file << (int) N[i] << " " << timing_data[i] << std::endl;
    }

    system("python plot_scaling.py MatMat");
    return 0;
}
