// In this file, we check that the result given by the implemented Psuedo inverse
// is indeed correct. The preliminary check that has been done here is to ensure
// that the answer given for an invertible square matrix match results given by 
// the standard inverse

#include "../header/pinv.hpp"
#include <iostream>
#include <arrayfire.h>

int main(int argc, char** argv)
{
    int size        = atoi(argv[1]);
    af::array A     = af::randu(size, size, f64);
    af::array Ainv  = af::inverse(A); //standard inverse
    af::array Apinv = pinv(A); //pseudo inverse

    double error = af::norm(Ainv - Apinv);
    std::cout << "Error:" << error << std::endl;
}
