// In this file, we check that the result given by the implemented Psuedo inverse
// is indeed correct. The preliminary check that has been done here is to ensure
// that the answer given for an invertible square matrix match results given by 
// the standard inverse

#include "MatrixData.hpp"
#include "MatrixFactorizer.hpp"

int main(int argc, char** argv)
{
    int size     = atoi(argv[1]);
    array A      = af::randu(size, size, f64);
    array A_inv  = af::inverse(A); //standard inverse

    // Creating an instance of MatrixData:
    MatrixData M(A);

    double error = af::norm(A_inv - MatrixFactorizer::getPseudoInverse(M));
    std::cout << "Error:" << error << std::endl;
}
