// In this file, we check that the result given by the implemented interpolation algorithm 
// is correct. We pass the matrix to the approximate to the getInterpolation function which 
// then provides an approximation A = U * S * V. We check that if the approximation given 
// reduces with increasing the rank of the approximation.

// This file tests that the algorithm is able to handle when passed with points in 2D.
// That is both x, and y coordinates are passed for both the source and the targets

#include "MatrixData.hpp"
#include "MatrixFactorizations/getInterpolation.hpp"
#include "utils/computeError.hpp"

// Inverse Quadric kernel:
// K(r) = 1 / (1 + r^2)
array interaction_kernel(const array &i, const array &j, const array &targets, const array &sources)
{   
    array x_targets = targets(af::span, 0);
    array x_sources = sources(af::span, 0);

    array y_targets = targets(af::span, 1);
    array y_sources = sources(af::span, 1);

    array x_diff = x_targets(i) - (x_sources.T())(j);
    array y_diff = y_targets(i) - (y_sources.T())(j);
    
    array r_squared = x_diff * x_diff + y_diff * y_diff;
    return(1 / (1 + r_squared));
}

int main(int argc, char** argv)
{
    // Printing backend information:
    af::info();
    cout << endl;

    int size = atoi(argv[1]);
    int rank = atoi(argv[2]);

    // Initializing the array which we need to approximate:
    // Location of points in 2D:
    // (p1(:, 0) is x-coords);(p1(:, 1) is y-coords) for p1
    // array p1 = -0.5 - af::randu(size, 2, f64); // r = 0.5 c = -1
    // array p2 =  0.5 + af::randu(size, 2, f64); // r = 0.5 c = 1

    // array p1_x = -0.5 - (0.5 + af::range(size).as(f64)) * 1 / size; //-0.5  - af::randu(size, 2, f64); // r = 0.5 c = -1
    // array p2_x =  1.5 + (0.5 + af::range(size).as(f64)) * 0.5 / size; //, 0.5  + af::randu(size, 2, f64); // r = 0.5 c = 1

    // array p1_y = -2.5 - (0.5 + af::range(size).as(f64)) * 0.5 / size; //-0.5  - af::randu(size, 2, f64); // r = 0.5 c = -1
    // array p2_y =  3.5 + (0.5 + af::range(size).as(f64)) * 1 / size; //, 0.5  + af::randu(size, 2, f64); // r = 0.5 c = 1

    // array p1 = af::join(1, p1_x, p1_y);
    // array p2 = af::join(1, p2_x, p2_y);

    // array std_x_p1 = (2 * af::randu(size, f64) - 1); 
    // array std_y_p1 = (2 * af::randu(size, f64) - 1);
    // array std_x_p2 = (2 * af::randu(size, f64) - 1);
    // array std_y_p2 = (2 * af::randu(size, f64) - 1);

    double std_x_p1d[] = { 0.2019070583020710, 
                          -0.4714255689668849,
                          -0.5478724624368818 
                         };

    double std_y_p1d[] = { 0.7228424431963103, 
                          -0.3557883453402335,
                          -0.1651224234375461 
                         };

    double std_x_p2d[] = { 0.2136248097130777, 
                           0.5227853548169850,
                          -0.3306246583502916
                         };

    double std_y_p2d[] = { 0.6825782935497307, 
                          -0.9390349387534234,
                           0.9165682775287667 
                         };

    array std_x_p1(3, 1, std_x_p1d);
    array std_y_p1(3, 1, std_y_p1d);

    array std_x_p2(3, 1, std_x_p2d);
    array std_y_p2(3, 1, std_y_p2d);

    af::print("std_x_p1", std_x_p1, 16);
    af::print("std_y_p1", std_y_p1, 16);
    af::print("std_x_p2", std_x_p2, 16);
    af::print("std_y_p2", std_y_p2, 16);

    double rx_p1, rx_p2, ry_p1, ry_p2;
    double cx_p1, cx_p2, cy_p1, cy_p2;

    rx_p1 = rx_p2 = ry_p1 = ry_p2 = 0.5;
    cx_p1 = cy_p1 = -1;
    cx_p2 = cy_p2 =  1;

    array p1_x = cx_p1 + rx_p1 * std_x_p1;
    array p2_x = cx_p2 + rx_p2 * std_x_p2;
    array p1_y = cy_p1 + ry_p1 * std_y_p1;
    array p2_y = cy_p2 + ry_p2 * std_y_p2;
 
    // array p1_x = -0.5 - (0.5 + af::range(size).as(f64)) * 1 / size; //-0.5  - af::randu(size, 2, f64); // r = 0.5 c = -1
    // array p2_x =  1.5 + (0.5 + af::range(size).as(f64)) * 0.5 / size; //, 0.5  + af::randu(size, 2, f64); // r = 0.5 c = 1

    // array p1_y = -2.5 - (0.5 + af::range(size).as(f64)) * 0.5 / size; //-0.5  - af::randu(size, 2, f64); // r = 0.5 c = -1
    // array p2_y =  3.5 + (0.5 + af::range(size).as(f64)) * 1 / size; //, 0.5  + af::randu(size, 2, f64); // r = 0.5 c = 1

    array p1 = af::join(1, p1_x, p1_y);
    array p2 = af::join(1, p2_x, p2_y);

    af::print("p1", p1, 16);
    af::print("p2", p2, 16);

    // p1 = af::sort(p1);
    // p2 = af::sort(p2);

    // double p1d[] = {-1.523423, -1.434653, -1.3124, -1.21241, -1.11241, -1.00142, -0.96356, -0.834633, -0.734663, -0.612523,
    //                 -1.523462, -1.43462346, -1.33465, -1.251255, -1.12346, -1.001245, -0.936236, -0.8256346, -0.72143124, -0.64444};

    // double p2d[] = {1.53452, 1.4234, 1.3346, 1.2346, 1.13267, 1.0675, 0.9124, 0.862346, 0.7643, 0.63253,
    //                 1.56346, 1.4643, 1.3133, 1.2364, 1.12346, 1.0214, 0.9632, 0.843643, 0.7735, 0.67544};

    // array p1(10, 2, p1d);
    // array p2(10, 2, p2d);

    // Creating an instance of MatrixData:
    MatrixData M(interaction_kernel, p1, p2);
    af_print(M.getArray());

    // Initializing the arrays U, S, V:
    array U, S, V;
    cout << "Using Chebyshev Nodes" << endl;
    MatrixFactorizer::getInterpolation(U, S, V, rank, "CHEBYSHEV", M);
    cout << "Printing shape of U, S and V" << endl;
    cout << "For U = " << U.dims(0) << "," << U.dims(1) << endl;
    cout << "For S = " << S.dims(0) << "," << S.dims(1) << endl;
    cout << "For V = " << V.dims(0) << "," << V.dims(1) << endl;
    
    af::print("U", U, 16);
    af::print("S", S, 16);
    af::print("V", V, 16);

    // Finding Z_approx:
    array Z_approx = af::matmul(U, S, V);
    printAllErrorNorms(Z_approx, M.getArray());

    // cout << "Using Legendre Nodes" << endl;
    // MatrixFactorizer::getInterpolation(U, S, V, rank, "LEGENDRE", M);
    // // Finding Z_approx:
    // Z_approx = af::matmul(U, S, V);
    // printAllErrorNorms(Z_approx, M.getArray());

    // cout << "Using Equispaced Nodes" << endl;
    // MatrixFactorizer::getInterpolation(U, S, V, rank, "EQUISPACED", M);
    // // Finding Z_approx:
    // Z_approx = af::matmul(U, S, V);
    // printAllErrorNorms(Z_approx, M.getArray());
}
