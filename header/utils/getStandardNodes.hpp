#ifndef __getStandardNodes_hpp__
#define __getStandardNodes_hpp__

#include <arrayfire.h>
#include <cstdlib>
#include <string>
using af::array;

// Returns the roots on the N-th Legendre polynomial
// These values are in the standard interval of [-1, 1]
void getLegendreNodes(const int N, array &nodes)
{
    // Holds the coefficients of the polynomial:
    array poly = af::join(0, af::constant(0, N, f64), af::constant(1, 1, f64));

    // Getting the scaled companion matrix:
    // This has been obtained following 
    // https://docs.scipy.org/doc/numpy-1.12.0/reference/generated/numpy.polynomial.legendre.legcompanion.html
    array companion = af::constant(0, N, N, f64);
    array scale     = 1 / af::sqrt(2 * af::range(N).as(f64) + 1);

    // Assigning the values to the sub-diagonal and super-diagonal:
    companion(af::seq(1, af::end, N+1)) =   af::range(N).as(f64)(af::seq(1, N-1))
                                          * scale(af::seq(0, N-2)) 
                                          * scale(af::seq(1, N-1));

    companion(af::seq(N, af::end, N+1)) = companion(af::seq(1, af::end, N+1));
    // Getting the eigen values of the companion matrix to get roots
    // Using SVD since ArrayFire doesn't give Eigenvalues:
    // Since the matrix is symmetric, we can use the singular values 
    // as it's equivalent to the eigenvalues.
    array U, S;
    af::svd(U, S, U, companion); // U isn't needed anyways

    // Taking the negatives of the singular values and appending:
    nodes = af::join(0, -1 * S(af::seq(0, af::end, 2)), S(af::seq(af::end-1, 0, -2)));
    nodes.eval();
}

// Returns the roots on the N-th Chebyshev polynomial of the first kind
// These values are in the standard interval of [-1, 1]
void getChebyshevNodes(const int N, array &nodes)
{
    nodes = -af::cos(af::Pi * (2 * af::range(N).as(f64) + 1) / (2 * N));
    nodes.eval();
}

// Returns the equispaced nodes in the standard interval of [-1, 1]
void getEquispacedNodes(const int N, array &nodes)
{
    nodes = -1 + (0.5 + af::range(N).as(f64)) * 2 / double(N); 
    nodes.eval();
}

void getStandardNodes(unsigned N_nodes, std::string nodes_type, array &standard_nodes)
{
    // Obtaining the standard Chebyshev nodes of the first kind:
    if(nodes_type == "CHEBYSHEV")
    {
        getChebyshevNodes(N_nodes, standard_nodes);
    }

    // Obtaining the standard Legendre nodes:
    else if(nodes_type == "LEGENDRE")
    {
        getLegendreNodes(N_nodes, standard_nodes);
    }

    // Obtaining the standard equispaced nodes:
    else if(nodes_type == "EQUISPACED")
    {
        getEquispacedNodes(N_nodes, standard_nodes);
    }

    else
    {
        cout << "Invalid choice for interpolation type" << endl;
        cout << "Please use either CHEBYSHEV, LEGENDRE or EQUISPACED" << endl;
        exit(1);
    }
}

#endif
