#ifndef __getChebyshev_hpp__
#define __getChebyshev_hpp__

#include "MatrixData.hpp"

// Returns the roots on the N-th Chebyshev polynomial
// These values are in the standard interval of [-1, 1]
array getChebyshevNodes(int N)
{
    return -af::cos(af::Pi*(2 * af::range(N, f64) + 1) / (2 * N));
}

// Maps the points from a domain having center c1 and radius r1
// To a domain having center c2 and radius r2
// Basically mapping from [c1 - r1 / 2, c1 + r1 / 2] --> [c2 - r2 / 2, c2 + r2 / 2]
void scalePoints(double c1, double r1, array x1, double c2, double r2, array& x2)
{   
    x2 = c2 + r2 * (x1 - c1) /r1;
}

// Returns the Chebyshev polynomials evaluated at x:
// T is of shape (x.elements(), N)
// where T[i, j] is the j-th Chebyshev polynomial evaluated at the i-th point x[i]
void getChebyshevPolynomials(int N, array x, array& T)
{
    // Initializing T:
    T = af::constant(0, x.elements(), N, f64);
    // The n-th Chebyshev polynomial is given by T_n(x) = cos(n * arccos(x))
    for(int i = 0; i < N, i++)
    {
        T.col(i) = af::cos(i * af::acos(x));
    }
}

// Returns the interpolation operator / L2L
void getChebyshevL2L(array x, array x_cheb, array& L2L)
{
    // Getting the Chebyshev polymonials:
    // Rank is the number of points in x_cheb
    rank = x_cheb.elements();

    // Declaring T_x and T_cheb:
    array T_x, T_cheb;
    // Getting T_x and T_cheb:
    getChebyshevPolynomials(rank, x, T_x);
    getChebyshevPolynomials(rank, x_cheb, T_cheb);

    L2L = (2 * af::matmul(T_x, T_cheb.T()) - 1) / rank;
}

namespace MatrixFactorizer
{
    void getChebyshev(array& U, array& S, array& V, int rank, MatrixData M) 
    {
        target_coords = M.getTargetCoordinates();
        source_coords = M.getSourceCoordinates();

        // Determining the center and radius of targets:
        c_targets = 0.5 * (af::max(target_coords) + af::min(target_coords));
        r_targets = 0.5 * (af::max(target_coords) - af::min(target_coords));

        // Determining the center and radius of sources:
        c_sources = 0.5 * (af::max(source_coords) + af::min(source_coords));
        r_sources = 0.5 * (af::max(source_coords) - af::min(source_coords));
    
        // Obtaining the standard Chebyshev nodes:
        cheb_nodes = getChebyshevNodes(rank);

        // Obtain the scaled Chebyshev nodes for the targets:
        array cheb_nodes_targets, cheb_nodes_sources;
        scalePoints(0, 1, cheb_nodes, c_targets, r_targets, cheb_nodes_targets);
        scalePoints(0, 1, cheb_nodes, c_sources, r_sources, cheb_nodes_sources);

        // Standard Locations of the coordinates:
        array standard_targets, standard_sources;
        scalePoints(c_targets, r_targets, target_coords, 0, 1, standard_targets);
        scalePoints(c_sources, r_sources, source_coords, 0, 1, cheb_nodes);
    }
}

#endif
