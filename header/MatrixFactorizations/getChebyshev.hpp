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
array getChebyshevPolynomials(int N, array x)
{
    // Initializing T:
    array T = af::constant(0, x.elements(), N, f64);
    // The n-th Chebyshev polynomial is given by T_n(x) = cos(n * arccos(x))
    for(int i = 0; i < N; i++)
    {
        T.col(i) = af::cos(i * af::acos(x));
    }

    return T;
}

// Returns the interpolation operator / L2L
array getChebyshevL2L(array x, array x_cheb)
{
    // Getting the Chebyshev polymonials:
    // Rank is the number of points in x_cheb
    int rank = x_cheb.elements();

    // Getting T_x and T_cheb:
    array T_x    = getChebyshevPolynomials(rank, x);
    array T_cheb = getChebyshevPolynomials(rank, x_cheb);

    array L2L = (2 * af::matmul(T_x, T_cheb.T()) - 1) / rank;
    return L2L;
}

array getM2L(array cheb_nodes_1, array cheb_nodes_2, MatrixData M)
{
    // Evaluating the Kernel function at the Chebyshev nodes:
    return M.buildArray(int(cheb_nodes_1.elements()), int(cheb_nodes_2.elements()),
                        cheb_nodes_1, cheb_nodes_2
                       );
}

namespace MatrixFactorizer
{
    void getChebyshev(array& U, array& S, array& V, int rank, MatrixData M) 
    {
        array target_coords = M.getTargetCoordinates();
        array source_coords = M.getSourceCoordinates();

        // Determining the center and radius of targets:
        double c_targets = 0.5 * (af::max<double>(target_coords) + af::min<double>(target_coords));
        double r_targets = 0.5 * (af::max<double>(target_coords) - af::min<double>(target_coords));

        // Determining the center and radius of sources:
        double c_sources = 0.5 * (af::max<double>(source_coords) + af::min<double>(source_coords));
        double r_sources = 0.5 * (af::max<double>(source_coords) - af::min<double>(source_coords));
    
        // Obtaining the standard Chebyshev nodes:
        array standard_cheb_nodes = getChebyshevNodes(rank);

        // Obtain the scaled Chebyshev nodes for the targets:
        array cheb_nodes_targets, cheb_nodes_sources;
        scalePoints(0, 1, standard_cheb_nodes, c_targets, r_targets, cheb_nodes_targets);
        scalePoints(0, 1, standard_cheb_nodes, c_sources, r_sources, cheb_nodes_sources);

        // Standard Locations of the coordinates:
        array standard_targets, standard_sources;
        scalePoints(c_targets, r_targets, target_coords, 0, 1, standard_targets);
        scalePoints(c_sources, r_sources, source_coords, 0, 1, standard_sources);
    
        U = getChebyshevL2L(standard_targets, standard_cheb_nodes);
        S = getM2L(cheb_nodes_targets, cheb_nodes_sources, M);
        V = getChebyshevL2L(standard_sources, standard_cheb_nodes).T();
    }
}

#endif
