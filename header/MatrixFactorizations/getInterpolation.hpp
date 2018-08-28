#ifndef __getInterpolation_hpp__
#define __getInterpolation_hpp__

#include <cstdlib>
#include "MatrixData.hpp"
#include "utils/determineCenterAndRadius.hpp"

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

// Maps the points from a domain having center c1 and radius r1
// To a domain having center c2 and radius r2
// Basically mapping from [c1 - r1 / 2, c1 + r1 / 2] --> [c2 - r2 / 2, c2 + r2 / 2]
void scalePoints(const double c1, const double r1, const array &x1,
                 const double c2, const double r2, array& x2)
{   
    x2 = c2 + r2 * (x1 - c1) /r1;
    x2.eval();
}

// Gives the interpolation operator / L2L for 1D:
void getL2L1D(array x, array x_nodes, array &L2L)
{
    // Size is the number of points in x
    // Rank is the number of points in x_nodes
    unsigned rank = x_nodes.elements();

    // Tiling x       to bring it to shape (size, rank)
    // Tiling x_nodes to bring it to shape (rank, rank)
    // We will be using x_nodes to get the denominator for 
    // the lagrange polynomials i.e. Π(x_i - x_j) where i != j
    // Similarly, we will be using x to get the numerator for
    // the lagrange polynomials i.e. num(P_i) = Π(x - x_j) where i =! j 
    x       = af::tile(x, 1, rank);
    x_nodes = af::tile(x_nodes, 1, rank);

    // Allowing broadcasting:
    af::gforSet(true);
    x       = x       - x_nodes(af::span, 0).T();
    x_nodes = x_nodes - x_nodes(af::span, 0).T();
    af::gforSet(false);

    // Performing Π(x_i - x_j) where i != j
    x_nodes = af::product(af::select(x_nodes == 0, 1, x_nodes), 1);

    // temp is used to get num(P_i) = Π(x - x_j) where i =! j
    array temp;

    // DO NOT USE GFOR HERE! THROWS WRONG RESULTS!
    for(int i = 0; i < x_nodes.dims(0); i++)
    {
        temp        = x;
        temp.col(i) = 1;
        L2L.col(i)  = af::product(temp, 1);
    }

    // Allowing broadcasting:
    af::gforSet(true);
    L2L = L2L / x_nodes.T();
    af::gforSet(false);

    L2L.eval();
}

// Gives the interpolation operator / L2L for 2D:
void getL2L2D(array &x, array &y, array &nodes, array &L2L)
{
    unsigned rank = nodes.elements();

    // Initializing:
    array L2L_x(x.elements(), rank, f64);
    array L2L_y(y.elements(), rank, f64);
    
    getL2L1D(x, nodes, L2L_x);
    getL2L1D(y, nodes, L2L_y);

    // Using a gfor loop for batched operations:
    gfor(af::seq i, rank * rank)
    {
        L2L(af::span, i) =  L2L_x(af::span, i % rank) 
                          * L2L_y(af::span, i / rank);
    }
}

// Gives the interpolation operator / L2L for 3D:
void getL2L3D(array &x, array &y, array &z, array &nodes, array &L2L)
{
    unsigned rank = nodes.elements();

    // Initializing:
    array L2L_x(x.elements(), rank, f64);
    array L2L_y(y.elements(), rank, f64);
    array L2L_z(z.elements(), rank, f64);

    getL2L1D(x, nodes, L2L_x);
    getL2L1D(y, nodes, L2L_y);
    getL2L1D(z, nodes, L2L_z);

    // Using a gfor loop for batched operations:
    gfor(af::seq i, rank * rank * rank)
    {
        L2L(af::span, i) =  L2L_x(af::span, i % rank) 
                          * L2L_y(af::span, (i / rank) % rank) 
                          * L2L_z(af::span, i / (rank * rank));
    }
}

void getM2L(const array &nodes_1, const array &nodes_2, MatrixData &M, array &M2L)
{
    // Evaluating the Kernel function at the Chebyshev nodes:
    M2L = M.buildArray(int(nodes_1.dims(0)), int(nodes_2.dims(0)),
                       nodes_1, nodes_2
                      );
    M2L.eval();
}

namespace MatrixFactorizer
{
    void getInterpolation(array &U, array &S, array &V, int n_nodes, string interpolation_type, MatrixData &M) 
    {
        array standard_nodes;
        // Obtaining the standard Chebyshev nodes of the first kind:
        if(interpolation_type == "CHEBYSHEV")
        {
            getChebyshevNodes(n_nodes, standard_nodes);
        }

        // Obtaining the standard Legendre nodes:
        else if(interpolation_type == "LEGENDRE")
        {
            getLegendreNodes(n_nodes, standard_nodes);
        }

        // Obtaining the standard equispaced nodes:
        else if(interpolation_type == "EQUISPACED")
        {
            getEquispacedNodes(n_nodes, standard_nodes);
        }

        else
        {
            cout << "Invalid choice for interpolation type" << endl;
            cout << "Please use either CHEBYSHEV, LEGENDRE or EQUISPACED" << endl;
            exit(1);
        }

        array target_coords = M.getTargetCoordinates();
        array source_coords = M.getSourceCoordinates();

        if(M.getDimensionality() == 1)
        {
            double c_targets, r_targets, c_sources, r_sources;
            // Determining the center and radius of targets:
            determineCenterAndRadius(target_coords, c_targets, r_targets);
            // Determining the center and radius of sources:
            determineCenterAndRadius(source_coords, c_sources, r_sources);

            // Obtain the scaled Chebyshev nodes for the targets:
            array nodes_targets, nodes_sources;

            scalePoints(0, 1, standard_nodes, c_targets, r_targets, nodes_targets);
            scalePoints(0, 1, standard_nodes, c_sources, r_sources, nodes_sources);

            // Standard Locations of the coordinates:
            array standard_targets, standard_sources;
            scalePoints(c_targets, r_targets, target_coords, 0, 1, standard_targets);
            scalePoints(c_sources, r_sources, source_coords, 0, 1, standard_sources);
        
            // Initializing U, S, V:
            U = af::constant(0, M.getNumRows(), n_nodes, f64);
            S = array(n_nodes, n_nodes, f64);
            V = af::constant(0, M.getNumColumns(), n_nodes, f64);

            getL2L1D(standard_targets, standard_nodes, U);
            getM2L(nodes_targets, nodes_sources, M, S);
            getL2L1D(standard_sources, standard_nodes, V);
        }

        else if(M.getDimensionality() == 2)
        {
            double c_x_targets, r_x_targets, c_x_sources, r_x_sources,
                   c_y_targets, r_y_targets, c_y_sources, r_y_sources;

            // Determining the centers and radii of targets:
            determineCenterAndRadius(target_coords(af::span, 0), c_x_targets, r_x_targets);
            determineCenterAndRadius(target_coords(af::span, 1), c_y_targets, r_y_targets);

            // Determining the centers and radii of sources:
            determineCenterAndRadius(source_coords(af::span, 0), c_x_sources, r_x_sources);
            determineCenterAndRadius(source_coords(af::span, 1), c_y_sources, r_y_sources);

            // Obtain the scaled Chebyshev nodes for the targets:
            array nodes_targets,   nodes_sources, 
                  nodes_targets_x, nodes_sources_x, 
                  nodes_targets_y, nodes_sources_y;

            scalePoints(0, 1, standard_nodes, c_x_targets, r_x_targets, nodes_targets_x);
            scalePoints(0, 1, standard_nodes, c_y_targets, r_y_targets, nodes_targets_y);
            scalePoints(0, 1, standard_nodes, c_x_sources, r_x_sources, nodes_sources_x);
            scalePoints(0, 1, standard_nodes, c_y_sources, r_y_sources, nodes_sources_y);

            nodes_targets_x = af::flat(af::tile(nodes_targets_x, 1, n_nodes));
            nodes_sources_x = af::flat(af::tile(nodes_sources_x, 1, n_nodes));
            nodes_targets_y = af::flat(af::tile(nodes_targets_y.T(), n_nodes));
            nodes_sources_y = af::flat(af::tile(nodes_sources_y.T(), n_nodes));

            // Joining along axis-1 so that the function can be passed to M2L:
            nodes_targets = af::join(1, nodes_targets_x, nodes_targets_y);
            nodes_sources = af::join(1, nodes_sources_x, nodes_sources_y);

            // Standard Locations of the coordinates:
            array standard_targets_x, standard_targets_y, standard_sources_x, standard_sources_y;
            scalePoints(c_x_targets, r_x_targets, target_coords(af::span, 0), 0, 1, standard_targets_x);
            scalePoints(c_y_targets, r_y_targets, target_coords(af::span, 1), 0, 1, standard_targets_y);
            scalePoints(c_x_sources, r_x_sources, source_coords(af::span, 0), 0, 1, standard_sources_x);
            scalePoints(c_y_sources, r_y_sources, source_coords(af::span, 1), 0, 1, standard_sources_y);

            U = array(standard_targets_x.dims(0), n_nodes * n_nodes, f64);
            V = array(standard_sources_x.dims(0), n_nodes * n_nodes, f64);

            getL2L2D(standard_targets_x, standard_targets_y, standard_nodes, U);
            getM2L(nodes_targets, nodes_sources, M, S);
            getL2L2D(standard_sources_x, standard_sources_y, standard_nodes, V);
        }

        else if(M.getDimensionality() == 3)
        {
            double c_x_targets, r_x_targets, c_x_sources, r_x_sources,
                   c_y_targets, r_y_targets, c_y_sources, r_y_sources,
                   c_z_targets, r_z_targets, c_z_sources, r_z_sources;

            // Determining the centers and radii of targets:
            determineCenterAndRadius(target_coords(af::span, 0), c_x_targets, r_x_targets);
            determineCenterAndRadius(target_coords(af::span, 1), c_y_targets, r_y_targets);
            determineCenterAndRadius(target_coords(af::span, 2), c_z_targets, r_z_targets);

            // Determining the centers and radii of sources:
            determineCenterAndRadius(source_coords(af::span, 0), c_x_sources, r_x_sources);
            determineCenterAndRadius(source_coords(af::span, 1), c_y_sources, r_y_sources);
            determineCenterAndRadius(source_coords(af::span, 2), c_z_sources, r_z_sources);

            // Obtain the scaled Chebyshev nodes for the targets:
            array nodes_targets, nodes_sources;
            array nodes_targets_x, nodes_targets_y, nodes_targets_z; 
            array nodes_sources_x, nodes_sources_y, nodes_sources_z;

            scalePoints(0, 1, standard_nodes, c_x_targets, r_x_targets, nodes_targets_x);
            scalePoints(0, 1, standard_nodes, c_y_targets, r_y_targets, nodes_targets_y);
            scalePoints(0, 1, standard_nodes, c_z_targets, r_z_targets, nodes_targets_z);

            scalePoints(0, 1, standard_nodes, c_x_sources, r_x_sources, nodes_sources_x);
            scalePoints(0, 1, standard_nodes, c_y_sources, r_y_sources, nodes_sources_y);
            scalePoints(0, 1, standard_nodes, c_z_sources, r_z_sources, nodes_sources_z);

            nodes_targets_x = af::flat(af::tile(nodes_targets_x, 1, n_nodes, n_nodes));
            nodes_sources_x = af::flat(af::tile(nodes_sources_x, 1, n_nodes, n_nodes));

            nodes_targets_y = af::flat(af::tile(nodes_targets_y.T(), n_nodes, 1, n_nodes));
            nodes_sources_y = af::flat(af::tile(nodes_sources_y.T(), n_nodes, 1, n_nodes));

            nodes_targets_z = af::flat(af::tile(af::moddims(nodes_targets_z, 1, 1, n_nodes), n_nodes, n_nodes, 1));
            nodes_sources_z = af::flat(af::tile(af::moddims(nodes_sources_z, 1, 1, n_nodes), n_nodes, n_nodes, 1));

            // Joining along axis-1 so that the function can be passed to M2L:
            nodes_targets = af::join(1, nodes_targets_x, nodes_targets_y, nodes_targets_z);
            nodes_sources = af::join(1, nodes_sources_x, nodes_sources_y, nodes_sources_z);

            // Standard Locations of the coordinates:
            array standard_targets_x, standard_targets_y, standard_targets_z;
            array standard_sources_x, standard_sources_y, standard_sources_z;

            scalePoints(c_x_targets, r_x_targets, target_coords(af::span, 0), 0, 1, standard_targets_x);
            scalePoints(c_y_targets, r_y_targets, target_coords(af::span, 1), 0, 1, standard_targets_y);
            scalePoints(c_z_targets, r_z_targets, target_coords(af::span, 2), 0, 1, standard_targets_z);
            
            scalePoints(c_x_sources, r_x_sources, source_coords(af::span, 0), 0, 1, standard_sources_x);
            scalePoints(c_y_sources, r_y_sources, source_coords(af::span, 1), 0, 1, standard_sources_y);
            scalePoints(c_z_sources, r_z_sources, source_coords(af::span, 2), 0, 1, standard_sources_z);

            U = array(standard_targets_x.dims(0), n_nodes * n_nodes * n_nodes, f64);
            V = array(standard_sources_x.dims(0), n_nodes * n_nodes * n_nodes, f64);

            getL2L3D(standard_targets_x, standard_targets_y, standard_targets_z, standard_nodes, U);
            getM2L(nodes_targets, nodes_sources, M, S);
            getL2L3D(standard_sources_x, standard_sources_y, standard_sources_z, standard_nodes, V);
        }

        else
        {
            cout << "Only dimension < 4 supported currently!" << endl;
            exit(1);
        }

        U.eval(); S.eval(); V.eval();
    }
}

#endif
