#ifndef __getInterpolation_hpp__
#define __getInterpolation_hpp__

#include <cstdlib>
#include "MatrixData.hpp"

// Returns the roots on the N-th Legendre polynomial
// These values are in the standard interval of [-1, 1]
array getLegendreNodes(int N)
{
    // Holds the coefficients of the polynomial:
    array poly = af::join(0, af::constant(0, N, f64), af::constant(1, 1, f64));

    // Getting the scaled companion matrix:
    // This has been obtained following https://docs.scipy.org/doc/numpy-1.12.0/reference/generated/numpy.polynomial.legendre.legcompanion.html
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
    return af::join(0, -1 * S(af::seq(0, af::end, 2)), S(af::seq(af::end-1, 0, -2)));
}

// Returns the roots on the N-th Chebyshev polynomial of the first kind
// These values are in the standard interval of [-1, 1]
array getChebyshevNodes(int N)
{
    return -af::cos(af::Pi * (2 * af::range(N).as(f64) + 1) / (2 * N));
}

// Returns the equispaced nodes in the standard interval of [-1, 1]
array getEquispacedNodes(int N)
{
    return (-1 + (0.5 + af::range(N).as(f64)) * 2 / double(N));
}


// Maps the points from a domain having center c1 and radius r1
// To a domain having center c2 and radius r2
// Basically mapping from [c1 - r1 / 2, c1 + r1 / 2] --> [c2 - r2 / 2, c2 + r2 / 2]
void scalePoints(double c1, double r1, array x1, double c2, double r2, array& x2)
{   
    x2 = c2 + r2 * (x1 - c1) /r1;
}

// Returns the interpolation operator / L2L
array getL2L(array x, array x_nodes)
{
    // Size is the number of points in x
    // Rank is the number of points in x_nodes
    int size = x.elements();
    int rank = x_nodes.elements();

    array L2L = af::constant(0, size, rank, f64);

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

    return L2L;
}

array getM2L(array nodes_1, array nodes_2, MatrixData M)
{
    // Evaluating the Kernel function at the Chebyshev nodes:
    return M.buildArray(int(nodes_1.dims(0)), int(nodes_2.dims(0)),
                        nodes_1, nodes_2
                       );
}

namespace MatrixFactorizer
{
    void getInterpolation(array& U, array& S, array& V, int rank, string interpolation_type, MatrixData M) 
    {
        array standard_nodes;
        // Obtaining the standard Chebyshev nodes of the first kind:
        if(interpolation_type == "CHEBYSHEV")
        {
            standard_nodes = getChebyshevNodes(rank);
        }

        // Obtaining the standard Legendre nodes:
        else if(interpolation_type == "LEGENDRE")
        {
            standard_nodes = getLegendreNodes(rank);
        }

        // Obtaining the standard equispaced nodes:
        else if(interpolation_type == "EQUISPACED")
        {
            standard_nodes = getEquispacedNodes(rank);
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
            // Determining the center and radius of targets:
            double c_targets = 0.5 * (af::max<double>(target_coords) + af::min<double>(target_coords));
            double r_targets = 0.5 * (af::max<double>(target_coords) - af::min<double>(target_coords));

            // Determining the center and radius of sources:
            double c_sources = 0.5 * (af::max<double>(source_coords) + af::min<double>(source_coords));
            double r_sources = 0.5 * (af::max<double>(source_coords) - af::min<double>(source_coords));

            // Obtain the scaled Chebyshev nodes for the targets:
            array nodes_targets, nodes_sources;
            scalePoints(0, 1, standard_nodes, c_targets, r_targets, nodes_targets);
            scalePoints(0, 1, standard_nodes, c_sources, r_sources, nodes_sources);

            // Standard Locations of the coordinates:
            array standard_targets, standard_sources;
            scalePoints(c_targets, r_targets, target_coords, 0, 1, standard_targets);
            scalePoints(c_sources, r_sources, source_coords, 0, 1, standard_sources);
        
            U = getL2L(standard_targets, standard_nodes);
            S = getM2L(nodes_targets, nodes_sources, M);
            V = getL2L(standard_sources, standard_nodes).T();
        }

        else if(M.getDimensionality() == 2)
        {
            // Determining the center and radius of targets:
            double c_x_targets = 0.5 * (  af::max<double>(target_coords(af::span, 0)) 
                                        + af::min<double>(target_coords(af::span, 0))
                                       );
            double r_x_targets = 0.5 * (  af::max<double>(target_coords(af::span, 0)) 
                                        - af::min<double>(target_coords(af::span, 0))
                                       );

            double c_y_targets = 0.5 * (  af::max<double>(target_coords(af::span, 1)) 
                                        + af::min<double>(target_coords(af::span, 1))
                                       );
            double r_y_targets = 0.5 * (  af::max<double>(target_coords(af::span, 1)) 
                                        - af::min<double>(target_coords(af::span, 1))
                                       );

            c_x_targets = c_y_targets = -1;
            r_x_targets = r_y_targets = 0.5;

            cout << "r = (" << r_x_targets << "," << r_y_targets << ")" << endl;
            cout << "c = (" << c_x_targets << "," << c_y_targets << ")" << endl;

            // Determining the center and radius of sources:
            double c_x_sources = 0.5 * (  af::max<double>(source_coords(af::span, 0)) 
                                        + af::min<double>(source_coords(af::span, 0))
                                       );
            double r_x_sources = 0.5 * (  af::max<double>(source_coords(af::span, 0)) 
                                        - af::min<double>(source_coords(af::span, 0))
                                       );

            double c_y_sources = 0.5 * (  af::max<double>(source_coords(af::span, 1)) 
                                        + af::min<double>(source_coords(af::span, 1))
                                       );
            double r_y_sources = 0.5 * (  af::max<double>(source_coords(af::span, 1)) 
                                        - af::min<double>(source_coords(af::span, 1))
                                       );

            c_x_sources = c_y_sources = 1;
            r_x_sources = r_y_sources = 0.5;

            cout << "For the sources:" << endl;
            cout << "r = (" << r_x_sources << "," << r_y_sources << ")" << endl;
            cout << "c = (" << c_x_sources << "," << c_y_sources << ")" << endl;
        
            // Obtain the scaled Chebyshev nodes for the targets:
            array nodes_targets, nodes_sources, nodes_targets_x, nodes_targets_y, nodes_sources_x, nodes_sources_y;
            scalePoints(0, 1, standard_nodes, c_x_targets, r_x_targets, nodes_targets_x);
            scalePoints(0, 1, standard_nodes, c_y_targets, r_y_targets, nodes_targets_y);
            scalePoints(0, 1, standard_nodes, c_x_sources, r_x_sources, nodes_sources_x);
            scalePoints(0, 1, standard_nodes, c_y_sources, r_y_sources, nodes_sources_y);

            nodes_targets_x = af::flat(af::tile(nodes_targets_x, 1, rank));
            nodes_sources_x = af::flat(af::tile(nodes_sources_x, 1, rank));
            nodes_targets_y = af::flat(af::tile(nodes_targets_y.T(), rank));
            nodes_sources_y = af::flat(af::tile(nodes_sources_y.T(), rank));

            // Joining along axis-1 so that the function can be passed to M2L:
            nodes_targets = af::join(1, nodes_targets_x, nodes_targets_y);
            nodes_sources = af::join(1, nodes_sources_x, nodes_sources_y);

            // Standard Locations of the coordinates:
            array standard_targets_x, standard_targets_y, standard_sources_x, standard_sources_y;
            scalePoints(c_x_targets, r_x_targets, target_coords(af::span, 0), 0, 1, standard_targets_x);
            scalePoints(c_y_targets, r_y_targets, target_coords(af::span, 1), 0, 1, standard_targets_y);
            scalePoints(c_x_sources, r_x_sources, source_coords(af::span, 0), 0, 1, standard_sources_x);
            scalePoints(c_y_sources, r_y_sources, source_coords(af::span, 1), 0, 1, standard_sources_y);

            array U_x = getL2L(standard_targets_x, standard_nodes);
            array U_y = getL2L(standard_targets_y, standard_nodes);
            
            af::print("U_x", U_x, 16);
            af::print("U_y", U_y, 16);

            U = af::constant(0, standard_targets_x.dims(0), rank * rank, f64);
            for(int i = 0; i < standard_targets_x.dims(0); i++)
            {
                for(int j = 0; j < rank * rank; j++)
                {
                    U(i, j) = U_x(i, j % rank) * U_x(i, j / rank);
                }
            }

            af::print("U", U, 16);

            S = getM2L(nodes_targets, nodes_sources, M);

            array V_x = getL2L(standard_sources_x, standard_nodes);
            array V_y = getL2L(standard_sources_y, standard_nodes);

            // af::print("V_x", V_x, 16);
            // af::print("V_y", V_y, 16);

            V = af::constant(0, standard_sources_x.dims(0), rank * rank, f64);
            for(int i = 0; i < standard_sources_x.dims(0); i++)
            {
                for(int j = 0; j < rank * rank; j++)
                {
                    V(i, j) = V_x(i, j % rank) * V_y(i, j / rank);
                }
            }

            // Taking transpose of V:
            V = V.T();
        }

        // else if(M.getDimensionality() == 3)
        // {
        //     // Determining the center and radius of targets:
        //     double c_x_targets = 0.5 * (  af::max<double>(target_coords(af::span, 0)) 
        //                                 + af::min<double>(target_coords(af::span, 0))
        //                                );
        //     double r_x_targets = 0.5 * (  af::max<double>(target_coords(af::span, 0)) 
        //                                 - af::min<double>(target_coords(af::span, 0))
        //                                );

        //     double c_y_targets = 0.5 * (  af::max<double>(target_coords(af::span, 1)) 
        //                                 + af::min<double>(target_coords(af::span, 1))
        //                                );
        //     double r_y_targets = 0.5 * (  af::max<double>(target_coords(af::span, 1)) 
        //                                 - af::min<double>(target_coords(af::span, 1))
        //                                );

        //     double c_z_targets = 0.5 * (  af::max<double>(target_coords(af::span, 2)) 
        //                                 + af::min<double>(target_coords(af::span, 2))
        //                                );
        //     double r_z_targets = 0.5 * (  af::max<double>(target_coords(af::span, 2)) 
        //                                 - af::min<double>(target_coords(af::span, 2))
        //                                );

        //     // Determining the center and radius of sources:
        //     double c_x_sources = 0.5 * (  af::max<double>(source_coords(af::span, 0)) 
        //                                 + af::min<double>(source_coords(af::span, 0))
        //                                );
        //     double r_x_sources = 0.5 * (  af::max<double>(source_coords(af::span, 0)) 
        //                                 - af::min<double>(source_coords(af::span, 0))
        //                                );

        //     double c_y_sources = 0.5 * (  af::max<double>(source_coords(af::span, 1)) 
        //                                 + af::min<double>(source_coords(af::span, 1))
        //                                );
        //     double r_y_sources = 0.5 * (  af::max<double>(source_coords(af::span, 1)) 
        //                                 - af::min<double>(source_coords(af::span, 1))
        //                                );

        //     double c_z_sources = 0.5 * (  af::max<double>(source_coords(af::span, 1)) 
        //                                 + af::min<double>(source_coords(af::span, 1))
        //                                );

        //     double r_z_sources = 0.5 * (  af::max<double>(source_coords(af::span, 2)) 
        //                                 - af::min<double>(source_coords(af::span, 2))
        //                                );

        //     // Obtain the scaled Chebyshev nodes for the targets:
        //     array nodes_targets, nodes_sources,
        //     array nodes_targets_x, nodes_targets_y, nodes_targets_z; 
        //     array nodes_sources_x, nodes_sources_y, nodes_sources_z;

        //     scalePoints(0, 1, standard_nodes, c_x_targets, r_x_targets, nodes_targets_x);
        //     scalePoints(0, 1, standard_nodes, c_y_targets, r_y_targets, nodes_targets_y);
        //     scalePoints(0, 1, standard_nodes, c_x_sources, r_x_sources, nodes_sources_x);
        //     scalePoints(0, 1, standard_nodes, c_y_sources, r_y_sources, nodes_sources_y);
        //     scalePoints(0, 1, standard_nodes, c_x_sources, r_x_sources, nodes_sources_x);
        //     scalePoints(0, 1, standard_nodes, c_y_sources, r_y_sources, nodes_sources_y);

        //     nodes_targets_x = af::flat(af::tile(nodes_targets_x, 1, rank));
        //     nodes_sources_x = af::flat(af::tile(nodes_sources_x, 1, rank));
        //     nodes_targets_y = af::flat(af::tile(nodes_targets_y.T(), rank));
        //     nodes_sources_y = af::flat(af::tile(nodes_sources_y.T(), rank));

        //     // Joining along axis-1 so that the function can be passed to M2L:
        //     nodes_targets = af::join(1, nodes_targets_x, nodes_targets_y);
        //     nodes_sources = af::join(1, nodes_sources_x, nodes_sources_y);

        //     // Standard Locations of the coordinates:
        //     array standard_targets_x, standard_targets_y, standard_sources_x, standard_sources_y;
        //     scalePoints(c_x_targets, r_x_targets, target_coords(af::span, 0), 0, 1, standard_targets_x);
        //     scalePoints(c_y_targets, r_y_targets, target_coords(af::span, 1), 0, 1, standard_targets_y);
        //     scalePoints(c_x_sources, r_x_sources, source_coords(af::span, 0), 0, 1, standard_sources_x);
        //     scalePoints(c_y_sources, r_y_sources, source_coords(af::span, 1), 0, 1, standard_sources_y);

        //     array U_x = getL2L(standard_targets_x, standard_nodes);
        //     array U_y = getL2L(standard_targets_y, standard_nodes);
         
        //     U = af::constant(0, standard_targets_x.dims(0), rank * rank, f64);
        //     for(int i = 0; i < standard_targets_x.dims(0); i++)
        //     {
        //         for(int j = 0; j < rank * rank; j++)
        //         {
        //             U(i, j) = U_x(i, j % rank) * U_x(i, j / rank);
        //         }
        //     }       

        //     S = getM2L(nodes_targets, nodes_sources, M);

        //     array V_x = getL2L(standard_sources_x, standard_nodes);
        //     array V_y = getL2L(standard_sources_y, standard_nodes);

        //     V = af::constant(0, standard_sources_x.dims(0), rank * rank, f64);
        //     for(int i = 0; i < standard_sources_x.dims(0); i++)
        //     {
        //         for(int j = 0; j < rank * rank; j++)
        //         {
        //             V(i, j) = V_x(i, j % rank) * V_y(i, j / rank);
        //         }
        //     }

        //     // Taking transpose of V:
        //     V = V.T();
        // }

        else
        {
            cout << "Only dimension < 4 supported currently!" << endl;
            exit(1);
        }
    }
}

#endif
