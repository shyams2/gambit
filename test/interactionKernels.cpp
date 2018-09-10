#include <arrayfire.h>
using af::array;

void computeRSquared2D(array &i, array &j, array &targets, array &sources, array &r_squared)
{
    array x_targets = targets(af::span, 0);
    array x_sources = sources(af::span, 0);

    array y_targets = targets(af::span, 1);
    array y_sources = sources(af::span, 1);

    array x_diff = x_targets(i) - (x_sources.T())(j);
    array y_diff = y_targets(i) - (y_sources.T())(j);
    
    // R^2 = (x_i - x_j)^2 + (y_i - y_j)^2;
    r_squared =   x_diff * x_diff + y_diff * y_diff;
    r_squared.eval();
}

void computeRSquared3D(array &i, array &j, array &targets, array &sources, array &r_squared)
{
    array x_targets = targets(af::span, 0);
    array x_sources = sources(af::span, 0);

    array y_targets = targets(af::span, 1);
    array y_sources = sources(af::span, 1);

    array z_targets = targets(af::span, 2);
    array z_sources = sources(af::span, 2);

    array x_diff = x_targets(i) - (x_sources.T())(j);
    array y_diff = y_targets(i) - (y_sources.T())(j);
    array z_diff = z_targets(i) - (z_sources.T())(j);
    
    // R^2 = (x_i - x_j)^2 + (y_i - y_j)^2 + (z_i - z_j)^2;
    r_squared   = x_diff * x_diff + y_diff * y_diff + z_diff * z_diff;
    r_squared.eval();
}

// Laplace Equation: -Δ u = 0
// ===============SINGLE LAYER LAPLACE EQUATION KERNEL 2D========================
// S(x, y) =  - (1 / 2π) ln(r); r != 0, 0 if r = 0
array laplaceSingleLayer2D(array &i, array &j, array &targets, array &sources)
{   
    array r_squared;
    computeRSquared2D(i, j, targets, sources, r_squared);
    array interaction = af::select(r_squared == 0, 0, (-1 / (4 * af::Pi)) * af::log(r_squared)); 

    interaction.eval();
    return interaction;
}

// ===============SINGLE LAYER LAPLACE EQUATION KERNEL 3D========================
// S(x, y) = (1 / 4πr); r != 0, 0 if r = 0
array laplaceSingleLayer3D(array &i, array &j, array &targets, array &sources)
{   
    array r_squared;
    computeRSquared3D(i, j, targets, sources, r_squared);
    array interaction = af::select(r_squared == 0, 0, 1 / (4 * af::Pi * af::sqrt(r_squared))); 

    interaction.eval();
    return interaction;
}

// ===============INVERSE QUADRIC KERNEL========================
// K(x, y) = 1 / (1 + r^2)
array inverseQuadricKernel(array &i, array &j, array &targets, array &sources)
{   
    array r           = targets(i) - (sources.T())(j);
    array interaction = 1 / (1 + r * r);

    interaction.eval();
    return interaction;
}

// ===============QUADRIC KERNEL========================
// K(x, y) = (1 + r^2)
array quadricKernel(array &i, array &j, array &targets, array &sources)
{   
    array r           = targets(i) - (sources.T())(j);
    array interaction = 1 + r * r;

    interaction.eval();
    return interaction;
}

// ===============GAUSSIAN KERNEL========================
// K(x, y) = exp(-r^2)
array gaussianKernel(array &i, array &j, array &targets, array &sources)
{
    array r           = targets(i) - (sources.T())(j);
    array interaction = af::exp(-r * r);

    interaction.eval();
    return interaction;    
}
