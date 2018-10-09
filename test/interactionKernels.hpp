#include <arrayfire.h>
using af::array;

void computeRSquared2D(const array &i, const array &j, const array &targets, const array &sources, array &r_squared)
{
    array x_targets = targets(af::span, 0);
    array x_sources = sources(af::span, 0);

    array y_targets = targets(af::span, 1);
    array y_sources = sources(af::span, 1);

    array x_diff = x_targets(i) - (x_sources.T())(j);
    array y_diff = y_targets(i) - (y_sources.T())(j);
    
    // R^2 = (x_i - x_j)^2 + (y_i - y_j)^2;
    r_squared = x_diff * x_diff + y_diff * y_diff;
    r_squared.eval();
}

void computeRSquared3D(const array &i, const array &j, const array &targets, const array &sources, array &r_squared)
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
    r_squared = x_diff * x_diff + y_diff * y_diff + z_diff * z_diff;
    r_squared.eval();
}

// Laplace Equation: -Δ u = 0
// ===============SINGLE LAYER LAPLACE EQUATION KERNEL========================
// In 2D:
// S(x, y) =  - (1 / 2π) ln(r); r != 0, 0 if r = 0
// In 3D:
// S(x, y) = (1 / 4πr); r != 0, 0 if r = 0
array laplaceSingleLayer(const array &i, const array &j, const array &targets, const array &sources)
{   
    array r_squared, interaction;

    if(targets.dims(1) == 2)
    {
        computeRSquared2D(i, j, targets, sources, r_squared);
        interaction = af::select(r_squared == 0, 0, (-1 / (4 * af::Pi)) * af::log(r_squared)); 
    }

    else // assuming dim = 3
    {
        computeRSquared3D(i, j, targets, sources, r_squared);
        interaction = af::select(r_squared == 0, 0, 1 / (4 * af::Pi * af::sqrt(r_squared))); 
    }

    interaction.eval();
    return interaction;
}

// ===============INVERSE QUADRIC KERNEL========================
// K(x, y) = 1 / (1 + r^2)
array inverseQuadricKernel(const array &i, const array &j, const array &targets, const array &sources)
{   
    array r, r_squared, interaction;
    
    if(targets.dims(0) == targets.elements())
    {
        r           = targets(i) - (sources.T())(j);
        interaction = 1 / (1 + r * r);
    }

    else if(targets.dims(1) == 2)
    {
        computeRSquared2D(i, j, targets, sources, r_squared);
        interaction = 1 / (1 + r_squared);
    }

    else // assuming dim = 3
    {
        computeRSquared3D(i, j, targets, sources, r_squared);
        interaction = 1 / (1 + r_squared);
    }

    interaction.eval();
    return interaction;
}

// ===============QUADRIC KERNEL========================
// K(x, y) = (1 + r^2)
array quadricKernel(const array &i, const array &j, const array &targets, const array &sources)
{   
    array r, r_squared, interaction;
    
    if(targets.dims(0) == targets.elements())
    {
        r           = targets(i) - (sources.T())(j);
        interaction = (1 + r * r);
    }

    else if(targets.dims(1) == 2)
    {
        computeRSquared2D(i, j, targets, sources, r_squared);
        interaction = (1 + r_squared);
    }

    else // assuming dim = 3
    {
        computeRSquared3D(i, j, targets, sources, r_squared);
        interaction = (1 + r_squared);
    }
 
    interaction.eval();
    return interaction;
}

// ===============GAUSSIAN KERNEL========================
// K(x, y) = exp(-r^2)
array gaussianKernel(const array &i, const array &j, const array &targets, const array &sources)
{
    array r, r_squared, interaction;
    
    if(targets.dims(0) == targets.elements())
    {
        r           = targets(i) - (sources.T())(j);
        interaction = af::exp(-r * r);
    }

    else if(targets.dims(1) == 2)
    {
        computeRSquared2D(i, j, targets, sources, r_squared);
        interaction = af::exp(-r_squared);
    }

    else // assuming dim = 3
    {
        computeRSquared3D(i, j, targets, sources, r_squared);
        interaction = af::exp(-r_squared);
    }

    interaction.eval();
    return interaction;    
}

// Stokes Equation(Incompressible Creeping Flows): -μ Δ u + ∇.p = 0, ∇.q = 0
// ===============SINGLE LAYER STOKES EQUATION KERNEL========================
// In 2D:
// S(x, y) = (1 / 4πμ) (-ln(r) I - (r ⊗ r) / r^2); r != 0, 0 if r = 0
// In 3D:
// S(x, y) = (1 / 8πμ) ( I / r - (r ⊗ r) / r^3); r != 0, 0 if r = 0
// Below, we've taken μ = 1:
array stokesSingleLayer(const array &i, const array &j, const array &targets, const array &sources)
{   
    array r, r_kron_r, interaction;

    if(targets.dims(1) == 2)
    {
        r_kron_r = array(targets.dims(0), sources.dims(0), 4, f64);
        
        array x_targets = targets(af::span, 0);
        array x_sources = sources(af::span, 0);

        array y_targets = targets(af::span, 1);
        array y_sources = sources(af::span, 1);

        array x_diff = x_targets(i) - (x_sources.T())(j);
        array y_diff = y_targets(i) - (y_sources.T())(j);

        r_kron_r(af::span, af::span, 0) = x_diff * x_diff;
        r_kron_r(af::span, af::span, 1) = x_diff * y_diff;
        r_kron_r(af::span, af::span, 2) = y_diff * x_diff;
        r_kron_r(af::span, af::span, 3) = y_diff * y_diff;

        // Important thing to note here: this is ||r||^2
        array r_squared = x_diff * x_diff + y_diff * y_diff;

        interaction = (-1 / (4 * af::Pi)) * 
                      (  0.5 * af::log(r_squared)
                       - r_kron_r / r_squared
                      );

        // Only choosing those elements which are not infty / NaNs:
        interaction = af::select(af::isNaN(interaction) || af::isInf(interaction), 0, interaction); 
    }

    else // assuming dim = 3
    {
        r_kron_r = array(targets.dims(0), sources.dims(0), 9, f64);
        
        array x_targets = targets(af::span, 0);
        array x_sources = sources(af::span, 0);

        array y_targets = targets(af::span, 1);
        array y_sources = sources(af::span, 1);

        array z_targets = targets(af::span, 2);
        array z_sources = sources(af::span, 2);

        array x_diff = x_targets(i) - (x_sources.T())(j);
        array y_diff = y_targets(i) - (y_sources.T())(j);
        array z_diff = z_targets(i) - (z_sources.T())(j);
    
        r_kron_r(af::span, af::span, 0) = x_diff * x_diff;
        r_kron_r(af::span, af::span, 1) = x_diff * y_diff;
        r_kron_r(af::span, af::span, 2) = x_diff * z_diff;
        r_kron_r(af::span, af::span, 3) = y_diff * x_diff;
        r_kron_r(af::span, af::span, 4) = y_diff * y_diff;
        r_kron_r(af::span, af::span, 5) = y_diff * z_diff;
        r_kron_r(af::span, af::span, 6) = z_diff * x_diff;
        r_kron_r(af::span, af::span, 7) = z_diff * y_diff;
        r_kron_r(af::span, af::span, 8) = z_diff * z_diff;

        // ||r|| = sqrt((x_i - x_j)^2 + (y_i - y_j)^2 + (z_i - z_j)^2);
        r = af::sqrt(x_diff * x_diff + y_diff * y_diff + z_diff * z_diff);

        interaction = (1 / (8 * af::Pi)) * 
                      (  1 / r
                       - r_kron_r / (af::pow(r, 3))
                      );

        // Only choosing those elements which are not infty / NaNs:
        interaction = af::select(af::isNaN(interaction) || af::isInf(interaction), 0, interaction); 
   }

    interaction.eval();
    return interaction;
}
