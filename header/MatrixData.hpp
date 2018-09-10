// Class that's used to create the matrix object that will then
// be passed to the matrix factorizing class. In this class, we
// have the following options for initialization:

// 1) Passing the matrix directly in the form of the array
// 2) Passing the size of the array and rank of array
// 3) Passing the size of array along with the minimum and maximum singular values
// 4) Pass the function that states how the entries of the matrix are generated.
//    These are particularly useful since there is no involved cost in storing the
//    matrix. It follows the convention used in a generic kernel where the index of
//    source and target points along with the distances arrays of source and targets are passed

#ifndef __MatrixData_hpp__
#define __MatrixData_hpp__

#include <iostream>
#include <vector>
#include <functional>
#include <stdlib.h>
#include <stdio.h>
#include <arrayfire.h>
#include <assert.h>
#include <cmath>
// Headers needed for file-writing with HighFive:
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>

using std::cout;
using std::endl;
using std::string;
using std::abs; // NOTE: using abs without std down casts to int
using af::array;

class MatrixData
{

private:
    // array being factorized
    // NOTE: Not necessarily allocated
    array A, *A_ptr;
    // Number of rows and columns:
    size_t n_rows, n_cols;
    // Dimensionality of the points considered:
    size_t n_dims;
    // The following get used when the kernel function is passed during instantiation:
    std::function<array(array, array, array&, array&)> matrixEntriesAF;
    // When double* is used instead:
    std::function<double(unsigned, unsigned, double*, size_t, double*, size_t, size_t)> matrixEntriesDouble;
    // If the data is passed in AF form:
    array sources_af, targets_af;
    // If the data is passed in double* form:
    double *sources, *targets;
    string array_datatype; // whether af::array or double*

public:

    // When array itself is passed when calling:
    MatrixData(array &input_array);

    // Builds a matrix with given numerical rank:
    // m    - number of rows
    // n    - number of colums
    // rank - rank of matrix
    MatrixData(size_t m, size_t n, size_t rank);

    // Builds a matrix with singular values ranging uniformly from singular_min to singular_max:
    // m            - number of rows
    // n            - number of colums
    // singular_min - Lowest singular value of matrix
    // singular_max - Maximum singular value of matrix
    MatrixData(size_t m, size_t n, double singular_min, double singular_max);

    // When the blueprint to create the matrix is provided: 
    // That is function to generate matrix is passed to contructor
    // Format to be followed is matrixEntries(i, j, sources, targets)
    // i       - index / indices to access from targets
    // j       - index / indices to access from sources
    // targets - 1D array of targets locations;
    // sources - 1D array of source locations;
    MatrixData(std::function<array(array, array, array&, array&)> matrixEntries,
               array &targets, array &sources
              );

    // Alternate form when the arguments to the function are double instead of af::array:
    MatrixData(std::function<double(unsigned, unsigned, double*, size_t, double*, size_t, size_t)> matrixEntries,
               double* targets, const size_t n_targets, double* sources, const size_t n_sources,
               const size_t n_dim
              );

    // The buildArray method of this class would give a matrix of the form:
    //    S O U R C E S
    // T  x x x x x x x
    // A  x x x x x x x
    // R  x x x x x x x
    // G  x x x x x x x
    // E  x x x x x x x
    // T  x x x x x x x
    // S  x x x x x x x
    // The entry (i, j) in this matrix would give the interaction between the
    // i-th target and the j-th source
    array buildArray();
    array buildArray(size_t n_rows, size_t n_cols, array targets, array sources);

    // Estimates the rank of the matrix encoded using SVD
    int estimateRank(double tolerance);

    // Gets the condition number κ using SVD by finding out S_max / S_min
    double getConditionNumber();

    // Gets the maximum singular value of the matrix using power iterations:
    double estimateSpectralNorm(double tolerance);

    // Returns the array A stored. In the case where function generating matrix
    // is specified, the matrix is built and returned
    array getArray();

    // Returns the entries of the requested rows:
    array getRow(unsigned i);

    // Returns the entries of the requested columns:
    array getCol(unsigned i);

    // Returns the number of rows and columns:
    size_t getNumRows();
    size_t getNumCols();
    // Gets the dimensionality:
    size_t getDimensionality();
    // Returns the source and target coordinates:
    array getTargetCoordinates();
    array getSourceCoordinates();
    // Dumps the array encoded to the file mentioned under the dataset name "array"
    void dumpArray(string file_name);
    // Loads back the data from the H5File onto the object:
    void loadArray(string file_name);
};

MatrixData::MatrixData(array& input_array)
{
    this->n_rows = input_array.dims(0);
    this->n_cols = input_array.dims(1);
    this->A_ptr  = &input_array;
}

MatrixData::MatrixData(size_t m, size_t n, size_t rank)
{
    this->n_rows = m;
    this->n_cols = n;

    array Q1, Q2, R, T;
    double min_m_n = std::min(m, n);

    // Obtaining the orthogonal matrices Q1 and Q2:
    af::qr(Q1, R, T, af::randn(m, n, f64));
    af::qr(Q2, R, T, af::randn(m, n, f64));

    array U = Q1(af::span, af::seq(min_m_n));
    array V = Q2(af::seq(min_m_n), af::span);
    array S = af::randu(min_m_n, f64);
    
    // Assigning singular value above rank index to 1e-16:
    S(af::seq(rank, af::end)) = 1e-16;
    // Making it a diagonal matrix:
    S       = af::diag(S, 0, false);
    // Getting A:
    this->A = af::matmul(U, S, V);
}

MatrixData::MatrixData(size_t m, size_t n, double singular_min, double singular_max)
{
    this->n_rows = m;
    this->n_cols = n;

    array Q1, Q2, R, T;
    double min_m_n = std::min(m, n);

    // Obtaining the orthogonal matrices Q1 and Q2:
    af::qr(Q1, R, T, af::randn(m, n, f64));
    af::qr(Q2, R, T, af::randn(m, n, f64));

    array U = Q1(af::span, af::seq(min_m_n));
    array V = Q2(af::seq(min_m_n), af::span);

    double step_size = (singular_max - singular_min) / (double) (min_m_n - 1);
    array S          = singular_min + af::range(min_m_n).as(U.type()) * step_size;

    // Making it a diagonal matrix:
    S       = af::diag(S, 0, false);
    this->A = af::matmul(U, S, V);
}

MatrixData::MatrixData(std::function<array(array, array, array&, array&)> matrixEntries,
                       array &targets, array &sources
                      )
{
    this->matrixEntriesAF = matrixEntries;
    this->n_rows          = targets.dims(0);
    this->n_cols          = sources.dims(0);

    this->targets_af = targets;
    this->sources_af = sources;

    this->n_dims  = 1;

    if(targets.elements() != targets.dims(0))
    {
        this->n_dims = targets.dims(1);
    }
}

// NOTE: Undeveloped Feature
MatrixData::MatrixData(std::function<double(unsigned, unsigned, double*, size_t, double*, size_t, size_t)> matrixEntries,
                       double* targets, const size_t n_targets, double* sources, const size_t n_sources,
                       const size_t n_dims
                      )
{
    this->matrixEntriesDouble = matrixEntries;
    this->n_rows              = n_targets;
    this->n_cols              = n_sources;

    this->targets = targets;
    this->sources = sources;

    this->n_dims  = n_dims;
}

array MatrixData::buildArray()
{   
    array array_to_return;
    // Allowing broadcasting:
    af::gforSet(true);
    array_to_return = this->matrixEntriesAF(af::range(this->n_rows),
                                            af::range(this->n_cols),
                                            this->targets_af, this->sources_af
                                           );
    af::gforSet(false);

    return array_to_return;
}

// Overloaded function when the new(interpolated) sources / target locations and 
// directly provide to function used to build the kernel operator's entries:
array MatrixData::buildArray(size_t n_rows, size_t n_cols, array targets, array sources)
{
    array array_to_return;
    // Allowing broadcasting:
    af::gforSet(true);
    array_to_return = this->matrixEntriesAF(af::range(n_rows),
                                            af::range(n_cols),
                                            targets, sources
                                           );
    af::gforSet(false);

    return array_to_return;
}

// TODO: Look into faster methods instead of SVD:
int MatrixData::estimateRank(double tolerance = 1e-12) 
{
    // Temporary array that exists only in the scope of this function:
    array temp;
    
    if((this->A).elements() == 0)
    {
        temp = MatrixData::buildArray();
    }

    else
    {
        temp = this->A;
    }

    int rank = 0;
    af::array U, S, Vh;
    af::svd(U, S, Vh, temp);
    for(rank = 0; rank < S.elements(); rank++)
    {
        if(S(rank).scalar<double>() < tolerance)
            break;
    }

    return rank;
}

// TODO: Look into faster methods instead of SVD:
double MatrixData::getConditionNumber() 
{
    // Temporary array that exists only in the scope of this function:
    array temp;
    
    if((this->A).elements() == 0)
    {
        temp = MatrixData::buildArray();
    }

    else
    {
        temp = this->A;
    }

    double kappa; // Condition number

    af::array U, S, Vh;
    af::svd(U, S, Vh, temp);

    kappa = (af::max(S) / af::min(S)).scalar<double>();
    return kappa;
}

double MatrixData::estimateSpectralNorm(double tolerance = 1e-12) 
{
    // Temporary array that exists only in the scope of this function:
    array temp;
    
    if((this->A).elements() == 0)
    {
        temp = MatrixData::buildArray();
    }

    else
    {
        temp = this->A;
    }

    af::array x, y;
    x  = af::randu(temp.dims(1), f64);
    // Normalizing:
    x /= af::norm(x);

    // Setting an initial high value:
    double eps = 1e14;
    af::array s_old, s_new;

    // For the first iteration, setting s_old = -1:
    s_old = af::constant(-1, 1, f64);

    while(eps > tolerance)
    {
        y  = af::matmul(temp, x);
        y /= af::norm(y);
        x  = af::matmul(temp.T(), y);
        x /= af::norm(x);

        // Estimating the maximum singular value:
        s_new  = af::matmul(y.T(), temp, x);
        eps    = af::abs(s_new - s_old).scalar<double>();
        s_old  = s_new;
    }

    return s_new.scalar<double>();
}

array MatrixData::getArray()
{
    if((this->A).elements() == 0)
    {
        return MatrixData::buildArray();
    }

    else
    {
        return this->A;
    }
}

array MatrixData::getRow(unsigned i)
{
    if((this->A).elements() == 0)
    {
        array array_to_return;
        // Allowing broadcasting:
        af::gforSet(true);
        array_to_return = this->matrixEntriesAF(af::constant(i, 1, u32),
                                                af::range(this->n_cols),
                                                this->targets_af, this->sources_af
                                               );
        af::gforSet(false);

        return array_to_return;
    }

    else
    {
        return (this->A).row(i);
    }
}

array MatrixData::getCol(unsigned i)
{
    if((this->A).elements() == 0)
    {
        array array_to_return;
        // Allowing broadcasting:
        af::gforSet(true);
        array_to_return = this->matrixEntriesAF(af::range(this->n_rows),
                                                af::constant(i, 1, u32),
                                                this->targets_af, this->sources_af
                                               );
        af::gforSet(false);

        return array_to_return;
    }

    else
    {
        return (this->A).col(i);
    }
}

size_t MatrixData::getNumRows()
{
    return this->n_rows;
}

size_t MatrixData::getNumCols()
{
    return this->n_cols;
}

size_t MatrixData::getDimensionality()
{
    return this->n_dims;
}

array MatrixData::getTargetCoordinates()
{
    return this->targets_af;
}

array MatrixData::getSourceCoordinates()
{
    return this->sources_af;
}

void MatrixData::dumpArray(string file_name)
{
    // NOTE: Currently data is being dumped in 1D form, which will then need to be reshaped:
    HighFive::File file(file_name, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);
    std::vector<size_t> dims(1);
    dims[0] = MatrixData::getNumRows() * MatrixData::getNumCols();

    // Create the dataset
    HighFive::DataSet dataset = file.createDataSet<double>("array", HighFive::DataSpace(dims));
    double *temp = new double[dims[0]];
    this->getArray().host(temp);
    // Write it
    dataset.write(temp);
    // Deleting the temporary variable used for transfer:
    delete[] temp;
}

void MatrixData::loadArray(string file_name)
{
    HighFive::File file(file_name, HighFive::File::ReadOnly);
    std::vector<double> temp;

    // We get the dataset
    HighFive::DataSet dataset = file.getDataSet("array");
    // We convert the hdf5 dataset to a single dimension vector
    dataset.read(temp);
    this->A = array(this->getNumRows(), this->getNumCols(), temp.data(), afHost);
}

#endif
