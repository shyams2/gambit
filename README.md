# gambit - A library consisting of fast low-rank factorizations and solvers:

The files of this repository are contained in:
* `header`: Here the factorization methods / utilities are declared.
* `test`: In this folder, the implementation of the various methods declared in `header` are checked. It also contains the timing runs for some standard functions like matvec, LU, QR, etc that I've tried out.

NOTE: Still in active development! Expect some adventures!

## Building the code:

At the moment, this library makes use of the [ArrayFire](www.github.com/arrayfire/arrayfire) library and the [Eigen3](http://eigen.tuxfamily.org/index.php?title=Main_Page) library. Additionally, we are also using the convenient [HighFive](https://github.com/BlueBrain/HighFive) library to perform convenient file writing in the HDF5 format. Build instructions for the same can be found on their pages. Make sure to set the environment variable `EIGEN_PATH` to the Eigen root folder since it's needed by the `Makefile`. It's also necessary to set the environment variable `HIGHFIVE_PATH` to the HighFive root directory.

## Note on naming conventions followed:

For the moment, I've followed these rules:

- Variables named with lower case with underscore. For example:`n_rows`, `n_columns`, `data_array`
- Functions are named following camelCase. For example:`findRank`, `getColumns`, `getRow`
- Classes are declared with the first letter being capitalized aswell. For example: `MatrixData`, `MatrixFactorizer`
- Constants are all defined using ALLCAPS. For example `PI`, `MAX_ROWS`
