#include <iostream>
#include <arrayfire.h>
#include <string>
#include <cstdlib>

using std::endl;
using std::cout;

void compute_error(af::array &Z_approx, af::array &Z, string &error_type)
{
    double abs_error, rel_error;

    // Calculating the error using L1-norm:
    if(error_type == "L1")
    {
        abs_error = af::sum<double>(af::abs(Z_approx - Z));
        rel_error = af::sum<double>(af::abs(Z_approx - Z)) / af::sum<double>(af::abs(Z));
    }

    // Calculating the error using L2-norm:
    if(error_type == "L2")
    {
        abs_error = af::norm(Z_approx - Z);
        rel_error = af::norm(Z_approx - Z) / af::norm(Z);
    }

    // Calculating the error using L∞-norm:
    if(error_type == "L-inf" || error_type == "L-max")
    {
        abs_error = af::max<double>(Z_approx - Z);
        rel_error = af::max<double>(Z_approx - Z) / af::max<double>(Z);
    }

    else
    {
        cout << "Invalid choice:Use L1, L2, L-inf" << endl;
        exit(1);
    }

    cout << "Absolute Error:" << abs_error << endl;
    cout << "Relative Error:" << rel_error << endl << endl;
}
