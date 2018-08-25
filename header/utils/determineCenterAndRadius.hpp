#include <iostream>
#include <arrayfire.h>
#include <string>
#include <cstdlib>

void determineBoundaries(const af::array&, double &max_val, double &min_val)
{
    max_val = af::max<double>(input_array);
    min_val = af::min<double>(input_array);
}

void determineCenterAndRadius(const af::array &input_array, double &center, double &radius)
{
    double max_val, min_val;
    determineBoundaries(input_array, max_val, min_val);
    c_targets = 0.5 * (max_val + min_val);
    r_targets = 0.5 * (max_val - min_val);
}
