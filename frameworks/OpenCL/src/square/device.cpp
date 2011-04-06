// Device kernels...

#include "square/app.hpp"

#define STRINGIFY(X) #X

std::string const program_text = STRINGIFY(


    // output = [ x^2 | x <- input]
    __kernel void square(__global const float* input,
                         __global float* output,
                         size_t const size)
    {
        size_t i = get_global_id(0);
        if (i < size) {
            output[i] = input[i] * input[i];
        }
    }


    // output = [ x^2 | x <- input]
    __kernel void square1(__global float* input,
                         __global float* output,
                         size_t const size)
    {
        size_t i = get_global_id(0);
        if (i < size) {
            output[i] = input[i] * input[i];
        }
    }


    );


std::string const & App::get_device_program_text()
{
    return program_text;
}
