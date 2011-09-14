This directroy contains examples of two types of hydrocodes, specifically Godunov and Two-step Lax-Wendroff stabilized with TVD, coded using CUDA by Jonathan Robey while at an internship at Los Alamos National Laboratory in the summer of 2011.

Compiling:
The codes were written using cuda 4.0.11 on a linux machine for a compute 2.0 device, but they appear to compile and run correctly on a compute 1.3 device.

The command 'make all' should be sufficient to compile assuming that the cuda libraries and header files are inthe include paths.

Running:
<exec> -i <arg_file>

init files can be generated using 'hydro_utils'
INITS: sod ibw sib expl corn cshk
In order to run one of these files, it is necessary to create an argument file, examples of these are provided in the input directory, and assume that the code will be run in the output directory.

