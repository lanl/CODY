MISH MPI
======

Description
-------

This is a C implementation of a 2-D Godunov hydrocode using MPI for mulitprocessing.

Usage
-----

The code below gives the format of the expected calls

````
mpirun -np *nproc* ./hydro *init*
````

The accepted values for init are given in the README.md file in the parent directory. *nproc* is the number of processes to use when running the code.

