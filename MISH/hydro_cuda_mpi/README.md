MISH MPI/CUDA
======

Description
-------

This is a CUDA implementation of a 2-D Godunov hydrocode which uses MPI for interprocess, and thus device, communication.

Usage
-----

The code below gives the format of the expected calls

````
mpirun -np *nproc* ./hydro *init*
````

The accepted values for *init* are given in the README.md file in the parent directory. *nproc* gives the number of processes, and thus devices to run the code on.

