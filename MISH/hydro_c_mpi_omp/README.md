MISH MPI/OMP
======

Description
-------

This is a C implementation of a 2-D Godunov hydrocode using MPI for interprocess communication and OpenMP for multithreading within that process.

Usage
-----

The code below gives the format of the expected calls

````
export OMP_NUM_THREADS=*nth*; mprirun -np *nproc* hydro *init*
````

The accepted values for *init* are given in the README.md file in the parent directory. *nproc* and *nth* give the number of processes and threads to use for the computation respectively.

