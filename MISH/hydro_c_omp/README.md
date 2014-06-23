MISH OMP
======

Description
-------

This is a C implementation of a 2-D Godunov hydrocode using OpenMP for Multithreading.

Usage
-----

The code below gives the format of the expected calls

````
export OMP_NUM_THREADS=*nth*; ./hydro *init*
````

The accepted values for *init* are given in the README.md file in the parent directory. *nth* is the maximum number of threads to use.

