//@HEADER
// ***************************************************
//
// HPCG: High Performance Conjugate Gradient Benchmark
//
// Contact:
// Michael A. Heroux ( maherou@sandia.gov)
// Jack Dongarra     (dongarra@eecs.utk.edu)
// Piotr Luszczek    (luszczek@eecs.utk.edu)
//
// ***************************************************
//@HEADER

/*!
 @file GenerateGeometry.cpp

 HPCG routine
 */

#pragma once

#include "Geometry.hpp"
#include "ComputeOptimalShapeXYZ.hpp"
#include "GenerateGeometry.hpp"

#include <cmath>
#include <cstdlib>

/*!
    Computes the factorization of the total number of processes into a
    3-dimensional process grid that is as close as possible to a cube. The
    quality of the factorization depends on the prime number structure of the
    total number of processes. It then stores this decompostion together with
    the parallel parameters of the run in the geometry data structure.

    @param[in]  size total number of MPI processes
    @param[in]  rank this process' rank among other MPI processes
    @param[in]  numThreads number of OpenMP threads in this process
    @param[in]  nx, ny, nz number of grid points for each local block in the x,
                y, and z dimensions, respectively
    @param[out] geom data structure that will store the above parameters and the
    factoring of total number of processes into three dimensions
*/
inline void
GenerateGeometry(
    int size,
    int rank,
    int numThreads,
    int nx,
    int ny,
    int nz,
    int stencilSize,
    Geometry *geom
) {
    using namespace std;

    int npx, npy, npz;
    //
    ComputeOptimalShapeXYZ(size, npx, npy, npz);
    // Now compute this process's indices in the 3D cube
    int ipz = rank / (npx * npy);
    int ipy = (rank - ipz * npx * npy) / npx;
    int ipx = rank % npx;
#if 0
    cout << "size = "<< size << endl
         << "nx  = " << nx << endl
         << "ny  = " << ny << endl
         << "nz  = " << nz << endl
         << "npx = " << npx << endl
         << "npy = " << npy << endl
         << "npz = " << npz << endl;
    }
    cout << "For rank = " << rank << endl
         << "ipx = " << ipx << endl
         << "ipy = " << ipy << endl
         << "ipz = " << ipz << endl;
    //     
    assert(size == npx * npy * npz);
#endif
    geom->size = size;
    geom->rank = rank;
    geom->numThreads = numThreads;
    geom->nx = nx;
    geom->ny = ny;
    geom->nz = nz;
    geom->npx = npx;
    geom->npy = npy;
    geom->npz = npz;
    geom->stencilSize = stencilSize;
    geom->ipx = ipx;
    geom->ipy = ipy;
    geom->ipz = ipz;
}
