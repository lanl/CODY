/**
 * Copyright (c) 2017      Los Alamos National Security, LLC
 *                         All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * LA-CC 10-123
 */

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
 @file Geometry.hpp

 HPCG data structure for problem geometry
 */

#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#include <ostream>

/*!
  This defines the type for integers that have local subdomain dimension.

  Define as "long long" when local problem dimension is > 2^31
*/
typedef int local_int_t;
//typedef long long local_int_t;

/*!
  This defines the type for integers that have global dimension

  Define as "long long" when global problem dimension is > 2^31
*/
//typedef int global_int_t;
typedef long long global_int_t;

// This macro should be defined if the global_int_t is not long long
// in order to stop complaints from non-C++11 compliant compilers.
//#define HPCG_NO_LONG_LONG

/*!
  This is a data structure to contain all processor geometry information
*/
struct Geometry_STRUCT {
    int size; //!< Number of MPI processes
    int rank; //!< This process' rank in the range [0 to size - 1]
    int numThreads; //!< This process' number of threads
    int nx;   //!< Number of x-direction grid points for each local subdomain
    int ny;   //!< Number of y-direction grid points for each local subdomain
    int nz;   //!< Number of z-direction grid points for each local subdomain
    int npx;  //!< Number of processors in x-direction
    int npy;  //!< Number of processors in y-direction
    int npz;  //!< Number of processors in z-direction
    int stencilSize; //!< Size of the stencil
    int ipx;  //!< Current rank's x location in the npx by npy by npz processor grid
    int ipy;  //!< Current rank's y location in the npx by npy by npz processor grid
    int ipz;  //!< Current rank's z location in the npx by npy by npz processor grid
};
typedef struct Geometry_STRUCT Geometry;

/*!
  Returns the rank of the MPI process that is assigned the global row index
  given as the input argument.

  @param[in] geom  The description of the problem's geometry.
  @param[in] index The global row index

  @return Returns the MPI rank of the process assigned the row
*/
inline int
ComputeRankOfMatrixRow(
    const Geometry & geom,
    global_int_t index
) {
    global_int_t gnx = geom.nx*geom.npx;
    global_int_t gny = geom.ny*geom.npy;

    global_int_t iz = index/(gny*gnx);
    global_int_t iy = (index-iz*gny*gnx)/gnx;
    global_int_t ix = index%gnx;
    global_int_t ipz = iz/geom.nz;
    global_int_t ipy = iy/geom.ny;
    global_int_t ipx = ix/geom.nx;
    int rank = ipx+ipy*geom.npx+ipz*geom.npy*geom.npx;
    return rank;
}

/**
 *
 */
inline global_int_t
getGlobalXYZ(
    const Geometry &geom
) {
    global_int_t res = global_int_t(geom.npx * geom.nx) *
                       global_int_t(geom.npy * geom.ny) *
                       global_int_t(geom.npz * geom.nz);

    return res;
}

////////////////////////////////////////////////////////////////////////////////
std::ostream &
operator<<(std::ostream &os, const Geometry &geom);

#endif // GEOMETRY_HPP
