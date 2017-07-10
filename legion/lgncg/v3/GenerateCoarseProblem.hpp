/**
 * Copyright (c) 2016-2017 Los Alamos National Security, LLC
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
    @file GenerateProblem.cpp

    HPCG routine
 */

#include "LegionStuff.hpp"
#include "LegionMatrices.hpp"

#include "GenerateGeometry.hpp"
#include "GenerateProblem.hpp"
#include "SetupHalo.hpp"

#include <cassert>

inline void
GenerateCoarseProblemTopLevel(
    LogicalSparseMatrix &Af,
    Context ctx,
    HighLevelRuntime *runtime
) {
    // Make local copies of geometry information.  Use global_int_t since the
    // RHS products in the calculations below may result in global range values.
    global_int_t nxf = Af.geom->nx;
    global_int_t nyf = Af.geom->ny;
    global_int_t nzf = Af.geom->nz;
    // Need fine grid dimensions to be divisible by 2
    assert(nxf % 2 == 0);
    assert(nyf % 2 == 0);
    assert(nzf % 2 == 0);
    //Coarse nx, ny, nz
    local_int_t nxc, nyc, nzc;
    nxc = nxf / 2;
    nyc = nyf / 2;
    nzc = nzf / 2;
    // This is the size of our subblock
    local_int_t localNumberOfRows = nxc * nyc * nzc;
    // If this assert fails, it most likely means that the local_int_t is set to
    // int and should be set to long long
    assert(localNumberOfRows > 0);
    // Construct the geometry and linear system
    Geometry *geomc = new Geometry();
    GenerateGeometry(
        Af.geom->size,
        Af.geom->rank,
        Af.geom->numThreads,
        nxc,
        nyc,
        nzc,
        Af.geom->stencilSize,
        geomc
    );
    //
    LogicalSparseMatrix *Ac = new LogicalSparseMatrix();
    Ac->allocate("Ac", *geomc, ctx, runtime);
    Ac->partition(geomc->size, ctx, runtime);
    Ac->geom = geomc;
    Af.Ac = Ac;
}

/*!
    Routine to construct a prolongation/restriction operator for a given fine
    grid matrix solution (as computed by a direct solver).

    @param[inout]  Af - The known system matrix, on output its coarse operator,
                   fine-to-coarse operator and auxiliary vectors will be defined.

    Note that the matrix Af is considered const because the attributes we are
    modifying are declared as mutable.

*/
inline void
GenerateCoarseProblem(
    const SparseMatrix &Af,
    Context ctx,
    HighLevelRuntime *lrt
) {
    const Geometry *AfGeom = Af.geom->data();
    assert(AfGeom);
    // Make local copies of geometry information.  Use global_int_t since the
    // RHS products in the calculations below may result in global range values.
    global_int_t nxf = AfGeom->nx;
    global_int_t nyf = AfGeom->ny;
    global_int_t nzf = AfGeom->nz;
    // Need fine grid dimensions to be divisible by 2
    assert(nxf % 2 == 0);
    assert(nyf % 2 == 0);
    assert(nzf % 2 == 0);
    //Coarse nx, ny, nz
    local_int_t nxc, nyc, nzc;
    nxc = nxf / 2;
    nyc = nyf / 2;
    nzc = nzf / 2;
    // This is the size of our subblock
    local_int_t localNumberOfRows = nxc * nyc * nzc;
    // If this assert fails, it most likely means that the local_int_t is set to
    // int and should be set to long long
    assert(localNumberOfRows > 0);
    // Construct the geometry and linear system
    GenerateGeometry(
        AfGeom->size,
        AfGeom->rank,
        AfGeom->numThreads,
        nxc,
        nyc,
        nzc,
        AfGeom->stencilSize,
        Af.Ac->geom->data()
    );
    //
    GenerateProblem(*Af.Ac, NULL, NULL, NULL, ctx, lrt);
    GetNeighborInfo(*Af.Ac);
}
