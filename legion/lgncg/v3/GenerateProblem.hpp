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
 @file GenerateProblem.hpp

 HPCG routine
 */

#pragma once

#include "hpcg.hpp"

#include "LegionStuff.hpp"
#include "LegionItems.hpp"
#include "LegionArrays.hpp"
#include "LegionMatrices.hpp"
#include "ReduceSum.hpp"

#include <map>
#include <sstream>
#include <cassert>

/**
 * Tally total number of non-zeros in simulation.
 */
static inline global_int_t
getTotalNumberOfNonZeros(
    SparseMatrix &A,
    global_int_t localNumberOfNonzeros,
    LegionRuntime::HighLevel::Context ctx,
    LegionRuntime::HighLevel::Runtime *runtime
) {
    Item< DynColl<global_int_t> > *dcars = A.dcAllRedSumGI;
    dcars->data()->localBuffer = localNumberOfNonzeros;
    //
    TaskLauncher tlLocalNZ(
        LOCAL_NONZEROS_TID,
        TaskArgument(NULL, 0)
    );
    //
    tlLocalNZ.add_region_requirement(
        RegionRequirement(
            dcars->logicalRegion,
            RO_E,
            dcars->logicalRegion
        )
    ).add_field(dcars->getFieldID());
    //
    Future future = runtime->execute_task(ctx, tlLocalNZ);
    //
    auto &dyncol = dcars->data()->dc;
    runtime->defer_dynamic_collective_arrival(ctx, dyncol, future);
    //
    dyncol = runtime->advance_dynamic_collective(ctx, dyncol);
    //
    Future fSum = runtime->get_dynamic_collective_result(ctx, dyncol);

    return fSum.get<global_int_t>();
}

/*!
    Reference version of GenerateProblem to generate the sparse matrix, right
    hand side, initial guess, and exact solution.

    @param[in]  A        The known system matrix
    @param[inout] b      The newly allocated and generated right hand side
    vector (if b!=0 on entry)
    @param[inout] x      The newly allocated solution vector with entries set to
    0.0 (if x!=0 on entry)
    @param[inout] xexact The newly allocated solution vector with entries set to
    the exact solution (if the xexact!=0 non-zero on entry)

    @see GenerateGeometry
*/
inline void
GenerateProblem(
    SparseMatrix &A,
    Array<floatType> *b,
    Array<floatType> *x,
    Array<floatType> *xexact,
    LegionRuntime::HighLevel::Context ctx,
    LegionRuntime::HighLevel::Runtime *runtime
) {
    using namespace std;
    const Geometry *const Ageom = A.geom->data();
    // Make local copies of geometry information.  Use global_int_t since the
    // RHS products in the calculations below may result in global range values.
    const global_int_t nx  = Ageom->nx;
    const global_int_t ny  = Ageom->ny;
    const global_int_t nz  = Ageom->nz;
    const global_int_t npx = Ageom->npx;
    const global_int_t npy = Ageom->npy;
    const global_int_t npz = Ageom->npz;
    const global_int_t ipx = Ageom->ipx;
    const global_int_t ipy = Ageom->ipy;
    const global_int_t ipz = Ageom->ipz;
    const global_int_t gnx = nx * npx;
    const global_int_t gny = ny * npy;
    const global_int_t gnz = nz * npz;
    // This is the size of our subblock.
    const local_int_t localNumberOfRows = nx * ny * nz;
    // If this assert fails, it most likely means that the local_int_t is set to
    // int and should be set to.  Throw an exception of the number of rows is
    // less than zero (can happen if int overflow)long long
    assert(localNumberOfRows > 0);
    // We are approximating a 27-point finite element/volume/difference 3D
    // stencil
    const local_int_t numberOfNonzerosPerRow = Ageom->stencilSize;

    // Total number of grid points in mesh
    const global_int_t totalNumberOfRows = ((global_int_t)localNumberOfRows)
                                         * ((global_int_t)Ageom->size);
    // If this assert fails, it most likely means that the global_int_t is set
    // to int and should be set to long long
    assert(totalNumberOfRows > 0);
    ////////////////////////////////////////////////////////////////////////////
    // Allocate arrays
    ////////////////////////////////////////////////////////////////////////////
    if (Ageom->rank == 0) {
        const size_t mn = localNumberOfRows * numberOfNonzerosPerRow;
        const size_t sparseMatMemInB = (
            sizeof(char)         * localNumberOfRows //nonzerosInRow
          + sizeof(global_int_t) * mn                //mtxIndG
          + sizeof(local_int_t)  * mn                //mtxIndL
          + sizeof(floatType)    * mn                //matrixValues
          + sizeof(floatType)    * localNumberOfRows //matrixDiagonal
          + sizeof(global_int_t) * localNumberOfRows //localToGlobalMap
        ) * Ageom->size;
        //
        const size_t vectorsMemInB = (
            (b      ? sizeof(floatType) * localNumberOfRows : 0)
          + (x      ? sizeof(floatType) * localNumberOfRows : 0)
          + (xexact ? sizeof(floatType) * localNumberOfRows : 0)
        ) * Ageom->size;
        //
        const size_t pMemInB = sparseMatMemInB + vectorsMemInB;
        const double pMemInMB = pMemInB / 1024.0 / 1024.0;
        cout << "--> Approximate Generate Problem Memory Footprint="
             << pMemInMB << " MB" << endl;
    }
    char *nonzerosInRow = A.nonzerosInRow->data();
    assert(nonzerosInRow);
    // Interpreted as 2D array
    Array2D<global_int_t> mtxIndG(
        localNumberOfRows, numberOfNonzerosPerRow, A.mtxIndG->data()
    );
    // Interpreted as 2D array
    Array2D<local_int_t> mtxIndL(
        localNumberOfRows, numberOfNonzerosPerRow, A.mtxIndL->data()
    );
    // Interpreted as 2D array
    Array2D<floatType> matrixValues(
        localNumberOfRows, numberOfNonzerosPerRow, A.matrixValues->data()
    );
    //
    floatType *matrixDiagonal = A.matrixDiagonal->data();
    //
    global_int_t *localToGlobalMap = A.localToGlobalMap->data();
    //
    floatType *bv      = nullptr;
    floatType *xv      = nullptr;
    floatType *xexactv = nullptr;
    if (b != 0) {
        bv = b->data(); // Only compute exact solution if requested
        assert(bv);
    }
    if (x != 0) {
        xv = x->data(); // Only compute exact solution if requested
        assert(xv);
    }
    if (xexact != 0) {
        xexactv = xexact->data(); // Only compute exact solution if requested
        assert(xexactv);
    }
    //
    global_int_t localNumberOfNonzeros = 0;
    //
    for (local_int_t iz=0; iz<nz; iz++) {
        global_int_t giz = ipz*nz+iz;
        for (local_int_t iy=0; iy<ny; iy++) {
            global_int_t giy = ipy*ny+iy;
            for (local_int_t ix=0; ix<nx; ix++) {
                global_int_t gix = ipx*nx+ix;
                local_int_t currentLocalRow = iz*nx*ny+iy*nx+ix;
                global_int_t currentGlobalRow = giz*gnx*gny+giy*gnx+gix;
                localToGlobalMap[currentLocalRow] = currentGlobalRow;
                char numberOfNonzerosInRow = 0;
                // Current index in current row
                global_int_t currentIndexG = 0;
                local_int_t currentNonZeroElemIndex = 0;
                for (int sz=-1; sz<=1; sz++) {
                    if (giz+sz>-1 && giz+sz<gnz) {
                        for (int sy=-1; sy<=1; sy++) {
                            if (giy+sy>-1 && giy+sy<gny) {
                                for (int sx=-1; sx<=1; sx++) {
                                    if (gix+sx>-1 && gix+sx<gnx) {
                                        global_int_t curcol = currentGlobalRow
                                                            + sz*gnx*gny
                                                            + sy*gnx+sx;
                                        if (curcol==currentGlobalRow) {
                                            matrixDiagonal[currentLocalRow] = 26.0;
                                            matrixValues(currentLocalRow, currentNonZeroElemIndex) = 26.0;
                                        } else {
                                            matrixValues(currentLocalRow, currentNonZeroElemIndex) = -1.0;
                                        }
                                        currentNonZeroElemIndex++;
                                        mtxIndG(currentLocalRow, currentIndexG++) = curcol;
                                        numberOfNonzerosInRow++;
                                    } // end x bounds test
                                } // end sx loop
                            } // end y bounds test
                        } // end sy loop
                    } // end z bounds test
                } // end sz loop
                nonzerosInRow[currentLocalRow] = numberOfNonzerosInRow;
                localNumberOfNonzeros += numberOfNonzerosInRow;
                if (b != 0) {
                    bv[currentLocalRow] = 26.0
                                        - ((double)(numberOfNonzerosInRow-1));
                }
                if (x != 0) {
                    xv[currentLocalRow] = 0.0;
                }
                if (xexact != 0) {
                    xexactv[currentLocalRow] = 1.0;
                }
            } // end ix loop
        } // end iy loop
    } // end iz loop

    SparseMatrixScalars *Asclrs   = A.sclrs->data();
    Asclrs->totalNumberOfRows     = totalNumberOfRows;
    Asclrs->localNumberOfRows     = localNumberOfRows;
    // Will be updated later to include external values in GetNeighborInfo.
    Asclrs->localNumberOfColumns  = localNumberOfRows;
    Asclrs->localNumberOfNonzeros = localNumberOfNonzeros;
    //
    Asclrs->totalNumberOfNonzeros = getTotalNumberOfNonZeros(
        A,
        localNumberOfNonzeros,
        ctx,
        runtime
    );
    // If this assert fails, it most likely means that the global_int_t is
    // set to int and should be set to long long This assert is usually the
    // first to fail as problem size increases beyond the 32-bit integer
    // range.  Throw an exception of the number of nonzeros is less than
    // zero (can happen if int overflow)
    assert(Asclrs->totalNumberOfNonzeros > 0);
}
