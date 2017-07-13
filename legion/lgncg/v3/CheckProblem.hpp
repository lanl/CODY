/**
 * Copyright (c)      2017 Los Alamos National Security, LLC
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
    @file CheckProblem.hpp

    HPCG routine
 */

#pragma once

#include "LegionStuff.hpp"
#include "LegionArrays.hpp"
#include "LegionMatrices.hpp"
#include "ProblemCommon.hpp"

#include <fstream>
#include "hpcg.hpp"
#include <cassert>


/*!
    Check the contents of the generated sparse matrix to see if values match
    expected contents.

    @param[in]  A        The known system matrix.

    @param[inout] b      The newly allocated and generated right hand side
                         vector (if b!=0 on entry).

    @param[inout] x      The newly allocated solution vector with entries set to
                         0.0 (if x!=0 on entry).

    @param[inout] xexact The newly allocated solution vector with entries set to
                         the exact solution (if the xexact!=0 non-zero on
                         entry).

    @see GenerateGeometry
*/

inline void
CheckProblem(
    SparseMatrix &A,
    Array<floatType> *b,
    Array<floatType> *x,
    Array<floatType> *xexact,
    Context ctx,
    Runtime *runtime
) {
    // Make local copies of geometry information.  Use global_int_t since the
    // RHS products in the calculations below may result in global range values.
    const Geometry *const Ageom = A.geom->data();
    const auto *const Asclrs = A.sclrs->data();
    //
    const global_int_t nx =  Ageom->nx;
    const global_int_t ny =  Ageom->ny;
    const global_int_t nz =  Ageom->nz;
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
    local_int_t localNumberOfRows = nx * ny * nz;
    // Total number of grid points in mesh.
    global_int_t totalNumberOfRows = ((global_int_t)localNumberOfRows)
                                   * ((global_int_t)Ageom->size);

    floatType *bv = 0;
    floatType *xv = 0;
    floatType *xexactv = 0;
    // Only compute exact solution if requested.
    if (b != 0) bv = b->data();
    if (x != 0) xv = x->data();
    if (xexact != 0) xexactv = xexact->data();

    const floatType *const AmatrixDiagonal = A.matrixDiagonal->data();
    assert(AmatrixDiagonal);
    //
    const local_int_t numberOfNonzerosPerRow = Ageom->stencilSize;
    // Interpreted as 2D array
    assert(A.matrixValues->data());
    Array2D<floatType> matrixValues(
        localNumberOfRows, numberOfNonzerosPerRow, A.matrixValues->data()
    );
    // Interpreted as 2D array
    assert(A.mtxIndG->data());
    Array2D<global_int_t> mtxIndG(
        localNumberOfRows, numberOfNonzerosPerRow, A.mtxIndG->data()
    );
    const char *const AnonzerosInRow = A.nonzerosInRow->data();
    assert(AnonzerosInRow);
    //
    const global_int_t *const AlocalToGlobalMap = A.localToGlobalMap->data();
    assert(AlocalToGlobalMap);
    //
    local_int_t localNumberOfNonzeros = 0;
    //
    for (local_int_t iz = 0; iz < nz; iz++) {
        global_int_t giz = ipz * nz + iz;
        for (local_int_t iy = 0; iy < ny; iy++) {
            global_int_t giy = ipy * ny + iy;
            for (local_int_t ix = 0; ix < nx; ix++) {
                global_int_t gix = ipx * nx + ix;
                local_int_t currentLocalRow = iz * nx * ny + iy * nx + ix;
                global_int_t currentGlobalRow = giz * gnx * gny + giy * gnx + gix;
                assert(AlocalToGlobalMap[currentLocalRow] == currentGlobalRow);
                char numberOfNonzerosInRow = 0;
                // Current index in current row.
                global_int_t currentIndexG = 0;
                local_int_t currentNonZeroElemIndex = 0;
                for (int sz = -1; sz <= 1; sz++) {
                    if (giz + sz >- 1 && giz + sz < gnz) {
                        for (int sy = -1; sy <= 1; sy++) {
                            if (giy + sy > -1 && giy + sy < gny) {
                                for (int sx = -1; sx <= 1; sx++) {
                                    if (gix + sx > -1 && gix + sx < gnx) {
                                        global_int_t curcol = currentGlobalRow + sz * gnx * gny + sy * gnx + sx;
                                        if (curcol == currentGlobalRow) {
                                            assert(AmatrixDiagonal[currentLocalRow] == 26.0);
                                            assert(matrixValues(currentLocalRow, currentNonZeroElemIndex) == 26.0);
                                        }
                                        else {
                                            assert(matrixValues(currentLocalRow, currentNonZeroElemIndex) == -1.0);
                                        }
                                        currentNonZeroElemIndex++;
                                        assert(mtxIndG(currentLocalRow, currentIndexG++) == curcol);
                                        numberOfNonzerosInRow++;
                                    } // end x bounds test
                                } // end sx loop
                            } // end y bounds test
                        } // end sy loop
                    } // end z bounds test
                } // end sz loop
                assert(AnonzerosInRow[currentLocalRow] == numberOfNonzerosInRow);
                localNumberOfNonzeros += numberOfNonzerosInRow;
                if (b!=0)      assert(bv[currentLocalRow] == 26.0 - ((floatType)(numberOfNonzerosInRow-1)));
                if (x!=0)      assert(xv[currentLocalRow] == 0.0);
                if (xexact!=0) assert(xexactv[currentLocalRow] == 1.0);
            } // end ix loop
        } // end iy loop
    } // end iz loop
    //
    global_int_t totalNumberOfNonzeros = getTotalNumberOfNonZeros(
        A,
        localNumberOfNonzeros,
        ctx,
        runtime
    );
    //
    assert(Asclrs->totalNumberOfRows == totalNumberOfRows);
    assert(Asclrs->totalNumberOfNonzeros == totalNumberOfNonzeros);
    assert(Asclrs->localNumberOfRows == localNumberOfRows);
    assert(Asclrs->localNumberOfNonzeros == localNumberOfNonzeros);
}
