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
      @file ComputeSYMGS.hpp

      HPCG routine
 */

#pragma once

#include "LegionStuff.hpp"
#include "LegionArrays.hpp"
#include "LegionMatrices.hpp"
#include "ExchangeHalo.hpp"

#include <cassert>

/**
 *
 */
struct ComputeSYMGSArgs {
    local_int_t localNumberOfColumns;
    local_int_t localNumberOfRows;
    int stencilSize;
};

/*!
    Computes one step of symmetric Gauss-Seidel:

    Assumption about the structure of matrix A:
    - Each row 'i' of the matrix has nonzero diagonal value whose address is
      matrixDiagonal[i]
    - Entries in row 'i' are ordered such that:
         - Lower triangular terms are stored before the diagonal element.
         - Upper triangular terms are stored after the diagonal element.
         - No other assumptions are made about entry ordering.

    Symmetric Gauss-Seidel notes:
    - We use the input vector x as the RHS and start with an initial guess for y
    of all zeros.
    - We perform one forward sweep.  x should be initially zero on the first GS
      sweep, but we do not attempt to exploit this fact.
    - We then perform one back sweep.
    - For simplicity we include the diagonal contribution in the for-j loop,
      then correct the sum after

    @param[in] A the known system matrix.

    @param[in] r the input vector.

    @param[inout] x On entry, x should contain relevant values, on exit x
                  contains the result of one symmetric GS sweep with r as the
                  RHS.

    @warning Early versions of this kernel (Version 1.1 and earlier) had the r
    and x arguments in reverse order, and out of sync with other kernels.

    @return returns 0 upon success and non-zero otherwise.

    @see ComputeSYMGS
*/
inline int
ComputeSYMGSKernel(
    Array<floatType>       &AmatrixValues,
    Array<local_int_t>     &AmtxIndL,
    Array<char>            &AnonzerosInRow,
    Array<floatType>       &AmatrixDiagonal,
    Array<floatType>       &r,
    Array<floatType>       &x,
    const ComputeSYMGSArgs &args
) {
    // Make sure x contain space for halo values.
    assert(x.length() == size_t(args.localNumberOfColumns));
    //
    const local_int_t nrow = args.localNumberOfRows;
    const local_int_t nnpr = args.stencilSize;
    //
    const floatType *const matrixDiagonal = AmatrixDiagonal.data();
    assert(matrixDiagonal);
    //
    const floatType *const rv = r.data();
    assert(rv);
    floatType *const xv = x.data();
    assert(xv);
    // Interpreted as 2D array
    Array2D<floatType> matrixValues(
        nrow, nnpr, AmatrixValues.data()
    );
    // Interpreted as 2D array
    Array2D<local_int_t> mtxIndL(
        nrow, nnpr, AmtxIndL.data()
    );
    const char *const nonzerosInRow = AnonzerosInRow.data();
    //
    for (local_int_t i = 0; i < nrow; i++) {
        const floatType *const currentValues = matrixValues(i);
        const local_int_t *const currentColIndices = mtxIndL(i);
        const int currentNumberOfNonzeros = nonzerosInRow[i];
        const floatType currentDiagonal = matrixDiagonal[i];
        double sum = rv[i]; // RHS value
        //
        for (int j = 0; j < currentNumberOfNonzeros; j++) {
            const local_int_t curCol = currentColIndices[j];
            sum -= currentValues[j] * xv[curCol];
        }
        // Remove diagonal contribution from previous loop.
        sum += xv[i] * currentDiagonal;
        //
        xv[i] = sum / currentDiagonal;
    }
    // Now the back sweep.
    for (local_int_t i = nrow - 1; i >= 0; i--) {
        const floatType *const currentValues = matrixValues(i);
        const local_int_t *const currentColIndices = mtxIndL(i);
        const int currentNumberOfNonzeros = nonzerosInRow[i];
        const floatType currentDiagonal = matrixDiagonal[i];
        double sum = rv[i]; // RHS value
        //
        for (int j = 0; j < currentNumberOfNonzeros; j++) {
            const local_int_t curCol = currentColIndices[j];
            sum -= currentValues[j] * xv[curCol];
        }
        // Remove diagonal contribution from previous loop.
        sum += xv[i] * currentDiagonal;
        xv[i] = sum / currentDiagonal;
    }
    //
    return 0;
}

/**
 *
 */
inline int
ComputeSYMGS(
    SparseMatrix &A,
    Array<floatType> &r,
    Array<floatType> &x,
    Context ctx,
    Runtime *lrt
) {
    ExchangeHalo(A, x, ctx, lrt);
    //
    const ComputeSYMGSArgs args = {
        .localNumberOfColumns = A.sclrs->data()->localNumberOfColumns,
        .localNumberOfRows    = A.sclrs->data()->localNumberOfRows,
        .stencilSize          = A.geom->data()->stencilSize
    };
    //
#ifdef LGNCG_TASKING
    //
    TaskLauncher tl(
        SYMGS_TID,
        TaskArgument(&args, sizeof(args))
    );
    //
    A.matrixValues->intent  (RO_E, tl, ctx, lrt);
    A.mtxIndL->intent       (RO_E, tl, ctx, lrt);
    A.nonzerosInRow->intent (RO_E, tl, ctx, lrt);
    A.matrixDiagonal->intent(RO_E, tl, ctx, lrt);
    //
    r.intent(RO_E, tl, ctx, lrt);
    x.intent(RW_E, tl, ctx, lrt);
    //
    lrt->execute_task(ctx, tl);
    //
    return 0;
#else
    return ComputeSYMGSKernel(
               *A.matrixValues,
               *A.mtxIndL,
               *A.nonzerosInRow,
               *A.matrixDiagonal,
               r,
               x,
               args
           );
#endif
}

/**
 *
 */
void
ComputeSYMGSTask(
    const Task *task,
    const std::vector<PhysicalRegion> &regions,
    Context ctx,
    Runtime *lrt
) {
    const auto *const args = (ComputeSYMGSArgs *)task->args;
    //
    int rid = 0;
    Array<floatType> matrixValues  (regions[rid++], ctx, lrt);
    Array<local_int_t> mtxIndL     (regions[rid++], ctx, lrt);
    Array<char> nonzerosInRow      (regions[rid++], ctx, lrt);
    Array<floatType> matrixDiagonal(regions[rid++], ctx, lrt);
    //
    Array<floatType> r(regions[rid++], ctx, lrt);
    Array<floatType> x(regions[rid++], ctx, lrt);
    //
    ComputeSYMGSKernel(
        matrixValues,
        mtxIndL,
        nonzerosInRow,
        matrixDiagonal,
        r,
        x,
        *args
    );
}

/**
 *
 */
inline void
registerSYMGSTasks(void)
{
#ifdef LGNCG_TASKING
    HighLevelRuntime::register_legion_task<ComputeSYMGSTask>(
        SYMGS_TID /* task id */,
        Processor::LOC_PROC /* proc kind  */,
        true /* single */,
        true /* index */,
        AUTO_GENERATE_ID,
        TaskConfigOptions(true /* leaf task */),
        "ComputeSYMGSTask"
    );
#endif
}
