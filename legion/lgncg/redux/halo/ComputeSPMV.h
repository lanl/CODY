/**
 * Copyright (c) 2016      Los Alamos National Security, LLC
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

#pragma once

#include <cassert>

/*!
    Routine to compute matrix vector product y = Ax where: Precondition: First
    call exchange_externals to get off-processor values of x

    This is the reference SPMV implementation.  It CANNOT be modified for the
    purposes of this benchmark.

    @param[in]  A the known system matrix
    @param[in]  x the known vector
    @param[out] y the On exit contains the result: Ax.

    @return returns 0 upon success and non-zero otherwise

    @see ComputeSPMV
*/
inline int
ComputeSPMV(
    const SparseMatrix &A,
    Vector &x,
    Vector &y
) {
    assert(x.localLength>=A.localNumberOfColumns); // Test vector lengths
    assert(y.localLength>=A.localNumberOfRows);

#if 0
    ExchangeHalo(A,x);
#endif
    const double *const xv = x.values;
    double *const yv = y.values;
    const local_int_t nrow = A.localNumberOfRows;
    //
    for (local_int_t i = 0; i < nrow; i++) {
        double sum = 0.0;
        const double *const cur_vals = A.matrixValues[i];
        const local_int_t *const cur_inds = A.mtxIndL[i];
        const int cur_nnz = A.nonzerosInRow[i];
        //
        for (int j = 0; j < cur_nnz; j++) {
            sum += cur_vals[j] * xv[cur_inds[j]];
        }
        yv[i] = sum;
    }
    return 0;
}
