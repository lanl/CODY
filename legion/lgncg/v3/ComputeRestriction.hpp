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
    @file ComputeRestriction.hpp

    HPCG routine
 */

#pragma once

#include "LegionStuff.hpp"
#include "LegionArrays.hpp"
#include "LegionMatrices.hpp"

/*!
    Routine to compute the coarse residual vector.

    @param[inout]  A - Sparse matrix object containing pointers to mgData->Axf,
                   the fine grid matrix-vector product and mgData->rc the coarse
                   residual vector.

    @param[in]    rf - Fine grid RHS.


    Note that the fine grid residual is never explicitly constructed.  We only
    compute it for the fine grid points that will be injected into corresponding
    coarse grid points.

    @return Returns zero on success and a non-zero value otherwise.
*/
inline int
ComputeRestriction(
    SparseMatrix &A,
    Array<floatType> &rf,
    Context,
    Runtime *
) {
    const floatType *const Axfv = A.mgData->Axf->data();
    const floatType *const rfv = rf.data();
    floatType *const rcv = A.mgData->rc->data();
    const local_int_t *const f2c = A.mgData->f2cOperator->data();
    const local_int_t nc = A.mgData->rc->length();

    for (local_int_t i = 0; i < nc; ++i) {
        rcv[i] = rfv[f2c[i]] - Axfv[f2c[i]];
    }
    //
    return 0;
}
