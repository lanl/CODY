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
    @file ComputeProlongation.hpp

    HPCG routine
 */

#pragma once

#include "LegionStuff.hpp"
#include "LegionArrays.hpp"
#include "LegionMatrices.hpp"

/*!
    Routine to compute the coarse residual vector.

    @param[in]  Af - Fine grid sparse matrix object containing pointers to
                current coarse grid correction and the f2c operator.

    @param[inout] xf - Fine grid solution vector, update with coarse grid
                  correction.

    Note that the fine grid residual is never explicitly constructed.  We only
    compute it for the fine grid points that will be injected into corresponding
    coarse grid points.

    @return Returns zero on success and a non-zero value otherwise.
*/
inline int
ComputeProlongation(
    SparseMatrix &Af,
    Array<floatType> &xf,
    Context,
    Runtime *
) {
    floatType *const xfv = xf.data();
    const floatType *const xcv = Af.mgData->xc->data();
    const local_int_t *const f2c = Af.mgData->f2cOperator->data();
    const local_int_t nc = Af.mgData->rc->length();

    for (local_int_t i = 0; i < nc; ++i) {
        xfv[f2c[i]] += xcv[i];
    }
    //
    return 0;
}
