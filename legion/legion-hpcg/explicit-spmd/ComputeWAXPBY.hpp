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
    @file ComputeWAXPBY.cpp

    HPCG routine
 */

#pragma once

#include "LegionArrays.hpp"

#include <cassert>


/*!
    Routine to compute the update of a vector with the sum of two
    scaled vectors where: w = alpha*x + beta*y

    This is the reference WAXPBY impmentation.  It CANNOT be modified for the
    purposes of this benchmark.

    @param[in] n the number of vector elements (on this processor)
    @param[in] alpha, beta the scalars applied to x and y respectively.
    @param[in] x, y the input vectors
    @param[out] w the output vector.

    @return returns 0 upon success and non-zero otherwise

    @see ComputeWAXPBY
*/

inline int
ComputeWAXPBY(
    const local_int_t n,
    const floatType alpha,
    const Array<floatType> &x,
    const floatType beta,
    const Array<floatType> &y,
    Array<floatType> &w,
    Context,
    Runtime *
) {
    // Test vector lengths
    assert(x.length() >= size_t(n));
    assert(y.length() >= size_t(n));

    const floatType *const xv = x.data();
    const floatType *const yv = y.data();
    floatType *const       wv = w.data();

    if (alpha == 1.0) {
        for (local_int_t i = 0; i < n; i++) wv[i] = xv[i] + beta * yv[i];
    }
    else if (beta == 1.0) {
        for (local_int_t i = 0; i < n; i++) wv[i] = alpha * xv[i] + yv[i];
    }
    else  {
        for (local_int_t i = 0; i < n; i++) wv[i] = alpha * xv[i] + beta * yv[i];
    }
    //
    return 0;
}
