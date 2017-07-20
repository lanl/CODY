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
 @file ComputeDotProduct_ref.cpp

 HPCG routine
 */

#pragma once

#include "LegionArrays.hpp"
#include "CollectiveOps.hpp"

#include "mytimer.hpp"

#include <cassert>

/*!
    Routine to compute the dot product of two vectors where:

    This is the reference dot-product implementation.  It _CANNOT_ be modified
    for the purposes of this benchmark.

    @param[in] n the number of vector elements (on this processor)
    @param[in] x, y the input vectors
    @param[in] result a pointer to scalar value, on exit will contain result.
    @param[out] timeAllreduce the time it took to perform the communication
    between processes

    @return returns 0 upon success and non-zero otherwise

    @see ComputeDotProduct
*/

#if 0
inline int
ComputeDotProductKernel(
    local_int_t n,
    Array<floatType> &x, Array<floatType> &y,
    floatType &result,
    double &timeAllreduce,
    Item< DynColl<floatType> > &dcReduceSum
) {
    assert(x.length() >= size_t(n));
    assert(y.length() >= size_t(n));
    //
    floatType local_result = 0.0;
    //
    const floatType *const xv = x.data();
    assert(xv);
    //
    const floatType *const yv = y.data();
    assert(yv);
    //
    if (yv == xv) {
        for (local_int_t i = 0; i < n; i++) local_result += xv[i] * xv[i];
    }
    else {
        for (local_int_t i = 0; i < n; i++) local_result += xv[i] * yv[i];
    }
    // Collect all partial sums.
    double t0 = mytimer();
    result = allReduce(local_result, dcReduceSum, ctx, runtime);
    timeAllreduce += mytimer() - t0;
    //
    return 0;
}
#endif

/**
 *
 */
inline int
ComputeDotProduct(
    local_int_t n,
    Array<floatType> &x,
    Array<floatType> &y,
    floatType &result,
    double &timeAllreduce,
    Item< DynColl<floatType> > &dcReduceSum,
    Context ctx,
    Runtime *runtime
) {
    assert(x.length() >= size_t(n));
    assert(y.length() >= size_t(n));
    //
    floatType local_result = 0.0;
    //
    const floatType *const xv = x.data();
    assert(xv);
    //
    const floatType *const yv = y.data();
    assert(yv);
    //
    if (yv == xv) {
        for (local_int_t i = 0; i < n; i++) local_result += xv[i] * xv[i];
    }
    else {
        for (local_int_t i = 0; i < n; i++) local_result += xv[i] * yv[i];
    }
    // Collect all partial sums.
    double t0 = mytimer();
    result = allReduce(local_result, dcReduceSum, ctx, runtime);
    timeAllreduce += mytimer() - t0;
    //
    return 0;
}

/**
 *
 */
void
ComputeDotProductTask(
    const Task *task,
    const std::vector<PhysicalRegion> &regions,
    Context ctx,
    Runtime *lrt
) {
}

inline void
registerDDotTasks(void)
{
#ifdef LGNCG_TASKING
    HighLevelRuntime::register_legion_task<ComputeWAXPBYTask>(
        WAXPBY_TID /* task id */,
        Processor::LOC_PROC /* proc kind  */,
        true /* single */,
        true /* index */,
        AUTO_GENERATE_ID,
        TaskConfigOptions(true /* leaf task */),
        "ComputeWAXPBYTask"
    );
#endif
}
