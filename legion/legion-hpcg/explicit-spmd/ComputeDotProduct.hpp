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

/**
 *
 */
struct ComputeDotProductArgs {
    local_int_t n;
};

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

inline int
ComputeDotProductKernel(
    Array<floatType> &x,
    Array<floatType> &y,
    const ComputeDotProductArgs &args,
    floatType &result
) {
    assert(x.length() >= size_t(args.n));
    assert(y.length() >= size_t(args.n));
    //
    floatType local_result = 0.0;
    //
    const floatType *const xv = x.data();
    assert(xv);
    //
    const floatType *const yv = y.data();
    assert(yv);
    //
    const local_int_t n = args.n;
    for (local_int_t i = 0; i < n; i++) local_result += xv[i] * yv[i];
    //
    result = local_result;
    //
    return 0;
}

/**
 *
 */
inline int
ComputeDotProduct(
    local_int_t n,
    Array<floatType> &x,
    Array<floatType> &y,
    Future &resultFuture,
    double &timeAllreduce, // FIXME
    Item< DynColl<floatType> > &dcReduceSum,
    Context ctx,
    Runtime *lrt
) {
    ComputeDotProductArgs args = {
        .n = n
    };
    //
    floatType localResult = 0.0;
    //
    int rc = 0;
#ifdef LGNCG_TASKING
    //
    TaskLauncher tl(
        DDOT_TID,
        TaskArgument(&args, sizeof(args))
    );
    //
    x.intent(RO_E, tl, ctx, lrt);
    y.intent(RO_E, tl, ctx, lrt);
    //
    auto f = lrt->execute_task(ctx, tl);
    // TODO pass future to allReduce not result.
    localResult = f.get_result<floatType>();
#else
    rc = ComputeDotProductKernel(x, y, args, localResult);
#endif
    double t0 = mytimer(); // FIXME
    resultFuture = allReduce(localResult, dcReduceSum, ctx, lrt);
    timeAllreduce += mytimer() - t0;
    //
    return rc;
}

/**
 *
 */
floatType
ComputeDotProductTask(
    const Task *task,
    const std::vector<PhysicalRegion> &regions,
    Context ctx,
    Runtime *lrt
) {
    const auto *const args = (ComputeDotProductArgs *)task->args;
    //
    Array<floatType> x(regions[0], ctx, lrt);
    Array<floatType> y(regions[1], ctx, lrt);
    //
    floatType localResult = 0.0;
    ComputeDotProductKernel(x, y, *args, localResult);
    //
    return localResult;
}

inline void
registerDDotTasks(void)
{
#ifdef LGNCG_TASKING
    HighLevelRuntime::register_legion_task<floatType, ComputeDotProductTask>(
        DDOT_TID /* task id */,
        Processor::LOC_PROC /* proc kind  */,
        true /* single */,
        true /* index */,
        AUTO_GENERATE_ID,
        TaskConfigOptions(true /* leaf task */),
        "ComputeDotProductTask"
    );
#endif
}
