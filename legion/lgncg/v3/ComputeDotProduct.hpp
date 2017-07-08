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

#include "LegionStuff.hpp"
#include "LegionArrays.hpp"
#include "mytimer.hpp"

#include <cassert>

/**
 * Tally total number of non-zeros in simulation.
 */
static inline floatType
allReduceSum(
    Item< DynColl<floatType> > &dc,
    LegionRuntime::HighLevel::Context ctx,
    LegionRuntime::HighLevel::Runtime *runtime
) {
    TaskLauncher tl(LOCAL_PARTIAL_SUM_TID, TaskArgument(NULL, 0));
    //
    tl.add_region_requirement(
        RegionRequirement(
            dc.logicalRegion,
            RO_E,
            dc.logicalRegion
        )
    ).add_field(dc.fid);
    //
    Future f = runtime->execute_task(ctx, tl);
    //
    DynamicCollective &dynCol = dc.data()->dc;
    runtime->defer_dynamic_collective_arrival(ctx, dynCol, f);
    dynCol  = runtime->advance_dynamic_collective(ctx, dynCol);
    //
    Future fSum = runtime->get_dynamic_collective_result(ctx, dynCol);
    //
    return fSum.get<floatType>();
}

/*!
    Routine to compute the dot product of two vectors where:

    This is the reference dot-product implementation.  It _CANNOT_ be modified
    for the purposes of this benchmark.

    @param[in] n the number of vector elements (on this processor)
    @param[in] x, y the input vectors
    @param[in] result a pointer to scalar value, on exit will contain result.
    @param[out] time_allreduce the time it took to perform the communication
    between processes

    @return returns 0 upon success and non-zero otherwise

    @see ComputeDotProduct
*/

inline int
ComputeDotProduct(
    local_int_t n,
    Array<floatType> &x,
    Array<floatType> &y,
    double &result,
    Item< DynColl<floatType> > &dcAllreduceSum,
    double &time_allreduce,
    LegionRuntime::HighLevel::Context ctx,
    LegionRuntime::HighLevel::Runtime *runtime
) {
    assert(x.length() >= size_t(n));
    assert(y.length() >= size_t(n));

    floatType &local_result = dcAllreduceSum.data()->localBuffer;
    local_result = 0.0;
    //
    floatType *xv = x.data();
    assert(xv);
    //
    floatType *yv = y.data();
    assert(yv);

    if (yv == xv) {
        for (local_int_t i = 0; i < n; i++) local_result += xv[i] * xv[i];
    }
    else {
        for (local_int_t i = 0; i < n; i++) local_result += xv[i] * yv[i];
    }
    // Collect all partial sums.
    floatType t0 = mytimer();
    result = allReduceSum(dcAllreduceSum, ctx, runtime);
    time_allreduce += mytimer() - t0;
    //
    return 0;
}
