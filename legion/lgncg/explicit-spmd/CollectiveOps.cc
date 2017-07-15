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

#include "CollectiveOps.hpp"

#include "LegionStuff.hpp"
#include "LegionItems.hpp"

#define MAX(x, y) x > y ? x : y
#define MIN(x, y) x < y ? x : y

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
const floatType FloatReduceSumAccumulate::identity = 0.0;

template<>
void
FloatReduceSumAccumulate::apply<true>(LHS &lhs, RHS rhs) {
    lhs += rhs;
}

template<>
void
FloatReduceSumAccumulate::apply<false>(LHS &lhs, RHS rhs) {
    assert(false);
}

template<>
void
FloatReduceSumAccumulate::fold<true>(RHS &rhs1, RHS rhs2) {
    rhs1 += rhs2;
}

template<>
void
FloatReduceSumAccumulate::fold<false>(RHS &rhs1, RHS rhs2) {
    assert(false);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
const floatType FloatReduceMaxAccumulate::identity = 0.0;

template<>
void
FloatReduceMaxAccumulate::apply<true>(LHS &lhs, RHS rhs) {
    lhs = MAX(lhs, rhs);
}

template<>
void
FloatReduceMaxAccumulate::apply<false>(LHS &lhs, RHS rhs) {
    assert(false);
}

template<>
void
FloatReduceMaxAccumulate::fold<true>(RHS &rhs1, RHS rhs2) {
    rhs1 = MAX(rhs1, rhs2);
}

template<>
void
FloatReduceMaxAccumulate::fold<false>(RHS &rhs1, RHS rhs2) {
    assert(false);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
const global_int_t IntReduceSumAccumulate::identity = 0;

template<>
void
IntReduceSumAccumulate::apply<true>(LHS &lhs, RHS rhs) {
    lhs += rhs;
}

template<>
void
IntReduceSumAccumulate::apply<false>(LHS &lhs, RHS rhs) {
    assert(false);
}

template<>
void
IntReduceSumAccumulate::fold<true>(RHS &rhs1, RHS rhs2) {
    rhs1 += rhs2;
}

template<>
void
IntReduceSumAccumulate::fold<false>(RHS &rhs1, RHS rhs2) {
    assert(false);
}

/**
 * The type of DynColl passed in changes the behavior of the all reduce.
 */
floatType
allReduce(
    floatType localResult,
    Item< DynColl<floatType> > &dc,
    Context ctx,
    Runtime *runtime
) {
    TaskLauncher tl(DYN_COLL_TASK_CONTRIB_FT_TID, TaskArgument(NULL, 0));
    //
    tl.add_region_requirement(
        RegionRequirement(
            dc.logicalRegion,
            RO_E,
            dc.logicalRegion
        )
    ).add_field(dc.fid);
    //
    dc.data()->localBuffer = localResult;
    //
    Future f = runtime->execute_task(ctx, tl);
    //
    DynamicCollective &dynCol = dc.data()->dc;
    //
    runtime->defer_dynamic_collective_arrival(ctx, dynCol, f);
    dynCol = runtime->advance_dynamic_collective(ctx, dynCol);
    //
    Future fSum = runtime->get_dynamic_collective_result(ctx, dynCol);
    //
    return fSum.get<floatType>();
}
