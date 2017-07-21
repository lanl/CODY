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

#include <typeinfo>

#define LGNCG_MAX(x, y) x > y ? x : y
#define LGNCG_MIN(x, y) x < y ? x : y

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
const floatType FloatReduceSumAccumulate::identity = 1.0;

template<>
void
FloatReduceSumAccumulate::apply<true>(LHS &lhs, RHS rhs) {
    lhs += rhs;
}

template<>
void
FloatReduceSumAccumulate::apply<false>(LHS &lhs, RHS rhs) {
    exit(1);
}

template<>
void
FloatReduceSumAccumulate::fold<true>(RHS &rhs1, RHS rhs2) {
    exit(1);
}

template<>
void
FloatReduceSumAccumulate::fold<false>(RHS &rhs1, RHS rhs2) {
    exit(1);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
const floatType FloatReduceMaxAccumulate::identity = 1.0;

template<>
void
FloatReduceMaxAccumulate::apply<true>(LHS &lhs, RHS rhs) {
    lhs = LGNCG_MAX(lhs, rhs);
}

template<>
void
FloatReduceMaxAccumulate::apply<false>(LHS &lhs, RHS rhs) {
    exit(1);
}

template<>
void
FloatReduceMaxAccumulate::fold<true>(RHS &rhs1, RHS rhs2) {
    exit(1);
}

template<>
void
FloatReduceMaxAccumulate::fold<false>(RHS &rhs1, RHS rhs2) {
    exit(1);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
const floatType FloatReduceMinAccumulate::identity = 1.0;

template<>
void
FloatReduceMinAccumulate::apply<true>(LHS &lhs, RHS rhs) {
    lhs = LGNCG_MIN(lhs, rhs);
}

template<>
void
FloatReduceMinAccumulate::apply<false>(LHS &lhs, RHS rhs) {
    exit(1);
}

template<>
void
FloatReduceMinAccumulate::fold<true>(RHS &rhs1, RHS rhs2) {
    exit(1);
}

template<>
void
FloatReduceMinAccumulate::fold<false>(RHS &rhs1, RHS rhs2) {
    exit(1);
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
    exit(1);
}

template<>
void
IntReduceSumAccumulate::fold<true>(RHS &rhs1, RHS rhs2) {
    rhs1 += rhs2;
}

template<>
void
IntReduceSumAccumulate::fold<false>(RHS &rhs1, RHS rhs2) {
    exit(1);
}

/**
 *
 */
floatType
dynCollTaskContribFT(
    const Task *task,
    const std::vector<PhysicalRegion> &regions,
    Context ctx,
    Runtime *runtime
) {
    Future f = task->futures[0];
    return f.get<floatType>();
}

/**
 *
 */
global_int_t
dynCollTaskContribGIT(
    const Task *task,
    const std::vector<PhysicalRegion> &regions,
    Context ctx,
    Runtime *runtime
) {
    Future f = task->futures[0];
    return f.get<global_int_t>();
}

/**
 *
 */
void
registerCollectiveOpsTasks(void)
{
    HighLevelRuntime::register_legion_task<global_int_t, dynCollTaskContribGIT>(
        DYN_COLL_TASK_CONTRIB_GIT_TID /* task id */,
        Processor::LOC_PROC /* proc kind  */,
        true /* single */,
        false /* index */,
        AUTO_GENERATE_ID,
        TaskConfigOptions(true /* leaf task */),
        "dynCollTaskContribGIT"
    );
    HighLevelRuntime::register_legion_task<floatType, dynCollTaskContribFT>(
        DYN_COLL_TASK_CONTRIB_FT_TID /* task id */,
        Processor::LOC_PROC /* proc kind  */,
        true /* single */,
        false /* index */,
        AUTO_GENERATE_ID,
        TaskConfigOptions(true /* leaf task */),
        "dynCollTaskContribFT"
    );
    HighLevelRuntime::register_reduction_op<FloatReduceSumAccumulate>(
        FLOAT_REDUCE_SUM_TID
    );
    HighLevelRuntime::register_reduction_op<FloatReduceMinAccumulate>(
        FLOAT_REDUCE_MIN_TID
    );
    HighLevelRuntime::register_reduction_op<FloatReduceMaxAccumulate>(
        FLOAT_REDUCE_MAX_TID
    );
    HighLevelRuntime::register_reduction_op<IntReduceSumAccumulate>(
        INT_REDUCE_SUM_TID
    );
}
