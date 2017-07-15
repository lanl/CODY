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

#pragma once

#include "Types.hpp"

#include "legion.h"

using namespace LegionRuntime::HighLevel;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template<typename TYPE>
struct DynColl {
    //
    int tid;
    //
    int64_t nArrivals = 0;
    //
    TYPE localBuffer;
    //
    DynamicCollective dc;

    /**
     *
     */
    DynColl(
        int tid,
        int64_t nArrivals
    ) : tid(tid)
      , nArrivals(nArrivals)
      , localBuffer(TYPE(0)) { }
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
class FloatReduceSumAccumulate {
public:
    typedef floatType LHS;
    typedef floatType RHS;
    static const floatType identity;

    template <bool EXCLUSIVE>
    static void apply(LHS &lhs, RHS rhs);

    template <bool EXCLUSIVE>
    static void fold(RHS &rhs1, RHS rhs2);
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
class FloatReduceMaxAccumulate {
public:
    typedef floatType LHS;
    typedef floatType RHS;
    static const floatType identity;

    template <bool EXCLUSIVE>
    static void apply(LHS &lhs, RHS rhs);

    template <bool EXCLUSIVE>
    static void fold(RHS &rhs1, RHS rhs2);
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
class IntReduceSumAccumulate {
public:
    typedef global_int_t LHS;
    typedef global_int_t RHS;
    static const global_int_t identity;

    template <bool EXCLUSIVE>
    static void apply(LHS &lhs, RHS rhs);

    template <bool EXCLUSIVE>
    static void fold(RHS &rhs1, RHS rhs2);
};

template <typename TYPE>
class Item;

/**
 * The type of DynColl passed in changes the behavior of the all reduce.
 */
floatType
allReduce(
    floatType localResult,
    Item< DynColl<floatType> > &dc,
    Context ctx,
    Runtime *runtime
);

/**
 * The type of DynColl passed in changes the behavior of the all reduce.
 */
global_int_t
allReduce(
    floatType localResult,
    Item< DynColl<global_int_t> > &dc,
    Context ctx,
    Runtime *runtime
);

/**
 *
 */
global_int_t
dynCollTaskContribGIT(
    const Task *task,
    const std::vector<PhysicalRegion> &regions,
    Context ctx,
    HighLevelRuntime *runtime
);

/**
 *
 */
floatType
dynCollTaskContribFT(
    const Task *task,
    const std::vector<PhysicalRegion> &regions,
    Context ctx, HighLevelRuntime *runtime
);

