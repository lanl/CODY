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
#include "LegionItems.hpp"

#include "legion.h"

#include <cfloat>

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
    {
        mInitLocalBuffer(tid, localBuffer);
    }

private:

    /**
     *
     */
    void
    mInitLocalBuffer(
        int tid,
        global_int_t &lb
    ) {
        switch(tid) {
            case INT_REDUCE_SUM_TID:
                lb = TYPE(0);
                break;
            default:
                assert(false);
        }
    }

    /**
     *
     */
    void
    mInitLocalBuffer(
        int tid,
        floatType &lb
    ) {
        switch (tid) {
            case FLOAT_REDUCE_SUM_TID:
                lb = TYPE(0);
                break;
            case FLOAT_REDUCE_MIN_TID:
                lb = TYPE(DBL_MAX);
                break;
            case FLOAT_REDUCE_MAX_TID:
                lb = TYPE(-DBL_MAX);
                break;
            default:
                assert(false);
        }
    }
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
class FloatReduceMinAccumulate {
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

/**
 * The type of DynColl passed in changes the behavior of the all reduce.
 */
template <typename TYPE>
TYPE
allReduce(
    floatType localResult,
    Item< DynColl<TYPE> > &dc,
    Context ctx,
    Runtime *runtime
) {
    int tid = 0;
    //
    if (typeid(TYPE) == typeid(floatType)) {
        tid = DYN_COLL_TASK_CONTRIB_FT_TID;
    }
    else if (typeid(TYPE) == typeid(global_int_t)) {
        tid = DYN_COLL_TASK_CONTRIB_GIT_TID;
    }
    else {
        assert(false);
    }
    //
    TaskLauncher tl(tid, TaskArgument(NULL, 0));
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
    Future fres = runtime->get_dynamic_collective_result(ctx, dynCol);
    return fres.get<TYPE>();
}
