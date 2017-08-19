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

/*!
    @file FutureMath.hpp
 */

#pragma once

#include "LegionStuff.hpp"

enum FutureMathOp {
    FMO_DIV,
    FMO_SQRT
};

/**
 *
 */
struct ComputeFutureArgs {
    FutureMathOp op;
};

/**
 *
 */
inline floatType
ComputeFutureKernel(
    Future *a,
    FutureMathOp op,
    Future *b
) {
    const floatType av = a->get_result<floatType>(disableWarnings);
    floatType bv = 0.0;
    //
    if (b) {
        bv = b->get_result<floatType>(disableWarnings);
    }
    //
    switch (op) {
        case FMO_DIV:  return  av / bv;
        case FMO_SQRT: return  sqrt(av);
        default: exit(1);
    }
    //
    return 0;
}

/**
 *
 */
inline Future
ComputeFuture(
    Future *a,
    FutureMathOp op,
    Future *b,
    Context ctx,
    Runtime *lrt
) {
#ifdef LGNCG_TASKING
    //
    ComputeFutureArgs args {
        .op = op
    };
    //
    TaskLauncher tl(
        FUTURE_MATH_TID,
        TaskArgument(&args, sizeof(args))
    );
    //
    tl.add_future(*a);
    if (b) tl.add_future(*b);
    //
    return lrt->execute_task(ctx, tl);
#else
    const floatType res = ComputeFutureKernel(a, op, b);
    return Future::from_value(lrt, res);
#endif
}

/**
 *
 */
floatType
ComputeFutureTask(
    const Task *task,
    const std::vector<PhysicalRegion> &regions,
    Context ctx,
    Runtime *lrt
) {
    const auto *const args = (ComputeFutureArgs *)task->args;
    //
    const int nf = task->futures.size();
    // Always at least one.
    Future af = task->futures[0];
    Future bf;
    if (nf == 2) {
        bf = task->futures[1];
    }
    //
    return ComputeFutureKernel(&af, args->op, nf == 2 ? &bf : NULL);
}

/**
 *
 */
inline void
registerFutureMathTasks(void)
{
#ifdef LGNCG_TASKING
    HighLevelRuntime::register_legion_task<floatType, ComputeFutureTask>(
        FUTURE_MATH_TID /* task id */,
        Processor::LOC_PROC /* proc kind  */,
        true /* single */,
        true /* index */,
        AUTO_GENERATE_ID,
        TaskConfigOptions(true /* leaf task */),
        "ComputeFutureTask"
    );
#endif
}
