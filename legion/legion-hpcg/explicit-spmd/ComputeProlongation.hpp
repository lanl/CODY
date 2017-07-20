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

/**
 *
 */
struct ComputeProlongationArgs {
    local_int_t nc;
};

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
ComputeProlongationKernel(
    Array<floatType>              &Afxc,
    Array<local_int_t>            &Aff2cOperator,
    Array<floatType>              &xf,
    const ComputeProlongationArgs &args
) {
    const floatType *const xcv = Afxc.data();
    assert(xcv);
    const local_int_t *const f2c = Aff2cOperator.data();
    assert(f2c);
    //
    floatType *const xfv = xf.data();
    assert(xfv);
    //
    const local_int_t nc = args.nc;

    for (local_int_t i = 0; i < nc; ++i) {
        xfv[f2c[i]] += xcv[i];
    }
    //
    return 0;
}

/**
 *
 */
inline int
ComputeProlongation(
    SparseMatrix &Af,
    Array<floatType> &xf,
    Context ctx,
    Runtime *lrt
) {
    const ComputeProlongationArgs args = {
        .nc = local_int_t(Af.mgData->rc->length())
    };
#ifdef LGNCG_TASKING
    //
    TaskLauncher tl(
        PROLONGATION_TID,
        TaskArgument(&args, sizeof(args))
    );
    //
    Af.mgData->xc->intent         (RO_E, tl, ctx, lrt);
    Af.mgData->f2cOperator->intent(RO_E, tl, ctx, lrt);
    //
    xf.intent(WO_E, tl, ctx, lrt);
    //
    auto f = lrt->execute_task(ctx, tl);
    f.wait(); // TODO RM
    return 0;
#else
    return ComputeProlongationKernel(
               *Af.mgData->xc,
               *Af.mgData->f2cOperator,
               xf,
               args
           );
#endif
}

/**
 *
 */
void
ComputeProlongationTask(
    const Task *task,
    const std::vector<PhysicalRegion> &regions,
    Context ctx,
    Runtime *lrt
) {
    const auto *const args = (ComputeProlongationArgs *)task->args;
    //
    int rid = 0;
    Array<floatType>   Afxc (regions[rid++], ctx, lrt);
    Array<local_int_t> Aff2c(regions[rid++], ctx, lrt);
    //
    Array<floatType> xf(regions[rid++], ctx, lrt);
    //
    ComputeProlongationKernel(
        Afxc,
        Aff2c,
        xf,
        *args
    );
}

/**
 *
 */
inline void
registerProlongationTasks(void)
{
#ifdef LGNCG_TASKING
    HighLevelRuntime::register_legion_task<ComputeProlongationTask>(
        PROLONGATION_TID /* task id */,
        Processor::LOC_PROC /* proc kind  */,
        true /* single */,
        true /* index */,
        AUTO_GENERATE_ID,
        TaskConfigOptions(true /* leaf task */),
        "ComputeProlongationTask"
    );
#endif
}
