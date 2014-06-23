/**
 * Copyright (c) 2014      Los Alamos National Security, LLC
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

#ifndef LGNCG_COMP_RESTRICTION_H_INCLUDED
#define LGNCG_COMP_RESTRICTION_H_INCLUDED

#include "vector.h"
#include "sparsemat.h"
#include "utils.h"
#include "cg-task-args.h"
#include "tids.h"

#include "legion.h"

/**
 * computes the coarse residual vector.
 */

namespace lgncg {

/**
 * responsible for setting up the task launch of the compute restriction
 * operation.
 */
static inline void
restriction(SparseMatrix &A,
            Vector &rf,
            LegionRuntime::HighLevel::Context &ctx,
            LegionRuntime::HighLevel::HighLevelRuntime *lrt)
{
#if 0
    using namespace LegionRuntime::HighLevel;
    int idx = 0;
    // no per-task arguments. CGTaskArgs has task-specific info.
    ArgumentMap argMap;
    // setup the global task args
    CGTaskArgs targs;
    targs.va = A.mgData->Axf;
    targs.vb = rf;
    targs.vc = A.mgData->rc;
    targs.vd = A.mgData->f2cOp;

    IndexLauncher il(LGNCG_RESTRICTION_TID, A.mgData->Axf.lDom(),
                     TaskArgument(&targs, sizeof(targs)), argMap);
    // A.mgData->Axf
    il.add_region_requirement(
        RegionRequirement(A.mgData->Axf.lp(), 0, READ_ONLY,
                          EXCLUSIVE, A.mgData->Axf.fid)
    );
    il.add_field(idx++, A.mgData->Axf.fid);
    // rf
    il.add_region_requirement(
        RegionRequirement(rf.lp(), 0, READ_ONLY, EXCLUSIVE, rf.lr)
    );
    il.add_field(idx++, rf.fid);
    // A.mgData->rc
    il.add_region_requirement(
        RegionRequirement(A.mgData->rc.lp(), 0, WRITE_DISCARD,
                          EXCLUSIVE, A.mgData->rc.lr)
    );
    il.add_field(idx++, A.mgData->rc.fid);
    // A.mgData->f2cOp
    // first push a partitioning scheme for the fine to coarse operator vector
    A.mgData->f2cOp.partition(rf.sgb().size(), ctx, lrt);
    il.add_region_requirement(
        RegionRequirement(A.mgData->f2cOp.lp(), 0, READ_ONLY,
                          EXCLUSIVE, A.mgData->f2cOp.lr)
    );
    il.add_field(idx++, A.mgData->f2cOp.fid);
    // execute the thing...
    (void)lrt->execute_index_space(ctx, il);
    // done with the f2c partition
    A.mgData->f2cOp.unpartition(ctx, lrt);
#endif
}

/**
 *
 */
inline void
restrictionTask(
    const LegionRuntime::HighLevel::Task *task,
    const std::vector<LegionRuntime::HighLevel::PhysicalRegion> &rgns,
    LegionRuntime::HighLevel::Context ctx,
    LegionRuntime::HighLevel::HighLevelRuntime *lrt
)
{
#if 0
    using namespace LegionRuntime::HighLevel;
    using namespace LegionRuntime::Accessor;
    using LegionRuntime::Arrays::Rect;
    // A.mgData->Axf, rf, A.mgData->rc, A.mgData->f2cOp
    assert(4 == rgns.size());
    size_t rid = 0;
    CGTaskArgs targs = *(CGTaskArgs *)task->args;
    // name the regions
    const PhysicalRegion &Axfpr = rgns[rid++];
    const PhysicalRegion &rfpr  = rgns[rid++];
    const PhysicalRegion &rcpr  = rgns[rid++];
    const PhysicalRegion &f2cpr = rgns[rid++];
    // convenience typedefs
    typedef RegionAccessor<AccessorType::Generic, double>  GDRA;
    typedef RegionAccessor<AccessorType::Generic, int64_t> GLRA;
    // vectors
    GDRA Axf = Axfpr.get_field_accessor(targs.va.fid).typeify<double>();
    GDRA rf  = rfpr.get_field_accessor(targs.vb.fid).typeify<double>();
    GDRA rc  = rcpr.get_field_accessor(targs.vc.fid).typeify<double>();
    GLRA f2c = f2cpr.get_field_accessor(targs.vd.fid).typeify<int64_t>();

    Rect<1> rec; ByteOffset boff[1];
    // Axf
    Rect<1> myGridBounds = targs.va.sgb()[getTaskID(task)];
    double *Axfp = Axf.raw_rect_ptr<1>(myGridBounds, rec, boff);
    bool offd = offsetsAreDense<1, double>(myGridBounds, boff);
    assert(offd);
    // rf
    myGridBounds = targs.vb.sgb()[getTaskID(task)];
    double *rfp = rf.raw_rect_ptr<1>(myGridBounds, rec, boff);
    offd = offsetsAreDense<1, double>(myGridBounds, boff);
    assert(offd);
    // rc
    myGridBounds = targs.vc.sgb()[getTaskID(task)];
    double *rcp = rc.raw_rect_ptr<1>(myGridBounds, rec, boff);
    offd = offsetsAreDense<1, double>(myGridBounds, boff);
    assert(offd);
    // f2c
    myGridBounds = targs.vd.sgb()[getTaskID(task)];
    int64_t *f2cp = f2c.raw_rect_ptr<1>(myGridBounds, rec, boff);
    offd = offsetsAreDense<1, int64_t>(myGridBounds, boff);
    assert(offd);
    // now, actually perform the computation
    int64_t nc = targs.vc.sgb()[getTaskID(task)].volume(); // length of rc
    for (int64_t i = 0; i < nc; ++i) {
        rcp[i] = rfp[f2cp[i]] - Axfp[f2cp[i]];
    }
#endif
}

}

#endif
