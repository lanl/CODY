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

#ifndef LGNCG_COMP_PROLONGATION_H_INCLUDED
#define LGNCG_COMP_PROLONGATION_H_INCLUDED

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
 * responsible for setting up the task launch of the compute prolongation
 * operation.
 */
static inline void
prolongation(SparseMatrix &Af,
             Vector &xf,
             LegionRuntime::HighLevel::Context &ctx,
             LegionRuntime::HighLevel::HighLevelRuntime *lrt)
{
    using namespace LegionRuntime::HighLevel;
    int idx = 0;
    // no per-task arguments. CGTaskArgs has task-specific info.
    ArgumentMap argMap;
    // setup the global task args
    CGTaskArgs targs;
    targs.va = xf;
    targs.vb = Af.mgData->xc;
    targs.vc = Af.mgData->f2cOp;
    for (int i = 0; i < Af.vals.lDom().get_volume(); ++i) {
        targs.va.sgb = xf.sgb()[i];
        targs.vb.sgb = Af.mgData->xc.sgb()[i];
        targs.vc.sgb = Af.mgData->f2cOp.sgb()[i];
        argMap.set_point(DomainPoint::from_point<1>(Point<1>(i)),
                         TaskArgument(&targs, sizeof(targs)));
    }
    IndexLauncher il(LGNCG_PROLONGATION_TID, xf.lDom(),
                     TaskArgument(NULL, 0), argMap);
    // xf
    il.add_region_requirement(
        RegionRequirement(xf.lp(), 0, READ_WRITE, EXCLUSIVE, xf.lr)
    );
    il.add_field(idx++, targs.va.fid);
    // Af.mgData->xc
    il.add_region_requirement(
        RegionRequirement(Af.mgData->xc.lp(), 0, READ_ONLY,
                          EXCLUSIVE, Af.mgData->xc.lr)
    );
    il.add_field(idx++, targs.vb.fid);
    // Af.mgData->f2cOp
    il.add_region_requirement(
        RegionRequirement(Af.mgData->f2cOp.lp(), 0, READ_ONLY,
                          EXCLUSIVE, Af.mgData->f2cOp.lr)
    );
    il.add_field(idx++, targs.vc.fid);
    // execute the thing...
    (void)lrt->execute_index_space(ctx, il);
}

/**
 *
 */
inline void
prolongationTask(
    const LegionRuntime::HighLevel::Task *task,
    const std::vector<LegionRuntime::HighLevel::PhysicalRegion> &rgns,
    LegionRuntime::HighLevel::Context ctx,
    LegionRuntime::HighLevel::HighLevelRuntime *lrt
)
{
    using namespace LegionRuntime::HighLevel;
    using namespace LegionRuntime::Accessor;
    using LegionRuntime::Arrays::Rect;
    assert(3 == rgns.size());
    size_t rid = 0;
    CGTaskArgs targs = *(CGTaskArgs *)task->local_args;
    // name the regions
    const PhysicalRegion &xfpr  = rgns[rid++];
    const PhysicalRegion &xcpr  = rgns[rid++];
    const PhysicalRegion &f2cpr = rgns[rid++];
    // convenience typedefs
    typedef RegionAccessor<AccessorType::Generic, double>  GDRA;
    typedef RegionAccessor<AccessorType::Generic, int64_t> GLRA;
    // vectors
    GDRA xf  = xfpr.get_field_accessor(targs.va.fid).typeify<double>();
    GDRA xc  = xcpr.get_field_accessor(targs.vb.fid).typeify<double>();
    GLRA f2c = f2cpr.get_field_accessor(targs.vc.fid).typeify<int64_t>();

    Rect<1> rec; ByteOffset boff[1];
    // xf
    Rect<1> myGridBounds = targs.va.sgb;
    double *xfp = xf.raw_rect_ptr<1>(myGridBounds, rec, boff);
    bool offd = offsetsAreDense<1, double>(myGridBounds, boff);
    assert(offd);
    // xc
    myGridBounds = targs.vb.sgb;
    const double *const xcp = xc.raw_rect_ptr<1>(myGridBounds, rec, boff);
    offd = offsetsAreDense<1, double>(myGridBounds, boff);
    assert(offd);
    // f2c
    myGridBounds = targs.vc.sgb;
    const int64_t *const f2cp = f2c.raw_rect_ptr<1>(myGridBounds, rec, boff);
    offd = offsetsAreDense<1, int64_t>(myGridBounds, boff);
    assert(offd);
    // now, actually perform the computation
    const int64_t nc = targs.vb.sgb.volume(); // length of xc
    //
    for (int64_t i = 0; i < nc; ++i) {
        xfp[f2cp[i]] += xcp[i];
    }
}

}

#endif
