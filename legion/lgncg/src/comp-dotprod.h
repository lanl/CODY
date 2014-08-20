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

#ifndef LGNCG_COMP_DOTPROD_H_INCLUDED
#define LGNCG_COMP_DOTPROD_H_INCLUDED

#include "tids.h"
#include "vector.h"
#include "utils.h"
#include "cg-task-args.h"

#include "legion.h"

namespace lgncg {

/**
 * responsible for setting up the task launch of the dot product task.
 */
static inline void
dotprod(Vector &x,
        Vector &y,
        double &result,
        LegionRuntime::HighLevel::Context &ctx,
        LegionRuntime::HighLevel::HighLevelRuntime *lrt)
{
    using namespace LegionRuntime::HighLevel;
    int idx = 0;
    // sanity - make sure that all launch domains are the same size
    assert(x.lDom().get_volume() == y.lDom().get_volume());
    // no per-task arguments. CGTaskArgs has task-specific info.
    ArgumentMap argMap;
    // setup the global task args
    CGTaskArgs targs;
    targs.va = x;
    targs.vb = y;
    for (int i = 0; i < x.lDom().get_volume(); ++i) {
        // FIXME optimization: don't really need to set sgb for all -- just for 1
        targs.va.sgb = x.sgb()[i];
        targs.vb.sgb = y.sgb()[i];
        argMap.set_point(DomainPoint::from_point<1>(Point<1>(i)),
                         TaskArgument(&targs, sizeof(targs)));
    }
    IndexLauncher il(LGNCG_DOTPROD_TID, x.lDom(),
                     TaskArgument(NULL, 0), argMap);
    // x's regions /////////////////////////////////////////////////////////////
    il.add_region_requirement(
        RegionRequirement(x.lp(), 0, READ_ONLY, EXCLUSIVE, x.lr)
    );
    il.add_field(idx++, x.fid);
    // y's regions /////////////////////////////////////////////////////////////
    il.add_region_requirement(
        RegionRequirement(y.lp(), 0, READ_ONLY, EXCLUSIVE, y.lr)
    );
    il.add_field(idx++, y.fid);
    // execute the thing...
    FutureMap fm = lrt->execute_index_space(ctx, il);
    // now sum up all the answers
    result = 0.0;
    typedef GenericPointInRectIterator<1> GPRI1D;
    for (GPRI1D pi(x.lDom().get_rect<1>()); pi; pi++) {
        Future f = fm.get_future(DomainPoint::from_point<1>(pi.p));
        result += f.get_result<double>();
    }
}

/**
 * computes the dot product of x' and y.
 * where: result = x' * y
 */
double
dotProdTask(const LegionRuntime::HighLevel::Task *task,
            const std::vector<LegionRuntime::HighLevel::PhysicalRegion> &rgns,
            LegionRuntime::HighLevel::Context ctx,
            LegionRuntime::HighLevel::HighLevelRuntime *lrt)
{
    using namespace LegionRuntime::HighLevel;
    using namespace LegionRuntime::Accessor;
    using LegionRuntime::Arrays::Rect;
    // x, y
    assert(2 == rgns.size());
    size_t rid = 0;
    CGTaskArgs targs = *(CGTaskArgs *)task->local_args;
#if 0 // nice debug
    printf("%d: sub-grid bounds: (%d) to (%d)\n",
            getTaskID(task), rect.lo.x[0], rect.hi.x[0]);
#endif
    // name the regions
    const PhysicalRegion &xpr = rgns[rid++];
    const PhysicalRegion &ypr = rgns[rid++];
    // convenience typedefs
    typedef RegionAccessor<AccessorType::Generic, double>  GDRA;
    // vectors
    GDRA x = xpr.get_field_accessor(targs.va.fid).typeify<double>();
    GDRA y = ypr.get_field_accessor(targs.vb.fid).typeify<double>();
    // x
    // this is the same for all vectors
    Rect<1> myGridBounds = targs.va.sgb;
    Rect<1> xsr; ByteOffset xOff[1];
    const double *const xp = x.raw_rect_ptr<1>(myGridBounds, xsr, xOff);
    bool offd = offsetsAreDense<1, double>(myGridBounds, xOff);
    assert(offd);
    // y
    Rect<1> ysr; ByteOffset yOff[1];
    const double *const yp = y.raw_rect_ptr<1>(myGridBounds, ysr, yOff);
    offd = offsetsAreDense<1, double>(myGridBounds, yOff);
    assert(offd);
    // now, actually perform the computation
    const int64_t lLen = myGridBounds.volume();
    double localRes = 0.0;
    for (int64_t i = 0; i < lLen; ++i) {
        localRes += xp[i] * yp[i];
    }
    return localRes;
}

} // end lgncg namespace

#endif
