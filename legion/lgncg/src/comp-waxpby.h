/**
 * Copyright (c) 2014      Los Alamos National Security, LLC
 *                         All rights reserved.
 */

#ifndef LGNCG_COMP_WAXPBY_H_INCLUDED
#define LGNCG_COMP_WAXPBY_H_INCLUDED

#include "tids.h"
#include "vector.h"
#include "sparsemat.h"
#include "utils.h"
#include "cg-task-args.h"

#include "legion.h"

namespace lgncg {

/**
 * responsible for setting up the task launch of the waxpby.
 */
static inline void
waxpby(double alpha,
       Vector &x,
       double beta,
       Vector &y,
       Vector &w,
       LegionRuntime::HighLevel::Context &ctx,
       LegionRuntime::HighLevel::HighLevelRuntime *lrt)
{
    using namespace LegionRuntime::HighLevel;
    int idx = 0;
    // sanity - make sure that all launch domains are the same size
    assert(x.lDom().get_volume() == y.lDom().get_volume() &&
           y.lDom().get_volume() == w.lDom().get_volume());
    // no per-task arguments. CGTaskArgs has task-specific info.
    ArgumentMap argMap;
    // setup the global task args
    CGTaskArgs targs;
    targs.alpha = alpha;
    targs.beta  = beta;
    targs.va = x;
    targs.vb = y;
    targs.vc = w;
    for (int i = 0; i < x.lDom().get_volume(); ++i) {
        // FIXME optimization: don't really need to set sgb for all -- just for 1
        targs.va.sgb = x.sgb()[i];
        targs.vb.sgb = y.sgb()[i];
        targs.vc.sgb = w.sgb()[i];
        argMap.set_point(DomainPoint::from_point<1>(Point<1>(i)),
                         TaskArgument(&targs, sizeof(targs)));
    }
    IndexLauncher il(LGNCG_WAXPBY_TID, x.lDom(), TaskArgument(NULL, 0), argMap);
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
    // w's regions /////////////////////////////////////////////////////////////
    il.add_region_requirement(
        RegionRequirement(w.lp(), 0, WRITE_DISCARD, EXCLUSIVE, w.lr)
    );
    il.add_field(idx++, w.fid);
    // execute the thing...
    (void)lrt->execute_index_space(ctx, il);
}

/**
 * computes the update of a vector with the sum of two scaled vectors
 * where: w = alpha * x + beta * y
 */
inline void
waxpbyTask(const LegionRuntime::HighLevel::Task *task,
           const std::vector<LegionRuntime::HighLevel::PhysicalRegion> &rgns,
           LegionRuntime::HighLevel::Context ctx,
           LegionRuntime::HighLevel::HighLevelRuntime *lrt)
{
    using namespace LegionRuntime::HighLevel;
    using namespace LegionRuntime::Accessor;
    using LegionRuntime::Arrays::Rect;
    // x, y, w
    assert(3 == rgns.size());
    size_t rid = 0;
    CGTaskArgs targs = *(CGTaskArgs *)task->local_args;
#if 0 // nice debug
    printf("%d: sub-grid bounds: (%d) to (%d)\n",
            getTaskID(task), rect.lo.x[0], rect.hi.x[0]);
#endif
    // name the regions
    const PhysicalRegion &xpr  = rgns[rid++];
    const PhysicalRegion &ypr  = rgns[rid++];
    const PhysicalRegion &wpr  = rgns[rid++];
    // convenience typedefs
    typedef RegionAccessor<AccessorType::Generic, double>  GDRA;
    // vectors
    GDRA x = xpr.get_field_accessor(targs.va.fid).typeify<double>();
    GDRA y = ypr.get_field_accessor(targs.vb.fid).typeify<double>();
    GDRA w = wpr.get_field_accessor(targs.vc.fid).typeify<double>();
    // x
    // this is the same for all vectors -- only do this once for x, y, and w
    Rect<1> myGridBounds = targs.va.sgb;
    Rect<1> xsr; ByteOffset xOff[1];
    double *xp = x.raw_rect_ptr<1>(myGridBounds, xsr, xOff);
    bool offd = offsetsAreDense<1, double>(myGridBounds, xOff);
    assert(offd);
    // y
    Rect<1> ysr; ByteOffset yOff[1];
    double *yp = y.raw_rect_ptr<1>(myGridBounds, ysr, yOff);
    offd = offsetsAreDense<1, double>(myGridBounds, yOff);
    assert(offd);
    // w
    Rect<1> wsr; ByteOffset wOff[1];
    double *wp = w.raw_rect_ptr<1>(myGridBounds, wsr, wOff);
    offd = offsetsAreDense<1, double>(myGridBounds, wOff);
    assert(offd);
    // now, actually perform the computation
    double alpha = targs.alpha;
    double beta  = targs.beta;
    int64_t lLen = myGridBounds.volume();
    for (int64_t i = 0; i < lLen; ++i) {
        wp[i] = (alpha * xp[i]) + (beta * yp[i]);
    }
}

} // end lgncg namespace

#endif
