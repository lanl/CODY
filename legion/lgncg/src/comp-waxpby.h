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

#ifndef LGNCG_COMP_WAXPBY_H_INCLUDED
#define LGNCG_COMP_WAXPBY_H_INCLUDED

#include "tids.h"
#include "vector.h"
#include "sparsemat.h"
#include "utils.h"
#include "cg-task-args.h"

#include "legion.h"

namespace {

struct waxpbyTaskArgs {
    double alpha;
    double beta;

    waxpbyTaskArgs(double alpha, double beta) : alpha(alpha), beta(beta) { ; }
};

}

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
    ArgumentMap argMap;
    waxpbyTaskArgs taskArgs(alpha, beta);
    IndexLauncher il(LGNCG_WAXPBY_TID, x.lDom(),
                     TaskArgument(&taskArgs, sizeof(taskArgs)), argMap);
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
    (void)ctx; (void)lrt;
    static const uint8_t xRID = 0;
    static const uint8_t yRID = 1;
    static const uint8_t wRID = 2;
    // x, y, w
    assert(3 == rgns.size());
    const waxpbyTaskArgs targs = *(waxpbyTaskArgs *)task->args;
    // name the regions
    const PhysicalRegion &xpr = rgns[xRID];
    const PhysicalRegion &ypr = rgns[yRID];
    const PhysicalRegion &wpr = rgns[wRID];
    // convenience typedefs
    typedef RegionAccessor<AccessorType::Generic, double>  GDRA;
    // vectors
    GDRA x = xpr.get_field_accessor(0).typeify<double>();
    const Domain xDom = lrt->get_index_space_domain(
        ctx, task->regions[xRID].region.get_index_space()
    );
    GDRA y = ypr.get_field_accessor(0).typeify<double>();
    GDRA w = wpr.get_field_accessor(0).typeify<double>();
    // this is the same for all vectors -- only do this once for x, y, and w
    Rect<1> myGridBounds = xDom.get_rect<1>();
    // x
    Rect<1> xsr; ByteOffset xOff[1];
    const double *const xp = x.raw_rect_ptr<1>(myGridBounds, xsr, xOff);
    bool offd = offsetsAreDense<1, double>(myGridBounds, xOff);
    assert(offd);
    // y
    Rect<1> ysr; ByteOffset yOff[1];
    const double *const yp = y.raw_rect_ptr<1>(myGridBounds, ysr, yOff);
    offd = offsetsAreDense<1, double>(myGridBounds, yOff);
    assert(offd);
    // w
    Rect<1> wsr; ByteOffset wOff[1];
    double *wp = w.raw_rect_ptr<1>(myGridBounds, wsr, wOff);
    offd = offsetsAreDense<1, double>(myGridBounds, wOff);
    assert(offd);
    // now, actually perform the computation
    const double alpha = targs.alpha;
    const double beta  = targs.beta;
    const int64_t lLen = myGridBounds.volume();
    for (int64_t i = 0; i < lLen; ++i) {
        wp[i] = (alpha * xp[i]) + (beta * yp[i]);
    }
}

} // end lgncg namespace

#endif

