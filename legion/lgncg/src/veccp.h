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

#ifndef LGNCG_VECCP_H_INCLUDED
#define LGNCG_VECCP_H_INCLUDED

#include "vector.h"
#include "utils.h"
#include "cg-task-args.h"

#include "legion.h"

namespace lgncg {

static inline void
veccp(const Vector &from,
      Vector &to,
      LegionRuntime::HighLevel::Context &ctx,
      LegionRuntime::HighLevel::HighLevelRuntime *lrt)
{
    using namespace LegionRuntime::HighLevel;
    int idx = 0;
    // sanity - make sure that both launch domains are the same size
    assert(from.lDom().get_volume() == to.lDom().get_volume());
    // and the # of sgbs
    assert(from.sgb().size() == to.sgb().size());
    // setup per-task arguments. we pass each task its sub-grid bounds
    CGTaskArgs targs;
    targs.va = from;
    targs.vb = to;
    ArgumentMap argMap;
    for (int i = 0; i < from.lDom().get_volume(); ++i) {
        targs.va.sgb = targs.vb.sgb = from.sgb()[i];
        argMap.set_point(DomainPoint::from_point<1>(Point<1>(i)),
                         TaskArgument(&targs, sizeof(targs)));
    }
    IndexLauncher il(LGNCG_VECCP_TID, from.lDom(),
                     TaskArgument(NULL, 0), argMap);
    // from
    il.add_region_requirement(
        RegionRequirement(from.lp(), 0, READ_ONLY, EXCLUSIVE, from.lr)
    );
    il.add_field(idx++, from.fid);
    // to
    il.add_region_requirement(
        RegionRequirement(to.lp(), 0, WRITE_DISCARD, EXCLUSIVE, to.lr)
    );
    il.add_field(idx++, to.fid);
    // execute the thing...
    (void)lrt->execute_index_space(ctx, il);
}

/**
 * veccp task
 */
inline void
veccpTask(const LegionRuntime::HighLevel::Task *task,
          const std::vector<LegionRuntime::HighLevel::PhysicalRegion> &rgns,
          LegionRuntime::HighLevel::Context ctx,
          LegionRuntime::HighLevel::HighLevelRuntime *lrt)
{
    using namespace LegionRuntime::HighLevel;
    using namespace LegionRuntime::Accessor;
    using LegionRuntime::Arrays::Rect;

    size_t rid = 0;
    CGTaskArgs targs = *(CGTaskArgs *)task->local_args;
#if 0 // nice debug
    printf("%d: sub-grid bounds: (%d) to (%d)\n",
            getTaskID(task), rect.lo.x[0], rect.hi.x[0]);
#endif
    // name the regions: from and to
    const PhysicalRegion &vf = rgns[rid++];
    const PhysicalRegion &vt = rgns[rid++];

    typedef RegionAccessor<AccessorType::Generic, double> GDRA;
    GDRA fm = vf.get_field_accessor(targs.va.fid).typeify<double>();
    GDRA to = vt.get_field_accessor(targs.vb.fid).typeify<double>();
#if 1
    Rect<1> fSubRect, tSubRect;
    ByteOffset fOff[1], tOff[1];
    const double *const fp = fm.raw_rect_ptr<1>(targs.va.sgb, fSubRect, fOff);
    double *tp = to.raw_rect_ptr<1>(targs.vb.sgb, tSubRect, tOff);
    // some sanity here...
    assert(fOff[0] == tOff[0]);
    assert(fSubRect == tSubRect);
    bool offd = offsetsAreDense<1, double>(fSubRect, fOff);
    assert(offd);
    const int64_t lLen = fSubRect.volume();
    // FIXME - use memmove
    for (int64_t i = 0; i < lLen; ++i) {
        tp[i] = fp[i];
    }
#else // slow reference implementation
    typedef GenericPointInRectIterator<1> GPRI1D;
    typedef DomainPoint DomPt;
    // start the real work
    for (GPRI1D pi(rect); pi; pi++) {
        to.write(DomPt::from_point<1>(pi.p),
                 fm.read(DomPt::from_point<1>(pi.p)));
    }
#endif
}

} // end lgncg namespace

#endif
