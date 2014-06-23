/**
 * Copyright (c) 2014      Los Alamos National Security, LLC
 *                         All rights reserved.
 */

#ifndef LGNCG_VEC_ZERO_H_INCLUDED
#define LGNCG_VEC_ZERO_H_INCLUDED

#include "vector.h"
#include "utils.h"
#include "cg-task-args.h"

#include "legion.h"

namespace lgncg {

static inline void
veczero(Vector &v,
        LegionRuntime::HighLevel::Context &ctx,
        LegionRuntime::HighLevel::HighLevelRuntime *lrt)
{
    using namespace LegionRuntime::HighLevel;
    int idx = 0;
    // setup per-task arguments. we pass each task its sub-grid bounds
    ArgumentMap argMap;
    for (int i = 0; i < v.lDom().get_volume(); ++i) {
        Rect<1> sgb = v.sgb()[i];
        argMap.set_point(DomainPoint::from_point<1>(Point<1>(i)),
                         TaskArgument(&sgb, sizeof(Rect<1>)));
    }
    IndexLauncher il(LGNCG_VEC_ZERO_TID, v.lDom(),
                     TaskArgument(&v, sizeof(Vector)), argMap);
    il.add_region_requirement(
        RegionRequirement(v.lp(), 0, WRITE_DISCARD, EXCLUSIVE, v.lr)
    );
    il.add_field(idx++, v.fid);
    // execute the thing...
    lrt->execute_index_space(ctx, il).wait_all_results();
}

/**
 * veczero task
 */
inline void
veczeroTask(const LegionRuntime::HighLevel::Task *task,
            const std::vector<LegionRuntime::HighLevel::PhysicalRegion> &rgns,
            LegionRuntime::HighLevel::Context ctx,
            LegionRuntime::HighLevel::HighLevelRuntime *lrt)
{
    using namespace LegionRuntime::HighLevel;
    using namespace LegionRuntime::Accessor;
    using LegionRuntime::Arrays::Rect;

    size_t rid = 0;
    // cache the arguments
    Vector vec = *(Vector *)task->args;
    Rect<1> rect = *(Rect<1> *)task->local_args;
    // name the region
    const PhysicalRegion &vpr = rgns[rid++];
    // convenience typedef
    typedef RegionAccessor<AccessorType::Generic, double> GDRA;
    GDRA v = vpr.get_field_accessor(vec.fid).typeify<double>();
    Rect<1> vSubRect;
    ByteOffset vOff[1];
    double *vp = v.raw_rect_ptr<1>(rect, vSubRect, vOff);
    bool offd = offsetsAreDense<1, double>(vSubRect, vOff);
    assert(offd);
    int64_t lLen = vSubRect.volume();
    for (int64_t i = 0; i < lLen; ++i) {
        vp[i] = 0.0;
    }
}

} // end lgncg namespace

#endif
