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
    ArgumentMap argMap;
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
    (void)ctx; (void)lrt;
    static const uint8_t fRID = 0;
    static const uint8_t tRID = 1;
    // name the regions: from and to
    const PhysicalRegion &vf = rgns[fRID];
    const PhysicalRegion &vt = rgns[tRID];
    typedef RegionAccessor<AccessorType::Generic, double> GDRA;
    GDRA fm = vf.get_field_accessor(0).typeify<double>();
    const Domain fDom = lrt->get_index_space_domain(
        ctx, task->regions[fRID].region.get_index_space()
    );
    GDRA to = vt.get_field_accessor(0).typeify<double>();
    for (GenericPointInRectIterator<1> itr(fDom.get_rect<1>()); itr; itr++) {
        to.write(
            DomainPoint::from_point<1>(itr.p),
            fm.read(DomainPoint::from_point<1>(itr.p))
        );
    }
}

} // end lgncg namespace

#endif
