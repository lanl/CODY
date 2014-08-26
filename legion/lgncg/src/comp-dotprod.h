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
#include "dotprod-accumulate.h"

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
    Future fResult = lrt->execute_index_space(ctx, il, LGNCG_DOTPROD_RED_ID);
    // now sum up all the answers
    result = fResult.get_result<double>();
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
    (void)ctx; (void)lrt;
    static const uint8_t xRID = 0;
    static const uint8_t yRID = 1;
    // x, y
    assert(2 == rgns.size());
#if 0 // nice debug
    printf("%d: sub-grid bounds: (%d) to (%d)\n",
            getTaskID(task), rect.lo.x[0], rect.hi.x[0]);
#endif
    // name the regions
    const PhysicalRegion &xpr = rgns[xRID];
    const PhysicalRegion &ypr = rgns[yRID];
    // convenience typedefs
    typedef RegionAccessor<AccessorType::Generic, double>  GDRA;
    // vectors
    GDRA x = xpr.get_field_accessor(0).typeify<double>();
    const Domain xDom = lrt->get_index_space_domain(
        ctx, task->regions[xRID].region.get_index_space()
    );
    GDRA y = ypr.get_field_accessor(0).typeify<double>();
    // now, actually perform the computation
    double localRes = 0.0;
    for (GenericPointInRectIterator<1> itr(xDom.get_rect<1>()); itr; itr++) {
        localRes += x.read(DomainPoint::from_point<1>(itr.p)) *
                    y.read(DomainPoint::from_point<1>(itr.p));
    }
    return localRes;
}

} // end lgncg namespace

#endif
