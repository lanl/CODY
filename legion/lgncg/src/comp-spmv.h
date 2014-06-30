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

#ifndef LGNCG_COMP_SPMV_H_INCLUDED
#define LGNCG_COMP_SPMV_H_INCLUDED

#include "vector.h"
#include "sparsemat.h"
#include "utils.h"
#include "cg-task-args.h"
#include "tids.h"

#include "legion.h"

/**
 * implements the SParse Matrix Vector multiply functionality.
 */

namespace lgncg {

/**
 * responsible for setting up the task launch of the sparse matrix vector
 * multiply.
 */
static inline void
spmv(const SparseMatrix &A,
     const Vector &x,
     Vector &y,
     LegionRuntime::HighLevel::Context &ctx,
     LegionRuntime::HighLevel::HighLevelRuntime *lrt)
{
    using namespace LegionRuntime::HighLevel;
    int idx = 0;
    // sanity - make sure that all launch domains are the same size
    assert(A.vals.lDom().get_volume() == x.lDom().get_volume() &&
           x.lDom().get_volume() == y.lDom().get_volume());
    // setup per-task args
    ArgumentMap argMap;
    CGTaskArgs targs;
    targs.sa = A;
    targs.va = x;
    targs.vb = y;
    for (int i = 0; i < A.vals.lDom().get_volume(); ++i) {
        targs.sa.vals.sgb = A.vals.sgb()[i];
        targs.sa.mIdxs.sgb = A.mIdxs.sgb()[i];
        targs.sa.nzir.sgb = A.nzir.sgb()[i];
        targs.va.sgb = x.sgb()[i];
        targs.vb.sgb = y.sgb()[i];
        argMap.set_point(DomainPoint::from_point<1>(Point<1>(i)),
                         TaskArgument(&targs, sizeof(targs)));
    }
    IndexLauncher il(LGNCG_SPMV_TID, A.vals.lDom(),
                     TaskArgument(NULL, 0), argMap);
    // A's regions /////////////////////////////////////////////////////////////
    // vals
    il.add_region_requirement(
        RegionRequirement(A.vals.lp(), 0, READ_ONLY, EXCLUSIVE, A.vals.lr)
    );
    il.add_field(idx++, A.vals.fid);
    // mIdxs
    il.add_region_requirement(
        RegionRequirement(A.mIdxs.lp(), 0, READ_ONLY, EXCLUSIVE, A.mIdxs.lr)
    );
    il.add_field(idx++, A.mIdxs.fid);
    // nzir
    il.add_region_requirement(
        RegionRequirement(A.nzir.lp(), 0, READ_ONLY, EXCLUSIVE, A.nzir.lr)
    );
    il.add_field(idx++, A.nzir.fid);
    // x's regions /////////////////////////////////////////////////////////////
    il.add_region_requirement(
        /* notice we are using the entire region here */
        //RegionRequirement(x.lr, 0, READ_ONLY, EXCLUSIVE, x.lr)
        // XXX what's the diff?
        RegionRequirement(x.lr, READ_ONLY, EXCLUSIVE, x.lr)
    );
    il.add_field(idx++, x.fid);
    // y's regions /////////////////////////////////////////////////////////////
    il.add_region_requirement(
        RegionRequirement(y.lp(), 0, WRITE_DISCARD, EXCLUSIVE, y.lr)
    );
    il.add_field(idx++, y.fid);
    // execute the thing...
    (void)lrt->execute_index_space(ctx, il);
}

/**
 * computes: y = A * x ; where y and x are vectors and A is a sparse matrix
 */
inline void
spmvTask(const LegionRuntime::HighLevel::Task *task,
          const std::vector<LegionRuntime::HighLevel::PhysicalRegion> &rgns,
          LegionRuntime::HighLevel::Context ctx,
          LegionRuntime::HighLevel::HighLevelRuntime *lrt)
{
    using namespace LegionRuntime::HighLevel;
    using namespace LegionRuntime::Accessor;
    using LegionRuntime::Arrays::Rect;

    // A (x3), x, b
    assert(5 == rgns.size());
    size_t rid = 0;
    CGTaskArgs targs = *(CGTaskArgs *)task->local_args;
#if 0 // nice debug
    printf("%d: sub-grid bounds: (%d) to (%d)\n",
            getTaskID(task), rect.lo.x[0], rect.hi.x[0]);
#endif
    // name the regions
    const PhysicalRegion &avpr = rgns[rid++];
    const PhysicalRegion &aipr = rgns[rid++];
    const PhysicalRegion &azpr = rgns[rid++];
    const PhysicalRegion &xpr  = rgns[rid++];
    const PhysicalRegion &ypr  = rgns[rid++];
    // convenience typedefs
    typedef RegionAccessor<AccessorType::Generic, double>  GDRA;
    typedef RegionAccessor<AccessorType::Generic, int64_t> GLRA;
    typedef RegionAccessor<AccessorType::Generic, uint8_t> GSRA;
    // sparse matrix
    GDRA av = avpr.get_field_accessor(targs.sa.vals.fid).typeify<double>();
    GLRA ai = aipr.get_field_accessor(targs.sa.mIdxs.fid).typeify<int64_t>();
    GSRA az = azpr.get_field_accessor(targs.sa.nzir.fid).typeify<uint8_t>();
    // vectors
    GDRA x = xpr.get_field_accessor(targs.va.fid).typeify<double>();
    GDRA y = ypr.get_field_accessor(targs.vb.fid).typeify<double>();

    Rect<1> avsr; ByteOffset avOff[1];
    Rect<1> myGridBounds = targs.sa.vals.sgb;
    // calculate nRows and nCols for the local subgrid
    assert(0 == myGridBounds.volume() % targs.sa.nCols);
    int64_t lNRows = myGridBounds.volume() / targs.sa.nCols;
    int64_t lNCols = targs.sa.nCols;
    double *avp = av.raw_rect_ptr<1>(myGridBounds, avsr, avOff);
    bool offd = offsetsAreDense<1, double>(myGridBounds, avOff);
    assert(offd);
    // remember that vals and mIdxs should be the same size
    Rect<1> aisr; ByteOffset aiOff[1];
    int64_t *aip = ai.raw_rect_ptr<1>(myGridBounds, aisr, aiOff);
    offd = offsetsAreDense<1, int64_t>(myGridBounds, aiOff);
    assert(offd);
    // nzir is smaller (by a stencil size factor).
    Rect<1> azsr; ByteOffset azOff[1];
    myGridBounds = targs.sa.nzir.sgb;
    uint8_t *azp = az.raw_rect_ptr<1>(myGridBounds, azsr, azOff);
    offd = offsetsAreDense<1, uint8_t>(myGridBounds, azOff);
    assert(offd);
    // x
    Rect<1> xsr; ByteOffset xOff[1];
    // notice that we aren't using the subgridBounds here -- need all of x
    myGridBounds = targs.va.bounds;
    double *xp = x.raw_rect_ptr<1>(myGridBounds, xsr, xOff);
    offd = offsetsAreDense<1, double>(myGridBounds, xOff);
    assert(offd);
    // y
    Rect<1> ysr; ByteOffset yOff[1];
    myGridBounds = targs.vb.sgb;
    double *yp = y.raw_rect_ptr<1>(myGridBounds, ysr, yOff);
    offd = offsetsAreDense<1, double>(myGridBounds, yOff);
    assert(offd);
    // now, actually perform the computation
    for (int64_t i = 0; i < lNRows; ++i) {
        double sum = 0.0;
        // get to base of next row of values
        const double *const cVals = (avp + (i * lNCols));
        // get to base of next row of "real" indices of values
        const int64_t *const cIndx = (aip + (i * lNCols));
        // capture how many non-zero values are in this particular row
        int64_t cnnz = azp[i];
        for (int64_t j = 0; j < cnnz; ++j) {
            sum += cVals[j] * xp[cIndx[j]];
        }
        yp[i] = sum;
    }
}

}

#endif
