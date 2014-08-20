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
 *
 * -- High Performance Conjugate Gradient Benchmark (HPCG)
 *    HPCG - 2.1 - January 31, 2014
 *
 *    Michael A. Heroux
 *    Scalable Algorithms Group, Computing Research Center
 *    Sandia National Laboratories, Albuquerque, NM
 *
 *    Piotr Luszczek
 *    Jack Dongarra
 *    University of Tennessee, Knoxville
 *    Innovative Computing Laboratory
 *    (C) Copyright 2013 All Rights Reserved
 *
 */

#ifndef LGNCG_SETUP_HALO_H_INCLUDED
#define LGNCG_SETUP_HALO_H_INCLUDED

#include "vector.h"
#include "sparsemat.h"
#include "utils.h"
#include "cg-task-args.h"
#include "tids.h"

#include "legion.h"

#include <map>

/**
 * Implements the halo setup. Well... Not really. Just updates indices for now.
 * This will implement halo setup at some point.
 */

namespace lgncg {

/**
 *
 */
static inline void
setupHalo(SparseMatrix &A,
          LegionRuntime::HighLevel::Context &ctx,
          LegionRuntime::HighLevel::HighLevelRuntime *lrt)
{
    using namespace LegionRuntime::HighLevel;
    int idx = 0;
    // setup per-task args
    ArgumentMap argMap;
    CGTaskArgs targs;
    targs.sa = A;
    for (int i = 0; i < A.vals.lDom().get_volume(); ++i) {
        targs.sa.nzir.sgb  = A.nzir.sgb()[i];
        targs.sa.mIdxs.sgb = A.mIdxs.sgb()[i];
        // notice the entire bounds
        targs.sa.g2g.sgb = A.g2g.bounds;
        argMap.set_point(DomainPoint::from_point<1>(Point<1>(i)),
                         TaskArgument(&targs, sizeof(targs)));
    }
    IndexLauncher il(LGNCG_SETUP_HALO_TID, A.nzir.lDom(),
                     TaskArgument(NULL, 0), argMap);
    // A's regions /////////////////////////////////////////////////////////////
    // # non 0s in row
    il.add_region_requirement(
        RegionRequirement(A.nzir.lp(), 0, READ_ONLY, EXCLUSIVE, A.nzir.lr)
    );
    il.add_field(idx++, A.nzir.fid);
    // mIdxs
    il.add_region_requirement(
        RegionRequirement(A.mIdxs.lp(), 0, READ_WRITE, EXCLUSIVE, A.mIdxs.lr)
    );
    il.add_field(idx++, A.mIdxs.fid);
    // global to global row lookup table (all of it)
    il.add_region_requirement(
        RegionRequirement(A.g2g.lr, 0, READ_ONLY, EXCLUSIVE, A.g2g.lr)
    );
    il.add_field(idx++, A.g2g.fid);
    // execute the thing...
    (void)lrt->execute_index_space(ctx, il);
}

/**
 *
 */
inline void
setupHaloTask(const LegionRuntime::HighLevel::Task *task,
              const std::vector<LegionRuntime::HighLevel::PhysicalRegion> &rgns,
              LegionRuntime::HighLevel::Context ctx,
              LegionRuntime::HighLevel::HighLevelRuntime *lrt)
{
    using namespace LegionRuntime::HighLevel;
    using namespace LegionRuntime::Accessor;
    using LegionRuntime::Arrays::Rect;
    // A (x3)
    assert(3 == rgns.size());
    size_t rid = 0;
    const CGTaskArgs targs = *(CGTaskArgs *)task->local_args;
    // name the regions
    const PhysicalRegion &azpr = rgns[rid++];
    const PhysicalRegion &aipr = rgns[rid++];
    const PhysicalRegion &ggpr = rgns[rid++];
    // convenience typedefs
    typedef RegionAccessor<AccessorType::Generic, uint8_t>  GSRA;
    typedef RegionAccessor<AccessorType::Generic, int64_t>  GLRA;
    typedef RegionAccessor<AccessorType::Generic, I64Tuple> GTRA;
    // sparse matrix
    GSRA az = azpr.get_field_accessor(targs.sa.nzir.fid).typeify<uint8_t>();
    GLRA ai = aipr.get_field_accessor(targs.sa.mIdxs.fid).typeify<int64_t>();
    GTRA gg = ggpr.get_field_accessor(targs.sa.g2g.fid).typeify<I64Tuple>();
    // primarily used for offset density tests
    Rect<1> sr; ByteOffset bOff[1];
    const Rect<1> mgb = targs.sa.nzir.sgb;
    // now start getting pointers to the data we are going to work with
    const uint8_t *const non0sInRow = az.raw_rect_ptr<1>(mgb, sr, bOff);
    bool offd = offsetsAreDense<1, uint8_t>(mgb, bOff);
    assert(offd);
    //
    int64_t *const mIdxs = ai.raw_rect_ptr<1>(targs.sa.mIdxs.sgb, sr, bOff);
    offd = offsetsAreDense<1, int64_t>(mgb, bOff);
    assert(offd);
    //
    const I64Tuple *const g2gMap = gg.raw_rect_ptr<1>(targs.sa.g2g.sgb, sr, bOff);
    offd = offsetsAreDense<1, I64Tuple>(mgb, bOff);
    assert(offd);
    ////////////////////////////////////////////////////////////////////////////
    // let the games begin...
    // stash local number of rows and columns
    const int64_t lNRows = mgb.volume();
    const int64_t lNCols = targs.sa.nCols;
    // now create a global to local map from the local to global map. we do this
    // primarily for fast lookups of local IDs given a global ID.
    // FIXME rename
    std::map<int64_t, int64_t> g2lMap;
    for (int64_t i = 0; i < targs.sa.g2g.sgb.volume(); ++i) {
        g2lMap[g2gMap[i].b] = g2gMap[i].a;
    }
    for (int64_t i = 0; i < lNRows; ++i) {
        const int64_t cNon0sInRow = non0sInRow[i];
        for (int64_t j = 0; j < cNon0sInRow; ++j) {
            int64_t *cIndxs = (mIdxs + (i * lNCols));
            const int64_t curIndex = cIndxs[j];
            // at [i][j]
            cIndxs[j] = g2lMap[curIndex];
        }
    }
}

} // end lgncg namespace

#endif
