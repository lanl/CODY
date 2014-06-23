/**
 * Copyright (c) 2014      Los Alamos National Security, LLC
 *                         All rights reserved.
 */

#ifndef LGNCG_MG_DATA_H_INCLUDED
#define LGNCG_MG_DATA_H_INCLUDED

#include "vector.h"
#include "sparsemat.h"

#include "legion.h"

namespace lgncg {

struct MGData {
    // call comp-symgs this many times prior to coarsening
    int64_t nPresmootherSteps;
    // call comp-symgs this many times after coarsening
    int64_t nPostsmootherSteps;
    // 1d array containing the fine operator local ids that will be injected
    // into coarse space.
    // length: number of rows in fine grid matrix
    Vector f2cOp;
    // coarse grid residual vector.
    // length: cnx * cny * cnz (coarse number of rows)
    Vector rc;
    // coarse grid solution vector
    // length: (coarse number of columns (not real number of columns -- e.g.
    // stencil size)). since we are dealing with square matrices, then that
    // corresponds to the coarse number of rows. so, cnx * cny * cnz
    Vector xc;
    // fine grid residual vector
    // length: fnx * fny * fnz (see: comments about xc -- similar thing here,
    // but for the fine grid matrix).
    Vector Axf;

    MGData(void)
    {
        nPresmootherSteps = 0;
        nPostsmootherSteps = 0;
    }

    MGData(int64_t nFineRows,
           int64_t nCoarseRows,
           LegionRuntime::HighLevel::Context &ctx,
           LegionRuntime::HighLevel::HighLevelRuntime *lrt)
    {
        nPresmootherSteps = 1;
        nPostsmootherSteps = 1;
        f2cOp.create<int64_t>(nFineRows, ctx, lrt);
        // about the size here: read comment above.
        Axf.create<double>(nFineRows, ctx, lrt);
        rc.create<double>(nCoarseRows, ctx, lrt);
        // about the size here: read comment above.
        xc.create<double>(nCoarseRows, ctx, lrt);
    }

    void
    free(LegionRuntime::HighLevel::Context &ctx,
         LegionRuntime::HighLevel::HighLevelRuntime *lrt)
    {
        f2cOp.free(ctx, lrt);
        rc.free(ctx, lrt);
        xc.free(ctx, lrt);
        Axf.free(ctx, lrt);
    }

    void
    partition(int64_t nParts,
              LegionRuntime::HighLevel::Context &ctx,
              LegionRuntime::HighLevel::HighLevelRuntime *lrt)
    {
        f2cOp.partition(nParts, ctx, lrt);
        rc.partition(nParts, ctx, lrt);
        xc.partition(nParts, ctx, lrt);
        Axf.partition(nParts, ctx, lrt);
    }
};

}

#endif
