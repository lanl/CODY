/**
 * Copyright (c) 2014      Los Alamos National Security, LLC
 *                         All rights reserved.
 */

#ifndef LGNCG_SPARSEMAT_H_INCLUDED
#define LGNCG_SPARSEMAT_H_INCLUDED

#include "vector.h"
#include "mg-data.h"

#include "legion.h"

namespace lgncg {

struct SparseMatrix {
    // total number of rows and columns
    int64_t nRows, nCols;
    // total number of non 0s
    int64_t tNon0;
    // current number of partitions
    int64_t nParts;
    // matrix values
    Vector vals;
    // matrix diagonal values
    Vector diag;
    // matrix idxs
    Vector mIdxs;
    // # non-0s in row
    Vector nzir;
    // FIXME - pointers here are probably not going to work... maybe???
    // coarse grid matrix. NULL indicates no next level.
    SparseMatrix *Ac;
    // multi-grid data. NULL indicates no MG data.
    MGData *mgData;

    SparseMatrix(void)
    {
        nRows = nCols = 0;
        nParts = 0;
        tNon0 = 0;
        Ac = NULL;
        mgData = NULL;
    }

    void
    create(int64_t nRows,
           int64_t nCols,
           LegionRuntime::HighLevel::Context &ctx,
           LegionRuntime::HighLevel::HighLevelRuntime *lrt)
    {
        using namespace LegionRuntime::HighLevel;
        using LegionRuntime::Arrays::Rect;

        this->tNon0 = 0;
        this->nParts = 0;
        this->nRows = nRows;
        this->nCols = nCols;
        this->vals.create<double>(nRows * nCols, ctx, lrt);
        this->diag.create<double>(nRows, ctx, lrt);
        this->mIdxs.create<int64_t>(nRows * nCols, ctx, lrt);
        this->nzir.create<uint8_t>(nRows, ctx, lrt);
    }

    void
    free(LegionRuntime::HighLevel::Context &ctx,
         LegionRuntime::HighLevel::HighLevelRuntime *lrt)
    {
        vals.free(ctx, lrt);
        diag.free(ctx, lrt);
        mIdxs.free(ctx, lrt);
        nzir.free(ctx, lrt);
        if (Ac) {
            Ac->free(ctx, lrt);
            delete Ac;
        }
        if (mgData) {
            mgData->free(ctx, lrt);
            delete mgData;
        }
    }

    void
    partition(int64_t nParts,
              LegionRuntime::HighLevel::Context &ctx,
              LegionRuntime::HighLevel::HighLevelRuntime *lrt)
    {
        this->nParts = nParts;
        vals.partition(nParts, ctx, lrt);
        diag.partition(nParts, ctx, lrt);
        mIdxs.partition(nParts, ctx, lrt);
        nzir.partition(nParts, ctx, lrt);
    }

    // TODO add unpartition
};

struct DSparseMatrix {
    // total number of rows and columns
    int64_t nRows, nCols;
    // total number of non 0s
    int64_t tNon0;
    // current number of partitions
    int64_t nParts;
    // matrix values
    DVector vals;
    // matrix diagonal values
    DVector diag;
    // matrix idxs
    DVector mIdxs;
    // # non-0s in row
    DVector nzir;

    void
    operator=(const SparseMatrix &rhs)
    {
        nRows = rhs.nRows;
        nCols = rhs.nCols;
        tNon0 = rhs.tNon0;
        nParts = rhs.nParts;
        vals = rhs.vals;
        diag = rhs.diag;
        mIdxs = rhs.mIdxs;
        nzir = rhs.nzir;
    }
};

} // end lgncg namespace

#endif
