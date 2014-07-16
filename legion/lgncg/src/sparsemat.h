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

#ifndef LGNCG_SPARSEMAT_H_INCLUDED
#define LGNCG_SPARSEMAT_H_INCLUDED

#include "geometry.h"
#include "vector.h"
#include "mg-data.h"

#include "legion.h"

namespace lgncg {

struct SparseMatrix {
    // geometry associated with this matrix
    Geometry geom;
    // total number of rows and columns
    int64_t nRows, nCols;
    // total number of non 0s
    int64_t tNon0;
    // matrix values
    Vector vals;
    // matrix diagonal values
    Vector diag;
    // matrix idxs
    Vector mIdxs;
    // # non-0s in row
    Vector nzir;
    // local to global lookup table. O(1) lookup when given a local ID
    // high 32 bits: global ID -- low 32 bits: local ID
    Vector l2g;
    // coarse grid matrix. NULL indicates no next level.
    SparseMatrix *Ac;
    // multi-grid data. NULL indicates no MG data.
    MGData *mgData;
    // current number of partitions XXX remove?
    int64_t nParts;
    /**
     *
     */
    SparseMatrix(void)
    {
        nRows = nCols = 0;
        nParts = 0;
        tNon0 = 0;
        Ac = NULL;
        mgData = NULL;
    }
    /**
     *
     */
    void
    create(const Geometry &geom,
           int64_t nRows,
           int64_t nCols,
           LegionRuntime::HighLevel::Context &ctx,
           LegionRuntime::HighLevel::HighLevelRuntime *lrt)
    {
        using namespace LegionRuntime::HighLevel;
        using LegionRuntime::Arrays::Rect;

        this->geom = geom;
        this->tNon0 = 0;
        this->nParts = 0;
        this->nRows = nRows;
        this->nCols = nCols;

        this->vals.create<double>  (nRows * nCols, ctx, lrt);
        this->mIdxs.create<int64_t>(nRows * nCols, ctx, lrt);

        this->diag.create<double> (nRows, ctx, lrt);
        this->nzir.create<uint8_t>(nRows, ctx, lrt);
        this->l2g.create<int64_t>(nRows, ctx, lrt);
    }
    /**
     *
     */
    void
    free(LegionRuntime::HighLevel::Context &ctx,
         LegionRuntime::HighLevel::HighLevelRuntime *lrt)
    {
        vals.free(ctx, lrt);
        diag.free(ctx, lrt);
        mIdxs.free(ctx, lrt);
        nzir.free(ctx, lrt);
        l2g.free(ctx, lrt);
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
        l2g.partition(nParts, ctx, lrt);
    }
    // TODO add unpartition
};

struct DSparseMatrix {
    // geometry associated with this matrix
    Geometry geom;
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
    // local to global lookup table. O(1) lookup when given a local ID
    // high 32 bits: global ID -- low 32 bits: local ID
    DVector l2g;

    void
    operator=(const SparseMatrix &rhs)
    {
        geom = rhs.geom;
        nRows = rhs.nRows;
        nCols = rhs.nCols;
        tNon0 = rhs.tNon0;
        nParts = rhs.nParts;
        vals = rhs.vals;
        diag = rhs.diag;
        mIdxs = rhs.mIdxs;
        nzir = rhs.nzir;
        l2g = rhs.l2g;
    }
};

} // end lgncg namespace

#endif
