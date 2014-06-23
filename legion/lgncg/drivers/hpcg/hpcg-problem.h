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

#ifndef HPCG_PROBLEM_H_INCLUDED
#define HPCG_PROBLEM_H_INCLUDED

#include "hpcg-geometry.h"

#include "lgncg.h"

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <inttypes.h>

// ideally, we would have some central TID registration station. for now, just
// assume lgncg won't use TID 0 and 1.

struct Problem {
    // problem geometry (at finest level)
    Geometry geom;
    // stencil size
    int64_t stencilSize;
    // the sparse matrix
    lgncg::SparseMatrix A;
    // vector x
    lgncg::Vector x;
    // vector b
    lgncg::Vector b;
    // generate problem task ID
    static int genProbTID;
    // populate fine to coarse operator task ID
    static int populatef2cTID;
    // number of mg levels
    int64_t nmgl;
    /**
     * registers any legion tasks that Problem may use.
     */
    static void
    init(void);
    /**
     * constructor
     */
    Problem(int64_t nx,
            int64_t ny,
            int64_t nz,
            int64_t stencilSize,
            int64_t numMGLevels,
            LegionRuntime::HighLevel::Context &ctx,
            LegionRuntime::HighLevel::HighLevelRuntime *lrt)
    {
        this->geom = Geometry(nx, ny, nz);
        this->stencilSize = stencilSize;
        this->nmgl = numMGLevels;
        // now init the matrix and vectors
        A.create(geom.nx * geom.ny * geom.nz, stencilSize, ctx, lrt);
        // length of vectors
        int64_t vLen = geom.nx * geom.ny * geom.nz;
        x.create<double>(vLen, ctx, lrt);
        b.create<double>(vLen, ctx, lrt);
    }
    /**
     * "destructor"
     */
    void
    free(LegionRuntime::HighLevel::Context &ctx,
         LegionRuntime::HighLevel::HighLevelRuntime *lrt)
    {
        A.free(ctx, lrt);
        x.free(ctx, lrt);
        b.free(ctx, lrt);
    }
    /**
     * wrapper for the task(s) responsible for setting the problem instance's
     * initial conditions.
     *
     * TODO: parallelize the routines here...
     */
    void
    setInitialConds(LegionRuntime::HighLevel::Context &ctx,
                    LegionRuntime::HighLevel::HighLevelRuntime *lrt)
    {
        using namespace LegionRuntime::HighLevel;
        using LegionRuntime::Arrays::Rect;
        // start timer
        double start = LegionRuntime::TimeStamp::get_current_time_in_micros();
        printf("setting initial conditions:\n");
        // generate problem at the finest level
        genProb(A, &x, &b, geom, stencilSize, ctx, lrt);
        // now do it for the rest
        lgncg::SparseMatrix *curLevMatPtr = &A;
        Geometry curFineGeom = geom;
        for (int64_t mgi = 1; mgi < this->nmgl; ++mgi) {
            Geometry curCoarseGeom;
            genCoarseProb(*curLevMatPtr, curFineGeom, curCoarseGeom,
                          stencilSize, ctx, lrt);
            // update for next round of coarse grid matrix generation
            curFineGeom = curCoarseGeom;
            curLevMatPtr = curLevMatPtr->Ac;
        }
        double stop = LegionRuntime::TimeStamp::get_current_time_in_micros();
        printf("  . done in: %7.3lf ms\n", (stop - start) * 1e-3);
    }
private:
    Problem(void) { ; }

    static void
    genProb(lgncg::SparseMatrix &A,
            lgncg::Vector *x,
            lgncg::Vector *b,
            const Geometry &geom,
            int64_t stencilSize,
            LegionRuntime::HighLevel::Context &ctx,
            LegionRuntime::HighLevel::HighLevelRuntime *lrt);

    static void
    genCoarseProb(lgncg::SparseMatrix &Af,
                  const Geometry &fineGeom,
                  Geometry &coarseGeom,
                  int64_t stencilSize,
                  LegionRuntime::HighLevel::Context &ctx,
                  LegionRuntime::HighLevel::HighLevelRuntime *lrt);

};

#endif
