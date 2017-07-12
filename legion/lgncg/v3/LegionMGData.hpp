/**
 * Copyright (c)      2017 Los Alamos National Security, LLC
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

//@HEADER
// ***************************************************
//
// HPCG: High Performance Conjugate Gradient Benchmark
//
// Contact:
// Michael A. Heroux ( maherou@sandia.gov)
// Jack Dongarra     (dongarra@eecs.utk.edu)
// Piotr Luszczek    (luszczek@eecs.utk.edu)
//
// ***************************************************
//@HEADER

/*!
      @file MGData.hpp

      HPCG data structure
 */

#pragma once

#include "LegionArrays.hpp"

#include <cassert>

class SparseMatrix;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
struct LogicalMGData : public LogicalMultiBase {
    // 1D array containing the fine operator local IDs that will be injected
    // into coarse space.
    LogicalArray<local_int_t> f2cOperator;
    // Coarse grid residual vector.
    LogicalArray<floatType>rc;
    // Coarse grid solution vector.
    LogicalArray<floatType> xc;
    // Fine grid residual vector.
    LogicalArray<floatType> Axf;

protected:

    /**
     * Order matters here. If you update this, also update unpack.
     */
    void
    mPopulateRegionList(void) {
        mLogicalItems = {
            &f2cOperator,
            &rc,
            &xc,
            &Axf
        };
    }

public:
    /**
     *
     */
    LogicalMGData(void) {
        mPopulateRegionList();
    }

    /**
     *
     */
    void
    allocate(
        const std::string &name,
        const Geometry &geom,
        LegionRuntime::HighLevel::Context ctx,
        LegionRuntime::HighLevel::HighLevelRuntime *lrt
    ) { /* Nothing to do. */ }

    /**
     *
     */
    void
    allocate(
        const std::string &name,
        SparseMatrix &A,
        LegionRuntime::HighLevel::Context ctx,
        LegionRuntime::HighLevel::HighLevelRuntime *lrt
    );

    /**
     *
     */
    void
    partition(
        int64_t nParts,
        LegionRuntime::HighLevel::Context ctx,
        LegionRuntime::HighLevel::HighLevelRuntime *lrt
    ) { /* Nothing to do. */ }

    /**
     *
     */
    void
    partition(
        SparseMatrix &A,
        Context ctx,
        HighLevelRuntime *lrt
    );
};


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
struct MGData : public PhysicalMultiBase {
    // Call ComputeSYMGS this many times prior to coarsening.
    int numberOfPresmootherSteps;
    // Call ComputeSYMGS this many times after coarsening.
    int numberOfPostsmootherSteps;
    //
    Array<local_int_t> *f2cOperator = nullptr;
    //
    Array<floatType> *rc = nullptr;
    //
    Array<floatType> *xc = nullptr;
    //
    Array<floatType> *Axf = nullptr;

    /**
     *
     */
    MGData(
        const std::vector<PhysicalRegion> &regions,
        size_t baseRID,
        Context ctx,
        HighLevelRuntime *runtime
    ) {
        numberOfPresmootherSteps  = 1;
        numberOfPostsmootherSteps = 1;
        mUnpack(regions, baseRID, IFLAG_NIL, ctx, runtime);
    }

    /**
     *
     */
    virtual
    ~MGData(void) {
        delete f2cOperator;
        delete rc;
        delete xc;
        delete Axf;
    }

protected:

    /**
     * MUST MATCH PACK ORDER IN mPopulateRegionList!
     */
    void
    mUnpack(
        const std::vector<PhysicalRegion> &regions,
        size_t baseRID,
        ItemFlags iFlags,
        Context ctx,
        HighLevelRuntime *rt
    ) {
        size_t cid = baseRID;
        // Populate members from physical regions.
        f2cOperator = new Array<local_int_t>(regions[cid++], ctx, rt);
        assert(f2cOperator->data());
        //
        rc = new Array<floatType>(regions[cid++], ctx, rt);
        assert(rc->data());
        //
        xc = new Array<floatType>(regions[cid++], ctx, rt);
        assert(xc->data());
        //
        Axf = new Array<floatType>(regions[cid++], ctx, rt);
        assert(Axf->data());
        // Calculate number of region entries for this structure.
        mNRegionEntries = cid - baseRID;
    }
};
