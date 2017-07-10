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

#pragma once

#include "LegionStuff.hpp"
#include "LegionArrays.hpp"
#include "LegionMatrices.hpp"

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
struct LogicalCGData : public LogicalMultiBase {
    //
    LogicalArray<floatType> r;  //!< residual vector
    //
    LogicalArray<floatType> z;  //!< preconditioned residual vector
    //
    LogicalArray<floatType> p;  //!< direction vector
    //
    LogicalArray<floatType> Ap; //!< Krylov vector

protected:

    /**
     * Order matters here. If you update this, also update unpack.
     */
    void
    mPopulateRegionList(void) {
        mLogicalItems = {&r, &z, &p, &Ap};
    }

public:
    /**
     *
     */
    LogicalCGData(void) {
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

    void
    allocate(
        const std::string &name,
        SparseMatrix &A,
        LegionRuntime::HighLevel::Context ctx,
        LegionRuntime::HighLevel::HighLevelRuntime *lrt
    ) {
        #define aalloca(sName, size, ctx, rtp)                                 \
        do {                                                                   \
            sName.allocate(name + "-" #sName, size, ctx, rtp);                 \
        } while(0)

        const local_int_t nrow = A.sclrs->data()->localNumberOfRows;
        const local_int_t ncol = A.sclrs->data()->localNumberOfColumns;
        //
        aalloca(r,  nrow, ctx, lrt);
        aalloca(z,  ncol, ctx, lrt);
        aalloca(p,  ncol, ctx, lrt);
        aalloca(Ap, nrow, ctx, lrt);

        #undef aalloca
    }

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
        LegionRuntime::HighLevel::Context ctx,
        LegionRuntime::HighLevel::HighLevelRuntime *lrt
    ) {
        const local_int_t nrow = A.sclrs->data()->localNumberOfRows;
        const local_int_t ncol = A.sclrs->data()->localNumberOfColumns;
        // First nrow items are 'local data'. After is remote data.
        std::vector<local_int_t> partLens;
        // First partition for local data.
        local_int_t totLen = nrow;
        partLens.push_back(nrow);
        // The rest are based on receive lengths.
        const int nNeighbors = A.sclrs->data()->numberOfRecvNeighbors;
        const local_int_t *const recvLength = A.recvLength->data();
        //
        for (int n = 0; n < nNeighbors; ++n) {
            const int recvl = recvLength[n];
            partLens.push_back(recvl);
            totLen += recvl;
        }
        assert(totLen == ncol);
        //
        z.partition(partLens, ctx, lrt);
        p.partition(partLens, ctx, lrt);
        //
        const int numParts = 1;
        r.partition( numParts, ctx, lrt);
        Ap.partition(numParts, ctx, lrt);
    }

    /**
     * Cleans up and returns all allocated resources.
     */
    void
    deallocate(
        LegionRuntime::HighLevel::Context ctx,
        LegionRuntime::HighLevel::HighLevelRuntime *lrt
    ) {
        for (auto *i : mLogicalItems) {
            i->deallocate(ctx, lrt);
        }
    }
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
struct CGData : public PhysicalMultiBase {
    //
    Array<floatType> *r = nullptr;
    //
    Array<floatType> *z = nullptr;
    //
    Array<floatType> *p = nullptr;
    //
    Array<floatType> *Ap = nullptr;

    /**
     *
     */
    CGData(void) = default;

    /**
     *
     */
    virtual
    ~CGData(void) {
        delete r;
        delete z;
        delete p;
        delete Ap;
    }

    /**
     *
     */
    CGData(
        const std::vector<PhysicalRegion> &regions,
        size_t baseRID,
        Context ctx,
        HighLevelRuntime *runtime
    ) {
        mUnpack(regions, baseRID, IFLAG_NIL, ctx, runtime);
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
        r = new Array<floatType>(regions[cid++], ctx, rt);
        assert(r->data());
        //
        z = new Array<floatType>(regions[cid++], ctx, rt);
        assert(z->data());
        //
        p = new Array<floatType>(regions[cid++], ctx, rt);
        assert(p->data());
        //
        Ap = new Array<floatType>(regions[cid++], ctx, rt);
        assert(Ap->data());
        // Calculate number of region entries for this structure.
        mNRegionEntries = cid - baseRID;
    }
};
