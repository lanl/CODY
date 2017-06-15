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
        const Geometry &geom,
        LegionRuntime::HighLevel::Context &ctx,
        LegionRuntime::HighLevel::HighLevelRuntime *lrt
    ) {
        const auto globalXYZ = getGlobalXYZ(geom);

        r.allocate(globalXYZ, ctx, lrt);
#if 0
#warning "Misallocation of structure..."
#endif
        z.allocate(globalXYZ, ctx, lrt);
#if 0
#warning "Misallocation of structure..."
#endif
        p.allocate(globalXYZ, ctx, lrt);
        Ap.allocate(globalXYZ, ctx, lrt);
    }

    /**
     *
     */
    void
    partition(
        int64_t nParts,
        LegionRuntime::HighLevel::Context &ctx,
        LegionRuntime::HighLevel::HighLevelRuntime *lrt
    ) {
        r.partition(nParts, ctx, lrt);
        z.partition(nParts, ctx, lrt);
        p.partition(nParts, ctx, lrt);
        Ap.partition(nParts, ctx, lrt);
    }

    /**
     * Cleans up and returns all allocated resources.
     */
    void
    deallocate(
        LegionRuntime::HighLevel::Context &ctx,
        LegionRuntime::HighLevel::HighLevelRuntime *lrt
    ) {
        r.deallocate(ctx, lrt);
        z.deallocate(ctx, lrt);
        p.deallocate(ctx, lrt);
        Ap.deallocate(ctx, lrt);
    }
};

struct CGData {
    Array<floatType> r;  //!< residual vector
    Array<floatType> z;  //!< preconditioned residual vector
    Array<floatType> p;  //!< direction vector
    Array<floatType> Ap; //!< Krylov vector
};
