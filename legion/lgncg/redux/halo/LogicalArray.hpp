/*
 * Copyright (c) 2014-2016 Los Alamos National Security, LLC
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
 */

#pragma once

#include "LegionStuff.hpp"

#include <vector>

/**
 * Partition vector item.
 */
struct PVecItem {
    // a list of sub-grid bounds. provides a task ID to sub-grid bounds mapping
    std::vector< Rect<1> > subgridBnds;
    // launch domain
    LegionRuntime::HighLevel::Domain lDom;
    // logical partition
    LegionRuntime::HighLevel::LogicalPartition lPart;

    /**
     * constructor
     */
    PVecItem(
        const std::vector< Rect<1> > &sgb,
        const LegionRuntime::HighLevel::Domain &lDom,
        const LegionRuntime::HighLevel::LogicalPartition &lp
    ) : subgridBnds(sgb), lDom(lDom), lPart(lp) { ; }

private:
    //
    PVecItem(void) { ; }
};

template<typename T>
struct LogicalArray {
    // Number of elements stored in the vector (the entire extent).
    int64_t length;
    // Field ID.
    LegionRuntime::HighLevel::FieldID fid;
    // The vector rectangle bounds.
    LegionRuntime::Arrays::Rect<1> bounds;
    // Logical region that represents array.
    LegionRuntime::HighLevel::LogicalRegion logicalRegion;
private:
    // Index space.
    LegionRuntime::HighLevel::IndexSpace mIndexSpace;
    // Field space.
    LegionRuntime::HighLevel::FieldSpace mFS;
    // Vector of partition items.
    std::vector<PVecItem> mPVec;
    // The following are used for vector equality tests. That is, equality in
    // the "are these vectors the same from legion's perspective."
    id_t mIndexSpaceID;
    //
    LegionRuntime::HighLevel::FieldSpaceID mFieldSpaceID;
    //
    LegionRuntime::HighLevel::RegionTreeID mRTreeID;

public:
    /**
     *
     */
    void
    allocate(
        int64_t nElems,
        LegionRuntime::HighLevel::Context &ctx,
        LegionRuntime::HighLevel::HighLevelRuntime *lrt
    ) {
        fid = 0;
        length = nElems;
        // calculate the size of the logicalRegion vec (inclusive)
        auto n = length - 1;
        // vec rect
        bounds = Rect<1>(Point<1>::ZEROES(), Point<1>(n));
        // vector domain
        Domain dom(Domain::from_rect<1>(bounds));
        // vec index space
        mIndexSpace = lrt->create_index_space(ctx, dom);
        // vec field space
        mFS = lrt->create_field_space(ctx);
        // vec field allocator
        FieldAllocator fa = lrt->create_field_allocator(ctx, mFS);
        // all elements are going to be of size T
        fa.allocate_field(sizeof(T), fid);
        // now create the logical region
        logicalRegion = lrt->create_logical_region(ctx, mIndexSpace, mFS);
        // stash some info for equality checks
        mIndexSpaceID = logicalRegion.get_index_space().get_id();
        mFieldSpaceID = logicalRegion.get_field_space().get_id();
        mRTreeID      = logicalRegion.get_tree_id();
        // at this point we don't have a logical partition...  TODO maybe we can
        // just check if we are partitioned or not and return the correct
        // handle..?
    }

    /**
     * Cleans up and returns all allocated resources.
     */
    void
    deallocate(
        LegionRuntime::HighLevel::Context &ctx,
        LegionRuntime::HighLevel::HighLevelRuntime *lrt
    ) {
        for (auto i = mPVec.begin(); i != mPVec.end(); i++) {
            lrt->destroy_logical_partition(ctx, i->lPart);
        }
        mPVec.clear();
        lrt->destroy_logical_region(ctx, logicalRegion);
        lrt->destroy_field_space(ctx, mFS);
        lrt->destroy_index_space(ctx, mIndexSpace);
    }

    /**
     * Returns whether or not two LogicalArrays are the same (as far as the
     * Legion RT is concerned).
     */
    static bool
    same(
        const LogicalArray &a,
        const LogicalArray &b
    ) {
        return a.mIndexSpaceID == b.mIndexSpaceID &&
               a.mFieldSpaceID == b.mFieldSpaceID &&
               a.mRTreeID      == b.mRTreeID;
    }

    /**
     * Returns current (latest) launch domain or the one at specified index.
     */
    LegionRuntime::HighLevel::Domain
    lDom(size_t index = -1) const
    {
        if (index == -1) {
            const PVecItem &psi = mPVec.back();
            return psi.lDom;
        }
        const PVecItem &psi = mPVec[index];
        return psi.lDom;
    }

    /**
     * Returns current (latest) logical partition or the one at specified index.
     */
    LegionRuntime::HighLevel::LogicalPartition
    logicalPartition(size_t index = -1) const
    {
        if (index == -1) {
            const PVecItem &psi = mPVec.back();
            return psi.lPart;
        }
        const PVecItem &psi = mPVec[index];
        return psi.lPart;
    }

    /**
     * returns current sub-grid bounds.
     */
    std::vector< LegionRuntime::Arrays::Rect<1> >
    sgb(size_t index = 0) const
    {
        const PVecItem &psi = mPVec[index];
        return psi.subgridBnds;
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
        using namespace LegionRuntime::HighLevel;
        using LegionRuntime::Arrays::Rect;

        // For now only allow even partitioning.
        assert(0 == length % nParts && "Uneven partitioning requested.");
        //
        int64_t inc = length / nParts; // the increment
        Rect<1> colorBounds(Point<1>(0), Point<1>(nParts - 1));
        Domain colorDomain = Domain::from_rect<1>(colorBounds);
        //          +
        //          |
        //          |
        //     (x1)-+-+
        //          | |
        //          | m / nSubregions
        //     (x0) + |
        int64_t x0 = 0, x1 = inc - 1;
        DomainColoring disjointColoring;
        // a list of sub-grid bounds.
        // provides a task ID to sub-grid bounds mapping.
        std::vector< Rect<1> > subgridBnds;
        for (int64_t color = 0; color < nParts; ++color) {
            Rect<1> subRect((Point<1>(x0)), (Point<1>(x1)));
            // cache the subgrid bounds
            subgridBnds.push_back(subRect);
#if 0 // nice debug
            printf("vec disjoint partition: (%d) to (%d)\n",
                    subRect.lo.x[0], subRect.hi.x[0]);
#endif
            disjointColoring[color] = Domain::from_rect<1>(subRect);
            x0 += inc;
            x1 += inc;
        }
        auto iPart = lrt->create_index_partition(
                         ctx, mIndexSpace,
                         colorDomain, disjointColoring,
                         true /* disjoint */
        );
        // logical partitions
        using LegionRuntime::HighLevel::LogicalPartition;
        auto lp = lrt->get_logical_partition(ctx, logicalRegion, iPart);
        // launch domain -- one task per color
        // launch domain
        LegionRuntime::HighLevel::Domain lDom = colorDomain;
        // add to the vector of partitions
        mPVec.push_back(PVecItem(subgridBnds, lDom, lp));
    }
};
