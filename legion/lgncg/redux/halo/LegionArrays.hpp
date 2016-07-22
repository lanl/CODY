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

#include "LegionItems.hpp"
#include "LegionStuff.hpp"

#include <vector>
#include <iomanip>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
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

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template<typename TYPE>
struct LogicalArray : public LogicalItem<TYPE> {
    std::vector<PVecItem> mPVec;
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
    launchDomain(size_t index = -1) const
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
        assert(0 == this->mLength % nParts && "Uneven partitioning requested.");
        //
        int64_t inc = this->mLength / nParts; // the increment
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
                         ctx, this->mIndexSpace,
                         colorDomain, disjointColoring,
                         true /* disjoint */
        );
        // logical partitions
        using LegionRuntime::HighLevel::LogicalPartition;
        auto lp = lrt->get_logical_partition(ctx, this->logicalRegion, iPart);
        // launch domain -- one task per color
        // launch domain
        LegionRuntime::HighLevel::Domain lDom = colorDomain;
        // add to the vector of partitions
        mPVec.push_back(PVecItem(subgridBnds, lDom, lp));
    }
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template<typename TYPE>
class PhysicalScalar {
protected:
    //
    size_t mLength = 0;
    //
    TYPE *mData = nullptr;
    //
    PhysicalScalar(void) = default;
public:
    //
    PhysicalScalar(
        const PhysicalRegion &physicalRegion,
        Context ctx,
        HighLevelRuntime *runtime
    ) {
        typedef RegionAccessor<AccessorType::Generic, TYPE>  GRA;
        GRA tAcc = physicalRegion.get_field_accessor(0).template typeify<TYPE>();
        //
        Domain tDom = runtime->get_index_space_domain(
            ctx, physicalRegion.get_logical_region().get_index_space()
        );
        Rect<1> subrect;
        ByteOffset inOffsets[1];
        auto subGridBounds = tDom.get_rect<1>();
        mLength = subGridBounds.volume();
        //
        mData = tAcc.template raw_rect_ptr<1>(
            subGridBounds, subrect, inOffsets
        );
        // Sanity.
        if (!mData || (subrect != subGridBounds) ||
            !offsetsAreDense<1, TYPE>(subGridBounds, inOffsets)) {
            // Signifies that something went south.
            mData = nullptr;
        }
        // It's all good...
    }

    /**
     *
     */
    TYPE *
    data(void) { return mData; }
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template<typename TYPE>
class PhysicalArray : public PhysicalScalar<TYPE> {
protected:
    //
    PhysicalArray(void) = default;
public:
    //
    PhysicalArray(
        const PhysicalRegion &physicalRegion,
        Context ctx,
        HighLevelRuntime *runtime
    ) : PhysicalScalar<TYPE>(physicalRegion, ctx, runtime) { }

    /**
     *
     */
    size_t
    length(void) { return this->mLength; }
};
