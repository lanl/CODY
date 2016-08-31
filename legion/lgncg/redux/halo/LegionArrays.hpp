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

#include <deque>
#include <vector>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template<typename TYPE>
struct LogicalArray : public LogicalItem<TYPE> {
public:
    /**
     *
     */
    LogicalArray(void) : LogicalItem<TYPE>() { }

    /**
     *
     */
    void
    allocate(
        int64_t nElems,
        Legion::Context &ctx,
        Legion::HighLevelRuntime *lrt
    ) {
        this->mAllocate(nElems, ctx, lrt);
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
     *
     */
    void
    partition(
        size_t nParts,
        Legion::Context &ctx,
        Legion::HighLevelRuntime *lrt
    ) {
        using namespace Legion;
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
        size_t x0 = 0, x1 = inc - 1;
        DomainColoring disjointColoring;
        // a list of sub-grid bounds.
        // provides a task ID to sub-grid bounds mapping.
        std::vector< Rect<1> > subGridBounds;
        for (int64_t color = 0; color < nParts; ++color) {
            Rect<1> subRect((Point<1>(x0)), (Point<1>(x1)));
            // cache the subgrid bounds
            subGridBounds.push_back(subRect);
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
        using Legion::LogicalPartition;
        this->logicalPartition = lrt->get_logical_partition(
                                     ctx, this->logicalRegion, iPart
                                 );
        // launch domain -- one task per color
        // launch domain
        this->launchDomain = colorDomain;
    }

    /**
     *
     */
    void
    partition(
        std::vector<size_t> partLens,
        Legion::Context &ctx,
        Legion::HighLevelRuntime *lrt
    ) {
        using namespace Legion;
        using LegionRuntime::Arrays::Rect;

        //
        const size_t nParts = partLens.size();
        Rect<1> colorBounds(Point<1>(0), Point<1>(nParts - 1));
        Domain colorDomain = Domain::from_rect<1>(colorBounds);
        //
        size_t parti = 0;
        size_t x0 = 0, x1 = partLens[parti] - 1;
        DomainColoring disjointColoring;
        // provides a task ID to sub-grid bounds mapping.
        std::vector< Rect<1> > subGridBounds;
        for (int64_t color = 0; color < nParts; ++color) {
            Rect<1> subRect((Point<1>(x0)), (Point<1>(x1)));
            // cache the subgrid bounds
            subGridBounds.push_back(subRect);
#if 0 // nice debug
            printf("vec len=%ld\n", (long)partLens[parti]);
            printf("vec disjoint partition: (%ld) to (%ld)\n",
                    (long)subRect.lo.x[0], (long)subRect.hi.x[0]);
#endif
            disjointColoring[color] = Domain::from_rect<1>(subRect);
            // slide window
            parti++;
            x0 += partLens[parti];
            x1 += partLens[parti];
        }
        auto iPart = lrt->create_index_partition(
                         ctx, this->mIndexSpace,
                         colorDomain, disjointColoring,
                         true /* disjoint */
        );
        // logical partitions
        using Legion::LogicalPartition;
        this->logicalPartition = lrt->get_logical_partition(
                                     ctx, this->logicalRegion, iPart
                                 );
        // launch domain -- one task per color
        // launch domain
        this->launchDomain = colorDomain;
    }
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template<typename TYPE>
class Array : public Item<TYPE> {
public:
    /**
     *
     */
    Array(void) = default;

    /**
     *
     */
    Array(
        const PhysicalRegion &physicalRegion,
        Context ctx,
        HighLevelRuntime *runtime
    ) : Item<TYPE>(physicalRegion, ctx, runtime) { }

    /**
     *
     */
    size_t
    length(void) { return this->mLength; }
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template<typename TYPE>
struct ArrayAllocator {
private:
    Context mCTX;
    //
    HighLevelRuntime *mRuntime = nullptr;

public:
    struct Logical {
        LogicalArray<TYPE> array;
    };

    struct Physical {
        Array<TYPE> array;
        PhysicalRegion region;
    };

    //
    Logical l;
    //
    Physical p;

    /**
     *
     */
    ArrayAllocator(void) = default;

    /**
     *
     */
    ~ArrayAllocator(void)
    {
        // SKG TODO
        // Do we want to unmap when we have allocated in a task? Will it flow to
        // the top-level task if we do???
        //l.array.unmapRegion(mCTX, mRuntime);
    }

    /**
     *
     */
    ArrayAllocator(
        uint64_t nItems,
        Legion::PrivilegeMode privMode,
        Legion::CoherenceProperty cohProp,
        Context ctx,
        HighLevelRuntime *runtime
    ) {
        // Allocate logical array.
        l.array.allocate(nItems, ctx, runtime);
        // Inline map and stash physical region.
        p.region = l.array.mapRegion(privMode, cohProp, ctx, runtime);
        // Instantiate physical array type from physical region.
        p.array = Array<TYPE>(p.region, ctx, runtime);
        // Cache for destructor.
        mCTX = ctx;
        mRuntime = runtime;
    }

    /**
     *
     */
    TYPE *
    data(void) { return p.array.data(); }

    /**
     *
     */
    void
    bindToLogicalRegion(LogicalRegion &lr) {
        lr = p.region.get_logical_region();
    }

    /**
     *
     */
    void
    unmapRegion(void) {
        l.array.unmapRegion(mCTX, mRuntime);
    }
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/**
 * Interprets 1D array as NxM 2D array.
 */
template<typename TYPE>
class Array2D {
protected:
    //
    size_t mNRows = 0;
    //
    size_t mNCols = 0;
    //
    TYPE *mBasePtr = nullptr;
public:
    /**
     *
     */
    Array2D(void) = default;

    /**
     *
     */
    ~Array2D(void) {
        mNRows = 0;
        mNCols = 0;
        // Don't free memory here because we don't know what allocated that
        // memory. Assume that it'll get cleaned up another way.
        mBasePtr = nullptr;
    }

    /**
     *
     */
    Array2D(
        size_t nRows,
        size_t nCols,
        TYPE *basePtr
    ) : mNRows(nRows)
      , mNCols(nCols)
      , mBasePtr(basePtr) { }

    /**
     *
     */
    const TYPE &
    operator()(size_t row, size_t col) const
    {
        return mBasePtr[(row*mNCols)+col];
    }

    /**
     *
     */
    TYPE &
    operator()(size_t row, size_t col)
    {
        return mBasePtr[(row*mNCols)+col];
    }
};

/**
 *
 */
template<
    Legion::PrivilegeMode PRIV_MODE,
    Legion::CoherenceProperty COH_PROP
>
static void
intent(
    Legion::IndexLauncher &launcher,
    const std::deque<LogicalItemBase> &targetArrays
) {
    for (auto &a : targetArrays) {
        launcher.add_region_requirement(
            RegionRequirement(
                a.logicalPartition,
                0,
                PRIV_MODE,
                COH_PROP,
                a.logicalRegion
            )
        ).add_field(a.fid);
    }
}
