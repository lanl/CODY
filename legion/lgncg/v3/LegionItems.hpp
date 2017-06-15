/*
 * Copyright (c) 2014-2017 Los Alamos National Security, LLC
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
#include "Geometry.hpp"

#include <cassert>
#include <deque>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
struct LogicalItemBase {
public:
    // Field ID.
    Legion::FieldID fid = 0;
    // Logical region that represents item.
    Legion::LogicalRegion logicalRegion;
    // Launch domain.
    Legion::Domain launchDomain;
    // Logical partition.
    Legion::LogicalPartition logicalPartition;

protected:
    /**
     *
     */
    std::deque<LogicalItemBase *>
    mRegionPack(void) {
        return std::deque<LogicalItemBase *>({this});
    }

public:
    /**
     *
     */
    void
    intent(
        Legion::PrivilegeMode privMode,
        Legion::CoherenceProperty cohProp,
        Legion::IndexLauncher &launcher
    ) {
        launcher.add_region_requirement(
            RegionRequirement(
                logicalPartition,
                0,
                privMode,
                cohProp,
                logicalRegion
            )
        ).add_field(fid);
    }
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/**
 * Base class for logical structures that contain multiple logical structures.
 */
struct LogicalMultiBase {
    // Launch domain
    LegionRuntime::HighLevel::Domain launchDomain;

protected:

    std::deque<LogicalItemBase *> mLogicalItems;

    /**
     *
     */
    LogicalMultiBase(void) = default;

    /**
     *
     */
    virtual void
    mPopulateRegionList(void) = 0;

    /**
     *
     */
    void
    mIntent(
        Legion::PrivilegeMode privMode,
        Legion::CoherenceProperty cohProp,
        const std::deque<LogicalItemBase *> &targetArrays,
        Legion::IndexLauncher &launcher
    ) {
        for (auto &a : targetArrays) {
            a->intent(privMode, cohProp, launcher);
        }
    }

public:
    /**
     *
     */
    virtual void
    allocate(
        const Geometry &geom,
        LegionRuntime::HighLevel::Context &ctx,
        LegionRuntime::HighLevel::HighLevelRuntime *lrt
    ) { ; }

    /**
     *
     */
    virtual void
    partition(
        int64_t nParts,
        LegionRuntime::HighLevel::Context &ctx,
        LegionRuntime::HighLevel::HighLevelRuntime *lrt
    ) = 0;

    /**
     * Cleans up and returns all allocated resources.
     */
    virtual void
    deallocate(
        LegionRuntime::HighLevel::Context &ctx,
        LegionRuntime::HighLevel::HighLevelRuntime *lrt
    ) = 0;

    /**
     *
     */
    virtual void
    intent(
        Legion::PrivilegeMode privMode,
        Legion::CoherenceProperty cohProp,
        Legion::IndexLauncher &launcher
    ) {
        mIntent(privMode, cohProp, mLogicalItems, launcher);
    }
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template<typename TYPE>
struct LogicalItem : public LogicalItemBase {
    /**
     *
     */
    LogicalItem(void) : LogicalItemBase() { }

    /**
     * Instantiate a LogicalItem from a LogicalRegion.
     */
    LogicalItem(
        const LogicalRegion &lr,
        Legion::Context &ctx,
        Legion::HighLevelRuntime *lrt
    ) : LogicalItemBase()
    {
        //fid
        logicalRegion = lr;
        //launchDomain
        //logicalPartition
        mBounds       = lrt->get_index_space_domain(
                            ctx, lr.get_index_space()
                        ).get_rect<1>();
        //
        mLength       = 1;
        mIndexSpace   = lr.get_index_space();
        mFS           = lr.get_field_space();
        mIndexSpaceID = mIndexSpace.get_id();
        mFieldSpaceID = mFS.get_id();
        mRTreeID      = lr.get_tree_id();
    }

private:
    // The vector rectangle bounds.
    LegionRuntime::Arrays::Rect<1> mBounds;

protected:
    // Number of elements stored in the item (the entire extent).
    int64_t mLength = 0;
    // Index space.
    Legion::IndexSpace mIndexSpace;
    // Field space.
    Legion::FieldSpace mFS;
    // The following are used for vector equality tests. That is, equality in
    // the "are these vectors the same from legion's perspective."
    Legion::IndexSpaceID mIndexSpaceID;
    //
    Legion::FieldSpaceID mFieldSpaceID;
    //
    Legion::RegionTreeID mRTreeID;

    /**
     *
     */
    void
    mAllocate(
        int64_t len,
        Legion::Context &ctx,
        Legion::HighLevelRuntime *lrt
    ) {
        mLength = len;
        // calculate the size of the logicalRegion vec (inclusive)
        size_t n = mLength - 1;
        // vec rect
        mBounds = Rect<1>(Point<1>::ZEROES(), Point<1>(n));
        // vector domain
        Domain dom(Domain::from_rect<1>(mBounds));
        // vec index space
        mIndexSpace = lrt->create_index_space(ctx, dom);
        // vec field space
        mFS = lrt->create_field_space(ctx);
        // vec field allocator
        FieldAllocator fa = lrt->create_field_allocator(ctx, mFS);
        // all elements are going to be of size T
        fa.allocate_field(sizeof(TYPE), fid);
        // now create the logical region
        logicalRegion = lrt->create_logical_region(ctx, mIndexSpace, mFS);
        // stash some info for equality checks
        mIndexSpaceID = logicalRegion.get_index_space().get_id();
        mFieldSpaceID = logicalRegion.get_field_space().get_id();
        mRTreeID      = logicalRegion.get_tree_id();
    }

public:
    //
    PhysicalRegion physicalRegion;

    /**
     *
     */
    void
    allocate(
        Legion::Context &ctx,
        Legion::HighLevelRuntime *lrt
    ) {
        mAllocate(1, ctx, lrt);
    }

    /**
     * Cleans up and returns all allocated resources.
     */
    void
    deallocate(
        Legion::Context &ctx,
        Legion::HighLevelRuntime *lrt
    ) {
        lrt->destroy_logical_region(ctx, logicalRegion);
        lrt->destroy_field_space(ctx, mFS);
        lrt->destroy_index_space(ctx, mIndexSpace);
    }

    /**
     * Returns whether or not two LogicalItems are the same (as far as the
     * Legion RT is concerned).
     */
    static bool
    same(
        const LogicalItem &a,
        const LogicalItem &b
    ) {
        return a.mIndexSpaceID == b.mIndexSpaceID &&
               a.mFieldSpaceID == b.mFieldSpaceID &&
               a.mRTreeID      == b.mRTreeID;
    }

    /**
     *
     */
    Legion::PhysicalRegion
    mapRegion(
        Legion::PrivilegeMode privMode,
        Legion::CoherenceProperty cohProp,
        Legion::Context &ctx,
        Legion::HighLevelRuntime *lrt
    ) {
        using namespace Legion;
        using namespace LegionRuntime::Accessor;
        using LegionRuntime::Arrays::Rect;
        RegionRequirement req(
            logicalRegion, privMode, cohProp, logicalRegion
        );
        req.add_field(fid);
        InlineLauncher inl(req);
        physicalRegion = lrt->map_region(ctx, inl);
        physicalRegion.wait_until_valid();
        //
        return physicalRegion;
    }

    /**
     *
     */
    void
    unmapRegion(
        Legion::Context &ctx,
        Legion::HighLevelRuntime *lrt
    ) {
        lrt->unmap_region(ctx, physicalRegion);
    }
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template<typename TYPE>
class Item {
protected:
    //
    size_t mLength = 0;
    //
    TYPE *mData = nullptr;
public:
    //
    LogicalRegion logicalRegion;

    /**
     *
     */
    Item(void) = default;

    /**
     *
     */
    Item(
        const PhysicalRegion &physicalRegion,
        Context ctx,
        HighLevelRuntime *runtime
    ) {
        // cache logical region
        logicalRegion = physicalRegion.get_logical_region();
        //
        using GRA = RegionAccessor<AccessorType::Generic, TYPE>;
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
struct PhysicalMultiBase {
protected:
    // Number of region entries.
    size_t mNRegionEntries = 0;

    /**
     * MUST MATCH PACK ORDER IN mPopulateRegionList!
     */
    virtual void
    mUnpack(
        const std::vector<PhysicalRegion> &regions,
        size_t baseRID,
        Context ctx,
        HighLevelRuntime *rt
    ) = 0;

    /**
     *
     */
    virtual void
    mVerifyUnpack(void) = 0;

public:
    /**
     *
     */
    size_t
    nRegionEntries(void) { return mNRegionEntries; }
};
