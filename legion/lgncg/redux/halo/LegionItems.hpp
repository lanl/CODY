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

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
struct LogicalItemBase {
    // Field ID.
    LegionRuntime::HighLevel::FieldID fid = 0;
    // Logical region that represents array.
    LegionRuntime::HighLevel::LogicalRegion logicalRegion;
    // launch domain
    LegionRuntime::HighLevel::Domain launchDomain;
    // logical partition
    LegionRuntime::HighLevel::LogicalPartition logicalPartition;
    // The vector rectangle bounds.
    LegionRuntime::Arrays::Rect<1> bounds;
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
protected:
    // Number of elements stored in the vector (the entire extent).
    int64_t mLength;
    // Index space.
    LegionRuntime::HighLevel::IndexSpace mIndexSpace;
    // Field space.
    LegionRuntime::HighLevel::FieldSpace mFS;
    // The following are used for vector equality tests. That is, equality in
    // the "are these vectors the same from legion's perspective."
    LegionRuntime::HighLevel::IndexSpaceID mIndexSpaceID;
    //
    LegionRuntime::HighLevel::FieldSpaceID mFieldSpaceID;
    //
    LegionRuntime::HighLevel::RegionTreeID mRTreeID;
    /**
     *
     */
    void
    mAllocate(
        int64_t len,
        LegionRuntime::HighLevel::Context &ctx,
        LegionRuntime::HighLevel::HighLevelRuntime *lrt
    ) {
        mLength = len;
        // calculate the size of the logicalRegion vec (inclusive)
        auto n = mLength - 1;
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
        fa.allocate_field(sizeof(TYPE), fid);
        // now create the logical region
        logicalRegion = lrt->create_logical_region(ctx, mIndexSpace, mFS);
        // stash some info for equality checks
        mIndexSpaceID = logicalRegion.get_index_space().get_id();
        mFieldSpaceID = logicalRegion.get_field_space().get_id();
        mRTreeID      = logicalRegion.get_tree_id();
    }
public:
    /**
     *
     */
    void
    allocate(
        LegionRuntime::HighLevel::Context &ctx,
        LegionRuntime::HighLevel::HighLevelRuntime *lrt
    ) {
        mAllocate(1, ctx, lrt);
    }

    /**
     * Cleans up and returns all allocated resources.
     */
    void
    deallocate(
        LegionRuntime::HighLevel::Context &ctx,
        LegionRuntime::HighLevel::HighLevelRuntime *lrt
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
    LegionRuntime::HighLevel::PhysicalRegion
    mapRegion(
        LegionRuntime::HighLevel::PrivilegeMode privMode,
        LegionRuntime::HighLevel::CoherenceProperty cohProp,
        LegionRuntime::HighLevel::Context &ctx,
        LegionRuntime::HighLevel::HighLevelRuntime *lrt
    ) const {
        using namespace LegionRuntime::HighLevel;
        using namespace LegionRuntime::Accessor;
        using LegionRuntime::Arrays::Rect;
        RegionRequirement req(
            logicalRegion, privMode, cohProp, logicalRegion
        );
        req.add_field(fid);
        InlineLauncher inl(req);
        PhysicalRegion reg = lrt->map_region(ctx, inl);
        reg.wait_until_valid();
        //
        return reg;
    }

    /**
     *
     */
    void
    unmapRegion(
        LegionRuntime::HighLevel::PhysicalRegion &mappedRegion,
        LegionRuntime::HighLevel::Context &ctx,
        LegionRuntime::HighLevel::HighLevelRuntime *lrt
    ) const {
        lrt->unmap_region(ctx, mappedRegion);
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
