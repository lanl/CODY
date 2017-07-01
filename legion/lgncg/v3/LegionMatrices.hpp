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

#pragma once

#include "LegionStuff.hpp"
#include "LegionItems.hpp"
#include "LegionArrays.hpp"

#include "hpcg.hpp"
#include "Geometry.hpp"

#include <vector>
#include <map>
#include <cassert>
#include <deque>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
struct SparseMatrixScalars {
    //Max number of non-zero elements in any row.
    local_int_t maxNonzerosPerRow = 0;
    //Total number of matrix rows across all processes.
    global_int_t totalNumberOfRows = 0;
    //Total number of matrix non-zeros across all processes.
    global_int_t totalNumberOfNonzeros = 0;
    //Number of rows local to this process.
    local_int_t localNumberOfRows = 0;
    //Number of columns local to this process.
    local_int_t localNumberOfColumns = 0;
    //Number of non-zeros local to this process.
    global_int_t localNumberOfNonzeros = 0;
    //Number of entries that are external to this process.
    local_int_t numberOfExternalValues = 0;
    //Number of neighboring processes that will be send local data.
    int numberOfSendNeighbors = 0;
    //Total number of entries to be sent.
    local_int_t totalToBeSent = 0;
    // Buffer used to collect local partial sums.
    floatType localPartialSum = 0;
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
struct DynCollBuffer {
    floatType floatTypeTarget;
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/**
 * Holds structures required for task synchronization.
 */
struct Synchronizers {
    // The PhaseBarriers that I own.
    PhaseBarriers mine;
    // Dense array of neighbor PhaseBarriers that will only have the first
    // nNeighbors - 1 entries populated, so BE CAREFUL ;). Wasteful, but done
    // this way for convenience. At most a task will have HPCG_STENCIL -1
    // neighbors.
    PhaseBarriers neighbors[HPCG_STENCIL - 1];
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
struct LogicalSparseMatrix : public LogicalMultiBase {
    using LogicalMultiBase::intent;
    //
    LogicalArray<Geometry> geoms;
    //
    LogicalArray<SparseMatrixScalars> sclrs;
    //
    LogicalArray<DynCollBuffer> dcBuffer;
    //
    LogicalArray<char> nonzerosInRow;
    //
    LogicalArray<global_int_t> mtxIndG;
    //
    LogicalArray<local_int_t> mtxIndL;
    //
    LogicalArray<floatType> matrixValues;
    //
    LogicalArray<floatType> matrixDiagonal;
    //
    LogicalArray<global_int_t> localToGlobalMap;
    // The SAME dynamic collective instance replicated because
    // IndexLauncher will be unhappy with different launch domains :-(
    LogicalArray<DynamicCollective> dcAllreduceSum;
    // Neighboring processes.
    LogicalArray<int> neighbors;
    // Number of items that will be sent on a per neighbor basis.
    LogicalArray<local_int_t> sendLength;
    // Number of items that will be received on a per neighbor basis.
    LogicalArray<local_int_t> recvLength;
    // Synchronization structures.
    LogicalArray<Synchronizers> synchronizers;
    //
    LogicalArray<BaseExtent> pullBEs;
    ////////////////////////////////////////////////////////////////////////////
    // IFLAG_W_GHOSTS structures.
    ////////////////////////////////////////////////////////////////////////////
    // Buffer that will be used  to pull data from during ExchangeHalo.
    LogicalArray<floatType> pullBuffer;

protected:
    // Number of shards used for SparseMatrix decomposition.
    int mSize = 0;
    //
    std::deque<LogicalItemBase *> mLogicalItemsGhost;

    /**
     * Order matters here. If you update this, also update unpack.
     */
    void
    mPopulateRegionList(void) {
        mLogicalItems = {&geoms,
                         &sclrs,
                         &dcBuffer,
                         &nonzerosInRow,
                         &mtxIndG,
                         &mtxIndL,
                         &matrixValues,
                         &matrixDiagonal,
                         &localToGlobalMap,
                         &dcAllreduceSum,
                         &neighbors,
                         &sendLength,
                         &recvLength,
                         &synchronizers,
                         &pullBEs
        };
        //
        mLogicalItemsGhost = {&pullBuffer};
    }

public:
    /**
     *
     */
    LogicalSparseMatrix(void) {
        // -1 signifies that regions have not yet been allocated.
        mSize = -1;
        mPopulateRegionList();
    }

    /**
     *
     */
    void
    intent(
        Legion::PrivilegeMode privMode,
        Legion::CoherenceProperty cohProp,
        ItemFlags iFlags,
        Legion::IndexLauncher &launcher,
        LegionRuntime::HighLevel::Context &ctx,
        LegionRuntime::HighLevel::HighLevelRuntime *lrt
    ) {
        intent(privMode, cohProp, launcher);
        //
        if (withGhosts(iFlags)) {
            for (auto &a : mLogicalItemsGhost) {
                // Share every sub-region with everyone. During task execution
                // each task will pick only the regions that are required.
                for (int color = 0; color < mSize; ++color) {
                    auto lsr = lrt->get_logical_subregion_by_color(
                                   ctx,
                                   a->logicalPartition,
                                   color
                    );
                    launcher.add_region_requirement(
                        RegionRequirement(
                            lsr,
                            0,
                            /* NOTE disregard for passed args here. */
                            READ_WRITE,
                            SIMULTANEOUS,
                            a->logicalRegion
                        )
                    ).add_field(a->fid)
                     .add_flags(NO_ACCESS_FLAG);
                }
            }
        }
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
        mSize = geom.size;
        const auto globalXYZ   = getGlobalXYZ(geom);
        const auto stencilSize = geom.stencilSize;

        geoms.allocate(mSize, ctx, lrt);
        sclrs.allocate(mSize, ctx, lrt);
        dcBuffer.allocate(mSize, ctx, lrt);
        //
        nonzerosInRow.allocate(globalXYZ, ctx, lrt);
        // Flattened to 1D from 2D.
        mtxIndG.allocate(globalXYZ * stencilSize, ctx, lrt);
        // Flattened to 1D from 2D.
        mtxIndL.allocate(globalXYZ * stencilSize, ctx, lrt);
        // Flattened to 1D from 2D.
        matrixValues.allocate(globalXYZ * stencilSize, ctx, lrt);
        // 2D thing in reference implementation, but not needed (1D suffices).
        matrixDiagonal.allocate(globalXYZ, ctx, lrt);
        //
        localToGlobalMap.allocate(globalXYZ, ctx, lrt);
        //
        dcAllreduceSum.allocate(mSize, ctx, lrt);
        //
        const int maxNumNeighbors = geom.stencilSize - 1;
        // Each task will have at most 26 neighbors.
        neighbors.allocate(mSize * maxNumNeighbors, ctx, lrt);
        // Each task will have at most 26 neighbors.
        sendLength.allocate(mSize * maxNumNeighbors, ctx, lrt);
        recvLength.allocate(mSize * maxNumNeighbors, ctx, lrt);
        //
        synchronizers.allocate(mSize, ctx, lrt);
        //
        pullBEs.allocate(mSize * maxNumNeighbors, ctx, lrt);
        ////////////////////////////////////////////////////////////////////////
        // IFLAG_W_GHOSTS structures.
        ////////////////////////////////////////////////////////////////////////
        // FIXME: A bit wasteful on storage.
        pullBuffer.allocate(globalXYZ, ctx, lrt);
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
        geoms.partition(nParts, ctx, lrt);
        sclrs.partition(nParts, ctx, lrt);
        dcBuffer.partition(nParts, ctx, lrt);
        nonzerosInRow.partition(nParts, ctx, lrt);
        mtxIndG.partition(nParts, ctx, lrt);
        mtxIndL.partition(nParts, ctx, lrt);
        matrixValues.partition(nParts, ctx, lrt);
        matrixDiagonal.partition(nParts, ctx, lrt);
        localToGlobalMap.partition(nParts, ctx, lrt);
        dcAllreduceSum.partition(nParts, ctx, lrt);
        neighbors.partition(nParts, ctx, lrt);
        sendLength.partition(nParts, ctx, lrt);
        recvLength.partition(nParts, ctx, lrt);
        synchronizers.partition(nParts, ctx, lrt);
        pullBEs.partition(nParts, ctx, lrt);
        ////////////////////////////////////////////////////////////////////////
        // IFLAG_W_GHOSTS structures.
        ////////////////////////////////////////////////////////////////////////
        pullBuffer.partition(nParts, ctx, lrt);
        //
        // For the DynamicCollectives we need partition info before population.
        mPopulateDynamicCollectives(nParts, ctx, lrt);
        // Just pick a structure that has a representative launch domain.
        launchDomain = geoms.launchDomain;
    }

    /**
     * Cleans up and returns all allocated resources.
     */
    void
    deallocate(
        LegionRuntime::HighLevel::Context &ctx,
        LegionRuntime::HighLevel::HighLevelRuntime *lrt
    ) {
        geoms.deallocate(ctx, lrt);
        sclrs.deallocate(ctx, lrt);
        dcBuffer.deallocate(ctx, lrt);
        mtxIndG.deallocate(ctx, lrt);
        mtxIndL.deallocate(ctx, lrt);
        nonzerosInRow.deallocate(ctx, lrt);
        matrixValues.deallocate(ctx, lrt);
        matrixDiagonal.deallocate(ctx, lrt);
        localToGlobalMap.deallocate(ctx, lrt);
        dcAllreduceSum.deallocate(ctx, lrt);
        neighbors.deallocate(ctx, lrt);
        sendLength.deallocate(ctx, lrt);
        recvLength.deallocate(ctx, lrt);
        synchronizers.deallocate(ctx, lrt);
        pullBEs.deallocate(ctx, lrt);
        ////////////////////////////////////////////////////////////////////////
        // IFLAG_W_GHOSTS structures.
        ////////////////////////////////////////////////////////////////////////
        pullBEs.deallocate(ctx, lrt);
        pullBuffer.deallocate(ctx, lrt);
    }

private:
    void
    mPopulateDynamicCollectives(
        int64_t nArrivals,
        LegionRuntime::HighLevel::Context &ctx,
        LegionRuntime::HighLevel::HighLevelRuntime *lrt
    ) {
        Array<DynamicCollective> dcs(
            dcAllreduceSum.mapRegion(RW_E, ctx, lrt), ctx, lrt
        );
        global_int_t dummy = 0;
        DynamicCollective dc = lrt->create_dynamic_collective(
            ctx,
            nArrivals /* Number of arrivals. */,
            INT_REDUCE_SUM_TID,
            &dummy,
            sizeof(dummy)
        );
        // Replicate
        DynamicCollective *dcsd = dcs.data();
        for (int64_t i = 0; i < nArrivals; ++i) {
            dcsd[i] = dc;
        }
        // Done, so unmap.
        dcAllreduceSum.unmapRegion(ctx, lrt);
    }
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
struct SparseMatrix : public PhysicalMultiBase {
    // Geometry info for this instance.
    Item<Geometry> *geom = nullptr;
    // Container for all scalar values.
    Item<SparseMatrixScalars> *sclrs = nullptr;
    //
    Item<DynCollBuffer> *dcBuffer = nullptr;
    //
    Array<char> *nonzerosInRow = nullptr;
    // Flattened to 1D from 2D.
    Array<global_int_t> *mtxIndG = nullptr;
    // Flattened to 1D from 2D.
    Array<local_int_t> *mtxIndL = nullptr;
    // Flattened to 1D from 2D.
    Array<floatType> *matrixValues = nullptr;
    //
    Array<floatType> *matrixDiagonal = nullptr;
    //
    Array<global_int_t> *localToGlobalMap = nullptr;
    //
    Item<DynamicCollective> *dcAllreduceSum = nullptr;
    //
    Array<int> *neighbors = nullptr;
    //
    Array<local_int_t> *sendLength = nullptr;
    //
    Array<local_int_t> *recvLength = nullptr;
    //
    Item<Synchronizers> *synchronizers = nullptr;
    // The bases and extents that I will be getting from my neighbors that lets
    // me know which contiguous set of points will make up values I need to
    // read.
    Array<BaseExtent> *pullBEs = nullptr;

    ////////////////////////////////////////////////////////////////////////////
    // Task-launch-local structures.
    ////////////////////////////////////////////////////////////////////////////
    Context lctx;
    //
    HighLevelRuntime *lrt = nullptr;
    // Global to local mapping. NOTE: only valid after a call to
    // PopulateGlobalToLocalMap.
    std::map< global_int_t, local_int_t > globalToLocalMap;
    // Only valid after a call to SetupHalo.
    local_int_t *elementsToSend = nullptr;
    // A mapping between neighbor IDs and their regions.
    std::map<int, PhysicalRegion> neighborToRegions;
    // A mapping between neighbor IDs and ghost Arrays.
    std::map< int, LogicalArray<floatType> > ghostArrays;
    // The Array that holds push values.
    Array<floatType> *pullBuffer = nullptr;

    /**
     *
     */
    int
    mSetupGhostStructures(
        const std::vector<PhysicalRegion> &regions,
        size_t baseRID,
        Context ctx,
        HighLevelRuntime *runtime
    ) {
        const int *const nd = neighbors->data();
        const SparseMatrixScalars *const sclrsd = sclrs->data();
        const int me = geom->data()->rank;
        // Setup my push Array.
        pullBuffer = new Array<floatType>(regions[baseRID + me], ctx, runtime);
        assert(pullBuffer->data());
        // Get neighbor regions.
        for (int n = 0; n < sclrsd->numberOfSendNeighbors; ++n) {
            const int nid = nd[n];
            neighborToRegions[nid] = regions[baseRID + nid];
        }
        // TODO Unmap unused regions.
        // Return number of regions that we have consumed.
        return geom->data()->size;
    }

    /**
     *
     */
    SparseMatrix(void) = default;

    /**
     *
     */
    SparseMatrix(
        const std::vector<PhysicalRegion> &regions,
        size_t baseRID,
        Context ctx,
        HighLevelRuntime *runtime
    ) {
        mDoConstruct(regions, baseRID, IFLAG_NIL, ctx, runtime);
    }

    /**
     *
     */
    SparseMatrix(
        const std::vector<PhysicalRegion> &regions,
        size_t baseRID,
        ItemFlags iFlags,
        Context ctx,
        HighLevelRuntime *runtime
    ) {
        mDoConstruct(regions, baseRID, iFlags, ctx, runtime);
    }

    /**
     *
     */
    virtual
    ~SparseMatrix(void) {
        delete geom;
        delete sclrs;
        delete dcBuffer;
        delete nonzerosInRow;
        delete mtxIndG;
        delete mtxIndL;
        delete matrixValues;
        delete matrixDiagonal;
        delete localToGlobalMap;
        delete dcAllreduceSum;
        delete neighbors;
        delete sendLength;
        delete recvLength;
        delete synchronizers;
        // Task-local allocation of non-region memory.
        if (elementsToSend) delete[] elementsToSend;
        if (withGhosts(mUnpackFlags)) {
            delete pullBuffer;
        }
    }

protected:

    /**
     *
     */
    void
    mDoConstruct(
        const std::vector<PhysicalRegion> &regions,
        size_t baseRID,
        ItemFlags iFlags,
        Context ctx,
        HighLevelRuntime *rt
    ) {
        lctx = ctx;
        lrt = rt;
        mUnpackFlags = iFlags;
        //
        mUnpack(regions, baseRID, iFlags, ctx, rt);
    }

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
        geom = new Item<Geometry>(regions[cid++], ctx, rt);
        assert(geom->data());
        //
        sclrs = new Item<SparseMatrixScalars>(regions[cid++], ctx, rt);
        assert(sclrs->data());
        //
        dcBuffer = new Item<DynCollBuffer>(regions[cid++], ctx, rt);
        assert(dcBuffer->data());
        //
        nonzerosInRow = new Array<char>(regions[cid++], ctx, rt);
        assert(nonzerosInRow->data());
        //
        mtxIndG = new Array<global_int_t>(regions[cid++], ctx, rt);
        assert(mtxIndG->data());
        //
        mtxIndL = new Array<local_int_t>(regions[cid++], ctx, rt);
        assert(mtxIndL->data());
        //
        matrixValues = new Array<floatType>(regions[cid++], ctx, rt);
        assert(matrixValues->data());
        //
        matrixDiagonal = new Array<floatType>(regions[cid++], ctx, rt);
        assert(matrixDiagonal->data());
        //
        localToGlobalMap = new Array<global_int_t>(regions[cid++], ctx, rt);
        assert(localToGlobalMap->data());
        //
        dcAllreduceSum = new Item<DynamicCollective>(regions[cid++], ctx, rt);
        assert(dcAllreduceSum->data());
        //
        neighbors = new Array<int>(regions[cid++], ctx, rt);
        assert(neighbors->data());
        //
        sendLength = new Array<local_int_t>(regions[cid++], ctx, rt);
        assert(sendLength->data());
        //
        recvLength = new Array<local_int_t>(regions[cid++], ctx, rt);
        assert(recvLength->data());
        //
        synchronizers = new Item<Synchronizers>(regions[cid++], ctx, rt);
        assert(synchronizers->data());
        //
        pullBEs = new Array<BaseExtent>(regions[cid++], ctx, rt);
        assert(pullBEs->data());
        if (withGhosts(iFlags)) {
            cid += mSetupGhostStructures(regions, cid, ctx, rt);
        }
        // Calculate number of region entries for this structure.
        mNRegionEntries = cid - baseRID;
    }
};

/**
 *
 */
inline global_int_t
localNonzerosTask(
    const Task *task,
    const std::vector<PhysicalRegion> &regions,
    Context ctx, HighLevelRuntime *runtime
) {
    Item<SparseMatrixScalars> sms(regions[0], ctx, runtime);
    return sms.data()->localNumberOfNonzeros;
}

/**
 *
 */
inline floatType
localPartialSumTask(
    const Task *task,
    const std::vector<PhysicalRegion> &regions,
    Context ctx, HighLevelRuntime *runtime
) {
    Item<SparseMatrixScalars> sms(regions[0], ctx, runtime);
    return sms.data()->localPartialSum;
}

/**
 *
 */
inline void
PopulateGlobalToLocalMap(
    SparseMatrix &A,
    LegionRuntime::HighLevel::Context ctx,
    LegionRuntime::HighLevel::Runtime *runtime
) {
    const Geometry *const Ageom = A.geom->data();
    // Make local copies of geometry information.  Use global_int_t since the
    // RHS products in the calculations below may result in global range values.
    const global_int_t nx  = Ageom->nx;
    const global_int_t ny  = Ageom->ny;
    const global_int_t nz  = Ageom->nz;
    const global_int_t npx = Ageom->npx;
    const global_int_t npy = Ageom->npy;
    const global_int_t ipx = Ageom->ipx;
    const global_int_t ipy = Ageom->ipy;
    const global_int_t ipz = Ageom->ipz;
    const global_int_t gnx = nx * npx;
    const global_int_t gny = ny * npy;
    //!< global-to-local mapping
    auto &globalToLocalMap = A.globalToLocalMap;
    //
    for (local_int_t iz = 0; iz < nz; iz++) {
        global_int_t giz = ipz*nz+iz;
        for (local_int_t iy = 0; iy < ny; iy++) {
            global_int_t giy = ipy*ny+iy;
            for (local_int_t ix = 0; ix < nx; ix++) {
                global_int_t gix = ipx*nx+ix;
                local_int_t currentLocalRow = iz*nx*ny+iy*nx+ix;
                global_int_t currentGlobalRow = giz*gnx*gny+giy*gnx+gix;
                globalToLocalMap[currentGlobalRow] = currentLocalRow;
            }
        }
    }
}

/**
 *
 */
inline void
SetupGhostArrays(
    SparseMatrix &A,
    Array<floatType> &x,
    LegionRuntime::HighLevel::Context &ctx,
    LegionRuntime::HighLevel::HighLevelRuntime *lrt
) {
    // Make sure that we aren't doing this again for something that already has
    // the ghosts setup.
    assert(!x.hasGhosts());
    //
    const SparseMatrixScalars *const Asclrs = A.sclrs->data();
    const int nNeighbors = Asclrs->numberOfSendNeighbors;

    for (int n = 0; n < nNeighbors; ++n) {
        //
        auto xis = x.logicalRegion.get_index_space();
        auto xip = lrt->get_index_partition(ctx, xis, 0 /* color */);
        auto xlp = lrt->get_logical_partition(ctx, x.logicalRegion, xip);
        LogicalRegion xSubReg = lrt->get_logical_subregion_by_color(
            ctx,
            xlp,
            DomainPoint::from_point<1>(n + 1) // First is private.
        );
        LogicalArray<floatType> dstArray(xSubReg, ctx, lrt);
        dstArray.setParentLogicalRegion(x.logicalRegion);
        // Cache in x.
        x.ghosts.push_back(dstArray);
    }
}
