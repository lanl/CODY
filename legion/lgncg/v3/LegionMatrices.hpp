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
#include <set>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
struct SparseMatrixScalars {
    // Max number of non-zero elements in any row.
    local_int_t maxNonzerosPerRow = 0;
    // Total number of matrix rows across all processes.
    global_int_t totalNumberOfRows = 0;
    // Total number of matrix non-zeros across all processes.
    global_int_t totalNumberOfNonzeros = 0;
    // Number of rows local to this process.
    local_int_t localNumberOfRows = 0;
    // Number of columns local to this process.
    local_int_t localNumberOfColumns = 0;
    // Number of non-zeros local to this process.
    global_int_t localNumberOfNonzeros = 0;
    // Number of entries that are external to this process.
    local_int_t numberOfExternalValues = 0;
    // Number of neighboring processes that will be sent local data.
    int numberOfSendNeighbors = 0;
    // Number of neighboring processes from which I'll receive data.
    int numberOfRecvNeighbors = 0;
    //Total number of entries to be sent.
    local_int_t totalToBeSent = 0;
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
    LogicalArray< DynColl<global_int_t> > dcAllRedSumGI;
    LogicalArray< DynColl<floatType> > dcAllRedSumFT;
    // Neighboring processes.
    LogicalArray<int> neighbors;
    // Number of items that will be sent on a per neighbor basis.
    LogicalArray<local_int_t> sendLength;
    // Number of items that will be received on a per neighbor basis.
    LogicalArray<local_int_t> recvLength;
    // Synchronization structures.
    LogicalArray<Synchronizers> synchronizers;
    ////////////////////////////////////////////////////////////////////////////
    // Vector index is for a given shard that is sharing pull region info.
    // Innermost vector is for neighboring regions that we are sharing.
    ////////////////////////////////////////////////////////////////////////////
    // TODO Free up resources properly in destructor.
    std::vector< std::vector< LogicalArray<floatType> *> > srcSharedRegions;
    // Similar structure the pullers will use to setup RegionRequirements.
    std::vector< std::vector< LogicalArray<floatType> *> > dstSharedRegions;

protected:
    // Number of shards used for SparseMatrix decomposition.
    int mSize = 0;
    //
    bool mSharedRegionsPopulated = false;

    /**
     * Order matters here. If you update this, also update unpack.
     */
    void
    mPopulateRegionList(void) {
        mLogicalItems = {&geoms,
                         &sclrs,
                         &nonzerosInRow,
                         &mtxIndG,
                         &mtxIndL,
                         &matrixValues,
                         &matrixDiagonal,
                         &localToGlobalMap,
                         &dcAllRedSumGI,
                         &dcAllRedSumFT,
                         &neighbors,
                         &sendLength,
                         &recvLength,
                         &synchronizers
        };
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
    mPopulateSharedRegions(
        Context ctx,
        HighLevelRuntime *lrt
    ) {
        using namespace std;
        // FIXME ugly.
        const int maxNumNeighbors = HPCG_STENCIL - 1;
        //
        Array<SparseMatrixScalars> aSparseMatrixScalars(
            sclrs.mapRegion(RO_E, ctx, lrt), ctx, lrt
        );
        SparseMatrixScalars *sclrsd = aSparseMatrixScalars.data();
        assert(sclrsd);
        //
        Array<int> aNeighbors(
            neighbors.mapRegion(RO_E, ctx, lrt), ctx, lrt
        );
        assert(aNeighbors.data());
        // For convenience we'll interpret this as a 2D array.
        Array2D<int> neighborsd(
            mSize, maxNumNeighbors, aNeighbors.data()
        );
        //
        Array<local_int_t> aSendLengths(
            sendLength.mapRegion(RO_E, ctx, lrt), ctx, lrt
        );
        assert(aSendLengths.data());
        Array2D<local_int_t> sendLengthsd(
            mSize, maxNumNeighbors, aSendLengths.data()
        );
        // Iterate over all shards and populate table that maps task IDs to
        // their index into the neighbor list.
        vector< map<int, int> > tidToNIdx(mSize);
        for (int shard = 0; shard < mSize; ++shard) {
            // Get total number of neighbors this shard has
            const SparseMatrixScalars &myScalars = sclrsd[shard];
            const int nNeighbors = myScalars.numberOfSendNeighbors;
            //
            for (int n = 0; n < nNeighbors; ++n) {
                int tid = neighborsd(shard, n);
                tidToNIdx[shard][tid] = n;
            }
        }
        // We are going to need mSize slots for the vectors.
        srcSharedRegions.resize(mSize);
        dstSharedRegions.resize(mSize);
        // Figure out how many receive neighbors each shard has.
        for (int shard = 0; shard < mSize; ++shard) {
            const int nRecvNeighbors = tidToNIdx[shard].size();
            dstSharedRegions[shard].resize(nRecvNeighbors);
        }
        //
        for (int shard = 0; shard < mSize; ++shard) {
            const int nNeighbors = sclrsd[shard].numberOfSendNeighbors;
            for (int n = 0; n < nNeighbors; ++n) {
                const int nid = neighborsd(shard, n);
                auto *sa = new LogicalArray<floatType>();
                string rName = "A-pullRegion-SourceRank=" + to_string(shard)
                             + "DestinationRank=" + to_string(nid);
                sa->allocate(
                    rName,
                    sendLengthsd(shard, n),
                    ctx,
                    lrt
                );
                srcSharedRegions[shard].push_back(sa);
                dstSharedRegions[nid][tidToNIdx[nid][shard]] = sa;
            }
        }
        //
        sclrs.unmapRegion(ctx, lrt);
        neighbors.unmapRegion(ctx, lrt);
        sendLength.unmapRegion(ctx, lrt);
        //
        mSharedRegionsPopulated = true;
    }

    /**
     *
     */
    void
    intent(
        Legion::PrivilegeMode privMode,
        Legion::CoherenceProperty cohProp,
        ItemFlags iFlags,
        int shard,
        Legion::TaskLauncher &launcher,
        LegionRuntime::HighLevel::Context ctx,
        LegionRuntime::HighLevel::HighLevelRuntime *lrt
    ) {
        LogicalMultiBase::intent(privMode, cohProp, shard, launcher, ctx, lrt);
        //
        if (withGhosts(iFlags)) {
            if (!mSharedRegionsPopulated) {
                mPopulateSharedRegions(ctx, lrt);
            }
            // FIXME ugly.
            const int maxNumNeighbors = HPCG_STENCIL - 1;
            //
            Array<SparseMatrixScalars> aSparseMatrixScalars(
                sclrs.mapRegion(RO_E, ctx, lrt), ctx, lrt
            );
            SparseMatrixScalars *sclrsd = aSparseMatrixScalars.data();
            assert(sclrsd);
            //
            Array<int> aNeighbors(
                neighbors.mapRegion(RO_E, ctx, lrt), ctx, lrt
            );
            assert(aNeighbors.data());
            // For convenience we'll interpret this as a 2D array.
            Array2D<int> neighborsd(
                mSize, maxNumNeighbors, aNeighbors.data()
            );
            const int nNeighbors = sclrsd[shard].numberOfSendNeighbors;
            // First nNeighbors regions are the ones I'm populating. That is,
            // I'm the source for the values and my neighbors pull from those.
            for (int n = 0; n < nNeighbors; ++n) {
                auto *ap = srcSharedRegions[shard][n];
                auto &lr = ap->logicalRegion;
                launcher.add_region_requirement(
                    RegionRequirement(
                        lr,
                        READ_WRITE,
                        SIMULTANEOUS,
                        lr
                    ).add_flags(NO_ACCESS_FLAG)
                ).add_field(ap->fid);
            }
            // Next nNeighbors regions are the ones I need for my computation.
            // That is, they are the once that other tasks have populated.
            const int nRecvNeighbors = dstSharedRegions[shard].size();
            for (int n = 0; n < nRecvNeighbors; ++n) {
                auto *ap = dstSharedRegions[shard][n];
                auto &lr = ap->logicalRegion;
                launcher.add_region_requirement(
                    RegionRequirement(
                        lr,
                        READ_ONLY,
                        SIMULTANEOUS,
                        lr
                    )
                ).add_field(ap->fid);
            }
            //
            sclrs.unmapRegion(ctx, lrt);
            neighbors.unmapRegion(ctx, lrt);
        }
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
    ) {
        #define aalloca(sName, size, ctx, rtp)                                 \
        do {                                                                   \
            sName.allocate(name + "-" #sName, size, ctx, rtp);                 \
        } while(0)

        mSize = geom.size;
        const auto globalXYZ   = getGlobalXYZ(geom);
        const auto stencilSize = geom.stencilSize;

        aalloca(geoms, mSize, ctx, lrt);
        aalloca(sclrs, mSize, ctx, lrt);
        //
        aalloca(nonzerosInRow, globalXYZ, ctx, lrt);
        // Flattened to 1D from 2D.
        aalloca(mtxIndG, globalXYZ * stencilSize, ctx, lrt);
        // Flattened to 1D from 2D.
        aalloca(mtxIndL, globalXYZ * stencilSize, ctx, lrt);
        // Flattened to 1D from 2D.
        aalloca(matrixValues, globalXYZ * stencilSize, ctx, lrt);
        // 2D thing in reference implementation, but not needed (1D suffices).
        aalloca(matrixDiagonal, globalXYZ, ctx, lrt);
        //
        aalloca(localToGlobalMap, globalXYZ, ctx, lrt);
        //
        aalloca(dcAllRedSumGI, mSize, ctx, lrt);
        aalloca(dcAllRedSumFT, mSize, ctx, lrt);
        //
        const int maxNumNeighbors = geom.stencilSize - 1;
        // Each task will have at most 26 neighbors.
        aalloca(neighbors,  mSize * maxNumNeighbors, ctx, lrt);
        aalloca(sendLength, mSize * maxNumNeighbors, ctx, lrt);
        aalloca(recvLength, mSize * maxNumNeighbors, ctx, lrt);
        //
        aalloca(synchronizers, mSize, ctx, lrt);

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
    ) {
        for (auto *i : mLogicalItems) {
            i->partition(nParts, ctx, lrt);
        }
        // For the DynamicCollectives we need partition info before population.
        const auto nArrivals = nParts;
        DynColl<global_int_t> dynColGI(INT_REDUCE_SUM_TID, nArrivals);
        mPopulateDynamicCollectives(dcAllRedSumGI, dynColGI, ctx, lrt);
        //
        DynColl<floatType> dynColFT(FLOAT_REDUCE_SUM_TID, nArrivals);
        mPopulateDynamicCollectives(dcAllRedSumFT, dynColFT, ctx, lrt);
        // Just pick a structure that has a representative launch domain.
        launchDomain = geoms.launchDomain;
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

private:

    /**
     *
     */
    template <typename TYPE>
    void
    mPopulateDynamicCollectives(
        LogicalArray< DynColl<TYPE> > &targetLogicalArray,
        DynColl<TYPE> &dynCol,
        LegionRuntime::HighLevel::Context ctx,
        LegionRuntime::HighLevel::HighLevelRuntime *lrt
    ) {
        Array< DynColl<TYPE> > dcs(
            targetLogicalArray.mapRegion(RW_E, ctx, lrt), ctx, lrt
        );
        //
        DynColl<TYPE> *dcsd = dcs.data();
        assert(dcsd);
        //
        dynCol.dc = lrt->create_dynamic_collective(
            ctx,
            dynCol.nArrivals /* Number of arrivals. */,
            dynCol.tid,
            &dynCol.localBuffer,
            sizeof(dynCol.localBuffer)
        );
        // Replicate
        for (int64_t i = 0; i < dynCol.nArrivals; ++i) {
            dcsd[i] = dynCol;
        }
        // Done, so unmap.
        targetLogicalArray.unmapRegion(ctx, lrt);
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
    Item< DynColl<global_int_t> > *dcAllRedSumGI = nullptr;
    //
    Item< DynColl<floatType> > *dcAllRedSumFT = nullptr;
    //
    Array<int> *neighbors = nullptr;
    //
    Array<local_int_t> *sendLength = nullptr;
    //
    Array<local_int_t> *recvLength = nullptr;
    //
    Item<Synchronizers> *synchronizers = nullptr;
    ////////////////////////////////////////////////////////////////////////////
    // Task-launch-specific structures.
    ////////////////////////////////////////////////////////////////////////////
    // Global to local mapping. NOTE: only valid after a call to
    // PopulateGlobalToLocalMap.
    std::map< global_int_t, local_int_t > globalToLocalMap;
    // Only valid after a call to SetupHalo.
    local_int_t *elementsToSend = nullptr;
    // A mapping between neighbor IDs and their regions.
    std::map<int, PhysicalRegion> nidToPullRegion;
    // Pull regions that I populate for consumption by other tasks.
    std::vector< Array<floatType> *> pullBuffers;

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
        // Number of regions consumed.
        int cid = baseRID;
        // Setup my push Arrays.
        for (int n = 0; n < sclrsd->numberOfSendNeighbors; ++n) {
            pullBuffers.push_back(
                new Array<floatType>(regions[cid++], ctx, runtime)
            );
        }
        // Get neighbor regions that I will pull from.
        for (int n = 0; n < sclrsd->numberOfRecvNeighbors; ++n) {
            const int nid = nd[n];
            nidToPullRegion[nid] = regions[cid++];
        }
        // Return number of regions that we have consumed.
        return cid - baseRID;
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
        delete nonzerosInRow;
        delete mtxIndG;
        delete mtxIndL;
        delete matrixValues;
        delete matrixDiagonal;
        delete localToGlobalMap;
        delete dcAllRedSumGI;
        delete dcAllRedSumFT;
        delete neighbors;
        delete sendLength;
        delete recvLength;
        delete synchronizers;
        // Task-local allocation of non-region memory.
        if (elementsToSend) delete[] elementsToSend;
        for (auto *i : pullBuffers) delete i;
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
        dcAllRedSumGI = new Item< DynColl<global_int_t> >(regions[cid++], ctx, rt);
        assert(dcAllRedSumGI->data());
        //
        dcAllRedSumFT = new Item< DynColl<floatType> >(regions[cid++], ctx, rt);
        assert(dcAllRedSumFT->data());
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
    Item< DynColl<global_int_t> > dc(regions[0], ctx, runtime);
    return dc.data()->localBuffer;
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
    Item< DynColl<floatType> > dc(regions[0], ctx, runtime);
    return dc.data()->localBuffer;
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
    LegionRuntime::HighLevel::Context ctx,
    LegionRuntime::HighLevel::HighLevelRuntime *lrt
) {
    // Make sure that we aren't doing this again for something that already has
    // the ghosts setup.
    assert(!x.hasGhosts());
    //
    const SparseMatrixScalars *const Asclrs = A.sclrs->data();
    assert(Asclrs);
    const int nNeighbors = Asclrs->numberOfSendNeighbors;

    for (int n = 0; n < nNeighbors; ++n) {
        auto xis = x.logicalRegion.get_index_space();
        auto xip = lrt->get_index_partition(ctx, xis, 0 /* color */);
        auto xlp = lrt->get_logical_partition(ctx, x.logicalRegion, xip);
        LogicalRegion xSubReg = lrt->get_logical_subregion_by_color(
            ctx,
            xlp,
            DomainPoint::from_point<1>(n + 1) // First is private.
        );
        auto *dst = new LogicalArray<floatType>(xSubReg, ctx, lrt);
        dst->setParentLogicalRegion(x.logicalRegion);
        // Cache in x.
        x.ghosts.push_back(dst);
    }
}
