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
#include <utility> // For pair

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
    //Number of neighboring processes that I'll get data from.
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
    //
    PhaseBarriers myPhaseBarriers;
    //
    std::map< int, std::vector<PhaseBarriers> > neighborPhaseBarriers;

    /**
     *
     */
    Synchronizers(void) = default;

    /**
     *
     */
    template <class Archive>
    void
    serialize(Archive &ar)
    {
        ar(myPhaseBarriers, neighborPhaseBarriers);
    }

    /**
     *
     */
    static void
    deserialize(
        char *rawData,
        size_t rawDataSizeInB,
        Synchronizers &outRes
    ) {
        // Deserialize argument data
        std::stringstream solvArgsSS;
        solvArgsSS.write(rawData, rawDataSizeInB);
        {   // Scoped to guarantee flushing, etc.
            cereal::BinaryInputArchive ia(solvArgsSS);
            ia(outRes);
        }
    }
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
struct LogicalSparseMatrix : public LogicalMultiBase {
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
    LogicalArray<DynamicCollective> dcAllreduceSum;
    // Neighboring processes.
    LogicalArray<int> neighbors;

protected:
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
                         &dcAllreduceSum,
                         &neighbors
        };
    }

public:
    /**
     *
     */
    LogicalSparseMatrix(void) {
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
        const int size         = geom.size;
        const auto globalXYZ   = getGlobalXYZ(geom);
        const auto stencilSize = geom.stencilSize;

        geoms.allocate(size, ctx, lrt);
        sclrs.allocate(size, ctx, lrt);
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
        dcAllreduceSum.allocate(size, ctx, lrt);
        //
        const int maxNumNeighbors = geom.stencilSize - 1;
        // Each task will have at most 26 neighbors.
        neighbors.allocate(size * maxNumNeighbors, ctx, lrt);
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
        nonzerosInRow.partition(nParts, ctx, lrt);
        mtxIndG.partition(nParts, ctx, lrt);
        mtxIndL.partition(nParts, ctx, lrt);
        matrixValues.partition(nParts, ctx, lrt);
        matrixDiagonal.partition(nParts, ctx, lrt);
        localToGlobalMap.partition(nParts, ctx, lrt);
        // For the DynamicCollectives we need partition info before population.
        dcAllreduceSum.partition(nParts, ctx, lrt);
        // Really out of place here, but it works, folks...
        mPopulateDynamicCollectives(nParts, ctx, lrt);
        neighbors.partition(nParts, ctx, lrt);
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
        mtxIndG.deallocate(ctx, lrt);
        mtxIndL.deallocate(ctx, lrt);
        nonzerosInRow.deallocate(ctx, lrt);
        matrixValues.deallocate(ctx, lrt);
        matrixDiagonal.deallocate(ctx, lrt);
        localToGlobalMap.deallocate(ctx, lrt);
        dcAllreduceSum.deallocate(ctx, lrt);
        neighbors.deallocate(ctx, lrt);
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
    // Global to local mapping. NOTE: only valid after a call to
    // PopulateGlobalToLocalMap.
    std::map< global_int_t, local_int_t > globalToLocalMap;

public:

    /**
     *
     */
    SparseMatrix(void) = default;

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
        delete dcAllreduceSum;
        delete neighbors;
    }

    /**
     *
     */
    SparseMatrix(
        const std::vector<PhysicalRegion> &regions,
        size_t baseRID,
        Context ctx,
        HighLevelRuntime *runtime
    ) {
        mUnpack(regions, baseRID, ctx, runtime);
        mVerifyUnpack();
    }

protected:
    /**
     * MUST MATCH PACK ORDER IN mPopulateRegionList!
     */
    void
    mUnpack(
        const std::vector<PhysicalRegion> &regions,
        size_t baseRID,
        Context ctx,
        HighLevelRuntime *rt
    ) {
        size_t cid = baseRID;
        // Populate members from physical regions.
        geom = new Item<Geometry>(regions[cid++], ctx, rt);
        //
        sclrs = new Item<SparseMatrixScalars>(regions[cid++], ctx, rt);
        //
        nonzerosInRow = new Array<char>(regions[cid++], ctx, rt);
        //
        mtxIndG = new Array<global_int_t>(regions[cid++], ctx, rt);
        //
        mtxIndL = new Array<local_int_t>(regions[cid++], ctx, rt);
        //
        matrixValues = new Array<floatType>(regions[cid++], ctx, rt);
        //
        matrixDiagonal = new Array<floatType>(regions[cid++], ctx, rt);
        //
        localToGlobalMap = new Array<global_int_t>(regions[cid++], ctx, rt);
        //
        dcAllreduceSum = new Item<DynamicCollective>(regions[cid++], ctx, rt);
        //
        neighbors = new Array<int>(regions[cid++], ctx, rt);
        // Calculate number of region entries for this structure.
        mNRegionEntries = cid - baseRID;
    }

    /**
     *
     */
    void
    mVerifyUnpack(void) {
        assert(geom->data());
        assert(sclrs->data());
        assert(nonzerosInRow->data());
        assert(mtxIndG->data());
        assert(mtxIndL->data());
        assert(matrixValues->data());
        assert(matrixDiagonal->data());
        assert(localToGlobalMap->data());
        assert(dcAllreduceSum);
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
