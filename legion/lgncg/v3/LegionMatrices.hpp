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
#include "MGData.hpp"

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
    //Max number of non-zero elements in any row
    local_int_t maxNonzerosPerRow = 0;
    //total number of matrix rows across all processes
    global_int_t totalNumberOfRows = 0;
    //total number of matrix nonzeros across all processes
    global_int_t totalNumberOfNonzeros = 0;
    //number of rows local to this process
    local_int_t localNumberOfRows = 0;
    //number of columns local to this process
    local_int_t localNumberOfColumns = 0;
    //number of nonzeros local to this process
    local_int_t localNumberOfNonzeros = 0;
    //number of entries that are external to this process
    local_int_t numberOfExternalValues = 0;
    //number of neighboring processes that will be send local data
    int numberOfSendNeighbors = 0;
    //number of neighboring processes that i'll get data from
    int numberOfRecvNeighbors = 0;
    //total number of entries to be sent
    local_int_t totalToBeSent = 0;
    // size of serialized buffer used for Synchronizers.
    size_t sizeofSynchronizersBuffer = 0;
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
                         &localToGlobalMap
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
        const auto size = geom.size;
        const auto globalXYZ = getGlobalXYZ(geom);
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
    }

};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
struct SparseMatrix {
    // Geometry info for this instance.
    Geometry *geom = nullptr;
    // Container for all scalar values.
    SparseMatrixScalars *sclrs = nullptr;
    //
    char *nonzerosInRow = nullptr;
    // Flattened to 1D from 2D.
    global_int_t *mtxIndG = nullptr;
    // Flattened to 1D from 2D.
    local_int_t *mtxIndL = nullptr;
    // Flattened to 1D from 2D.
    floatType *matrixValues = nullptr;
    //
    floatType *matrixDiagonal = nullptr;
    //
    global_int_t *localToGlobalMap = nullptr;

protected:
    // Number of region entries.
    size_t mNRegionEntries = 0;

public:

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
        mUnpack(regions, baseRID, ctx, runtime);
        mVerifyUnpack();
    }

    /**
     *
     */
    size_t
    nRegionEntries(void) { return mNRegionEntries; }
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
        geom = Item<Geometry>(regions[cid++], ctx, rt).data();
        //
        sclrs = Item<SparseMatrixScalars>(regions[cid++], ctx, rt).data();
        //
        nonzerosInRow = Array<char>(regions[cid++], ctx, rt).data();
        //
        mtxIndG = Array<global_int_t>(regions[cid++], ctx, rt).data();
        //
        mtxIndL = Array<local_int_t>(regions[cid++], ctx, rt).data();
        //
        matrixValues = Array<floatType>(regions[cid++], ctx, rt).data();
        //
        matrixDiagonal = Array<floatType>(regions[cid++], ctx, rt).data();
        //
        localToGlobalMap = Array<global_int_t>(regions[cid++], ctx, rt).data();
        // Calculate number of region entries for this structure.
        mNRegionEntries = cid - baseRID;
    }

    /**
     *
     */
    void
    mVerifyUnpack(void) {
        assert(geom);
        assert(sclrs);
        assert(nonzerosInRow);
        assert(mtxIndG);
        assert(mtxIndL);
        assert(matrixValues);
        assert(matrixDiagonal);
        assert(localToGlobalMap);
    }
};
