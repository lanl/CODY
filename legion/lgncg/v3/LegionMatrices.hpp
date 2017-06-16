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
    }
};
