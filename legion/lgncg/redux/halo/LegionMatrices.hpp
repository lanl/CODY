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
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
struct LogicalSparseMatrix {
    // launch domain
    LegionRuntime::HighLevel::Domain launchDomain;
    //
    LogicalArray<Geometry> geometries;
    //
    LogicalArray<SparseMatrixScalars> localData;
    // The number of nonzeros in a row will always be 27 or fewer
    LogicalArray<LogicalRegion> lrNonzerosInRow;
    // Logical regions that point to regions of matrix indices as global values
    LogicalArray<LogicalRegion> lrMtxIndG;
    // Logical regions that point to matrix indices as local values
    LogicalArray<LogicalRegion> lrMtxIndL;
    // Logical regions that point to values of matrix entries
    LogicalArray<LogicalRegion> lrMatrixValues;
    // Logical regions that point to values of matrix diagonal entries
    LogicalArray<LogicalRegion> lrMatrixDiagonal;
    // Logical regions that point to local-to-global mapping arrays
    LogicalArray<LogicalRegion> lrLocalToGlobalMap;
    // Logical regions that point to serialized global-to-local mapping data
    LogicalArray<LogicalRegion> lrGlobalToLocalMap;
    // Logical regions that point to array containing elements to send to
    // neighboring processes
    LogicalArray<LogicalRegion> lrElementsToSend;
    // Logical regions that point to array of neighboring processes
    LogicalArray<LogicalRegion> lrNeighbors;
    //Logical regions that point to lengths of messages received from
    //neighboring processes SKG TODO RM? may not be needed
    LogicalArray<LogicalRegion> lrReceiveLength;
    //Logical regions that point to lengths of messages sent to neighboring
    //processes
    LogicalArray<LogicalRegion> lrSendLength;
    // Holds serialization data for Synchronizers.
    LogicalArray<char> synchronizersData;
    //
    std::vector<PhaseBarriers> ownerPhaseBarriers;
    // For each color (shard), keep track of its neighbor PhaseBarriers; The
    // vector index is the target task ID that does not own the PhaseBarriers,
    // but rather uses them for synchronization. The mapping is between neighbor
    // IDs (PhaseBarrier owners) and the PhaseBarriers that a particular ID is
    // sharing.
    std::vector< std::map< int, std::vector<PhaseBarriers> > > neighborPhaseBarriers;
    // Coarse grid matrix
    mutable struct SparseMatrix_STRUCT *Ac = nullptr;
    // Pointer to the coarse level data for this fine matrix
    mutable MGData *mgData;
    // Given the characteristics of this matrix, provide the required length of
    // vectors that will include halo cells.
    size_t requiredVectorLen = 0;
    // Given the characteristics of this matrix, provide the required
    // partitioning information of vectors that will interact with this matrix.
    std::vector<size_t> targetVectorPartLens;

    /**
     *
     */
    LogicalSparseMatrix(void) = default;

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
        //
                geometries.allocate(size, ctx, lrt);
                 localData.allocate(size, ctx, lrt);
           lrNonzerosInRow.allocate(size, ctx, lrt);
                 lrMtxIndG.allocate(size, ctx, lrt);
                 lrMtxIndL.allocate(size, ctx, lrt);
            lrMatrixValues.allocate(size, ctx, lrt);
          lrMatrixDiagonal.allocate(size, ctx, lrt);
        lrLocalToGlobalMap.allocate(size, ctx, lrt);
        lrGlobalToLocalMap.allocate(size, ctx, lrt);
          lrElementsToSend.allocate(size, ctx, lrt);
               lrNeighbors.allocate(size, ctx, lrt);
           lrReceiveLength.allocate(size, ctx, lrt);
              lrSendLength.allocate(size, ctx, lrt);
           // Dummy
         synchronizersData.allocate(size, ctx, lrt);
    }

    /**
     *
     */
    void
    allocateSynchronizers(
        size_t n,
        LegionRuntime::HighLevel::Context &ctx,
        LegionRuntime::HighLevel::HighLevelRuntime *lrt
    ) {
        // Deallocate dummy
        synchronizersData.deallocate(ctx, lrt);
        synchronizersData.allocate(n, ctx, lrt);
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
                geometries.partition(nParts, ctx, lrt);
                 localData.partition(nParts, ctx, lrt);
           lrNonzerosInRow.partition(nParts, ctx, lrt);
                 lrMtxIndG.partition(nParts, ctx, lrt);
                 lrMtxIndL.partition(nParts, ctx, lrt);
            lrMatrixValues.partition(nParts, ctx, lrt);
          lrMatrixDiagonal.partition(nParts, ctx, lrt);
        lrLocalToGlobalMap.partition(nParts, ctx, lrt);
        lrGlobalToLocalMap.partition(nParts, ctx, lrt);
          lrElementsToSend.partition(nParts, ctx, lrt);
               lrNeighbors.partition(nParts, ctx, lrt);
           lrReceiveLength.partition(nParts, ctx, lrt);
              lrSendLength.partition(nParts, ctx, lrt);
         synchronizersData.partition(nParts, ctx, lrt);
        // just pick a structure that has a representative launch domain.
        launchDomain = geometries.launchDomain;
    }

    /**
     *
     */
    void
    partitionSynchronizers(
        std::vector<size_t> partLens,
        LegionRuntime::HighLevel::Context &ctx,
        LegionRuntime::HighLevel::HighLevelRuntime *lrt
    ) {
        synchronizersData.partition(partLens, ctx, lrt);
    }

    /**
     * Cleans up and returns all allocated resources.
     */
    void
    deallocate(
        LegionRuntime::HighLevel::Context &ctx,
        LegionRuntime::HighLevel::HighLevelRuntime *lrt
    ) {
                geometries.deallocate(ctx, lrt);
                 localData.deallocate(ctx, lrt);
           lrNonzerosInRow.deallocate(ctx, lrt);
                 lrMtxIndG.deallocate(ctx, lrt);
                 lrMtxIndL.deallocate(ctx, lrt);
            lrMatrixValues.deallocate(ctx, lrt);
          lrMatrixDiagonal.deallocate(ctx, lrt);
        lrLocalToGlobalMap.deallocate(ctx, lrt);
        lrGlobalToLocalMap.deallocate(ctx, lrt);
          lrElementsToSend.deallocate(ctx, lrt);
               lrNeighbors.deallocate(ctx, lrt);
           lrReceiveLength.deallocate(ctx, lrt);
              lrSendLength.deallocate(ctx, lrt);
         synchronizersData.deallocate(ctx, lrt);
    }
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
class SparseMatrix {
protected:
    int cNItemsToUnpack = 0;
    // Physical Item Container
    struct PIC {
        Item<Geometry> geom;
        Item<SparseMatrixScalars> localData;
        Item<LogicalRegion> nonzerosInRow;
        Item<LogicalRegion> mtxIndG;
        Item<LogicalRegion> mtxIndL;
        Item<LogicalRegion> matrixValues;
        Item<LogicalRegion> matrixDiagonal;
        Item<LogicalRegion> localToGlobalMap;
        Item<LogicalRegion> globalToLocalMap;
        Item<LogicalRegion> elementsToSend;
        Item<LogicalRegion> neighbors;
        Item<LogicalRegion> receiveLength;
        Item<LogicalRegion> sendLength;
        Array<char>         synchronizersData;
    };
public:
    //
    local_int_t maxNonzerosPerRow;
    //
    global_int_t totalNumberOfRows = 0;
    //
    global_int_t totalNumberOfNonzeros = 0;
    //
    local_int_t localNumberOfRows = 0;
    //
    local_int_t localNumberOfColumns = 0;
    //
    local_int_t localNumberOfNonzeros = 0;
    //
    local_int_t numberOfExternalValues = 0;
    //
    int numberOfSendNeighbors;
    //
    int numberOfRecvNeighbors;
    //
    local_int_t totalToBeSent;
    //
    Geometry *geom = nullptr;
    //
    SparseMatrixScalars *localData = nullptr;
    //
    char *nonzerosInRow = nullptr;
    //
    global_int_t *mtxIndG = nullptr;
    //
    local_int_t *mtxIndL = nullptr;
    //
    floatType *matrixValues = nullptr;
    //
    floatType *matrixDiagonal = nullptr;
    //
    global_int_t *localToGlobalMap = nullptr;
    //
    local_int_t *elementsToSend = nullptr;
    //
    int *neighbors = nullptr;
    //
    local_int_t *receiveLength = nullptr;
    //
    local_int_t *sendLength = nullptr;
    //
    Synchronizers *synchronizers = nullptr;
    //
    std::map< global_int_t, local_int_t > globalToLocalMap;
    //
    PIC pic;

    /**
     *
     */
    SparseMatrix(void) = default;

    /**
     *
     */
    SparseMatrix(
        Geometry *geom,
        SparseMatrixScalars *localData,
        char *nonzerosInRow,
        global_int_t *mtxIndG,
        local_int_t *mtxIndL,
        floatType *matrixValues,
        floatType *matrixDiagonal,
        global_int_t *localToGlobalMap,
        local_int_t *elementsToSend,
        int *neighbors,
        local_int_t *receiveLength,
        local_int_t *sendLength,
        Synchronizers *synchronizers
    ) {
        this->geom = geom;
        this->localData = localData;
        this->nonzerosInRow = nonzerosInRow;
        this->mtxIndG = mtxIndG;
        this->mtxIndL = mtxIndL;
        this->matrixValues = matrixValues;
        this->matrixDiagonal = matrixDiagonal;
        this->localToGlobalMap = localToGlobalMap;
        this->elementsToSend = elementsToSend;
        this->neighbors = neighbors;
        this->receiveLength = receiveLength;
        this->sendLength = sendLength;
        this->synchronizers = synchronizers;
    }

    /**
     *
     */
    SparseMatrix(
        const std::vector<LegionRuntime::HighLevel::PhysicalRegion> &regions,
        size_t baseRegionID,
        LegionRuntime::HighLevel::Context ctx,
        LegionRuntime::HighLevel::Runtime *runtime
    ) {
        size_t crid = baseRegionID;
        //
        pic.geom = Item<Geometry>(regions[crid++], ctx, runtime);
        geom = pic.geom.data();
        assert(geom);
        //
        pic.localData = Item<SparseMatrixScalars>(
            regions[crid++], ctx, runtime
        );
        localData = pic.localData.data();
        assert(localData);
        //
        pic.nonzerosInRow = Item<LogicalRegion>(
            regions[crid++], ctx, runtime
        );
        //
        pic.mtxIndG = Item<LogicalRegion>(
            regions[crid++], ctx, runtime
        );
        //
        pic.mtxIndL = Item<LogicalRegion>(
            regions[crid++], ctx, runtime
        );
        //
        pic.matrixValues = Item<LogicalRegion>(
            regions[crid++], ctx, runtime
        );
        //
        pic.matrixDiagonal = Item<LogicalRegion>(
            regions[crid++], ctx, runtime
        );
        //
        pic.localToGlobalMap = Item<LogicalRegion>(
            regions[crid++], ctx, runtime
        );
        //
        pic.globalToLocalMap = Item<LogicalRegion>(
            regions[crid++], ctx, runtime
        );
        //
        pic.elementsToSend = Item<LogicalRegion>(
            regions[crid++], ctx, runtime
        );
        //
        pic.neighbors = Item<LogicalRegion>(
            regions[crid++], ctx, runtime
        );
        //
        pic.receiveLength = Item<LogicalRegion>(
            regions[crid++], ctx, runtime
        );
        //
        pic.sendLength = Item<LogicalRegion>(
            regions[crid++], ctx, runtime
        );
        pic.synchronizersData = Array<char>(
            regions[crid++], ctx, runtime
        );
        //
        cNItemsToUnpack = crid - baseRegionID;
    }

    /**
     *
     */
    int
    nRegionEntries(void) { return cNItemsToUnpack; }

    /**
     *
     */
    void
    refreshMembers(const SparseMatrix &sm) {
        geom                   = sm.geom;
        localData              = sm.localData;
        //
        maxNonzerosPerRow      = localData->maxNonzerosPerRow;
        totalNumberOfRows      = localData->totalNumberOfRows;
        totalNumberOfNonzeros  = localData->totalNumberOfNonzeros;
        localNumberOfRows      = localData->localNumberOfRows;
        localNumberOfColumns   = localData->localNumberOfColumns;
        localNumberOfNonzeros  = localData->localNumberOfNonzeros;
        numberOfExternalValues = localData->numberOfExternalValues;
        numberOfSendNeighbors  = localData->numberOfSendNeighbors;
        numberOfRecvNeighbors  = localData->numberOfRecvNeighbors;
        totalToBeSent          = localData->totalToBeSent;
        //
        nonzerosInRow          = sm.nonzerosInRow;
        mtxIndG                = sm.mtxIndG;
        mtxIndL                = sm.mtxIndL;
        matrixValues           = sm.matrixValues;
        matrixDiagonal         = sm.matrixDiagonal;
        localToGlobalMap       = sm.localToGlobalMap;
        elementsToSend         = sm.elementsToSend;
        neighbors              = sm.neighbors;
        receiveLength          = sm.receiveLength;
        sendLength             = sm.sendLength;
        synchronizers          = sm.synchronizers;
    }
};
