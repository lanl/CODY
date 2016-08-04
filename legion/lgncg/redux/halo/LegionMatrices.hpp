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
struct SparseMatrixLocalData {
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
};
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
struct LogicalSparseMatrix {
public:
    // launch domain
    LegionRuntime::HighLevel::Domain launchDomain;
public:
    //
    LogicalArray<Geometry> geometries;
    //
    LogicalArray<SparseMatrixLocalData> localData;
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
    // Coarse grid matrix
    mutable struct SparseMatrix_STRUCT *Ac;
    // Pointer to the coarse level data for this fine matrix
    mutable MGData *mgData;
    //number of entries that are external to this process
    local_int_t numberOfExternalValues;
    //number of neighboring processes that will be send local data
    int numberOfSendNeighbors;
    //total number of entries to be sent
    local_int_t totalToBeSent;
    //elements to send to neighboring processes
    local_int_t *elementsToSend;
    //neighboring processes
    int *neighbors;
    //lenghts of messages received from neighboring processes
    local_int_t *receiveLength;
    //lenghts of messages sent to neighboring processes
    local_int_t *sendLength;
    //send buffer for non-blocking sends
    floatType *sendBuffer;

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
        // just pick a structure that has a representative launch domain.
        launchDomain = geometries.launchDomain;
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
    }
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
class SparseMatrix {
protected:
    static constexpr int cNItemsToUnpack = 9;
    // Physical Item Container
    struct PIC {
        Item<Geometry> geom;
        Item<SparseMatrixLocalData> localData;
        Item<LogicalRegion> nonzerosInRow;
        Item<LogicalRegion> mtxIndG;
        Item<LogicalRegion> mtxIndL;
        Item<LogicalRegion> matrixValues;
        Item<LogicalRegion> matrixDiagonal;
        Item<LogicalRegion> localToGlobalMap;
        Item<LogicalRegion> globalToLocalMap;
    };
public:
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
    //
    Geometry *geom = nullptr;
    //
    SparseMatrixLocalData *localData = nullptr;
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
        SparseMatrixLocalData *localData,
        char *nonzerosInRow,
        global_int_t *mtxIndG,
        local_int_t *mtxIndL,
        floatType *matrixValues,
        floatType *matrixDiagonal,
        global_int_t *localToGlobalMap,
        const std::map< global_int_t, local_int_t > &globalToLocalMap
    ) {
        this->geom = geom;
        this->localData = localData;
        this->nonzerosInRow = nonzerosInRow;
        this->mtxIndG = mtxIndG;
        this->mtxIndL = mtxIndL;
        this->matrixValues = matrixValues;
        this->matrixDiagonal = matrixDiagonal;
        this->localToGlobalMap = localToGlobalMap;
        this->globalToLocalMap = globalToLocalMap;
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
        size_t curRID = baseRegionID;
        //
        pic.geom = Item<Geometry>(regions[curRID++], ctx, runtime);
        geom = pic.geom.data();
        assert(geom);
        //
        pic.localData = Item<SparseMatrixLocalData>(
            regions[curRID++], ctx, runtime
        );
        localData = pic.localData.data();
        assert(localData);
        //
        pic.nonzerosInRow = Item<LogicalRegion>(
            regions[curRID++], ctx, runtime
        );
        //
        pic.mtxIndG = Item<LogicalRegion>(regions[curRID++], ctx, runtime);
        //
        pic.mtxIndL = Item<LogicalRegion>(regions[curRID++], ctx, runtime);
        //
        pic.matrixValues = Item<LogicalRegion>(regions[curRID++], ctx, runtime);
        //
        pic.matrixDiagonal = Item<LogicalRegion>(
            regions[curRID++], ctx, runtime
        );
        //
        pic.localToGlobalMap = Item<LogicalRegion>(
            regions[curRID++], ctx, runtime
        );
        //
        pic.globalToLocalMap = Item<LogicalRegion>(
            regions[curRID++], ctx, runtime
        );
    }

    /**
     *
     */
    static int
    nRegionEntries(void) { return cNItemsToUnpack; }

    /**
     *
     */
    void
    refreshMembers(const SparseMatrix &sm) {
        geom                  = sm.geom;
        localData             = sm.localData;
        totalNumberOfRows     = localData->totalNumberOfRows;
        totalNumberOfNonzeros = localData->totalNumberOfNonzeros;
        localNumberOfRows     = localData->localNumberOfRows;
        localNumberOfColumns  = localData->localNumberOfColumns;
        localNumberOfNonzeros = localData->localNumberOfNonzeros;
        nonzerosInRow         = sm.nonzerosInRow;
        mtxIndG               = sm.mtxIndG;
        mtxIndL               = sm.mtxIndL;
        matrixValues          = sm.matrixValues;
        matrixDiagonal        = sm.matrixDiagonal;
        localToGlobalMap      = sm.localToGlobalMap;
        // TODO consider something else here (memory usage)
        globalToLocalMap      = sm.globalToLocalMap;
    }
};
