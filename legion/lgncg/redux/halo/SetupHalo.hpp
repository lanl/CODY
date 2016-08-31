/**
 * Copyright (c) 2016      Los Alamos National Security, LLC
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
 *
 * LA-CC 10-123
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

#include "LegionMatrices.hpp"

#include "hpcg.hpp"
#include "mytimer.hpp"

#include <vector>
#include <map>
#include <set>
#include <cassert>

/*!
    Reference version of SetupHalo that prepares system matrix data structure
    and creates data necessary for communication of boundary values of this
    process.

    @param[inout] A    The known system matrix

    @see ExchangeHalo
*/
inline void
SetupHalo(
    SparseMatrix &A,
    LegionRuntime::HighLevel::Context ctx,
    LegionRuntime::HighLevel::Runtime *runtime
) {
    using namespace std;
    // Extract Matrix pieces
    local_int_t localNumberOfRows = A.localNumberOfRows;
    char *nonzerosInRow = A.nonzerosInRow;
    // These have already been allocated, so just interpret at 2D array
    Array2D<global_int_t> mtxIndG(
        localNumberOfRows, A.maxNonzerosPerRow , A.mtxIndG
    );
    Array2D<local_int_t> mtxIndL(
        localNumberOfRows, A.maxNonzerosPerRow , A.mtxIndL
    );

    // Scan global IDs of the nonzeros in the matrix.  Determine if the column
    // ID matches a row ID.  If not: 1) We call the ComputeRankOfMatrixRow
    // function, which tells us the rank of the processor owning the row ID.  We
    // need to receive this value of the x vector during the halo exchange.  2)
    // We record our row ID since we know that the other processor will need
    // this value from us, due to symmetry.

    std::map< int, std::set< global_int_t> > sendList, receiveList;
    typedef std::map< int, std::set< global_int_t> >::iterator map_iter;
    typedef std::set<global_int_t>::iterator set_iter;
    std::map< local_int_t, local_int_t > externalToLocalMap;

    for (local_int_t i = 0; i < localNumberOfRows; i++) {
        global_int_t currentGlobalRow = A.localToGlobalMap[i];
        for (int j = 0; j<nonzerosInRow[i]; j++) {
            global_int_t curIndex = mtxIndG(i, j);
            int rankIdOfColumnEntry = ComputeRankOfMatrixRow(
                                          *(A.geom), curIndex
                                      );
            // If column index is not a row index, then it comes from another
            // processor
            if (A.geom->rank != rankIdOfColumnEntry) {
                receiveList[rankIdOfColumnEntry].insert(curIndex);
                // Matrix symmetry means we know the neighbor process wants my
                // value
                sendList[rankIdOfColumnEntry].insert(currentGlobalRow);
            }
        }
    }

    // Count number of matrix entries to send and receive
    local_int_t totalToBeSent = 0;
    for (map_iter curNeighbor = sendList.begin();
         curNeighbor != sendList.end(); ++curNeighbor) {
        totalToBeSent += (curNeighbor->second).size();
    }
    local_int_t totalToBeReceived = 0;
    for (map_iter curNeighbor = receiveList.begin();
         curNeighbor != receiveList.end(); ++curNeighbor) {
        totalToBeReceived += (curNeighbor->second).size();
    }

#ifdef HPCG_DETAILED_DEBUG
    // These are all attributes that should be true, due to symmetry
    HPCG_fout << "totalToBeSent = " << totalToBeSent
              << " totalToBeReceived = " << totalToBeReceived << endl;
    // Number of sent entry should equal number of received
    assert(totalToBeSent == totalToBeReceived);
    // Number of send-to neighbors should equal number of receive-from
    assert(sendList.size() == receiveList.size());
    // Each receive-from neighbor should be a send-to neighbor, and send the
    // same number of entries
    for (map_iter curNeighbor = receiveList.begin();
         curNeighbor != receiveList.end(); ++curNeighbor) {
        assert(sendList.find(curNeighbor->first)!=sendList.end());
        assert(sendList[curNeighbor->first].size() ==
               receiveList[curNeighbor->first].size());
    }
#endif

    // Build the arrays and lists needed by the ExchangeHalo function.
    //double * sendBuffer = new double[totalToBeSent];
    //
    ArrayAllocator<local_int_t> *aaElementsToSend = nullptr;
    local_int_t *elementsToSend = nullptr;
    if (totalToBeSent != 0) {
        aaElementsToSend = new ArrayAllocator<local_int_t>(
            totalToBeSent, WO_E, ctx, runtime
        );
        elementsToSend = aaElementsToSend->data();
        assert(elementsToSend);
    }
    //
    ArrayAllocator<int> *aaNeighbors = nullptr;
    int *neighbors = nullptr;
    if (sendList.size() != 0) {
        aaNeighbors = new ArrayAllocator<int>(
            sendList.size(), WO_E, ctx, runtime
        );
        neighbors = aaNeighbors->data();
        assert(neighbors);
    }
    //
    ArrayAllocator<local_int_t> *aaReceiveLength = nullptr;
    local_int_t *receiveLength = nullptr;
    if (receiveList.size() != 0) {
        aaReceiveLength = new ArrayAllocator<local_int_t>(
            receiveList.size(), WO_E, ctx, runtime
        );
        receiveLength = aaReceiveLength->data();
        assert(receiveLength);
    }
    //
    ArrayAllocator<local_int_t> *aaSendLength = nullptr;
    local_int_t *sendLength = nullptr;
    if (sendList.size() != 0) {
        aaSendLength = new ArrayAllocator<local_int_t>(
            sendList.size(), WO_E, ctx, runtime
        );
        sendLength = aaSendLength->data();
        assert(sendLength);
    }
    int neighborCount = 0;
    local_int_t receiveEntryCount = 0;
    local_int_t sendEntryCount = 0;
    for (map_iter curNeighbor = receiveList.begin();
         curNeighbor != receiveList.end(); ++curNeighbor, ++neighborCount) {
        // rank of current neighbor we are processing
        int neighborId = curNeighbor->first;
        // store rank ID of current neighbor
        neighbors[neighborCount] = neighborId;
        receiveLength[neighborCount] = receiveList[neighborId].size();
        // Get count if sends/receives
        sendLength[neighborCount] = sendList[neighborId].size();
        for (set_iter i = receiveList[neighborId].begin();
             i != receiveList[neighborId].end(); ++i, ++receiveEntryCount) {
            // The remote columns are indexed at end of internals
            externalToLocalMap[*i] = localNumberOfRows + receiveEntryCount;
        }
        for (set_iter i = sendList[neighborId].begin();
             i != sendList[neighborId].end(); ++i, ++sendEntryCount) {
            // store local ids of entry to send
            elementsToSend[sendEntryCount] = A.globalToLocalMap[*i];
        }
    }
    // Convert matrix indices to local IDs
    for (local_int_t i = 0; i< localNumberOfRows; i++) {
        for (int j = 0; j < nonzerosInRow[i]; j++) {
            global_int_t curIndex = mtxIndG(i, j);
            int rankIdOfColumnEntry = ComputeRankOfMatrixRow(*(A.geom), curIndex);
            // My column index, so convert to local index
            if (A.geom->rank==rankIdOfColumnEntry) {
                mtxIndL(i, j) = A.globalToLocalMap[curIndex];
            }
            // If column index is not a row index, then it comes from another
            // processor
            else {
                mtxIndL(i, j) = externalToLocalMap[curIndex];
            }
        }
    }
    // Store contents in our matrix struct
    // Convenience pointer to A.localData
    auto *AlD = A.localData;
    //
    AlD->numberOfExternalValues = externalToLocalMap.size();
    AlD->localNumberOfColumns   = AlD->localNumberOfRows
                                + AlD->numberOfExternalValues;
    AlD->numberOfSendNeighbors  = sendList.size();
    AlD->numberOfRecvNeighbors  = receiveList.size();
    AlD->totalToBeSent          = totalToBeSent;
    //
    if (aaElementsToSend) {
        aaElementsToSend->bindToLogicalRegion(*(A.pic.elementsToSend.data()));
        delete aaElementsToSend;
    }
    if (aaNeighbors) {
        aaNeighbors->bindToLogicalRegion(*(A.pic.neighbors.data()));
        delete aaNeighbors;
    }
    if (aaReceiveLength) {
         aaReceiveLength->bindToLogicalRegion(*(A.pic.receiveLength.data()));
         delete aaReceiveLength;
    }
    if (aaSendLength) {
        aaSendLength->bindToLogicalRegion(*(A.pic.sendLength.data()));
        delete aaSendLength;
    }

#if 0
    cout << " For rank " << A.geom->rank
            << " of " << A.geom->size
            << ", number of neighbors = " << AlD->numberOfSendNeighbors << endl;
    for (int i = 0; i < AlD->numberOfSendNeighbors; i++) {
        cout << "     rank " << A.geom->rank
             << " neighbor " << neighbors[i]
             << " send/recv length = " << sendLength[i]
             << "/" << receiveLength[i] << endl;
        for (local_int_t j = 0; j<sendLength[i]; ++j)
            cout << "       rank " << A.geom->rank
                 << " elementsToSend[" << j << "] = "
                 << elementsToSend[j] << endl;
    }
#endif
}

/**
 *
 */
inline void
populateSynchronizers(
    LogicalSparseMatrix &A,
    const Geometry &geom,
    SparseMatrixScalars *smScalars,
    LegionRuntime::HighLevel::Context ctx,
    LegionRuntime::HighLevel::Runtime *lrt
) {
    using namespace std;
    //
    cout << "*** Populating Synchronizers ..." << endl;
    const double startTime = mytimer();
    const int nShards = geom.size;
    //
    Array<LogicalRegion> alrSynchronizers(
        A.lrSynchronizers.mapRegion(RW_E, ctx, lrt), ctx, lrt
    );
    LogicalRegion *lrSynchronizers = alrSynchronizers.data();
    assert(lrSynchronizers);
    //
    for (int s = 0; s < nShards; ++s) {
        Synchronizers syncs = {
            .myPhaseBarriers = A.ownerPhaseBarriers[s],
            .neighborPhaseBarriers = A.neighborPhaseBarriers[s]
        };
        std::stringstream ss;
        {   // Scoped to guarantee flushing, etc.
            cereal::BinaryOutputArchive oa(ss);
            oa(syncs);
        }
        // Get size of serialized buffer.
        ss.seekp(0, ios::end);
        auto regionSizeInB = ss.tellp();
        // Allocate region-based backing store for serialized data.
        ArrayAllocator<char> lrStore(
            regionSizeInB, WO_E, ctx, lrt
        );
        char *lrStoreP = lrStore.data(); assert(lrStoreP);
        // Write back into logical region.
        memmove(lrStoreP, ss.str().c_str(), regionSizeInB);
        // Bind to shard's logical region.
        lrStore.bindToLogicalRegion(lrSynchronizers[s]);
        // Update shard's metadata associated with this buffer.
        smScalars[s].sizeofSynchronizersBuffer = regionSizeInB;
        // TODO Done, so unmap?
        lrStore.unmapRegion();
    }
    // TODO Done, so unmap?
    A.lrSynchronizers.unmapRegion(ctx, lrt);
    //
    const double endTime = mytimer();
    double time = endTime - startTime;
    cout << "--> Time=" << time << "s" << endl;
}

/**
 *
 */
inline void
SetupHaloTopLevel(
    LogicalSparseMatrix &A,
    const Geometry &geom,
    LegionRuntime::HighLevel::Context ctx,
    LegionRuntime::HighLevel::Runtime *lrt
) {
    using namespace std;
    //
    cout << "*** Setting Up Structures for SPMD Exchanges..." << endl;
    const double startTime = mytimer();
    const int nShards = geom.size;
    // Extract required info from logical structures.
    //
    cout << "--> Memory for SparseMatrixScalars="
         << (sizeof(SparseMatrixScalars) * nShards) / 1024.0
         << "kB" << endl;
    Array<SparseMatrixScalars> aSparseMatrixScalars(
        A.localData.mapRegion(RW_E, ctx, lrt), ctx, lrt
    );
    SparseMatrixScalars *smScalars = aSparseMatrixScalars.data();
    assert(smScalars);
    // Calculate the entire extent (including ghost cells) for a given halo'd
    // vector. This is simply the sum of all localNumberOfColumns.
    for (int s = 0; s < nShards; ++s) {
        const auto locLen = smScalars[s].localNumberOfColumns;
        A.requiredVectorLen += locLen;
        A.targetVectorPartLens.push_back(locLen);
    }
    cout << "--> Halo'd Vector Total Length=" << A.requiredVectorLen << endl;
    //
    Array<LogicalRegion> alrNeighbors(
        A.lrNeighbors.mapRegion(RO_E, ctx, lrt), ctx, lrt
    );
    LogicalRegion *lrpNeighbors = alrNeighbors.data();
    assert(lrpNeighbors);
    // Determine total number of PhaseBarriers required for synchronization.
    // Each shard will create two PhaseBarriers that it owns. A ready
    // PhaseBarrier to notify its consumers and a done PhaseBarrier to receive
    // notifications from its consumers.
    // 2x for ready/done pairs
    const int totpb = 2 * nShards;
    cout << "--> Total Number of PhaseBarriers Created="
         << totpb << endl;
    //
    A.ownerPhaseBarriers.resize(nShards);
    A.neighborPhaseBarriers.resize(nShards);
    // Iterate over all shards
    for (int shard = 0; shard < nShards; ++shard) {
        // Get total number of neighbors this shard has
        LogicalItem<LogicalRegion> lrNeighbors(lrpNeighbors[shard], ctx, lrt);
        Array<int> aNeighbors(lrNeighbors.mapRegion(RO_E, ctx, lrt), ctx, lrt);
        int *neighbors = aNeighbors.data(); assert(neighbors);
        const int nNeighbors = smScalars[shard].numberOfSendNeighbors;
#if 0
        cout << "Rank " << shard << " Has "
             << nNeighbors << " Send Neighbors " << endl;
        for (int i = 0; i < nNeighbors;  ++i) {
            cout << neighbors[i] << " ";
        }
        cout << endl;
#endif
        // Create my PhaseBarriers that I will then share with my neighbors.
        PhaseBarriers pbs = {
            .ready = lrt->create_phase_barrier(ctx, 1),
            .done  = lrt->create_phase_barrier(ctx, nNeighbors)
        };
        cout << "<-- task " << shard << " " << pbs.done << endl;
        A.ownerPhaseBarriers[shard] = pbs;
        // Share my PhaseBarriers with my neighbors
        for (int n = 0; n < nNeighbors; ++n) {
            auto &neighborMap = A.neighborPhaseBarriers[neighbors[n]];
            neighborMap[shard].push_back(pbs);
        }
        // Done with this data, so unmap.
        lrNeighbors.unmapRegion(ctx, lrt);
    }
    // No longer needed, so unmap.
    A.lrNeighbors.unmapRegion(ctx, lrt);
    //
    const double endTime = mytimer();
    double time = endTime - startTime;
    cout << "--> Time=" << time << "s" << endl;
    //
    populateSynchronizers(A, geom, smScalars, ctx, lrt);
    //
    A.localData.unmapRegion(ctx, lrt);
}
