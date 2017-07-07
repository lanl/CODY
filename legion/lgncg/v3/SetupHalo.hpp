/**
 * Copyright (c) 2016-2017 Los Alamos National Security, LLC
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

#include "mytimer.hpp"

#include "LegionStuff.hpp"
#include "LegionItems.hpp"
#include "LegionArrays.hpp"
#include "LegionMatrices.hpp"

#include <map>
#include <set>
#include <vector>
#include <cassert>

inline void
GetNeighborInfo(
    SparseMatrix &A
) {
    using namespace std;
    // Extract Matrix pieces
    SparseMatrixScalars *Asclrs = A.sclrs->data();
    Geometry *Ageom = A.geom->data();
    //
    const local_int_t numberOfNonzerosPerRow = Ageom->stencilSize;
    //
    local_int_t localNumberOfRows = Asclrs->localNumberOfRows;
    //
    char *nonzerosInRow = A.nonzerosInRow->data();
    // Interpreted as 2D array
    Array2D<global_int_t> mtxIndG(
        localNumberOfRows, numberOfNonzerosPerRow, A.mtxIndG->data()
    );
    //
    global_int_t *AlocalToGlobalMap = A.localToGlobalMap->data();

    std::map< int, std::set< global_int_t> > sendList, receiveList;
    typedef std::map< int, std::set< global_int_t> >::iterator map_iter;
    typedef std::set<global_int_t>::iterator set_iter;
    std::map< local_int_t, local_int_t > externalToLocalMap;

    for (local_int_t i = 0; i < localNumberOfRows; i++) {
        global_int_t currentGlobalRow = AlocalToGlobalMap[i];
        for (int j = 0; j < nonzerosInRow[i]; j++) {
            global_int_t curIndex = mtxIndG(i, j);
            int rankIdOfColumnEntry = ComputeRankOfMatrixRow(*(Ageom), curIndex);
            // If column index is not a row index, then it comes from another
            // processor
            if (Ageom->rank!=rankIdOfColumnEntry) {
                receiveList[rankIdOfColumnEntry].insert(curIndex);
                // Matrix symmetry means we know the neighbor process wants my
                // value
                sendList[rankIdOfColumnEntry].insert(currentGlobalRow);
            }
        }
    }
    // Count number of matrix entries to send and receive.
    local_int_t totalToBeSent = 0;
    for (map_iter curNeighbor = sendList.begin();
         curNeighbor != sendList.end(); ++curNeighbor)
    {
        totalToBeSent += (curNeighbor->second).size();
    }
    // Build the arrays and lists needed by the ExchangeHalo function.
    local_int_t *elementsToSend = new local_int_t[totalToBeSent];
    int *neighbors = new int[sendList.size()];
    local_int_t *receiveLength = new local_int_t[receiveList.size()];
    local_int_t *sendLength = new local_int_t[sendList.size()];
    int neighborCount = 0;
    local_int_t receiveEntryCount = 0;
    for (map_iter curNeighbor = receiveList.begin();
         curNeighbor != receiveList.end();
         ++curNeighbor, ++neighborCount)
    {
        // rank of current neighbor we are processing
        int neighborId = curNeighbor->first;
        // store rank ID of current neighbor
        neighbors[neighborCount] = neighborId;
        receiveLength[neighborCount] = receiveList[neighborId].size();
        // Get count of sends/receives
        sendLength[neighborCount] = sendList[neighborId].size();
        for (set_iter i = receiveList[neighborId].begin();
             i != receiveList[neighborId].end();
             ++i, ++receiveEntryCount)
        {
            // The remote columns are indexed at end of internals
            externalToLocalMap[*i] = localNumberOfRows + receiveEntryCount;
        }
    }
    // Store contents in our matrix struct.
    Asclrs->numberOfExternalValues = externalToLocalMap.size();
    //
    Asclrs->localNumberOfColumns = Asclrs->localNumberOfRows
                                 + Asclrs->numberOfExternalValues;
    //
    Asclrs->numberOfSendNeighbors = sendList.size();
    Asclrs->totalToBeSent = totalToBeSent;
    //
    for (int i = 0; i < Asclrs->numberOfSendNeighbors; ++i) {
        A.neighbors->data()[i]  = neighbors[i];
        A.sendLength->data()[i] = sendLength[i];
        A.recvLength->data()[i] = receiveLength[i];
    }
    //
    delete[] elementsToSend;
    delete[] neighbors;
    delete[] receiveLength;
    delete[] sendLength;
}

/*!
  Reference version of SetupHalo that prepares system matrix data structure and
  creates data necessary for communication of boundary values of this process.

  @param[inout] A    The known system matrix

  @see ExchangeHalo
*/
inline void
SetupHalo(
    SparseMatrix &A,
    LegionRuntime::HighLevel::Context ctx,
    LegionRuntime::HighLevel::Runtime *lrt
) {
    using namespace std;
    //
    PopulateGlobalToLocalMap(A, ctx, lrt);
    // Extract Matrix pieces
    SparseMatrixScalars *Asclrs = A.sclrs->data();
    Geometry *Ageom = A.geom->data();
    //
    const local_int_t numberOfNonzerosPerRow = Ageom->stencilSize;
    //
    local_int_t localNumberOfRows = Asclrs->localNumberOfRows;
    //
    char *nonzerosInRow = A.nonzerosInRow->data();
    // Interpreted as 2D array
    Array2D<global_int_t> mtxIndG(
        localNumberOfRows, numberOfNonzerosPerRow, A.mtxIndG->data()
    );
    // Interpreted as 2D array
    Array2D<local_int_t> mtxIndL(
        localNumberOfRows, numberOfNonzerosPerRow, A.mtxIndL->data()
    );
    //
    global_int_t *AlocalToGlobalMap = A.localToGlobalMap->data();

    std::map< int, std::set< global_int_t> > sendList, receiveList;
    typedef std::map< int, std::set< global_int_t> >::iterator map_iter;
    typedef std::set<global_int_t>::iterator set_iter;
    std::map< local_int_t, local_int_t > externalToLocalMap;

    for (local_int_t i = 0; i < localNumberOfRows; i++) {
        global_int_t currentGlobalRow = AlocalToGlobalMap[i];
        for (int j = 0; j < nonzerosInRow[i]; j++) {
            global_int_t curIndex = mtxIndG(i, j);
            int rankIdOfColumnEntry = ComputeRankOfMatrixRow(*(Ageom), curIndex);
            // If column index is not a row index, then it comes from another
            // processor
            if (Ageom->rank!=rankIdOfColumnEntry) {
                receiveList[rankIdOfColumnEntry].insert(curIndex);
                // Matrix symmetry means we know the neighbor process wants my
                // value
                sendList[rankIdOfColumnEntry].insert(currentGlobalRow);
            }
        }
    }
    // Count number of matrix entries to send and receive.
    local_int_t totalToBeSent = 0;
    for (map_iter curNeighbor = sendList.begin();
         curNeighbor != sendList.end(); ++curNeighbor)
    {
        totalToBeSent += (curNeighbor->second).size();
    }
    local_int_t totalToBeReceived = 0;
    for (map_iter curNeighbor = receiveList.begin();
         curNeighbor != receiveList.end(); ++curNeighbor)
    {
        totalToBeReceived += (curNeighbor->second).size();
    }
    // Build the arrays and lists needed by the ExchangeHalo function.
    local_int_t *elementsToSend = new local_int_t[totalToBeSent];
    int *neighbors = new int[sendList.size()];
    local_int_t *receiveLength = new local_int_t[receiveList.size()];
    local_int_t *sendLength = new local_int_t[sendList.size()];
    int neighborCount = 0;
    local_int_t receiveEntryCount = 0;
    local_int_t sendEntryCount = 0;
    for (map_iter curNeighbor = receiveList.begin();
         curNeighbor != receiveList.end();
         ++curNeighbor, ++neighborCount)
    {
        // rank of current neighbor we are processing
        int neighborId = curNeighbor->first;
        // store rank ID of current neighbor
        neighbors[neighborCount] = neighborId;
        receiveLength[neighborCount] = receiveList[neighborId].size();
        // Get count of sends/receives
        sendLength[neighborCount] = sendList[neighborId].size();
        for (set_iter i = receiveList[neighborId].begin();
             i != receiveList[neighborId].end();
             ++i, ++receiveEntryCount)
        {
            // The remote columns are indexed at end of internals
            externalToLocalMap[*i] = localNumberOfRows + receiveEntryCount;
        }
        for (set_iter i = sendList[neighborId].begin();
             i != sendList[neighborId].end();
             ++i, ++sendEntryCount)
        {
            // Store local ids of entry to send.
            elementsToSend[sendEntryCount] = A.globalToLocalMap[*i];
        }
    }
    //
    for (local_int_t i = 0; i< localNumberOfRows; i++) {
        for (int j = 0; j<nonzerosInRow[i]; j++) {
            global_int_t curIndex = mtxIndG(i, j);
            int rankIdOfColumnEntry = ComputeRankOfMatrixRow(*(Ageom), curIndex);
            // My column index, so convert to local index
            if (Ageom->rank == rankIdOfColumnEntry) {
                mtxIndL(i, j) = A.globalToLocalMap[curIndex];
            }
            // If column index is not a row index, then it comes from another processor
            else {
                mtxIndL(i, j) = externalToLocalMap[curIndex];
            }
        }
    }
    // Store contents in our matrix struct.
    A.elementsToSend = elementsToSend;
#if 0 // Debug
    {
        const int me = Ageom->rank;
        local_int_t *neighbors = A.neighbors->data();
        char fName[128];
        sprintf(fName, "tx-list-%d.txt", Ageom->rank);
        FILE *f = fopen(fName, "w+");
        assert(f);
        for (int n = 0; n < Asclrs->numberOfSendNeighbors; ++n) {
            const local_int_t sendl = A.sendLength->data()[n];
            fprintf(
                f, "Task %d Sending Following Items to Task %d\n",
                me, neighbors[n]
            );
            for (int i = 0; i < sendl; ++i) {
                fprintf(f, " elementsToSend[%d] = %d\n", i, elementsToSend[i]);
            }
            fprintf(f, "\n");
        }
        fclose(f);
    }
#endif

    {
        BaseExtent *ApullBEs = A.pullBEs->data();
        auto &nidToPullRegion = A.nidToPullRegion;

        const int nNeighbors = Asclrs->numberOfSendNeighbors;
        for (int n = 0; n < nNeighbors; ++n) {
            const int nid = neighbors[n];
            auto it = nidToPullRegion.find(nid);
            // Make sure we found it.
            assert(it != nidToPullRegion.end());
            // Grab Logical Region from Physical.
            LogicalRegion lr = it->second.get_logical_region();
            // Create Array structure from LogicalRegion.
            LogicalArray<floatType> la(lr, ctx, lrt);
            // Extract my piece of the pull buffer.
            la.partition(ApullBEs[n], ctx, lrt);
            LogicalRegion subReg = lrt->get_logical_subregion_by_color(
                ctx,
                la.logicalPartition,
                DomainPoint::from_point<1>(0) // Only one partition.
            );
            // Stash sub-region Array.
            auto *lsr = new LogicalArray<floatType>(subReg, ctx, lrt);
            lsr->setParentLogicalRegion(lr);
            A.nidToRemotePullBuffer[nid] = lsr;
        }
    }

    delete[] neighbors;
    delete[] receiveLength;
    delete[] sendLength;
    // delete[] elementsToSend; Don't delete. Stored in sparse matrix.
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
         << (sizeof(SparseMatrixScalars) * nShards) / 1024.0 / 1024.0
         << " MB" << endl;
    Array<SparseMatrixScalars> aSparseMatrixScalars(
        A.sclrs.mapRegion(RW_E, ctx, lrt), ctx, lrt
    );
    SparseMatrixScalars *smScalars = aSparseMatrixScalars.data();
    assert(smScalars);
    //
    const int maxNumNeighbors = geom.stencilSize - 1;
    cout << "--> Memory for Neighbors="
         << (sizeof(int) * nShards * maxNumNeighbors) / 1024.0 / 1024.0
         << " MB" << endl;
    //
    Array<int> aNeighbors(
        A.neighbors.mapRegion(RO_E, ctx, lrt), ctx, lrt
    );
    // For convenience we'll interpret this as a 2D array.
    Array2D<int> neighbors(geom.size, maxNumNeighbors, aNeighbors.data());
    //
    cout << "--> Memory for Synchronizers="
         << (sizeof(Synchronizers) * nShards) / 1024.0 / 1024.0
         << " MB" << endl;
    Array<Synchronizers> aSynchronizers(
            A.synchronizers.mapRegion(RW_E, ctx, lrt), ctx, lrt
    );
    Synchronizers *synchronizers = aSynchronizers.data();
    assert(synchronizers);
    //
    cout << "--> Memory for Send/Recv Lengths="
         << (sizeof(local_int_t) * nShards * maxNumNeighbors) / 1024.0 / 1024.0
         << " MB" << endl;
    Array<local_int_t> aSendLengths(
            A.sendLength.mapRegion(RO_E, ctx, lrt), ctx, lrt
    );
    Array2D<local_int_t> ASendLengths(
        geom.size, maxNumNeighbors, aSendLengths.data()
    );
    //
    cout << "--> Memory for Bases/Extents="
         << (sizeof(BaseExtent) * nShards * maxNumNeighbors) / 1024.0 / 1024.0
         << " MB" << endl;
    Array<BaseExtent> aPullBEs(
            A.pullBEs.mapRegion(RW_E, ctx, lrt), ctx, lrt
    );
    Array2D<BaseExtent> pullBEs(
        geom.size, maxNumNeighbors, aPullBEs.data()
    );
    // Iterate over all shards and populate table that maps task IDs to their
    // index into the neighbor list.
    vector< map<int, int> > tidToNIdx(nShards);
    for (int shard = 0; shard < nShards; ++shard) {
        // Get total number of neighbors this shard has
        const SparseMatrixScalars &myScalars = smScalars[shard];
        const int nNeighbors = myScalars.numberOfSendNeighbors;
        //
        for (int n = 0; n < nNeighbors; ++n) {
            int tid = neighbors(shard, n);
            tidToNIdx[shard][tid] = n;
        }
    }
    // Iterate over all the shards and populate the synchronization structures.
    for (int shard = 0; shard < nShards; ++shard) {
        // Get total number of neighbors this shard has
        const SparseMatrixScalars &myScalars = smScalars[shard];
        const int nNeighbors = myScalars.numberOfSendNeighbors;
#if 1 // Debug
        cout << "Rank " << shard << " Has "
             << nNeighbors << " Pull Neighbors: " << endl;
        for (int n = 0; n < nNeighbors; ++n) {
            cout << neighbors(shard, n) << " ";
        }
        cout << endl;
#endif // Debug
        // Create PhaseBarriers for shard.
        PhaseBarriers pbs = {
            // Means I am ready for neighboring tasks to PULL values.
            .ready = lrt->create_phase_barrier(ctx, 1),
            // Means All pulls are complete.
            .done  = lrt->create_phase_barrier(ctx, nNeighbors)
        };
        //
        Synchronizers &mySync = synchronizers[shard];
        // This one is mine.
        mySync.mine = pbs;
        // Share mine with my all neighbors.
        for (int n = 0; n < nNeighbors; ++n) {
            // Get nth neighbor ID.
            const int myn = neighbors(shard, n);
            // Share my Synchronizers with that neighbor.
            synchronizers[myn].neighbors[tidToNIdx[myn][shard]] = pbs;
        }
        // Share base/extent (partitioning data) info with my neighbors.
        int curBase = 0;
        for (int n = 0; n < nNeighbors; ++n) {
            // Get nth neighbor ID.
            const int myn = neighbors(shard, n);
            const local_int_t sendl = ASendLengths(shard, n);
            pullBEs(myn, tidToNIdx[myn][shard]) = BaseExtent(curBase, sendl);
            curBase += sendl;
        }
    }
#if 0 // Debug Synchronizers
    for (int shard = 0; shard < nShards; ++shard) {
        const SparseMatrixScalars &myScalars = smScalars[shard];
        const int nNeighbors = myScalars.numberOfSendNeighbors;
        const Synchronizers &mySync = synchronizers[shard];
        cout << "task=" << shard << " sync data: " << endl;
        cout << "--------> mine.ready=" << mySync.mine.ready << endl;
        cout << "--------> mine.done =" << mySync.mine.done << endl;
        for (int n = 0; n < nNeighbors; ++n) {
            cout << "--------> .ready=" << mySync.neighbors[n].ready << endl;
            cout << "--------> .done =" << mySync.neighbors[n].done<< endl;
        }
    }
#endif
#if 0 // Debug BaseExtent
    for (int shard = 0; shard < nShards; ++shard) {
        const SparseMatrixScalars &myScalars = smScalars[shard];
        const int nNeighbors = myScalars.numberOfSendNeighbors;
        cout << "task=" << shard << " base/extent data: " << endl;
        for (int n = 0; n < nNeighbors; ++n) {
            BaseExtent be = pullBEs(shard, n);
            cout << "nid=" << neighbors(shard, n) << endl;
            cout << "--------> .base=" << be.base << endl;
            cout << "--------> .xtnt=" << be.extent << endl;
        }
    }
#endif
    // Cleanup and reporting.
    A.sclrs.unmapRegion(ctx, lrt);
    A.neighbors.unmapRegion(ctx, lrt);
    A.synchronizers.unmapRegion(ctx, lrt);
    A.sendLength.unmapRegion(ctx, lrt);
    A.pullBEs.unmapRegion(ctx, lrt);
    //
    const double initEnd = mytimer();
    const double initTime = initEnd - startTime;
    cout << "--> Time=" << initTime << "s" << endl;
}
