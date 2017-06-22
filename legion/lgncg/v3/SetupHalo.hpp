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
#include <cassert>

/*!
  Reference version of SetupHalo that prepares system matrix data structure and
  creates data necessary for communication of boundary values of this process.

  @param[inout] A    The known system matrix

  @see ExchangeHalo
*/
void
SetupHalo(
    SparseMatrix & A
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
        for (int j=0; j<nonzerosInRow[i]; j++) {
            global_int_t curIndex = mtxIndG(i, j);
            int rankIdOfColumnEntry = ComputeRankOfMatrixRow(*(Ageom), curIndex);
            // If column index is not a row index, then it comes from another processor
            if (Ageom->rank!=rankIdOfColumnEntry) {
                receiveList[rankIdOfColumnEntry].insert(curIndex);
                // Matrix symmetry means we know the neighbor process wants my value
                sendList[rankIdOfColumnEntry].insert(currentGlobalRow);
            }
        }
    }
    // Count number of matrix entries to send and receive
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
#ifdef HPCG_DETAILED_DEBUG
    // These are all attributes that should be true, due to symmetry
    HPCG_fout << "totalToBeSent = " << totalToBeSent
              << " totalToBeReceived = " << totalToBeReceived << endl;
    // Number of sent entry should equal number of received
    assert(totalToBeSent == totalToBeReceived);
    // Number of send-to neighbors should equal number of receive-from
    assert(sendList.size()==receiveList.size());
    // Each receive-from neighbor should be a send-to neighbor, and send the
    // same number of entries
    for (map_iter curNeighbor = receiveList.begin();
         curNeighbor != receiveList.end(); ++curNeighbor)
    {
        assert(sendList.find(curNeighbor->first) != sendList.end());
        assert(sendList[curNeighbor->first].size() ==
               receiveList[curNeighbor->first].size());
    }
#endif
    // Build the arrays and lists needed by the ExchangeHalo function.
    double *sendBuffer = new double[totalToBeSent];
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
    for (local_int_t i=0; i< localNumberOfRows; i++) {
        for (int j=0; j<nonzerosInRow[i]; j++) {
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
    Asclrs->numberOfExternalValues = externalToLocalMap.size();
    //
#if 0 // FIXME
    Asclrs->localNumberOfColumns = Asclrs->localNumberOfRows
                                 + Asclrs->numberOfExternalValues;
#endif
    //
    Asclrs->numberOfSendNeighbors = sendList.size();
    Asclrs->totalToBeSent = totalToBeSent;
    //
    for (int i = 0; i < Asclrs->numberOfSendNeighbors; ++i) {
        A.neighbors->data()[i] = neighbors[i];
    }
#if 0
    A.elementsToSend = elementsToSend;
    A.receiveLength = receiveLength;
    A.sendLength = sendLength;
    A.sendBuffer = sendBuffer;
#endif

#ifdef HPCG_DETAILED_DEBUG
    HPCG_fout << " For rank " << A.geom->rank << " of " << A.geom->size
              << ", number of neighbors = " << A.numberOfSendNeighbors << endl;
    for (int i = 0; i < A.numberOfSendNeighbors; i++) {
        HPCG_fout << "     rank " << A.geom->rank
                  << " neighbor " << neighbors[i]
                  << " send/recv length = " << sendLength[i]
                  << "/" << receiveLength[i] << endl;
        for (local_int_t j = 0; j<sendLength[i]; ++j)
            HPCG_fout << "       rank " << A.geom->rank
                      << " elementsToSend[" << j << "] = " << elementsToSend[j]
                      << endl;
    }
#endif
}
