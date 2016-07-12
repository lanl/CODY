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
#include "LegionArrays.hpp"

#include "Geometry.hpp"
#include "MGData.hpp"

#include <vector>
#include <iomanip>
#include <map>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template<typename T>
struct LogicalSparseMatrix {
    char *title;
    //
    Geometry *mGeom;
    //total number of matrix rows across all processes
    global_int_t totalNumberOfRows;
    //total number of matrix nonzeros across all processes
    global_int_t totalNumberOfNonzeros;
    //number of rows local to this process
    local_int_t localNumberOfRows;
    //number of columns local to this process
    local_int_t localNumberOfColumns;
    //number of nonzeros local to this process
    local_int_t localNumberOfNonzeros;
    //The number of nonzeros in a row will always be 27 or fewer
    char *nonzerosInRow;
    //matrix indices as global values
    global_int_t **mtxIndG;
    //matrix indices as local values
    local_int_t **mtxIndL;
    //values of matrix entries
    double **matrixValues;
    //values of matrix diagonal entries
    double **matrixDiagonal;
    //global-to-local mapping
    std::map<global_int_t, local_int_t> globalToLocalMap;
    //local-to-global mapping
    std::vector<global_int_t> localToGlobalMap;
    //
    mutable bool isDotProductOptimized;
    //
    mutable bool isSpmvOptimized;
    //
    mutable bool isMgOptimized;
    //
    mutable bool isWaxpbyOptimized;
    /*!
      This is for storing optimized data structres created in OptimizeProblem and
      used inside optimized ComputeSPMV().
      */
    // Coarse grid matrix
    mutable struct SparseMatrix_STRUCT *Ac;
    // Pointer to the coarse level data for this fine matrix
    mutable MGData *mgData;
    // pointer that can be used to store implementation-specific data
    void *optimizationData;
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
    double *sendBuffer;
    /**
     *
     */
    LogicalSparseMatrix(void) = default;

    /**
     *
     */
    LogicalSparseMatrix(Geometry *geom)
    {
        title = 0;
        mGeom = geom;
        totalNumberOfRows = 0;
        totalNumberOfNonzeros = 0;
        localNumberOfRows = 0;
        localNumberOfColumns = 0;
        localNumberOfNonzeros = 0;
        nonzerosInRow = 0;
        mtxIndG = 0;
        mtxIndL = 0;
        matrixValues = 0;
        matrixDiagonal = 0;
        // Optimization is ON by default. The code that switches it OFF is in
        // the functions that are meant to be optimized.
        isDotProductOptimized = true;
        isSpmvOptimized = true;
        isMgOptimized = true;
        isWaxpbyOptimized = true;
        //
        numberOfExternalValues = 0;
        numberOfSendNeighbors = 0;
        totalToBeSent = 0;
        elementsToSend = 0;
        neighbors = 0;
        receiveLength = 0;
        sendLength = 0;
        sendBuffer = 0;
        //
        // Fine-to-coarse grid transfer initially not defined.
        mgData = 0;
        Ac = 0;
    }
};
