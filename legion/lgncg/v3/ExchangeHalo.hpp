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

#include "LegionMatrices.hpp"
#include "LegionArrays.hpp"
#include "Geometry.hpp"

#include <cstdlib>

/*!
    Communicates data that is at the border of the part of the domain assigned to
    this processor.

    @param[in]    A The known system matrix
    @param[inout] x On entry: the local vector entries followed by entries to be
    communicated; on exit: the vector with non-local entries updated by other
    processors
 */
inline void
ExchangeHalo(
      const SparseMatrix &A,
      Array<floatType> &x
) {
    using namespace std;
    // Extract Matrix pieces
    const SparseMatrixScalars *Asclrs = A.sclrs->data();
    const Geometry *Ageom = A.geom->data();
    const local_int_t localNumberOfRows = Asclrs->localNumberOfRows;
    const int num_neighbors = Asclrs->numberOfSendNeighbors;
    //local_int_t * receiveLength = A.receiveLength;
    local_int_t *sendLength = A.sendLength->data();
    int *neighbors = A.neighbors->data();
    //double *sendBuffer = A.sendBuffer;
    const local_int_t totalToBeSent = Asclrs->totalToBeSent;
    // Non-region memory populated during SetupHalo().
    local_int_t *elementsToSend = A.elementsToSend;
    assert(elementsToSend);

    double *const xv = x.data();
    assert(xv);

    // Number of Legion tasks.
    const int size = Ageom->size;
    // My task ID.
    const int rank = Ageom->rank;
#if 0
    //
    // Externals are at end of locals.
    //
    double * x_external = (double *) xv + localNumberOfRows;

    // Post receives first
    // TODO: Thread this loop
    for (int i = 0; i < num_neighbors; i++) {
      local_int_t n_recv = receiveLength[i];
      MPI_Irecv(x_external, n_recv, MPI_DOUBLE, neighbors[i], MPI_MY_TAG, MPI_COMM_WORLD, request+i);
      x_external += n_recv;
    }


    //
    // Fill up send buffer
    //

    // TODO: Thread this loop
    for (local_int_t i=0; i<totalToBeSent; i++) sendBuffer[i] = xv[elementsToSend[i]];

    //
    // Send to each neighbor
    //

    // TODO: Thread this loop
    for (int i = 0; i < num_neighbors; i++) {
      local_int_t n_send = sendLength[i];
      MPI_Send(sendBuffer, n_send, MPI_DOUBLE, neighbors[i], MPI_MY_TAG, MPI_COMM_WORLD);
      sendBuffer += n_send;
    }

    //
    // Complete the reads issued above
    //

    MPI_Status status;
    // TODO: Thread this loop
    for (int i = 0; i < num_neighbors; i++) {
      if ( MPI_Wait(request+i, &status) ) {
        std::exit(-1); // TODO: have better error exit
      }
    }

    delete [] request;

    return;
#endif
}
