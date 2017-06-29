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
    SparseMatrix &A,
    Array<floatType> &x,
    LegionRuntime::HighLevel::Context ctx,
    LegionRuntime::HighLevel::Runtime *lrt
) {
    using namespace std;
    // Extract Matrix pieces
    const SparseMatrixScalars *Asclrs = A.sclrs->data();
    Synchronizers *syncs = A.synchronizers->data();
    PhaseBarriers &myPBs = syncs->mine;
    const int nNeighbors = Asclrs->numberOfSendNeighbors;
    //local_int_t * receiveLength = A.receiveLength;
    //local_int_t *sendLength = A.sendLength->data();
    const int *const neighbors = A.neighbors->data();
    //double *sendBuffer = A.sendBuffer;
    const local_int_t totalToBeSent = Asclrs->totalToBeSent;
    // Non-region memory populated during SetupHalo().
    local_int_t *elementsToSend = A.elementsToSend;
    assert(elementsToSend);

    floatType *const xv = x.data();
    assert(xv);

    floatType *pushBuffer = A.pushBuffer->data();
    assert(pushBuffer);

    // Fill up push buffer.
    for (local_int_t i = 0; i < totalToBeSent; i++) {
        pushBuffer[i] = xv[elementsToSend[i]];
    }
    // Local copy done.
    myPBs.ready = lrt->advance_phase_barrier(ctx, myPBs.ready);

    for (int n = 0; n < nNeighbors; ++n) {
        //
        const int nid = neighbors[n];
        // Source
        auto srcIt = A.ghostArrays.find(nid);
        assert(srcIt != A.ghostArrays.end());
        LogicalArray<floatType> &srcArray = srcIt->second.first;
        //
        RegionRequirement srcrr(
            srcArray.logicalRegion,
            READ_ONLY,
            EXCLUSIVE,
            srcIt->second.second
        );
        srcrr.add_field(srcArray.fid);
        int srcVol = lrt->get_index_space_domain(ctx, srcArray.logicalRegion.get_index_space()).get_volume();
        //
        auto xis = x.logicalRegion.get_index_space();
        auto xip = lrt->get_index_partition(ctx, xis, 0);
        auto xlp = lrt->get_logical_partition(ctx, x.logicalRegion, xip);
        LogicalRegion xSubReg = lrt->get_logical_subregion_by_color(
            ctx,
            xlp,
            DomainPoint::from_point<1>(n + 1) // First is private.
        );
        RegionRequirement dstrr(
            xSubReg,
            WRITE_DISCARD,
            EXCLUSIVE,
            x.logicalRegion
        );
        dstrr.add_field(0); // FIXME

        int dstVol = lrt->get_index_space_domain(ctx, xSubReg.get_index_space()).get_volume();

        assert(srcVol == dstVol);

        CopyLauncher icl;
        icl.add_copy_requirements(srcrr, dstrr);
        lrt->issue_copy_operation(ctx, icl);
    }
}
