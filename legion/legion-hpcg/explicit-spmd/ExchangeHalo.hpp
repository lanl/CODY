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

#include "Geometry.hpp"
#include "LegionArrays.hpp"
#include "VectorOps.hpp"
#include "LegionMatrices.hpp"
#include "RegionToRegionCopy.hpp"

#include <cstdlib>

/*
 * PhaseBarriers debug: -level barrier=2 -logfile barriers.log
 */

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
    Context ctx,
    Runtime *lrt
) {
    using namespace std;
    // Extract Matrix pieces
    const SparseMatrixScalars *const Asclrs = A.sclrs->data();
    const int nNeighbors = Asclrs->numberOfSendNeighbors;
    // Nothing to do.
    if (nNeighbors == 0) return;
    // Else we have neighbors and data to move around.
    const int *const neighbors = A.neighbors->data();
    // Non-region memory populated during SetupHalo().
    const local_int_t *const elementsToSend = A.elementsToSend;
    assert(elementsToSend);
    // Setup ghost regions if not already there.
    if (!x.hasGhosts()) {
        assert(false && "x does not have ghost regions setup.");
    }
    //
    const floatType *const xv = x.data();
    assert(xv);
    //
    Synchronizers *syncs = A.synchronizers->data();
    PhaseBarriers &myPBs = syncs->mine;
    //
    myPBs.done.wait();
    myPBs.done = lrt->advance_phase_barrier(ctx, myPBs.done);
    // Fill up pull buffers (the buffers that neighboring task will pull from).
    const local_int_t *const sendLengthsd = A.sendLength->data();
    assert(sendLengthsd);
    //
    for (int n = 0, txidx = 0; n < nNeighbors; ++n) {
        floatType *const pbd = A.pullBuffers[n]->data();
        assert(pbd);
        //
        for (int i = 0; i < sendLengthsd[n]; ++i) {
            pbd[i] = xv[elementsToSend[txidx++]];
        }
    }
    myPBs.ready.arrive(1);
    myPBs.ready = lrt->advance_phase_barrier(ctx, myPBs.ready);
    //
    for (int n = 0; n < nNeighbors; ++n) {
        //
        const int nid = neighbors[n];
        // Source
        auto srcIt = A.nidToPullRegion.find(nid);
        assert(srcIt != A.nidToPullRegion.end());
        auto srclr = srcIt->second.get_logical_region();
        //
        RegionRequirement srcrr(
            srclr, RO_E, srclr
        );
        // Only ever one field for all of our structures.
        static const int srcFid = 0;
        srcrr.add_field(srcFid);
        // Destination.
        LogicalArray<floatType> *dstArray = x.ghosts[n];
        assert(dstArray->hasParentLogicalRegion());
        //
        RegionRequirement dstrr(
            dstArray->logicalRegion,
            WO_E,
            dstArray->getParentLogicalRegion()
        );
        dstrr.add_field(dstArray->fid);
        //
        TaskLauncher tl(
            REGION_TO_REGION_COPY_TID,
            TaskArgument(NULL, 0)
        );
        tl.add_region_requirement(srcrr);
        tl.add_region_requirement(dstrr);
        //
        syncs->neighbors[n].ready = lrt->advance_phase_barrier(
            ctx, syncs->neighbors[n].ready
        );
        tl.add_wait_barrier(syncs->neighbors[n].ready);
        //
        tl.add_arrival_barrier(syncs->neighbors[n].done);
        syncs->neighbors[n].done = lrt->advance_phase_barrier(
            ctx, syncs->neighbors[n].done
        );
        //
        lrt->execute_task(ctx, tl);
    }
}

/**
 *
 */
void
ExchangeHaloTask(
    const Task *task,
    const std::vector<PhysicalRegion> &regions,
    Context ctx,
    Runtime *lrt
) {
}

/**
 *
 */
inline void
registerExchangeHaloTasks(void)
{
    HighLevelRuntime::register_legion_task<ExchangeHaloTask>(
        EXCHANGE_HALO_TID /* task id */,
        Processor::LOC_PROC /* proc kind  */,
        true /* single */,
        false /* index */,
        AUTO_GENERATE_ID,
        TaskConfigOptions(false /* leaf task */),
        "ExchangeHaloTask"
    );
}
