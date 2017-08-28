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

#if 0
#define LGNCG_DO_TASKY_EXCHANGE
#endif

#ifdef LGNCG_DO_TASKY_EXCHANGE
/**
 *
 */
struct ExchangeHaloArgs {
    int nTxNeighbors;
    int nRxNeighbors;
};

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
    // Extract Matrix pieces
    const SparseMatrixScalars *const Asclrs = A.sclrs->data();
    const int nTxNeighbors = Asclrs->numberOfSendNeighbors;
    const int nRxNeighbors = Asclrs->numberOfRecvNeighbors;
    const int *const neighbors = A.neighbors->data();
    // Nothing to do.
    if (nTxNeighbors == 0) return;
    // Make sure that x's ghosts are already setup.
    if (!x.hasGhosts()) {
        assert(false && "x does not have ghost regions setup.");
    }
    ExchangeHaloArgs args {
        .nTxNeighbors = nTxNeighbors,
        .nRxNeighbors = nRxNeighbors
    };
    TaskLauncher tl(
        EXCHANGE_HALO_TID,
        TaskArgument(&args, sizeof(args))
    );
    LogicalRegion xPrivateLR;
    {
        auto xis = x.logicalRegion.get_index_space();
        auto xip = lrt->get_index_partition(ctx, xis, 0 /* color */);
        auto xlp = lrt->get_logical_partition(ctx, x.logicalRegion, xip);
        xPrivateLR = lrt->get_logical_subregion_by_color(
            ctx,
            xlp,
            DomainPoint::from_point<1>(0) // First is private.
        );
    }
    // x (private partition).
    RegionRequirement xrr(
        xPrivateLR, RO_E, xPrivateLR
    );
    tl.add_region_requirement(xrr).add_field(x.fid);
    // Matrix pieces.
    A.neighbors->intent(RO_E, tl, ctx, lrt);
    A.elementsToSend->intent(RO_E, tl, ctx, lrt);
    A.sendLength->intent(RO_E, tl, ctx, lrt);
    A.synchronizers->intent(RW_E, tl, ctx, lrt);
    // Pull Buffers.
    for (int n = 0; n < nTxNeighbors; ++n) {
        A.pullBuffers[n]->intent(RW_E, tl, ctx, lrt);
    }
    // Pull regions in neighbor order.
    for (int n = 0; n < nRxNeighbors; ++n) {
        const int nid = neighbors[n];
        auto srcIt = A.nidToPullRegion.find(nid);
        assert(srcIt != A.nidToPullRegion.end());
        auto srclr = srcIt->second.get_logical_region();
        //
        RegionRequirement srcrr(
            srclr, RO_E, srclr
        );
        static const int fid = 0;
        tl.add_region_requirement(srcrr).add_field(fid);
    }
    for (int n = 0; n < nTxNeighbors; ++n) {
        LogicalArray<floatType> *dstArray = x.ghosts[n];
        assert(dstArray->hasParentLogicalRegion());
        RegionRequirement dstrr(
            dstArray->logicalRegion,
            WO_E,
            dstArray->getParentLogicalRegion()
        );
        tl.add_region_requirement(dstrr).add_field(dstArray->fid);
    }
    //
    lrt->execute_task(ctx, tl);
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
    const auto *const args = (ExchangeHaloArgs *)task->args;
    const int nTxNeighbors = args->nTxNeighbors;
    const int nRxNeighbors = args->nRxNeighbors;
    int rid = 0;
    // x
    Array<floatType> x(regions[rid++], ctx, lrt);
    const floatType *const xv = x.data();
    assert(xv);
    // Matrix pieces.
    Array<int> Aneighbors(regions[rid++], ctx, lrt);
    const int *const neighbors = Aneighbors.data();
    assert(neighbors);
    //
    Array<local_int_t> AelementsToSend(regions[rid++], ctx, lrt);
    const local_int_t *const elementsToSend = AelementsToSend.data();
    assert(elementsToSend);
    //
    Array<local_int_t> AsendLength(regions[rid++], ctx, lrt);
    const local_int_t *const sendLengthsd = AsendLength.data();
    assert(sendLengthsd);
    //
    Array<Synchronizers> Asyncs(regions[rid++], ctx, lrt);
    Synchronizers *syncs = Asyncs.data();
    assert(syncs);
    PhaseBarriers &myPBs = syncs->mine;
    //
    std::vector<Array<floatType> *> ApullBuffers;
    for (int n = 0; n < nTxNeighbors; ++n) {
        auto *item = new Array<floatType>(regions[rid++], ctx, lrt);
        assert(item->data());
        ApullBuffers.push_back(item);
    }
    std::vector<LogicalRegion> srclrs;
    for (int n = 0; n < nRxNeighbors; ++n) {
        srclrs.push_back(regions[rid++].get_logical_region());
    }
    std::vector<LogicalRegion> dstlrs;
    for (int n = 0; n < nRxNeighbors; ++n) {
        dstlrs.push_back(regions[rid++].get_logical_region());
    }
    //
    myPBs.done.wait();
    myPBs.done = lrt->advance_phase_barrier(ctx, myPBs.done);
    // Fill up pull buffers (the buffers that neighboring task will pull from).
    for (int n = 0, txidx = 0; n < nTxNeighbors; ++n) {
        floatType *const pbd = ApullBuffers[n]->data();
        //
        for (int i = 0; i < sendLengthsd[n]; ++i) {
            pbd[i] = xv[elementsToSend[txidx++]];
        }
    }
    myPBs.ready.arrive(1);
    myPBs.ready = lrt->advance_phase_barrier(ctx, myPBs.ready);
    //
    for (int n = 0; n < nTxNeighbors; ++n) {
        static const int fid = 0;
        // Source region requirement.
        RegionRequirement srcrr(
            srclrs[n], RO_E, srclrs[n]
        );
        srcrr.add_field(fid);
        // Destination region requirement.
        //
        RegionRequirement dstrr(
            dstlrs[n], WO_E, dstlrs[n]
        );
        dstrr.add_field(fid);
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
    // Cleanup
    for (int n = 0; n < nTxNeighbors; ++n) {
        delete ApullBuffers[n];
    }
}
#endif

/**
 *
 */
inline void
registerExchangeHaloTasks(void)
{
#ifdef LGNCG_DO_TASKY_EXCHANGE
    HighLevelRuntime::register_legion_task<ExchangeHaloTask>(
        EXCHANGE_HALO_TID /* task id */,
        Processor::LOC_PROC /* proc kind  */,
        true /* single */,
        false /* index */,
        AUTO_GENERATE_ID,
        TaskConfigOptions(false /* leaf task */),
        "ExchangeHaloTask"
    );
#endif
}

#ifndef LGNCG_DO_TASKY_EXCHANGE
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
    const local_int_t *const elementsToSend = A.elementsToSend->data();
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
#endif
