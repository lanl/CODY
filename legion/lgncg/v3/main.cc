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

#include "LegionStuff.hpp"
#include "LegionArrays.hpp"
#include "LegionMatrices.hpp"
#include "LegionCGData.hpp"
#include "LegionMGData.hpp"
#include "IOOps.hpp"

#include "hpcg.hpp"
#include "mytimer.hpp"
#include "Geometry.hpp"
#include "CheckAspectRatio.hpp"
#include "GenerateGeometry.hpp"
#include "GenerateProblem.hpp"
#include "GenerateCoarseProblem.hpp"
#include "SetupHalo.hpp"
#include "VectorOps.hpp"
#include "CG.hpp"

#include <iostream>
#include <cstdlib>
#include <iomanip>

using namespace std;

/**
 * Generate Problem Task ///////////////////////////////////////////////////////
 */
void
genProblemTask(
    const Task *task,
    const vector<PhysicalRegion> &regions,
    Context ctx, HighLevelRuntime *runtime
) {
    const int taskID = getTaskID(task);
    const HPCG_Params params = *(HPCG_Params *)task->args;
    // Number of processes, my process ID
    const int size = params.commSize, rank = taskID;
    //
    local_int_t nx, ny,nz;
    nx = (local_int_t)params.nx;
    ny = (local_int_t)params.ny;
    nz = (local_int_t)params.nz;
    // Used to check return codes on function calls
    int ierr = 0;
    ierr = CheckAspectRatio(0.125, nx, ny, nz, "local problem", rank == 0);
    if (ierr) exit(ierr);
    ////////////////////////////////////////////////////////////////////////////
    // Problem setup Phase /////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    size_t rid = 0;
    //
    SparseMatrix A(regions, rid, ctx, runtime);
    rid += A.nRegionEntries();
    //
    GenerateGeometry(
        size,
        rank,
        params.numThreads,
        nx, ny, nz,
        params.stencilSize,
        A.geom->data()
    );
    ierr = CheckAspectRatio(
        0.125,
        A.geom->data()->npx,
        A.geom->data()->npy,
        A.geom->data()->npz,
        "process grid",
        rank == 0
    );
    if (ierr) exit(ierr);
    //
    SparseMatrix *curLevelMatrix = &A;
    for (int level = 1; level < NUM_MG_LEVELS; ++level) {
        curLevelMatrix->Ac = new SparseMatrix(regions, rid, ctx, runtime);
        rid += curLevelMatrix->Ac->nRegionEntries();
        curLevelMatrix = curLevelMatrix->Ac;
    }
    //
    Array<floatType> b     (regions[rid++], ctx, runtime);
    Array<floatType> x     (regions[rid++], ctx, runtime);
    Array<floatType> xexact(regions[rid++], ctx, runtime);
    //
    const int levelZero = 0;
    GenerateProblem(A, &b, &x, &xexact, levelZero, ctx, runtime);
    GetNeighborInfo(A);
    //
    curLevelMatrix = &A;
    for (int level = 1; level < NUM_MG_LEVELS; ++level) {
        GenerateCoarseProblem(*curLevelMatrix, level, ctx, runtime);
        curLevelMatrix = curLevelMatrix->Ac;
    }
}

/**
 *
 */
static void
generateInitGeometry(
    int nShards,
    const HPCG_Params &params,
    Geometry &globalGeom
) {
    // Number of processes, my process ID
    const int size = params.commSize, rank = 0;
    //
    local_int_t nx, ny,nz;
    nx = (local_int_t)params.nx;
    ny = (local_int_t)params.ny;
    nz = (local_int_t)params.nz;
    // Generate geometry so we can get things like npx, npy, npz, etc.
    GenerateGeometry(
        size,
        rank,
        params.numThreads,
        nx, ny, nz,
        params.stencilSize,
        &globalGeom
    );
}

/**
 *
 */
static void
createLogicalStructures(
    LogicalSparseMatrix     &A,
    LogicalArray<floatType> &x,
    LogicalArray<floatType> &y,
    LogicalArray<floatType> &xexact,
    const Geometry          &geom,
    Context ctx,
    HighLevelRuntime *runtime
) {
    cout << "*** Creating Logical Structures..." << endl;
    const double initStart = mytimer();
    // Application structures
    // First calculate global XYZ for the problem.
    global_int_t globalXYZ = getGlobalXYZ(geom);
    //
    A.allocate("A", geom, ctx, runtime);
    A.partition(geom.size, ctx, runtime);
    A.geom = new Geometry(geom);
    x.allocate("x", globalXYZ, ctx, runtime);
    x.partition(geom.size, ctx, runtime);
    y.allocate("y", globalXYZ, ctx, runtime);
    y.partition(geom.size, ctx, runtime);
    xexact.allocate("xexact", globalXYZ, ctx, runtime);
    xexact.partition(geom.size, ctx, runtime);
    //
    cout << "*** Creating Logical MG Structures..." << endl;
    LogicalSparseMatrix *curLevelMatrix = &A;
    for (int level = 1; level < NUM_MG_LEVELS; ++level) {
        GenerateCoarseProblemTopLevel(*curLevelMatrix, level, ctx, runtime);
        curLevelMatrix = curLevelMatrix->Ac;
    }
    //
    const double initEnd = mytimer();
    const double initTime = initEnd - initStart;
    cout << "--> Time=" << initTime << "s" << endl;
}

/**
 *
 */
static void
destroyLogicalStructures(
    LogicalSparseMatrix     &A,
    LogicalArray<floatType> &x,
    LogicalArray<floatType> &y,
    LogicalArray<floatType> &xexact,
    Context ctx,
    HighLevelRuntime *lrt
) {
    cout << "*** Destroying Logical Structures..." << endl;
    const double start = mytimer();
    //
    LogicalSparseMatrix *curLevelMatrix = &A;
    for (int level = 0; level < NUM_MG_LEVELS; ++level) {
        curLevelMatrix->deallocate(ctx, lrt);
        curLevelMatrix = curLevelMatrix->Ac;
    }
    //
    x.deallocate(ctx, lrt);
    y.deallocate(ctx, lrt);
    xexact.deallocate(ctx, lrt);
    //
    const double end = mytimer();
    const double totalTime = end - start;
    cout << "--> Time=" << totalTime << "s" << endl;
}

/**
 * Main Task ///////////////////////////////////////////////////////////////////
 * First task that gets spawned. Responsible for setup, etc.
 */
void
mainTask(
    const Task *,
    const vector<PhysicalRegion> &,
    Context ctx, HighLevelRuntime *runtime
) {
    // Ask the mapper how many shards we can have.
    const size_t nShards = getNumProcs();
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    // At this point we need to know some run parameters so we can allocate and
    // partition the logical data structures. We'll use this info to calculate
    // global values, etc. NOTE: not all members will contain valid data after
    // generateInitGeometry returns (e.g., rank, ipx, ipy, ipz).
    cout << "*** Number of Shards (~ NUMPE)=" << nShards << endl;;
    //
    // We only care about passing nShards for this bit. rank doesn't make sense
    // in this context.
    const SPMDMeta meta = {.rank = 0, .nRanks = int(nShards)};
    HPCG_Params params;
    //
    HPCG_Init(params, meta);
    //
    Geometry initGeom;
    generateInitGeometry(nShards, params, initGeom);
    cout << "*** Problem Information:"   << endl;
    cout << "--> size=" << initGeom.size << endl;
    cout << "--> npx="  << initGeom.npx  << endl;
    cout << "--> npy="  << initGeom.npy  << endl;
    cout << "--> npz="  << initGeom.npz  << endl;
    cout << "--> nx="   << initGeom.nx   << endl;
    cout << "--> ny="   << initGeom.ny   << endl;
    cout << "--> nz="   << initGeom.nz   << endl;
    ////////////////////////////////////////////////////////////////////////////
    cout << "*** Starting Initialization..." << endl;;
    // Application structures.
    LogicalSparseMatrix A;
    LogicalArray<floatType> b, x, xexact;
    //
    createLogicalStructures(
        A,
        b,
        x,
        xexact,
        initGeom,
        ctx,
        runtime
    );
    {
        cout << "*** Launching Initialization Tasks..." << endl;;
        const double initStart = mytimer();
        //
        MustEpochLauncher mel;
        //
        for (int shard = 0; shard < initGeom.size; ++shard) {
            TaskLauncher launcher(
                GEN_PROB_TID,
                TaskArgument(&params, sizeof(params))
            );
            // Add all matrix levels.
            LogicalSparseMatrix *curLevelMatrix = &A;
            for (int level = 0; level < NUM_MG_LEVELS; ++level) {
                curLevelMatrix->intent(RW_E, shard, launcher, ctx, runtime);
                curLevelMatrix = curLevelMatrix->Ac;
            }
            //
            b.intent(     RW_E, shard, launcher, ctx, runtime);
            x.intent(     RW_E, shard, launcher, ctx, runtime);
            xexact.intent(RW_E, shard, launcher, ctx, runtime);
            //
            mel.add_single_task(DomainPoint::from_point<1>(shard), launcher);
        }
        //
        cout << "*** Waiting for Initialization Tasks" << endl;
        //
        FutureMap fm = runtime->execute_must_epoch(ctx, mel);
        fm.wait_all_results();
        //
        const double initEnd = mytimer();
        double initTime = initEnd - initStart;
        cout << "--> Time=" << initTime << "s" << endl;
    }
    // Now that we have all the setup information stored in LogicalRegions,
    // perform the top-level setup required for inter-task synchronization using
    // PhaseBarriers.
    {
        LogicalSparseMatrix *curLevelMatrix = &A;
        for (int level = 0; level < NUM_MG_LEVELS; ++level) {
            SetupHaloTopLevel(*curLevelMatrix, level, ctx, runtime);
            curLevelMatrix = curLevelMatrix->Ac;
        }
    }

    ////////////////////////////////////////////////////////////////////////////
    // Launch the tasks to begin the solve.
    ////////////////////////////////////////////////////////////////////////////
    {
        cout << "*** Starting Solve..." << endl;
        const double start = mytimer();
        //
        MustEpochLauncher mel;
        //
        for (int shard = 0; shard < initGeom.size; ++shard) {
            TaskLauncher launcher(
                START_SOLVE_TID,
                TaskArgument(&params, sizeof(params))
            );
            const ItemFlags aif = IFLAG_W_GHOSTS;
            // Add all matrix levels.
            LogicalSparseMatrix *curLevelMatrix = &A;
            for (int level = 0; level < NUM_MG_LEVELS; ++level) {
                curLevelMatrix->intent(RW_E, aif, shard, launcher, ctx, runtime);
                curLevelMatrix = curLevelMatrix->Ac;
            }
            b.intent(     RW_E,         shard, launcher, ctx, runtime);
            x.intent(     RW_E,         shard, launcher, ctx, runtime);
            xexact.intent(RW_E,         shard, launcher, ctx, runtime);
            //
            mel.add_single_task(DomainPoint::from_point<1>(shard), launcher);
        }
        //
        FutureMap fm = runtime->execute_must_epoch(ctx, mel);
        fm.wait_all_results();
        //
        const double end = mytimer();
        const double totalTime = end - start;
        cout << "--> Time=" << totalTime << "s" << endl;
    }
    //
    cout << "*** Cleaning Up..." << endl;
    destroyLogicalStructures(
        A,
        b,
        x,
        xexact,
        ctx,
        runtime
    );
}

/**
 *
 */
static void
allocateMGData(
    SparseMatrix &A,
    int level,
    Context ctx,
    HighLevelRuntime *lrt
) {
    string levels = to_string(level);
    string matrixName = level == 0 ? "A" : "A-L" + levels;
    //
    LogicalMGData lMGData;
    lMGData.allocate(matrixName, A, ctx, lrt);
    lMGData.partition(A, ctx, lrt);
    //
    std::vector<PhysicalRegion> mgRegions;
    mgRegions.push_back(lMGData.f2cOperator.mapRegion(RW_E, ctx, lrt));
    mgRegions.push_back(         lMGData.rc.mapRegion(RW_E, ctx, lrt));
    mgRegions.push_back(         lMGData.xc.mapRegion(RW_E, ctx, lrt));
    mgRegions.push_back(        lMGData.Axf.mapRegion(RW_E, ctx, lrt));
    //
    const int mgDataBaseRID = 0;
    A.mgData = new MGData(mgRegions, mgDataBaseRID, ctx, lrt);
}

/**
 *
 */
static void
destroySolveLocalStructures(
    SparseMatrix &A,
    CGData &cgData,
    Context ctx,
    HighLevelRuntime *lrt
) {
    SparseMatrix *curLevelMatrix = &A;
    for (int level = 1; level < NUM_MG_LEVELS; ++level) {
        // These were mapped inline in startSolveTask, so explicitly unmap.
        curLevelMatrix->mgData->unmapRegions(ctx, lrt);
        curLevelMatrix = curLevelMatrix->Ac;
    }
    //
    cgData.unmapRegions(ctx, lrt);
}

/**
 *
 */
void
startSolveTask(
    const Task *task,
    const vector<PhysicalRegion> &regions,
    Context ctx,
    HighLevelRuntime *lrt
) {
    const HPCG_Params params = *(HPCG_Params *)task->args;
    //
    size_t rid = 0;
    const ItemFlags aif = IFLAG_W_GHOSTS;
    //
    SparseMatrix A(regions, rid, aif, ctx, lrt);
    rid += A.nRegionEntries();
    //
    SparseMatrix *curLevelMatrix = &A;
    for (int level = 1; level < NUM_MG_LEVELS; ++level) {
        curLevelMatrix->Ac = new SparseMatrix(regions, rid, aif, ctx, lrt);
        rid += curLevelMatrix->Ac->nRegionEntries();
        curLevelMatrix = curLevelMatrix->Ac;
    }
    //
    Array<floatType> b     (regions[rid++], ctx, lrt);
    Array<floatType> x     (regions[rid++], ctx, lrt);
    Array<floatType> xexact(regions[rid++], ctx, lrt);
    ////////////////////////////////////////////////////////////////////////////
    // Private data for this task.
    ////////////////////////////////////////////////////////////////////////////
    // CGData
    LogicalCGData lCGData;
    lCGData.allocate("cgdata", A, ctx, lrt);
    lCGData.partition(A, ctx, lrt);
    // Map CG data locally.
    vector<PhysicalRegion> cgRegions;
    const int nCGDataRegions = 4;
    cgRegions.reserve(nCGDataRegions);
    //
    cgRegions.push_back( lCGData.r.mapRegion(RW_E, ctx, lrt));
    cgRegions.push_back( lCGData.z.mapRegion(RW_E, ctx, lrt));
    cgRegions.push_back( lCGData.p.mapRegion(RW_E, ctx, lrt));
    cgRegions.push_back(lCGData.Ap.mapRegion(RW_E, ctx, lrt));
    //
    const int cgDataBaseRID = 0;
    CGData data(cgRegions, cgDataBaseRID, ctx, lrt);
    // MGData
    curLevelMatrix = &A;
    for (int level = 1; level < NUM_MG_LEVELS; ++level) {
        allocateMGData(*curLevelMatrix, level - 1, ctx, lrt);
        f2cOperatorPopulate(*curLevelMatrix, ctx, lrt);
        curLevelMatrix = curLevelMatrix->Ac;
    }
    // Sanity
    assert(getTaskID(task) == A.geom->data()->rank);
    ////////////////////////////////////////////////////////////////////////////
    // Setup halo information for all levels before we begin.
    ////////////////////////////////////////////////////////////////////////////
    curLevelMatrix = &A;
    for (int level = 0; level < NUM_MG_LEVELS; ++level) {
        SetupHalo(*curLevelMatrix, ctx, lrt);
        curLevelMatrix = curLevelMatrix->Ac;
    }

    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    // Start of benchmark.
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////

    const int refMaxIters  = 50;
    const int optMaxIters  = 10 * refMaxIters;
    //
    // int numberOfCalls = 10; FIXME
    int numberOfCalls = 2;
    // Check if QuickPath option is enabled.  If the running time is set to
    // zero, we minimize all paths through the program.
    bool quickPath = (params.runningTime == 0);
    //QuickPath means we do on one call of each block of repetitive code
    if (quickPath) numberOfCalls = 1;
    //
    int niters          = 0;
    double normr        = 0.0;
    double normr0       = 0.0;
    int errCount        = 0;
    double optWorstTime = 0.0;
    // Set tolerance to zero to make all runs do maxIters iterations
    double tolerance    = 0.0;

    const bool doMG = true;
    vector<double> optTimes(9, 0.0);
    //
    for (int i = 0; i < numberOfCalls; ++i) {
        ZeroVector(x, ctx, lrt); // Start x at all zeros.
        double lastCummulativeTime = optTimes[0];
        int ierr = CG(A,
                      data,
                      b,
                      x,
                      optMaxIters,
                      tolerance,
                      niters,
                      normr,
                      normr0,
                      &optTimes[0],
                      doMG,
                      ctx,
                      lrt
                   );
        if (ierr) ++errCount; // count the number of errors in CG
#if 0
        // the number of failures to reduce residual
        if (normr / normr0 > refTolerance) ++toleranceFailures;
        // pick the largest number of iterations to guarantee convergence
        if (niters > optNiters) optNiters = niters;
#endif
        double current_time = optTimes[0] - lastCummulativeTime;
        if (current_time > optWorstTime) optWorstTime = current_time;
    }
    // Cleanup task-local strucutres allocated for solve.
    destroySolveLocalStructures(A, data, ctx, lrt);
}

/**
 * Main ////////////////////////////////////////////////////////////////////////
 * Responsible for RT initialization.
 */
int
main(int argc, char **argv)
{
    LegionInit();
    return Runtime::start(argc, argv);
}
