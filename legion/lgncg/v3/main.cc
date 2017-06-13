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

#include "hpcg.hpp"
#include "mytimer.hpp"
#include "Geometry.hpp"
#include "CheckAspectRatio.hpp"
#include "GenerateGeometry.hpp"
#if 0
#include "GenerateProblem.hpp"
#include "SetupHalo.hpp"
#endif

#include "LegionStuff.hpp"
#include "LegionArrays.hpp"
#include "LegionMatrices.hpp"

#include <iostream>
#include <cstdlib>
#include <iomanip>

#include <unistd.h>

using namespace std;

LegionRuntime::Logger::Category Logger("LGNCG");

using namespace LegionRuntime::HighLevel;

/**
 * Generate Problem Task ///////////////////////////////////////////////////////
 */
void
genProblemTask(
    const Task *task,
    const std::vector<PhysicalRegion> &regions,
    Context ctx, HighLevelRuntime *runtime
) {
    const int nShards = int(getNumProcs());
    const int taskID = getTaskID(task);
    //
    HPCG_Params params;
    //
    SPMDMeta spmdMeta = {.rank = taskID, .nRanks = nShards};
    HPCG_Init(params, spmdMeta);
    // Check if QuickPath option is enabled.  If the running time is set to
    // zero, we minimize all paths through the program
    bool quickPath = (params.runningTime == 0);
    // Number of processes, my process ID
    int size = params.comm_size, rank = params.comm_rank;
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
    // Construct the geometry and linear system
    size_t rid = 0;
    //
    SparseMatrix A(regions, rid, ctx, runtime);
    rid += A.nRegionEntries();
    Array<floatType> b     (regions[rid++], ctx, runtime);
    Array<floatType> x     (regions[rid++], ctx, runtime);
    Array<floatType> xexact(regions[rid++], ctx, runtime);
    //
    Geometry *geom = A.geom;
    GenerateGeometry(size, rank, params.numThreads,
                     nx, ny, nz, params.stencilSize, geom);
    //
    ierr = CheckAspectRatio(0.125, geom->npx, geom->npy, geom->npz,
                            "process grid", rank == 0);
    if (ierr) exit(ierr);
    // Use this array for collecting timing information
    std::vector<double> times(10, 0.0);
    //
    double setup_time = mytimer();
    //
#if 0
    GenerateProblem(A, &b, &x, &xexact, ctx, runtime);
#endif
}

/**
 *
 */
static void
generateInitGeometry(
    int nShards,
    Geometry &globalGeom
) {
    const int stencilSize = 27;
    // We only care about passing nShards for this bit. rank doesn't make sense
    // in this context.
    const SPMDMeta meta = {.rank = 0, .nRanks = nShards};
    HPCG_Params params;
    HPCG_Init(params, meta);
    // Number of processes, my process ID
    const int size = params.comm_size, rank = params.comm_rank;
    //
    local_int_t nx, ny,nz;
    nx = (local_int_t)params.nx;
    ny = (local_int_t)params.ny;
    nz = (local_int_t)params.nz;
    // Generate geometry so we can get things like npx, npy, npz, etc.
    GenerateGeometry(size, rank, params.numThreads,
                     nx, ny, nz, stencilSize, &globalGeom);
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
    Context ctx, HighLevelRuntime *runtime
) {
    cout << "*** Creating Logical Structures..." << endl;
    const double initStart = mytimer();
    // Application structures
    // First calculate global XYZ for the problem.
    global_int_t globalXYZ = getGlobalXYZ(geom);
    //
    A.allocate(geom, ctx, runtime);
    A.partition(geom.size, ctx, runtime);
    x.allocate(globalXYZ, ctx, runtime);
    x.partition(geom.size, ctx, runtime);
    y.allocate(globalXYZ, ctx, runtime);
    y.partition(geom.size, ctx, runtime);
    xexact.allocate(globalXYZ, ctx, runtime);
    xexact.partition(geom.size, ctx, runtime);
    const double initEnd = mytimer();
    const double initTime = initEnd - initStart;
    cout << "--> Time=" << initTime << "s" << endl;
}

/**
 *
 */
static void
destroyLogicalStructures(
    LogicalSparseMatrix &A,
    LogicalArray<floatType> &x,
    LogicalArray<floatType> &y,
    LogicalArray<floatType> &xexact,
    Context ctx, HighLevelRuntime *runtime
) {
    cout << "*** Destroying Logical Structures..." << endl;
    const double initStart = mytimer();
    A.deallocate(ctx, runtime);
    x.deallocate(ctx, runtime);
    y.deallocate(ctx, runtime);
    xexact.deallocate(ctx, runtime);
    const double initEnd = mytimer();
    const double initTime = initEnd - initStart;
    cout << "--> Time=" << initTime << "s" << endl;
}

/**
 * Main Task ///////////////////////////////////////////////////////////////////
 * First task that gets spawned. Responsible for setup, etc.
 */
void
mainTask(
    const Task *,
    const std::vector<PhysicalRegion> &,
    Context ctx, HighLevelRuntime *runtime
) {
    const size_t nShards = getNumProcs();
    // ask the mapper how many shards we can have
    // TODO FIXME
    assert(nShards > 1 && "Run with at least 2 shards (e.g., try -ll:cpu 2)");
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    // At this point we need to know some run parameters so we can allocate and
    // partition the logical data structures. We'll use this info to calculate
    // global values, etc. NOTE: not all members will contain valid data after
    // generateInitGeometry returns (e.g., rank, ipx, ipy, ipz).
    cout << "*** Number of Shards (~ NUMPE)=" << nShards << endl;;
    Geometry initGeom;
    generateInitGeometry(nShards, initGeom);
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
        //
        cout << "*** Launching Initialization Tasks..." << endl;;
        const double initStart = mytimer();
        IndexLauncher launcher(
            GEN_PROB_TID,
            A.launchDomain,
            TaskArgument(nullptr, 0),
            ArgumentMap()
        );
        A.intent(WO_E, launcher);
        b.intent(WO_E, launcher);
        x.intent(WO_E, launcher);
        xexact.intent(WO_E, launcher);
        //
        auto futureMap = runtime->execute_index_space(ctx, launcher);
        cout << "*** Waiting for Initialization Tasks" << endl;
        futureMap.wait_all_results();
        const double initEnd = mytimer();
        double initTime = initEnd - initStart;
        cout << "--> Time=" << initTime << "s" << endl;
    }
#if 0
    // Now that we have all the setup information stored in LogicalRegions,
    // perform the top-level setup required for inter-task synchronization using
    // PhaseBarriers.
    SetupHaloTopLevel(A, initGeom, ctx, runtime);

    // SKG TODO RM For testing
    LogicalArray<floatType> testV;
    testV.allocate(A.requiredVectorLen, ctx, runtime);
    testV.partition(A.targetVectorPartLens, ctx, runtime);

    ////////////////////////////////////////////////////////////////////////////
    // Launch the tasks to begin the solve.
    ////////////////////////////////////////////////////////////////////////////
    //
    IndexLauncher launcher(
        START_SOLVE,
        A.launchDomain,
        TaskArgument(nullptr, 0),
        ArgumentMap()
    );
    // First give access to all logical subregions. Will be first nShards
    // regions.
    for (int color = 0; color < nShards; ++color) {
        auto lsr = runtime->get_logical_subregion_by_color(
                       ctx,
                       testV.logicalPartition,
                       color
        );
        launcher.add_region_requirement(
            RegionRequirement(
                lsr,
                0,
                RW_S,
                testV.logicalRegion
            )
        ).add_field(testV.fid)
         .add_flags(NO_ACCESS_FLAG);
    }
    // TODO FIXME - add real access flags on a per-item basis.
    // The application structures will always be at the end of the logical
    // subregions, which will always be of length nShards.
    intent<RW_S>(
        launcher,
        {A.geometries, A.localData, A.lrNonzerosInRow,
         A.lrMtxIndG, A.lrMtxIndL, A.lrMatrixValues,
         A.lrMatrixDiagonal, A.lrLocalToGlobalMap,
         A.lrGlobalToLocalMap, A.lrElementsToSend,
         A.lrNeighbors, A.lrReceiveLength, A.lrSendLength,
         A.synchronizersData, b}
    );
    //
    MustEpochLauncher mel;
    mel.add_index_task(launcher);
    FutureMap fm = runtime->execute_must_epoch(ctx, mel);
    fm.wait_all_results();
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
#endif
}

/**
 *
 */
void
startSolveTask(
    const Task *task,
    const std::vector<PhysicalRegion> &regions,
    Context ctx, HighLevelRuntime *lrt
) {
#if 0
    const int nShards = getNumProcs();
    const int taskID = getTaskID(task);
    const int nSubRegionReqs = nShards;
    // First n elements will be sub-region requirements for x.
    SparseMatrix A(regions, nSubRegionReqs, ctx, lrt);
    //
    Array<floatType> x(regions[taskID], ctx, lrt);
    //
    const int yrid = nSubRegionReqs + A.nRegionEntries();
    Array<floatType> y(regions[yrid], ctx, lrt);

#if (DEBUG_PHASE_BARRIER_SERIALIZATION == 1) // Serialization debug
    {
        size_t dataSize = A.localData->sizeofSynchronizersBuffer;
        Synchronizers syncs;
        char *synchronizersDataP = A.pic.synchronizersData.data();
        assert(synchronizersDataP);
        Synchronizers::deserialize(synchronizersDataP, dataSize, syncs);
        cout << "--> task " << taskID << " "
             << syncs.myPhaseBarriers.done << endl;
        sleep(taskID + 1);
        if (taskID == 0) printf("--> ALL DATA=");
        for (size_t i = 0; i < dataSize; ++i) {
            printf("%02hhX", (unsigned char)synchronizersDataP[i]);
            fflush(stdout);
        }
        if (taskID == nShards - 1) printf("\n");
    }
#endif

    // Test access
    LogicalRegion *lrTest = A.pic.matrixValues.data();
    assert(lrTest);
    sleep(taskID + 1);
    cout << taskID << ": -----> " << (*lrTest) << endl;
    // Now "dereference" the LogicalRegion.
    LogicalItem<floatType> lifTest(
       *lrTest, ctx, lrt
    );
    PhysicalRegion prTest = lifTest.mapRegion(RO_E, ctx, lrt);
#endif
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