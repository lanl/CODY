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

#include "hpcg.hpp"
#include "mytimer.hpp"
#include "Geometry.hpp"
#include "VectorOps.hpp"
#include "CheckAspectRatio.hpp"
#include "GenerateGeometry.hpp"
#include "GenerateProblem.hpp"
#include "GenerateCoarseProblem.hpp"
#include "SetupHalo.hpp"
#include "CG.hpp"
#include "TestNorms.hpp"
#include "CheckProblem.hpp"
#include "ComputeResidual.hpp"

#include <iostream>
#include <cstdlib>

using namespace std;

/**
 * Generate Problem Task ///////////////////////////////////////////////////////
 */
void
genProblemTask(
    const Task *task,
    const vector<PhysicalRegion> &regions,
    Context ctx,
    HighLevelRuntime *runtime
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
    cout << "--> Time=" << initTime << " s" << endl;
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
    cout << "--> Time=" << totalTime << " s" << endl;
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
    cout << endl;
    cout << "*****************************************************" << endl;
    cout << "*** Run Statistics..." << endl;
    cout << "*****************************************************" << endl;
    cout << endl;
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    // At this point we need to know some run parameters so we can allocate and
    // partition the logical data structures. We'll use this info to calculate
    // global values, etc. NOTE: not all members will contain valid data after
    // generateInitGeometry returns (e.g., rank, ipx, ipy, ipz).
    cout << "*** Number of Shards=" << nShards << endl;
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
    cout << "--> nmg="  << NUM_MG_LEVELS << endl;
    ////////////////////////////////////////////////////////////////////////////
    cout << "*** Starting Initialization..." << endl;;
    // Application structures.
    LogicalSparseMatrix A;
    LogicalArray<floatType> b, x, xexact;
    //
    createLogicalStructures(
        A, b, x, xexact, initGeom, ctx, runtime
    );
    // Time to initialize problem before start of benchmark (phase 1).
    const double initStart = mytimer();
    {
        cout << "*** Launching Initialization Tasks..." << endl;;
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
        FutureMap fm = runtime->execute_must_epoch(ctx, mel);
        fm.wait_all_results();
        //
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
    // Capture phase 1 initialization time to pass to benchmark tasks.
    params.phase1InitTime = mytimer() - initStart;
    cout << "--> Time=" << params.phase1InitTime << " s" << endl;

    ////////////////////////////////////////////////////////////////////////////
    // Launch the tasks to begin the benchmark.
    ////////////////////////////////////////////////////////////////////////////
    {
        cout << endl;
        cout << "*****************************************************" << endl;
        cout << "*** Starting Benchmark..." << endl;
        cout << "*****************************************************" << endl;
        //
        const double start = mytimer();
        //
        MustEpochLauncher mel;
        //
        for (int shard = 0; shard < initGeom.size; ++shard) {
            TaskLauncher launcher(
                START_BENCHMARK_TID,
                TaskArgument(&params, sizeof(params))
            );
            const ItemFlags aif = IFLAG_W_GHOSTS;
            // Add all matrix levels.
            LogicalSparseMatrix *curLevelMatrix = &A;
            for (int level = 0; level < NUM_MG_LEVELS; ++level) {
                curLevelMatrix->intent(RW_E, aif, shard, launcher, ctx, runtime);
                curLevelMatrix = curLevelMatrix->Ac;
            }
            b.intent(     RW_E, shard, launcher, ctx, runtime);
            x.intent(     RW_E, shard, launcher, ctx, runtime);
            xexact.intent(RW_E, shard, launcher, ctx, runtime);
            //
            mel.add_single_task(DomainPoint::from_point<1>(shard), launcher);
        }
        //
        FutureMap fm = runtime->execute_must_epoch(ctx, mel);
        fm.wait_all_results();
        //
        const double totalTime = mytimer() - start;
        //
        cout << "*****************************************************" << endl;
        cout << "*** Benchmark Complete..." << endl;
        cout << "*****************************************************" << endl;
        //
        cout << "--> Time=" << totalTime << " s" << endl;
    }
    //
    cout << "*** Cleaning Up..." << endl;
    //
    destroyLogicalStructures(
        A, b, x, xexact, ctx, runtime
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
    const string levels = to_string(level);
    const string matrixName = level == 0 ? "A" : "A-L" + levels;
    // TODO deallocate.
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
        // These were mapped inline in startBenchmarkTask, so explicitly unmap.
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
startBenchmarkTask(
    const Task *task,
    const vector<PhysicalRegion> &regions,
    Context ctx,
    HighLevelRuntime *lrt
) {
    //
    double setup_time = mytimer();
    // Number of levels including first.
    const int numberOfMgLevels = NUM_MG_LEVELS;
    // Use this array for collecting timing information.
    std::vector<double> times(10, 0.0);
    //
    const HPCG_Params params = *(HPCG_Params *)task->args;
    // Check if QuickPath option is enabled.  If the running time is set to
    // zero, we minimize all paths through the program.
    const bool quickPath = (params.runningTime == 0);
    //
    size_t rid = 0;
    const ItemFlags aif = IFLAG_W_GHOSTS;
    //
    SparseMatrix A(regions, rid, aif, ctx, lrt);
    rid += A.nRegionEntries();
    //
    SparseMatrix *curLevelMatrix = &A;
    for (int level = 1; level < numberOfMgLevels; ++level) {
        curLevelMatrix->Ac = new SparseMatrix(regions, rid, aif, ctx, lrt);
        rid += curLevelMatrix->Ac->nRegionEntries();
        curLevelMatrix = curLevelMatrix->Ac;
    }
    //
    Array<floatType> b     (regions[rid++], ctx, lrt);
    Array<floatType> x     (regions[rid++], ctx, lrt);
    Array<floatType> xexact(regions[rid++], ctx, lrt);
    //
#ifdef LGNCG_TASKING
    lrt->unmap_all_regions(ctx); // TODO Is this needed?
#endif
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
    for (int level = 1; level < numberOfMgLevels; ++level) {
        allocateMGData(*curLevelMatrix, level - 1, ctx, lrt);
        f2cOperatorPopulate(*curLevelMatrix, ctx, lrt);
        curLevelMatrix = curLevelMatrix->Ac;
    }
    // Setup halo information for all levels before we begin.
    curLevelMatrix = &A;
    for (int level = 0; level < numberOfMgLevels; ++level) {
        SetupHalo(*curLevelMatrix, ctx, lrt);
        curLevelMatrix = curLevelMatrix->Ac;
    }
    // Capture total time of setup.
    setup_time = mytimer() - setup_time;
    setup_time += params.phase1InitTime;
    // Save it for reporting.
    times[9] = setup_time;
    //
    const int rank = A.geom->data()->rank;
    //
    if (rank == 0) {
        bool taskingEnabled = false;
#ifdef LGNCG_TASKING
        taskingEnabled = true;
#endif
        cout << "--> Implementation=Legion" << endl;
        cout << "--> Tasking="
             << (taskingEnabled ? "enabled" : "disabled")
             << endl;
        cout << "--> Total problem setup time in main (s) = "
             << setup_time << endl;
    }

    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    // Start of benchmark.
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////

    // Sanity
    assert(getTaskID(task) == rank);
    // Used to check return codes on function calls.
    int ierr = 0;
    int numberOfCalls = 10;
    //QuickPath means we do on one call of each block of repetitive code.
    if (quickPath) numberOfCalls = 1;
    //
    const auto *const Asclrs = A.sclrs->data();

    ////////////////////////////////////////////////////////////////////////////
    // Problem Sanity Phase
    ////////////////////////////////////////////////////////////////////////////

    {
        SparseMatrix *curLevelMatrix = &A;
        Array<floatType> *curb = &b;
        Array<floatType> *curx = &x;
        Array<floatType> *curxexact = &xexact;
        //
        for (int level = 0; level< numberOfMgLevels; ++level) {
            CheckProblem(*curLevelMatrix, curb, curx, curxexact, ctx, lrt);
            // Make the nextcoarse grid the next level.
            curLevelMatrix = curLevelMatrix->Ac;
            // No vectors after the top level.
            curb = NULL;
            curx = NULL;
            curxexact = NULL;
        }
    }

    ////////////////////////////////////////////////////////////////////////////
    // Reference CG Timing Phase                                              //
    ////////////////////////////////////////////////////////////////////////////

    int niters          = 0;
    int totalNiters_ref = 0;
    floatType normr     = 0.0;
    floatType normr0    = 0.0;
    int refMaxIters     = 50;
    // Only need to run the residual reduction analysis once
    numberOfCalls = 1;
    // Compute the residual reduction for the natural ordering and reference
    // kernels.
    std::vector<double> ref_times(9, 0.0);
    // Set tolerance to zero to make all runs do maxIters iterations.
    double tolerance = 0.0;
    int err_count = 0;
    for (int i = 0; i < numberOfCalls; ++i) {
        ZeroVector(x, ctx, lrt);
        ierr = CG(A, data, b, x, refMaxIters, tolerance, niters,
                  normr, normr0, &ref_times[0], true, ctx, lrt
               );
        // Count the number of errors in CG.
        if (ierr) ++err_count;
        totalNiters_ref += niters;
    }
    if (rank == 0 && err_count) {
        cerr << err_count << " error(s) in call(s) to reference CG." << endl;
    }
    //
    double refTolerance = normr / normr0;

    ////////////////////////////////////////////////////////////////////////////
    // Optimized CG Setup Phase                                               //
    ////////////////////////////////////////////////////////////////////////////

    niters                 = 0;
    normr                  = 0.0;
    normr0                 = 0.0;
    err_count              = 0;
    int tolerance_failures = 0;
    //
    int optMaxIters = 10 * refMaxIters;
    int optNiters = refMaxIters;
    floatType opt_worst_time = 0.0;
    //
    std::vector<double> opt_times(9, 0.0);
    // Compute the residual reduction and residual count for the user ordering
    // and optimized kernels.
    for (int i = 0; i < numberOfCalls; ++i) {
        // Start x at all zeros.
        ZeroVector(x, ctx, lrt);
        floatType last_cummulative_time = opt_times[0];
        ierr = CG(A, data, b, x, optMaxIters, refTolerance, niters,
                  normr, normr0, &opt_times[0], true, ctx, lrt
               );
        // Count the number of errors in CG.
        if (ierr) ++err_count;
        // The number of failures to reduce residual.
        if (normr / normr0 > refTolerance) ++tolerance_failures;
        // Pick the largest number of iterations to guarantee convergence.
        if (niters > optNiters) optNiters = niters;
        //
        floatType current_time = opt_times[0] - last_cummulative_time;
        if (current_time > opt_worst_time) opt_worst_time = current_time;
    }
    //
    Future localOptWorstF = Future::from_value(lrt, opt_worst_time);
    opt_worst_time = allReduce(
        localOptWorstF, *A.dcAllRedMaxFT, ctx, lrt
    ).get<floatType>();
    //
    if (rank == 0 && err_count) {
        cerr << err_count << " error(s) in call(s) to optimized CG." << endl;
    }
    if (tolerance_failures) {
        if (rank == 0) {
            cerr << "Failed to reduce the residual "
                 << tolerance_failures << " times." << endl;
        }
    }

    ////////////////////////////////////////////////////////////////////////////
    // Optimized CG Timing Phase                                              //
    ////////////////////////////////////////////////////////////////////////////

    // Here we finally run the benchmark phase The variable total_runtime is the
    // target benchmark execution time in seconds.

    double total_runtime = params.runningTime;
    // Run at least once, account for rounding.
    int numberOfCgSets = int(total_runtime / opt_worst_time) + 1;
    //
    if (rank == 0) {
        cout << "Projected running time: " << total_runtime << " seconds" << endl;
        cout << "Number of CG sets: " << numberOfCgSets << endl;
    }

    /* This is the timed run for a specified amount of time. */
    optMaxIters = optNiters;
    // Force optMaxIters iterations
    double optTolerance = 0.0;
    TestNormsData testnormsData;
    testnormsData.samples = numberOfCgSets;
    testnormsData.values = new double[numberOfCgSets];

    const double optTimeStart = mytimer();
    for (int i = 0; i < numberOfCgSets; ++i) {
        // Zero out x.
        ZeroVector(x, ctx, lrt);
        ierr = CG(A, data, b, x, optMaxIters, optTolerance, niters,
                  normr, normr0, &times[0], true, ctx, lrt
               );
        if (ierr) {
            cerr << "Error in call to CG: " << ierr << ".\n" << endl;
        }
        if (rank==0) {
            cout << "Call [" << i << "] Scaled Residual ["
                 << normr / normr0 << "]" << endl;
        }
        // Record scaled residual from this run.
        testnormsData.values[i] = normr / normr0;
    }
    //
    if (rank == 0) {
        const double optTimeEnd = mytimer();
        const double runTime = optTimeEnd - optTimeStart;
        //
        const double aveRuntime = runTime / double(numberOfCgSets);
        cout << numberOfCgSets << " CG set complete in " << runTime
             << " s" << endl;
        //
        cout << endl << "--> Average Run Time for CG="
             << aveRuntime << " s" << endl << endl;
    }
    // Compute difference between known exact solution and computed solution All
    // processors are needed here.
    floatType residual = 0.0;
    ierr = ComputeResidual(
        Asclrs->localNumberOfRows,
        x, xexact, residual, *A.dcAllRedMaxFT, ctx, lrt
    );
    if (ierr) {
        cerr << "Error in call to compute_residual: " << ierr << ".\n" << endl;
    }
    if (rank == 0) {
        cout << "Difference between computed and exact = "
             << residual << ".\n" << endl;
    }
#if 0
    // Test Norm Results.
    ierr = TestNorms(testnormsData);

    ////////////////////////////////////////////////////////////////////////////
    // Report Results                                                         //
    ////////////////////////////////////////////////////////////////////////////
    // Report results to YAML file.
    ReportResults(
        A,
        numberOfMgLevels,
        numberOfCgSets,
        refMaxIters,
        optMaxIters,
        &times[0],
        testCGData,
        testSymmetryData,
        testnormsData,
        global_failure,
        quickPath,
        ctx,
        lrt
    );
#endif
    ////////////////////////////////////////////////////////////////////////////
    // Cleanup task-local strucutres allocated for solve.
    ////////////////////////////////////////////////////////////////////////////
    destroySolveLocalStructures(A, data, ctx, lrt);
    lCGData.deallocate(ctx, lrt);
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
