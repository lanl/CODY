/**
 * Copyright (c) 2016      Los Alamos National Security, LLC
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

#include "LegionStuff.hpp"
#include "LegionArrays.hpp"
#include "LegionMatrices.hpp"
#include "CGMapper.h"

#include <iostream>
#include <cstdlib>
#include <deque>
#include <iomanip>

using namespace std;

LegionRuntime::Logger::Category Logger("LGNCG");

/**
 * SPMD Main Task //////////////////////////////////////////////////////////////
 */
void
spmdInitTask(
    const Task *task,
    const std::vector<PhysicalRegion> &regions,
    Context ctx, HighLevelRuntime *runtime
) {
    static const int ridParams = 0;
    static const int ridGeom   = 1;
    //
    const int nShards = *(int *)task->args;
    const int taskID = getTaskID(task);

    // XXX SKG LogicalRegion TEST
    PhysicalScalar<LogicalRegion> testPR(regions[2], ctx, runtime);
    LogicalRegion *TEST = testPR.data();
    assert(TEST);

    LogicalArray<double> foo;
    foo.allocate(1+taskID, ctx, runtime);
    auto pr = foo.mapRegion(ctx, runtime);

    (*TEST) = foo.logicalRegion;

    RegionAccessor<AccessorType::Generic, double> pracc =
    pr.get_field_accessor(foo.fid).typeify<double>();

    for (GenericPointInRectIterator<1> pir(foo.bounds); pir; pir++) {
        cout << "<--"<< (taskID + 1) << endl;
        pracc.write(DomainPoint::from_point<1>(pir.p), (double)(taskID + 1));
    }

    runtime->unmap_region(ctx, pr);

    cout << "SUCCESS!!!!" << endl;

    return;
    // END XXX SKG LogicalRegion TEST

    //
    PhysicalScalar<HPCG_Params> psHPCGParams(regions[ridParams], ctx, runtime);
    HPCG_Params *params = psHPCGParams.data();
    assert(params);
    //
    SPMDMeta spmdMeta = {.rank = taskID, .nRanks = nShards};
    HPCG_Init(*params, spmdMeta);
    // Check if QuickPath option is enabled.  If the running time is set to
    // zero, we minimize all paths through the program
    bool quickPath = (params->runningTime == 0);
    // Number of processes, my process ID
    int size = params->comm_size, rank = params->comm_rank;
    //
    local_int_t nx, ny,nz;
    nx = (local_int_t)params->nx;
    ny = (local_int_t)params->ny;
    nz = (local_int_t)params->nz;
    // Used to check return codes on function calls
    int ierr = 0;
    ierr = CheckAspectRatio(0.125, nx, ny, nz, "local problem", rank == 0);
    if (ierr) exit(ierr);
    ////////////////////////////////////////////////////////////////////////////
    // Problem setup Phase /////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    // Construct the geometry and linear system
    PhysicalScalar<Geometry> psGeometry(regions[ridGeom], ctx, runtime);
    Geometry *geom = psGeometry.data();
    assert(geom);
    GenerateGeometry(size, rank, params->numThreads, nx, ny, nz, geom);
    //
    ierr = CheckAspectRatio(0.125, geom->npx, geom->npy, geom->npz,
                            "process grid", rank == 0);
    if (ierr) exit(ierr);
    // Use this array for collecting timing information
    std::vector<double> times(10, 0.0);
    //
    double setup_time = mytimer();
}

/**
 *
 */
static void
generateInitGeometry(
    int nShards,
    Geometry &globalGeom
) {
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
    GenerateGeometry(size, rank, params.numThreads, nx, ny, nz, &globalGeom);
}

/**
 *
 */
static global_int_t
getGlobalXYZ(
    const Geometry &geom
) {
    return geom.npx * geom.nx *
           geom.npy * geom.ny *
           geom.npz * geom.nz;
}

/**
 *
 */
static void
createLogicalStructures(
    LogicalArray<HPCG_Params>   &hpcgParams,
    LogicalArray<Geometry>      &geometries,
    LogicalSparseMatrix<double> &A,
    LogicalArray<double>        &x,
    LogicalArray<double>        &y,
    const Geometry              &geom,
    Context ctx, HighLevelRuntime *runtime
) {
    cout << "*** Creating Logical Structures..." << endl;;
    const double initStart = mytimer();
    // Metadata - one for each SPMD task.
    hpcgParams.allocate(geom.size, ctx, runtime);
    hpcgParams.partition(geom.size, ctx, runtime);
    geometries.allocate(geom.size, ctx, runtime);
    geometries.partition(geom.size, ctx, runtime);
    // Application structures
    // First calculate global XYZ for the problem.
    global_int_t globalXYZ = getGlobalXYZ(geom);
    //A TODO
    x.allocate(globalXYZ, ctx, runtime);
    y.allocate(globalXYZ, ctx, runtime);
    const double initEnd = mytimer();
    const double initTime = initEnd - initStart;
    cout << "    Done in: " << initTime << " s" << endl;;
}

/**
 *
 */
static void
destroyLogicalStructures(
    LogicalArray<HPCG_Params>   &hpcgParams,
    LogicalArray<Geometry>      &geometries,
    LogicalSparseMatrix<double> &A,
    LogicalArray<double>        &x,
    LogicalArray<double>        &y,
    Context ctx, HighLevelRuntime *runtime
) {
    cout << "*** Destroying Logical Structures..." << endl;;
    const double initStart = mytimer();
    hpcgParams.deallocate(ctx, runtime);
    geometries.deallocate(ctx, runtime);
    // A TODO
    x.deallocate(ctx, runtime);
    y.deallocate(ctx, runtime);
    const double initEnd = mytimer();
    const double initTime = initEnd - initStart;
    cout << "    Done in: " << initTime << " s" << endl;;
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
    // ask the mapper how many shards we can have
    const int nShards = runtime->get_tunable_value(
                            ctx, CGMapper::TID_NUM_SHARDS
                        );
    cout << "*** Number of Shards (~ NUMPE): " << nShards << endl;;
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    // At this point we need to know some run parameters so we can allocate and
    // partition the logical data structures. We'll use this info to calculate
    // global values, etc. NOTE: not all members will contain valid data after
    // generateInitGeometry returns (e.g., rank, ipx, ipy, ipz).
    Geometry initGeom;
    generateInitGeometry(nShards, initGeom);
    cout << "*** Problem Information:"   << endl;
    cout << "    size=" << initGeom.size << endl;
    cout << "    npx="  << initGeom.npx  << endl;
    cout << "    npy="  << initGeom.npy  << endl;
    cout << "    npz="  << initGeom.npz  << endl;
    cout << "    nx="   << initGeom.nx   << endl;
    cout << "    ny="   << initGeom.ny   << endl;
    cout << "    nz="   << initGeom.nz   << endl;
    ////////////////////////////////////////////////////////////////////////////
    cout << "*** Starting Initialization..." << endl;;
    // Legion array holding HPCG parameters that impact run.
    LogicalArray<HPCG_Params> hpcgParams;
    // Legion array holding problem geometries.
    LogicalArray<Geometry> geometries;
    // Application structures.
    LogicalSparseMatrix<double> A;
    LogicalArray<double> x, y;
    //
    createLogicalStructures(
        hpcgParams,
        geometries,
        A,
        x,
        y,
        initGeom,
        ctx,
        runtime
    );
    // XXX SKG LogicalRegion TEST
    LogicalArray<LogicalRegion> testLR;
    testLR.allocate(initGeom.size, ctx, runtime);
    testLR.partition(initGeom.size, ctx, runtime);
    // END XXX SKG LogicalRegion TEST
    cout << "*** Launching Initialization Tasks..." << endl;;
    const double initStart = mytimer();
    IndexLauncher launcher(
        SPMD_INIT_TID,
        hpcgParams.launchDomain(),
        TaskArgument(&nShards, sizeof(nShards)),
        ArgumentMap()
    );
    //
    launcher.add_region_requirement(
        RegionRequirement(
            hpcgParams.logicalPartition(),
            0,
            WRITE_DISCARD,
            EXCLUSIVE,
            hpcgParams.logicalRegion
        )
    ).add_field(hpcgParams.fid);
    //
    launcher.add_region_requirement(
        RegionRequirement(
            geometries.logicalPartition(),
            0,
            WRITE_DISCARD,
            EXCLUSIVE,
            geometries.logicalRegion
        )
    ).add_field(geometries.fid);
    // XXX SKG LogicalRegion TEST
    launcher.add_region_requirement(
        RegionRequirement(
            testLR.logicalPartition(),
            0,
            WRITE_DISCARD,
            EXCLUSIVE,
            testLR.logicalRegion
        )
    ).add_field(testLR.fid);
    // END XXX SKG LogicalRegion TEST
    //
    auto futureMap = runtime->execute_index_space(ctx, launcher);
    cout << "*** Waiting for Initialization Tasks" << endl;
    futureMap.wait_all_results();
    const double initEnd = mytimer();
    const double initTime = initEnd - initStart;
    cout << "*** Initialization Time (s): " << initTime << endl;

    // XXX SKG LogicalRegion TEST
    auto prs = testLR.mapRegion(ctx, runtime);
    RegionAccessor<AccessorType::Generic, LogicalRegion> pracc =
        prs.get_field_accessor(testLR.fid).typeify<LogicalRegion>();

    for (GenericPointInRectIterator<1> pir(testLR.bounds); pir; pir++) {
        auto lr = pracc.read(DomainPoint::from_point<1>(pir.p));
        RegionRequirement req(
            lr, READ_ONLY, EXCLUSIVE, lr
        );
        req.add_field(0);
        InlineLauncher inl(req);
        PhysicalRegion reg = runtime->map_region(ctx, inl);
        reg.wait_until_valid();
        Domain tDom = runtime->get_index_space_domain(
            ctx, reg.get_logical_region().get_index_space()
        );
        RegionAccessor<AccessorType::Generic, double> blacc =
            reg.get_field_accessor(0).typeify<double>();
        for (GenericPointInRectIterator<1> piri(tDom.get_rect<1>()); piri; piri++) {
            auto blah = blacc.read(DomainPoint::from_point<1>(piri.p));
            cout << "-->" << blah;
        }
        runtime->unmap_region(ctx, reg);
        cout << endl;
    }
    // END XXX SKG LogicalRegion TEST
    //
    cout << "*** Cleaning Up..." << endl;
    destroyLogicalStructures(
        hpcgParams,
        geometries,
        A,
        x,
        y,
        ctx,
        runtime
    );
    // XXX SKG LogicalRegion TEST
    testLR.deallocate(ctx, runtime);
    // END XXX SKG LogicalRegion TEST
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
