/**
 * Copyright (c) 2014      Los Alamos National Security, LLC
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

#include "hpcg-problem.h"
#include "lgncg.h"

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <inttypes.h>

namespace {
// convenience namespace aliases
namespace lrthl = LegionRuntime::HighLevel;

// task ids
enum TID {
    MAIN_TASK_ID = 0,
    CHECK_CG_TASK_ID
};

const char DRIVER_NAME[] = "lgn-hpcg";
const int DRIVER_VER     = 0;
const int DRIVER_SUBVER  = 6;

// FIXME HACK
#ifndef PRId64
#define PRId64 "lli"
#endif

struct DriverParams {
    int64_t nx, ny, nz;
    int64_t nSubrgns;
    int64_t maxIters;
    int64_t nMGLevels;
    double tolerance;
    // flag indicating whether or not we are going to use a preconditioner
    bool doPreconditioning;
    // approximating a 27-point finite element/volume/difference 3D stencil
    static const int64_t stencilSize = 27;
    /**
     * default constructor.
     */
    DriverParams(void) {
        nx = 16;
        ny = 16;
        nz = 16;
        // FIXME - change back to 2 once MG is working...
        nSubrgns  = 1;
        maxIters  = 128;
        nMGLevels = 4;
        tolerance = 0.001;
        // SKG FIXME - make doPreconditioning true by default
        doPreconditioning = false;
    }
};

/**
 * takes command line arguments and populates the driver params.
 */
void
getParamsFromArgs(DriverParams &params)
{
    using namespace LegionRuntime::HighLevel;
    // process any user-supplied arguments
    const InputArgs &cArgs = HighLevelRuntime::get_input_args();
    for (int i = 1; i < cArgs.argc; ++i) {
        if (!strcmp(cArgs.argv[i], "-s")) {
            params.nSubrgns = (int64_t)strtol(cArgs.argv[++i], NULL, 10);
        }
        if (!strcmp(cArgs.argv[i], "-i")) {
            params.maxIters = (int64_t)strtol(cArgs.argv[++i], NULL, 10);
        }
        if (!strcmp(cArgs.argv[i], "-nx")) {
            params.nx = (int64_t)strtol(cArgs.argv[++i], NULL, 10);
        }
        if (!strcmp(cArgs.argv[i], "-ny")) {
            params.ny = (int64_t)strtol(cArgs.argv[++i], NULL, 10);
        }
        if (!strcmp(cArgs.argv[i], "-nz")) {
            params.nz = (int64_t)strtol(cArgs.argv[++i], NULL, 10);
        }
        if (!strcmp(cArgs.argv[i], "-nmg")) {
            params.nMGLevels = (int64_t)strtol(cArgs.argv[++i], NULL, 10);
        }
        if (!strcmp(cArgs.argv[i], "-no-precond")) {
            params.doPreconditioning = false;
        }
        if (!strcmp(cArgs.argv[i], "-h") || !strcmp(cArgs.argv[i], "-help")) {
            std::cout << "usage: " << DRIVER_NAME
                      << " [-nx X] [-ny Y] [-nz Z] [-i MAX_ITERS]"
                         " [-s N_SUBR] [-nmg N_MGL] [-no-precond] [-h | --help]"
                      << std::endl;
            exit(EXIT_SUCCESS);
        }
    }
    assert(params.nx * params.ny * params.nz > 0);
    assert(params.nSubrgns > 0);
}

/**
 * let's the user know the setup and other interesting things.
 */
void
echoBanner(const DriverParams &params,
           const Problem &problem)
{
    printf("%s %d.%d\n", DRIVER_NAME, DRIVER_VER, DRIVER_SUBVER);
    printf("Legion Setup Summary:\n"
           "  Number of Subregions: %"PRId64"\n", params.nSubrgns);
    printf("Global Problem Dimensions:\n"
           "  Global nx: %"PRId64"\n"
           "  Global ny: %"PRId64"\n"
           "  Global nz: %"PRId64"\n",
           params.nx, params.ny, params.nz);
    printf("********** Problem Summary ***********:\n");
    printf("Linear System Information:\n"
           "  Number of Equations: %"PRId64"\n"
           "  Number of Nonzero Terms: %"PRId64"\n",
           problem.A.nRows, problem.A.tNon0);
}

/**
 * interface for launching tasks that verify that we are close enough to the
 * known, correct solution.
 */
void
checkcg(lgncg::Vector &x,
        LegionRuntime::HighLevel::Context &ctx,
        LegionRuntime::HighLevel::HighLevelRuntime *lrt)
{
    using namespace LegionRuntime::HighLevel;
    using namespace lgncg;

    int idx = 0;
    // no per-task arguments.
    ArgumentMap argMap;
    DVector targs;
    targs = x;
    for (int i = 0; i < x.lDom().get_volume(); ++i) {
        targs.sgb = x.sgb()[i];
        argMap.set_point(DomainPoint::from_point<1>(Point<1>(i)),
                         TaskArgument(&targs, sizeof(targs)));
    }
    IndexLauncher il(CHECK_CG_TASK_ID, x.lDom(),
                     TaskArgument(NULL, 0), argMap);
    // x
    il.add_region_requirement(
        RegionRequirement(x.lp(), 0, READ_ONLY, EXCLUSIVE, x.lr)
    );
    il.add_field(idx++, x.fid);
    // execute the thing...
    FutureMap fm = lrt->execute_index_space(ctx, il);
    fm.wait_all_results();
}

/**
 * makes sure we are getting close enough to the known, correct answer.
 */
void
checkCGTask(const lrthl::Task *task,
            const std::vector<lrthl::PhysicalRegion> &rgns,
            lrthl::Context ctx,
            lrthl::HighLevelRuntime *lrt)
{
    using namespace LegionRuntime::HighLevel;
    using namespace LegionRuntime::Accessor;
    using namespace lgncg;
    using LegionRuntime::Arrays::Rect;
    // silence warnings
    (void)ctx; (void)lrt;
    // we are expecting exactly one region in this routine
    assert(1 == rgns.size());
    size_t rid = 0;
    const int tid = task->index_point.point_data[0];
    // get the args
    const lgncg::DVector &vx = *(lgncg::DVector *)task->local_args;
    // name the physical regions
    const PhysicalRegion &xpr = rgns[rid++];

    typedef RegionAccessor<AccessorType::Generic, double> GDRA;
    GDRA x = xpr.get_field_accessor(vx.fid).typeify<double>();
    typedef GenericPointInRectIterator<1> GPRI1D;
    typedef DomainPoint DomPt;
    // this is what we expect the final answer to be
    static const double expectedAnswer = 1.0;
    static const double epsilon = 0.1;
    for (GPRI1D xpi(vx.sgb); xpi; xpi++) {
        double val = x.read(DomPt::from_point<1>(xpi.p));
        if (fabs(expectedAnswer - val) > epsilon) {
            printf("%d WARNING: check failure. "
                   "expected: %lf actual: %lf\n",
                   tid, expectedAnswer, val);
        }
    }
}

/* ////////////////////////////////////////////////////////////////////////// */
/* ////////////////////////////////////////////////////////////////////////// */
/**
 * main task implementation. consider this the real main. that is, this is the
 * first task that legion will start. from here all other things happen.
 */
void
mainTask(const lrthl::Task *task,
         const std::vector<lrthl::PhysicalRegion> &rgns,
         lrthl::Context ctx,
         lrthl::HighLevelRuntime *lrt)
{
    (void)rgns; (void)task;
    DriverParams params;
    getParamsFromArgs(params);
    // now that we have our params populated, construct the problem
    Problem problem(params.nx, params.ny,
                    params.nz, params.stencilSize,
                    params.nMGLevels, ctx, lrt);
    // now that we have all problem-related info, let the user know the setup
    echoBanner(params, problem);
    // sets initial conditions for all grid levels
    problem.setInitialConds(ctx, lrt);
    double start = LegionRuntime::TimeStamp::get_current_time_in_micros();
    printf("starting solve...\n");
    // solve the thing
    lgncg::cgSolv(problem.A,
                  problem.b,
                  params.tolerance,
                  params.maxIters,
                  problem.x,
                  params.nSubrgns,
                  params.doPreconditioning,
                  ctx,
                  lrt);
    double stop = LegionRuntime::TimeStamp::get_current_time_in_micros();
    printf("  . done in: %7.3lf ms\n", (stop - start) * 1e-3);
    start = LegionRuntime::TimeStamp::get_current_time_in_micros();
    printf("starting cg checks...\n");
    checkcg(problem.x, ctx, lrt);
    stop = LegionRuntime::TimeStamp::get_current_time_in_micros();
    printf("  . done in: %7.3lf ms\n", (stop - start) * 1e-3);
    // cleanup
    problem.free(ctx, lrt);
}

/**
 * registers all driver-related tasks.
 */
void
init(void)
{
    using namespace LegionRuntime::HighLevel;
    // no output buffering
    setbuf(stdout, NULL);
    // let legion know what the top-level task ID is
    HighLevelRuntime::set_top_level_task_id(MAIN_TASK_ID);
    // register the main task
    HighLevelRuntime::register_legion_task<mainTask>(
        MAIN_TASK_ID /* task id */,
        Processor::LOC_PROC /* proc kind  */,
        true /* single */,
        false /* index */,
        AUTO_GENERATE_ID,
        TaskConfigOptions(false /* leaf task */),
        "main-task"
    );
    HighLevelRuntime::register_legion_task<checkCGTask>(
        CHECK_CG_TASK_ID /* task id */,
        Processor::LOC_PROC /* proc kind  */,
        false /* single */,
        true /* index */,
        AUTO_GENERATE_ID,
        TaskConfigOptions(true /* leaf task */),
        "check-cg-task"
    );
    // now do any Problem-related init
    Problem::init();
    // also init the solver
    lgncg::init();
}

/**
 * pre-HighLevelRuntime::start initialization.
 */
void
preStartInit(void)
{
    // register driver tasks and all other tasks required for this calculation
    init();
}
} // end namespace

/* ////////////////////////////////////////////////////////////////////////// */
/* ////////////////////////////////////////////////////////////////////////// */
int
main(int argc,
     char **argv)
{
    using namespace LegionRuntime::HighLevel;

    preStartInit();
    // and... go!
    HighLevelRuntime::start(argc, argv);
}
