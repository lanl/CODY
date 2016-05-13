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

/**
 * A simple example illustrating a halo exchange using control replication. 
 */

#include <iostream>
#include <cstdlib>

#include "hpcg.hpp"
#include "Geometry.hpp"
#include "CheckAspectRatio.hpp"

#include "LegionStuff.hpp"
#include "CGMapper.h"

using namespace std;

LegionRuntime::Logger::Category Logger("LGNCG");

/**
 * SPMD Main Task //////////////////////////////////////////////////////////////
 *
 */
void
spmdMainTask(
    const Task *task,
    const std::vector<PhysicalRegion> &regions,
    Context ctx, HighLevelRuntime *runtime
) {
    HPCG_Params params;
    HPCG_Init(params);
    // Check if QuickPath option is enabled.  If the running time is set to
    // zero, we minimize all paths through the program
    bool quickPath = (params.runningTime == 0);
    // Number of MPI processes, My process ID
    // TODO FIXME get from mapper
    int size = params.comm_size, rank = params.comm_rank;
    //
    local_int_t nx,ny,nz;
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
}

/**
 * Main Task ///////////////////////////////////////////////////////////////////
 * Responsible for setting up the SPMD launch.
 */
void
mainTask(
    const Task *task,
    const std::vector<PhysicalRegion> &regions,
    Context ctx, HighLevelRuntime *runtime
) {
    // step 1: ask the mapper how many shards we should even have
    //int nShards = runtime->get_tunable_value(ctx, CGMapper::TID_NUM_SHARDS);
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
