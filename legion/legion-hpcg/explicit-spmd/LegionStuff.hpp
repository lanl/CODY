/**
 * Copyright (c) 2016-2017 Los Alamos National Security, LLC
 *                         All rights reserved.
 *
 * Copyright (c) 2016      Stanford University
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

#pragma once

#include "TaskTIDs.hpp"
#include "Types.hpp"

#include "legion.h"

#include <vector>
#include <map>

#include <unistd.h>

using namespace LegionRuntime::HighLevel;
using namespace LegionRuntime::Accessor;
using namespace LegionRuntime::HighLevel;

#define RW_E READ_WRITE, EXCLUSIVE
#define RO_E READ_ONLY , EXCLUSIVE
#define WO_E WRITE_ONLY, EXCLUSIVE

#define RW_S READ_WRITE, SIMULTANEOUS
#define RO_S READ_ONLY , SIMULTANEOUS
#define WO_S WRITE_ONLY, SIMULTANEOUS

////////////////////////////////////////////////////////////////////////////////
// Task forward declarations.
////////////////////////////////////////////////////////////////////////////////
void
mainTask(
    const Task *task,
    const std::vector<PhysicalRegion> &regions,
    Context ctx, HighLevelRuntime *runtime
);

void
genProblemTask(
    const Task *task,
    const std::vector<PhysicalRegion> &regions,
    Context ctx, HighLevelRuntime *runtime
);

void
startBenchmarkTask(
    const Task *task,
    const std::vector<PhysicalRegion> &regions,
    Context ctx, HighLevelRuntime *runtime
);

global_int_t
dynCollTaskContribGIT(
    const Task *task,
    const std::vector<PhysicalRegion> &regions,
    Context ctx, HighLevelRuntime *runtime
);

floatType
dynCollTaskContribFT(
    const Task *task,
    const std::vector<PhysicalRegion> &regions,
    Context ctx, HighLevelRuntime *runtime
);

void
regionToRegionCopyTask(
    const Task *task,
    const std::vector<PhysicalRegion> &regions,
    Context ctx, HighLevelRuntime *runtime
);

void
registerCollectiveOpsTasks(void);

void
registerVectorOpTasks(void);

////////////////////////////////////////////////////////////////////////////////
// Task Registration
////////////////////////////////////////////////////////////////////////////////
inline void
registerTasks(void)
{
    TaskVariantRegistrar tvr(MAIN_TID, "mainTask");
    tvr.add_constraint(ProcessorConstraint(Processor::LOC_PROC));
    Runtime::preregister_task_variant<mainTask>(tvr, "mainTask");
    Runtime::set_top_level_task_id(MAIN_TID);
    //
    HighLevelRuntime::register_legion_task<genProblemTask>(
        GEN_PROB_TID /* task id */,
        Processor::LOC_PROC /* proc kind  */,
        true /* single */,
        true /* index */,
        AUTO_GENERATE_ID,
        TaskConfigOptions(false /* leaf task */),
        "genProblemTask"
    );
    HighLevelRuntime::register_legion_task<startBenchmarkTask>(
        START_BENCHMARK_TID /* task id */,
        Processor::LOC_PROC /* proc kind  */,
        true /* single */,
        true /* index */,
        AUTO_GENERATE_ID,
        TaskConfigOptions(false /* leaf task */),
        "startBenchmarkTask"
    );
    HighLevelRuntime::register_legion_task<regionToRegionCopyTask>(
        REGION_TO_REGION_COPY_TID /* task id */,
        Processor::LOC_PROC /* proc kind  */,
        true /* single */,
        true /* index */,
        AUTO_GENERATE_ID,
        TaskConfigOptions(true /* leaf task */),
        "regionToRegionCopyTask"
    );
    //
    registerCollectiveOpsTasks();
    //
    registerVectorOpTasks();
}

////////////////////////////////////////////////////////////////////////////////
inline void
updateMappers(
    Machine machine,
    HighLevelRuntime *runtime,
    const std::set<Processor> &local_procs
) {
#if 0 // Disable for now.
    for (const auto &p : local_procs) {
        runtime->replace_default_mapper(new CGMapper(machine, runtime, p), p);
    }
#endif
}

////////////////////////////////////////////////////////////////////////////////
inline void
LegionInit(void)
{
    registerTasks();
    HighLevelRuntime::set_registration_callback(updateMappers);
}

/**
 * Courtesy of some other legion code.
 */
template <unsigned DIM, typename T>
inline bool
offsetsAreDense(
    const Rect<DIM> &bounds,
    const LegionRuntime::Accessor::ByteOffset *offset
) {
#ifdef LGNCG_ASSUME_DENSE_OFFSETS
    return true;
#else
    off_t exp_offset = sizeof(T);
    for (unsigned i = 0; i < DIM; i++) {
        bool found = false;
        for (unsigned j = 0; j < DIM; j++)
            if (offset[j].offset == exp_offset) {
                found = true;
                exp_offset *= (bounds.hi[j] - bounds.lo[j] + 1);
                break;
            }
        if (!found) return false;
    }
    return true;
#endif
}

/**
 * Courtesy of some other legion code.
 */
inline bool
offsetMismatch(
    int i,
    const LegionRuntime::Accessor::ByteOffset *off1,
    const LegionRuntime::Accessor::ByteOffset *off2
) {
#ifdef LGNCG_ASSUME_OFFSETS_MATCH
    return false;
#else
    while (i-- > 0) {
        if ((off1++)->offset != (off2++)->offset) return true;
    }
    return false;
#endif
}

/**
 * Convenience routine to get a task's ID.
 */
inline int
getTaskID(
    const Task *task
) {
    return task->index_point.point_data[0];
}

/**
 *
 */
inline size_t
getNumProcs(void)
{
    size_t nProc = 0;
    std::set<Processor> allProcs;
    Realm::Machine::get_machine().get_all_processors(allProcs);
    for (auto &p : allProcs) {
        if (p.kind() == Processor::LOC_PROC) nProc++;
    }
    return nProc;
}
