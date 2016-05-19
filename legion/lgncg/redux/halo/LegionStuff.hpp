/**
 * Copyright (c) 2016      Los Alamos National Security, LLC
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

#include "CGMapper.h"
#include "legion.h"

#include <vector>

using namespace LegionRuntime::HighLevel;
using namespace LegionRuntime::HighLevel;
using namespace LegionRuntime::Accessor;

////////////////////////////////////////////////////////////////////////////////
// Task IDs
////////////////////////////////////////////////////////////////////////////////
enum {
    MAIN_TID = 0,
    SPMD_INIT_TID
};

////////////////////////////////////////////////////////////////////////////////
// SPMD metadata.
////////////////////////////////////////////////////////////////////////////////
struct SPMDMeta {
    // Number of participants in SPMD computation.
    int nRanks;
};

////////////////////////////////////////////////////////////////////////////////
// SPMD context information.
////////////////////////////////////////////////////////////////////////////////
struct SPMDContext {
    // Unique identifier for SPMD work.
    int rank;
    // Number of participants in SPMD.
    int nRanks;
};

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
spmdInitTask(
    const Task *task,
    const std::vector<PhysicalRegion> &regions,
    Context ctx, HighLevelRuntime *runtime
);

////////////////////////////////////////////////////////////////////////////////
// Task Registration
////////////////////////////////////////////////////////////////////////////////
// TODO: be explicit about leaf tasks in registration.
inline void
registerTasks(void) {
    {
    TaskVariantRegistrar tvr(MAIN_TID, "mainTask");
    tvr.add_constraint(ProcessorConstraint(Processor::LOC_PROC));
    Runtime::preregister_task_variant<mainTask>(tvr, "mainTask");
    Runtime::set_top_level_task_id(MAIN_TID);
    }
    //
    HighLevelRuntime::register_legion_task<spmdInitTask>(
        SPMD_INIT_TID /* task id */,
        Processor::LOC_PROC /* proc kind  */,
        true /* single */,
        true /* index */,
        AUTO_GENERATE_ID,
        TaskConfigOptions(true /* leaf task */),
        "spmdInitTask"
    );
}

////////////////////////////////////////////////////////////////////////////////
static void
updateMappers(
    Machine machine,
    HighLevelRuntime *runtime,
    const std::set<Processor> &local_procs
) {
    for (const auto &p : local_procs) {
        runtime->replace_default_mapper(new CGMapper(machine, runtime, p), p);
    }
}

////////////////////////////////////////////////////////////////////////////////
inline void
LegionInit(void)
{
    registerTasks();
    Runtime::set_registration_callback(updateMappers);
}

////////////////////////////////////////////////////////////////////////////////
class FieldHelperBase {
protected:
  // this is shared by all FieldHelper<T>'s
  static FieldID sNextID;
};

////////////////////////////////////////////////////////////////////////////////
template <typename T>
class FieldHelper : protected FieldHelperBase {
protected:
    const char *name;
    FieldID fid;

public:
    static const FieldID ASSIGN_STATIC_ID = AUTO_GENERATE_ID - 1;
    ////////////////////////////////////////////////////////////////////////////
    FieldHelper(
        const char *mName,
        FieldID mFID = ASSIGN_STATIC_ID
    ) : name(mName),
        fid(mFID)
    {
        if(fid == ASSIGN_STATIC_ID) fid = sNextID++;
    }

    ////////////////////////////////////////////////////////////////////////////
    ~FieldHelper(void) = default;

    ////////////////////////////////////////////////////////////////////////////
    operator
    FieldID(void) const
    {
        assert(fid != AUTO_GENERATE_ID);
        return fid;
    }

    ////////////////////////////////////////////////////////////////////////////
    void
    allocate(
        Runtime *runtime,
        Context ctx,
        FieldSpace fs
    ) {
        FieldAllocator fa = runtime->create_field_allocator(ctx, fs);
        fid = fa.allocate_field(sizeof(T), fid);
        runtime->attach_name(fs, fid, name);
    }

    ////////////////////////////////////////////////////////////////////////////
    template <typename AT>
    RegionAccessor<AT, T>
    accessor(const PhysicalRegion &pr)
    {
        assert(fid != AUTO_GENERATE_ID);
        return pr.get_field_accessor(
                  fid
               ).template typeify<T>().template convert<AT>();
    }

    ////////////////////////////////////////////////////////////////////////////
    template <typename AT>
    RegionAccessor<AT, T>
    foldAccessor(const PhysicalRegion &pr)
    {
        assert(fid != AUTO_GENERATE_ID);
        std::vector<FieldID> fields;
        pr.get_fields(fields);
        assert((fields.size() == 1) && (fields[0] == fid));
        return pr.get_accessor().template typeify<T>().template convert<AT>();
    }
};

/**
 * courtesy of some other legion code.
 */
template <unsigned DIM, typename T>
inline bool
offsetsAreDense(const Rect<DIM> &bounds,
                const LegionRuntime::Accessor::ByteOffset *offset)
{
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
}

/**
 * courtesy of some other legion code.
 */
inline bool
offsetMismatch(int i,
               const LegionRuntime::Accessor::ByteOffset *off1,
               const LegionRuntime::Accessor::ByteOffset *off2)
{
    while (i-- > 0) {
        if ((off1++)->offset != (off2++)->offset) return true;
    }
    return false;
}

/**
 * convenience routine to get a task's ID
 */
inline int
getTaskID(
    const LegionRuntime::HighLevel::Task *task
) {
    return task->index_point.point_data[0];
}
