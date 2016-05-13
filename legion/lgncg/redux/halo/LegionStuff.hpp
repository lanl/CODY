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

#include <vector>

#include "TaskIDs.hpp"

#include "legion.h"

using namespace LegionRuntime::HighLevel;
using namespace LegionRuntime::HighLevel;
using namespace LegionRuntime::Accessor;

////////////////////////////////////////////////////////////////////////////////
// Task forward declarations.
////////////////////////////////////////////////////////////////////////////////
void
mainTask(
    const Task *task,
    const std::vector<PhysicalRegion> &regions,
    Context ctx, HighLevelRuntime *runtime
);

////////////////////////////////////////////////////////////////////////////////
inline void
registerTasks(void) {
    TaskVariantRegistrar tvr(MAIN_TID, "mainTask");
    tvr.add_constraint(ProcessorConstraint(Processor::LOC_PROC));
    Runtime::preregister_task_variant<mainTask>(tvr, "mainTask");
    Runtime::set_top_level_task_id(MAIN_TID);
    //
}

////////////////////////////////////////////////////////////////////////////////
inline void
LegionInit(void)
{
    registerTasks();
}

////////////////////////////////////////////////////////////////////////////////
class FieldHelperBase {
protected:
  // this is shared by all FieldHelper<T>'s
  static FieldID sNextID;
};
//

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
        const char *_name,
        FieldID _fid = ASSIGN_STATIC_ID
    ) : name(_name),
        fid(_fid)
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
    accessor(const PhysicalRegion& pr)
    {
        assert(fid != AUTO_GENERATE_ID);
        return pr.get_field_accessor(
                  fid
               ).template typeify<T>().template convert<AT>();
    }

    ////////////////////////////////////////////////////////////////////////////
    template <typename AT>
    RegionAccessor<AT, T>
    foldAccessor(const PhysicalRegion& pr)
    {
        assert(fid != AUTO_GENERATE_ID);
        std::vector<FieldID> fields;
        pr.get_fields(fields);
        assert((fields.size() == 1) && (fields[0] == fid));
        return pr.get_accessor().template typeify<T>().template convert<AT>();
    }
};
