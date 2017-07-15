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

/*!
    @file ProblemCommon.hpp

    HPCG routine
 */

#pragma once

#include "LegionStuff.hpp"
#include "LegionArrays.hpp"
#include "LegionMatrices.hpp"
#include "CollectiveOps.hpp"

/**
 * Tally total number of non-zeros in simulation.
 * TODO bring into coll ops.
 */
static inline global_int_t
getTotalNumberOfNonZeros(
    SparseMatrix &A,
    global_int_t localNumberOfNonzeros,
    Context ctx,
    Runtime *runtime
) {
    Item< DynColl<global_int_t> > *dcars = A.dcAllRedSumGI;
    dcars->data()->localBuffer = localNumberOfNonzeros;
    //
    TaskLauncher tlLocalNZ(
        LOCAL_NONZEROS_TID,
        TaskArgument(NULL, 0)
    );
    //
    tlLocalNZ.add_region_requirement(
        RegionRequirement(
            dcars->logicalRegion,
            RO_E,
            dcars->logicalRegion
        )
    ).add_field(dcars->getFieldID());
    //
    Future future = runtime->execute_task(ctx, tlLocalNZ);
    //
    auto &dyncol = dcars->data()->dc;
    runtime->defer_dynamic_collective_arrival(ctx, dyncol, future);
    //
    dyncol = runtime->advance_dynamic_collective(ctx, dyncol);
    //
    Future fSum = runtime->get_dynamic_collective_result(ctx, dyncol);

    return fSum.get<global_int_t>();
}
