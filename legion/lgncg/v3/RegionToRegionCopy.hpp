/**
 * Copyright (c)      2017 Los Alamos National Security, LLC
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

#pragma once

#include "LegionStuff.hpp"
#include "LegionArrays.hpp"
#include "VectorOps.hpp"

/**
 *
 */
inline void
regionToRegionCopyTask(
    const Task *task,
    const std::vector<PhysicalRegion> &regions,
    Context ctx,
    HighLevelRuntime *lrt
) {
    int rid = 0;
    Array<floatType> src(regions[rid++], ctx, lrt);
    Array<floatType> dst(regions[rid++], ctx, lrt);

#if 1 // Debug
    std::pair<int , int> tArgs = *(std::pair<int, int> *)task->args;
    const int parentColor = tArgs.first;
    const int targetNeighbor = tArgs.second;
    std::string pcs = std::to_string(parentColor);
    std::string tns = std::to_string(targetNeighbor);
    std::string fName = "r2r-t" + pcs + "-pull-from-n" + tns + ".txt";
    PrintVector(src, fName, ctx, lrt);
#endif

    assert(src.length() == dst.length());

    const floatType *const sv = src.data();
    assert(sv);

    floatType *const dv = dst.data();
    assert(dv);

    for (size_t i = 0; i < src.length(); ++i) {
        dv[i] = sv[i];
    }
}
