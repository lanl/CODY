/*
 * Copyright (c) 2014-2016 Los Alamos National Security, LLC
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

#pragma once

#include "LegionStuff.hpp"
#include "LegionItems.hpp"
#include "LegionArrays.hpp"

#include "hpcg.hpp"
#include "Geometry.hpp"
#include "MGData.hpp"

#include <vector>
#include <map>
#include <cassert>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
struct SparseMatrixScalars {
    //Max number of non-zero elements in any row
    local_int_t maxNonzerosPerRow = 0;
    //total number of matrix rows across all processes
    global_int_t totalNumberOfRows = 0;
    //total number of matrix nonzeros across all processes
    global_int_t totalNumberOfNonzeros = 0;
    //number of rows local to this process
    local_int_t localNumberOfRows = 0;
    //number of columns local to this process
    local_int_t localNumberOfColumns = 0;
    //number of nonzeros local to this process
    local_int_t localNumberOfNonzeros = 0;
    //number of entries that are external to this process
    local_int_t numberOfExternalValues = 0;
    //number of neighboring processes that will be send local data
    int numberOfSendNeighbors = 0;
    //number of neighboring processes that i'll get data from
    int numberOfRecvNeighbors = 0;
    //total number of entries to be sent
    local_int_t totalToBeSent = 0;
    // size of serialized buffer used for Synchronizers.
    size_t sizeofSynchronizersBuffer = 0;
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/**
 * Holds structures required for task synchronization.
 */
struct Synchronizers {
    //
    PhaseBarriers myPhaseBarriers;
    //
    std::map< int, std::vector<PhaseBarriers> > neighborPhaseBarriers;

    /**
     *
     */
    Synchronizers(void) = default;

    /**
     *
     */
    template <class Archive>
    void
    serialize(Archive &ar)
    {
        ar(myPhaseBarriers, neighborPhaseBarriers);
    }

    /**
     *
     */
    static void
    deserialize(
        char *rawData,
        size_t rawDataSizeInB,
        Synchronizers &outRes
    ) {
        // Deserialize argument data
        std::stringstream solvArgsSS;
        solvArgsSS.write(rawData, rawDataSizeInB);
        {   // Scoped to guarantee flushing, etc.
            cereal::BinaryInputArchive ia(solvArgsSS);
            ia(outRes);
        }
    }
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
struct LogicalSparseMatrix {
public:
    // launch domain
    LegionRuntime::HighLevel::Domain launchDomain;
    //
    LogicalArray<Geometry> geometries;
    //
    LogicalArray<SparseMatrixScalars> localData;
    //
    LogicalArray2D<floatType> matrixValues;

    /**
     *
     */
    LogicalSparseMatrix(void) = default;

    /**
     *
     */
    void
    allocate(
        const Geometry &geom,
        LegionRuntime::HighLevel::Context &ctx,
        LegionRuntime::HighLevel::HighLevelRuntime *lrt
    ) {
        const auto size = geom.size;
        const auto globalXYZ = getGlobalXYZ(geom);
        const auto stencilSize = geom.stencilSize;

        geometries.allocate(size, ctx, lrt);
        localData.allocate(size, ctx, lrt);
        matrixValues.allocate(globalXYZ, stencilSize, ctx, lrt);
    }

    /**
     *
     */
    void
    partition(
        int64_t nParts,
        LegionRuntime::HighLevel::Context &ctx,
        LegionRuntime::HighLevel::HighLevelRuntime *lrt
    ) {
                geometries.partition(nParts, ctx, lrt);
                 localData.partition(nParts, ctx, lrt);
              matrixValues.partition(nParts, ctx, lrt);
        // just pick a structure that has a representative launch domain.
        launchDomain = geometries.launchDomain;
    }

    /**
     * Cleans up and returns all allocated resources.
     */
    void
    deallocate(
        LegionRuntime::HighLevel::Context &ctx,
        LegionRuntime::HighLevel::HighLevelRuntime *lrt
    ) {
                geometries.deallocate(ctx, lrt);
                 localData.deallocate(ctx, lrt);
    }
};
