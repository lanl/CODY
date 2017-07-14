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

/*!
    @file hpcg.hpp

    HPCG data structures and functions
 */

#pragma once

#include "LegionStuff.hpp"

#define HPCG_STENCIL  27
#define NUM_MG_LEVELS 4

struct HPCG_Params {
    int commSize ; //!< Total number of shards.
    int numThreads; //!< This process' number of threads.
    int nx; //!< Number of x-direction grid points for each local subdomain.
    int ny; //!< Number of y-direction grid points for each local subdomain.
    int nz; //!< Number of z-direction grid points for each local subdomain.
    //!< Number of seconds to run the timed portion of the benchmark.
    int runningTime;
    int stencilSize; //!< Size of the stencil
};

/**
 *
 */
inline void
emit(const HPCG_Params &params) {
    using namespace std;
    //
    cout << "runningTime: " << params.runningTime << endl;
    cout << "commSize: "    << params.commSize << endl;
    cout << "numThreads: "  << params.numThreads << endl;
    cout << "nx: "          << params.nx << endl;
    cout << "ny: "          << params.ny << endl;
    cout << "nz: "          << params.nz << endl;
}

extern int HPCG_Init(HPCG_Params &params, const SPMDMeta &spmdMeta);
extern int HPCG_Finalize(void);
