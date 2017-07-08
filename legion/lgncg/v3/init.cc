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

#include <cstdio>
#include <cstring>

#include "hpcg.hpp"
#include "ReadHpcgDat.hpp"

#include "LegionStuff.hpp"

static int
startswith(
    const char *s,
    const char *prefix
) {
    size_t n = strlen( prefix );
    if (strncmp( s, prefix, n ))
        return 0;
    return 1;
}

int
HPCG_Init(
    HPCG_Params &params,
    const SPMDMeta &spmdMeta
) {
    int iparams[4];
    char cparams[4][6] = {"--nx=", "--ny=", "--nz=", "--rt="};
    // Initialize iparams
    for (int i = 0; i < 4; ++i) iparams[i] = 0;
    // process any user-supplied arguments
    const InputArgs &cArgs = HighLevelRuntime::get_input_args();
    for (int i = 1; i < cArgs.argc; ++i) {
        for (int j = 0; j < 4; ++j) {
            if (startswith(cArgs.argv[i], cparams[j])) {
                if (sscanf(cArgs.argv[i] + strlen(cparams[j]),
                    "%d", iparams+j) != 1 || iparams[j] < 10) {
                    iparams[j] = 0;
                }
            }
        }
    }
    // Check if --rt was specified on the command line
    // Assume runtime was not specified and will be read from the hpcg.dat file
    int *rt = iparams + 3;
    // If --rt was specified, we already have the runtime, so don't read it from
    // file
    if (iparams[3]) rt = 0;
    // no geometry arguments on the command line
    if (!iparams[0] && !iparams[1] && !iparams[2]) {
        ReadHpcgDat(iparams, rt);
    }
    // Check for small or unspecified nx, ny, nz values If any dimension is less
    // than 16, make it the max over the other two dimensions, or 16, whichever
    // is largest
    for (int i = 0; i < 3; ++i) {
        if (iparams[i] < 16)
            for (int j = 1; j <= 2; ++j)
                if (iparams[(i+j) % 3] > iparams[i])
                    iparams[i] = iparams[(i+j) % 3];
        if (iparams[i] < 16)
            iparams[i] = 16;
    }
    //
    params.nx = iparams[0];
    params.ny = iparams[1];
    params.nz = iparams[2];
    //
    params.runningTime = iparams[3];
    //
    params.commSize = spmdMeta.nRanks;
    //
    params.numThreads = 1;
    //
    params.stencilSize = HPCG_STENCIL;
    //
    return 0;
}
