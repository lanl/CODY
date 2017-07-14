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

#include "hpcg.hpp"
#include "Geometry.hpp"
#include "LegionArrays.hpp"

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

/**
 *
 */
inline void
emit(const Geometry &geom) {
    using namespace std;
    cout << "Size: "       << geom.size << endl;
    cout << "Rank: "       << geom.rank << endl;
    cout << "numThreads: " << geom.numThreads << endl;
    cout << "nx: "         << geom.nx << endl;
    cout << "ny: "         << geom.ny << endl;
    cout << "nz: "         << geom.nz << endl;
    cout << "stencilSize: "<< geom.stencilSize<< endl;
    cout << "npx: "        << geom.npx << endl;
    cout << "npy: "        << geom.npy << endl;
    cout << "npz: "        << geom.npz << endl;
    cout << "ipx: "        << geom.ipx << endl;
    cout << "ipy: "        << geom.ipy << endl;
    cout << "ipz: "        << geom.ipz << endl;
}
