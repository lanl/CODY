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

#include "hpcg.hpp"
#include "Geometry.hpp"

/**
 *
 */
std::ostream &
operator<<(std::ostream &os, const HPCG_Params &params) {
    using namespace std;
    //
    os << "runningTime: " << params.runningTime << endl;
    os << "comm_size: "   << params.comm_size << endl;
    os << "comm_rank: "   << params.comm_rank << endl;
    os << "numThreads: "  << params.numThreads << endl;
    os << "nx: "          << params.nx << endl;
    os << "ny: "          << params.ny << endl;
    os << "nz: "          << params.nz << endl;
    return os;
}

std::ostream &
operator<<(std::ostream &os, const Geometry &geom) {
    using namespace std;
    os << "Size: "       << geom.size << endl;
    os << "Rank: "       << geom.rank << endl;
    os << "numThreads: " << geom.numThreads << endl;
    os << "nx: "         << geom.nx << endl;
    os << "ny: "         << geom.ny << endl;
    os << "nz: "         << geom.nz << endl;
    os << "npx: "        << geom.npx << endl;
    os << "npy: "        << geom.npy << endl;
    os << "npz: "        << geom.npz << endl;
    os << "ipx: "        << geom.ipx << endl;
    os << "ipy: "        << geom.ipy << endl;
    os << "ipz: "        << geom.ipz << endl;
    return os;
}
