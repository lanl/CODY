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
      @file LegionMGData.cc

      HPCG data structure
 */

#include "LegionMGData.hpp"

#include "LegionArrays.hpp"
#include "LegionMatrices.hpp"

#include <cassert>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

    void
LogicalMGData::allocate(
    const std::string &name,
    SparseMatrix &A,
    LegionRuntime::HighLevel::Context ctx,
    LegionRuntime::HighLevel::HighLevelRuntime *lrt
) {
    #define aalloca(sName, size, ctx, rtp)                                     \
    do {                                                                       \
        sName.allocate(name + "-" #sName, size, ctx, rtp);                     \
    } while(0)

    assert(A.Ac);

    auto *Afsclrs = A.sclrs->data();
    auto *Acsclrs = A.Ac->sclrs->data();
    //
    const local_int_t nrowf = Afsclrs->localNumberOfRows;
    const local_int_t ncolf = Afsclrs->localNumberOfColumns;
    //
    const local_int_t nrowc = Acsclrs->localNumberOfRows;
    const local_int_t ncolc = Acsclrs->localNumberOfColumns;
    //
    aalloca(f2cOperator,  nrowf, ctx, lrt);
    aalloca(rc,           nrowc, ctx, lrt);
    aalloca(xc,           ncolc, ctx, lrt);
    aalloca(Axf,          ncolf, ctx, lrt);

    #undef aalloca
}
