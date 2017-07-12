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


/**
 *
 */
void
LogicalMGData::partition(
    SparseMatrix &A,
    Context ctx,
    HighLevelRuntime *lrt
) {
    assert(A.Ac);
    //
    auto *Asclrsf = A.sclrs->data();
    auto *Asclrsc = A.Ac->sclrs->data();
    //
    const local_int_t nrowf = Asclrsf->localNumberOfRows;
    const local_int_t ncolf = Asclrsf->localNumberOfColumns;
    const local_int_t nrowc = Asclrsc->localNumberOfRows;
    const local_int_t ncolc = Asclrsc->localNumberOfColumns;
    // First nrow items are 'local data'. After is remote data.
    std::vector<local_int_t> partLensf, partLensc;
    // First partition for local data.
    local_int_t totLenf = nrowf, totLenc = nrowc;
    partLensf.push_back(nrowf);
    partLensc.push_back(nrowc);
    // The rest are based on receive lengths.
    const int nNeighbors = A.sclrs->data()->numberOfRecvNeighbors;
    assert(nNeighbors == Asclrsc->numberOfSendNeighbors);
    //
    const local_int_t *const recvLengthf = A.recvLength->data();
    const local_int_t *const recvLengthc = A.Ac->recvLength->data();
    //
    for (int n = 0; n < nNeighbors; ++n) {
        const int recvlf = recvLengthf[n];
        const int recvlc = recvLengthc[n];
        //
        partLensf.push_back(recvlf);
        partLensc.push_back(recvlc);
        //
        totLenf += recvlf;
        totLenc += recvlc;
    }
    assert(totLenf == ncolf);
    assert(totLenc == ncolc);
    //
    f2cOperator.partition(1, ctx, lrt);
    rc.partition(1, ctx, lrt);
    //
    xc.partition(partLensc, ctx, lrt);
    Axf.partition(partLensf, ctx, lrt);
}
