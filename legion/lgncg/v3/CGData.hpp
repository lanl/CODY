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
 @file CGData.hpp

 HPCG data structure
 */

#pragma once

#include "LegionMatrices.hpp"
#include "LegionArrays.hpp"

struct LogicalCGData {
    // Launch domain
    LegionRuntime::HighLevel::Domain launchDomain;
    //
    LogicalArray<floatType> r;  //!< residual vector
    LogicalArray<floatType> z;  //!< preconditioned residual vector
    LogicalArray<floatType> p;  //!< direction vector
    LogicalArray<floatType> Ap; //!< Krylov vector

protected:
    std::deque<LogicalItemBase *> mLogicalItems;

    /**
     * Order matters here. If you update this, also update unpack.
     */
    void
    mPopulateRegionList(void) {
        mLogicalItems = {&r, &z, &p, &Ap};
    }
};

struct CGData {
    Array<floatType> r;  //!< residual vector
    Array<floatType> z;  //!< preconditioned residual vector
    Array<floatType> p;  //!< direction vector
    Array<floatType> Ap; //!< Krylov vector
};

/*!
 Constructor for the data structure of CG vectors.

    @param[in]  A    the data structure that describes the problem matrix and
                     its structure
    @param[out] data the data structure for CG vectors that will be allocated to
                     get it ready for use in CG iterations
 */
inline void
InitializeSparseCGData(
    SparseMatrix &A,
    CGData &data)
{
  local_int_t nrow = A.localNumberOfRows;
  local_int_t ncol = A.localNumberOfColumns;
  InitializeVector(data.r, nrow);
  InitializeVector(data.z, ncol);
  InitializeVector(data.p, ncol);
  InitializeVector(data.Ap, nrow);
  return;
}

/*!
 Destructor for the CG vectors data.

 @param[inout] data the CG vectors data structure whose storage is deallocated
 */
inline void DeleteCGData(CGData & data) {

  DeleteVector (data.r);
  DeleteVector (data.z);
  DeleteVector (data.p);
  DeleteVector (data.Ap);
  return;
}
