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
    @file ComputeMG.hpp

    HPCG routine
 */

#include <cassert>
#include <iostream>

/*!

    @param[in] A the known system matrix.

    @param[in] r the input vector.

    @param[inout] x On exit contains the result of the multigrid V-cycle with r
    as the RHS, x is the approximation to Ax = r.

    @return returns 0 upon success and non-zero otherwise.

    @see ComputeMG
*/
inline int
ComputeMG(
    const SparseMatrix & A,
    const Vector & r,
    Vector & x
) {
    assert(x.localLength==A.localNumberOfColumns); // Make sure x contain space for halo values

    ZeroVector(x); // initialize x to zero

    int ierr = 0;
    if (A.mgData!=0) { // Go to next coarse level if defined
      int numberOfPresmootherSteps = A.mgData->numberOfPresmootherSteps;
      for (int i=0; i< numberOfPresmootherSteps; ++i) ierr += ComputeSYMGS_ref(A, r, x);
      if (ierr!=0) return ierr;
      ierr = ComputeSPMV_ref(A, x, *A.mgData->Axf); if (ierr!=0) return ierr;
      // Perform restriction operation using simple injection
      ierr = ComputeRestriction_ref(A, r);  if (ierr!=0) return ierr;
      ierr = ComputeMG_ref(*A.Ac,*A.mgData->rc, *A.mgData->xc);  if (ierr!=0) return ierr;
      ierr = ComputeProlongation_ref(A, x);  if (ierr!=0) return ierr;
      int numberOfPostsmootherSteps = A.mgData->numberOfPostsmootherSteps;
      for (int i=0; i< numberOfPostsmootherSteps; ++i) ierr += ComputeSYMGS_ref(A, r, x);
      if (ierr!=0) return ierr;
    }
    else {
      ierr = ComputeSYMGS_ref(A, r, x);
      if (ierr!=0) return ierr;
    }
    return 0;
}

