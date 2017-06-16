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
 @file Vector.hpp

 HPCG data structures for dense vectors
 */

#pragma once

#include <cassert>
#include <cstdlib>
#include <iomanip>

#include "Geometry.hpp"
#include "hpcg.hpp"
#include "LegionArrays.hpp"

/*!
  Fill the input vector with zero values.

  @param[inout] v - On entrance v is initialized, on exit all its values are zero.
 */
inline void
ZeroVector(Array<floatType> &v)
{
    local_int_t localLength = v.length();
    double *vv = v.data();
    for (int i=0; i<localLength; ++i) vv[i] = 0.0;
    return;
}

/*!
  Copy input vector to output vector.

  @param[in] v Input vector
  @param[in] w Output vector
 */
inline void
CopyVector(
    Array<floatType> &v,
    Array<floatType> &w
) {
    local_int_t localLength = v.length();
    assert(w.length() >= localLength);
    floatType *vv = v.data();
    double *wv = w.data();
    for (int i = 0; i < localLength; ++i) wv[i] = vv[i];
}

/**
 *
 */
void
PrintVector(Array<floatType> &v)
{
    using namespace std;
    const decltype(v.length()) w = 2;
    const decltype(v.length()) brk = 80;

    floatType *vd = v.data();

    for (decltype(v.length()) i = 0; i < v.length(); ++i) {
        if (0 == i % (brk / w)) cout << endl;
        cout << fixed << setw(w - 1) << setprecision(w - 1)
             << setfill('0') << vd[i] << " ";
    }
}

#if 0
/*!
  Multiply (scale) a specific vector entry by a given value.

  @param[inout] v Vector to be modified
  @param[in] index Local index of entry to scale
  @param[in] value Value to scale by
 */
inline void ScaleVectorValue(Vector & v, local_int_t index, double value) {
  assert(index>=0 && index < v.localLength);
  double * vv = v.values;
  vv[index] *= value;
  return;
}
/*!
  Fill the input vector with pseudo-random values.

  @param[in] v
 */
inline void FillRandomVector(Vector & v) {
  local_int_t localLength = v.localLength;
  double * vv = v.values;
  for (int i=0; i<localLength; ++i) vv[i] = rand() / (double)(RAND_MAX) + 1.0;
  return;
}


/*!
  Deallocates the members of the data structure of the known system matrix provided they are not 0.

  @param[in] A the known system matrix
 */
inline void DeleteVector(Vector & v) {

  delete [] v.values;
  v.localLength = 0;
  return;
}
#endif
