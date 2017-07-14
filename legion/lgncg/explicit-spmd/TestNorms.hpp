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
      @file TestNorms.hpp

      HPCG routine
 */

#pragma once

#include <cmath>

struct TestNormsData {
    double *values;  //!< Sample values.
    double mean;     //!< Mean of all sampes.
    double variance; //!< Variance of mean.
    int samples;     //!< Number of samples.
    bool pass;       //!< Pass/fail indicator.
};

/*!
    Computes the mean and standard deviation of the array of norm results.

    @param[in] testnormsData data structure with the results of norm test.

    @return Returns 0 upon success or non-zero otherwise
*/
inline int
TestNorms(
    TestNormsData & testnormsData
) {
    double mean_delta = 0.0;
    for (int i = 0; i < testnormsData.samples; ++i) {
        mean_delta += (testnormsData.values[i] - testnormsData.values[0]);
    }
    double mean = testnormsData.values[0]
                + mean_delta / (double)testnormsData.samples;
    testnormsData.mean = mean;

    // Compute variance.
    double sumdiff = 0.0;
    for (int i = 0; i < testnormsData.samples; ++i) {
        sumdiff += (testnormsData.values[i] - mean)
                 * (testnormsData.values[i] - mean);
    }
    testnormsData.variance = sumdiff / (double)testnormsData.samples;
    // Determine if variation is sufficiently small to declare success.
    testnormsData.pass = (testnormsData.variance < 1.0e-6);
    //
    return 0;
}
