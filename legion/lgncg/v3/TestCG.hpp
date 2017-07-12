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
    @file TestCG.hpp

    HPCG data structure
 */

#pragma once

#include "LegionStuff.hpp"
#include "LegionArrays.hpp"
#include "LegionMatrices.hpp"
#include "LegionCGData.hpp"
#include "VectorOps.hpp"

#include "hpcg.hpp"

struct TestCGData {
    //!< Number of successful tests.
    int count_pass;
    //!< number of unsuccessful tests.
    int count_fail;
    //!< Expected number of test CG iterations without preconditioning with
    //   diagonally dominant matrix (~12).
    int expected_niters_no_prec;
    //!< Expected number of test CG iterations with preconditioning and with
    //   diagonally dominant matrix (~1-2).
    int expected_niters_prec;
    //!< Maximum number of test CG iterations without preconditioner.
    int niters_max_no_prec;
    //!< Maximum number of test CG iterations without preconditioner.
    int niters_max_prec;
    //!< Residual norm achieved during test CG iterations.
    double normr;
};

/*!
    Test the correctness of the Preconditined CG implementation by using a
    system matrix with a dominant diagonal.

    @param[in]    geom The description of the problem's geometry.

    @param[in]    A    The known system matrix.

    @param[in]    data the data structure with all necessary CG vectors preallocated.

    @param[in]    b    The known right hand side vector.

    @param[inout] x    On entry: the initial guess; on exit: the new approximate
                       solution.

    @param[out]   testcg_data the data structure with the results of the test
                  including pass/fail information.

    @return Returns zero on success and a non-zero value otherwise.

    @see CG()
 */
inline int
TestCG(
    SparseMatrix &A,
    CGData &data,
    Array<floatType> &b,
    Array<floatType> &x,
    TestCGData &testcg_data,
    Context ctx,
    Runtime *lrt
) {
    using namespace std;
    // Use this array for collecting timing information.
    std::vector<double> times(8, 0.0);
    //
    const auto *const Asclrs = A.sclrs->data();
    const local_int_t nrow = Asclrs->localNumberOfRows; 
    // Temporary storage for holding original diagonal and RHS.
    LogicalArray<floatType> origDiagAl, exaggeratedDiagAl, origBl;
    origDiagAl.allocate(       "origDiagA",        nrow, ctx, lrt);
    exaggeratedDiagAl.allocate("exaggeratedDiagA", nrow, ctx, lrt);
    origBl.allocate(           "origB",            nrow, ctx, lrt);
    //
    Array<floatType> origDiagA(
        origDiagAl.mapRegion(RW_E, ctx, lrt),
        ctx,
        lrt
    );
    Array<floatType> exaggeratedDiagA(
        exaggeratedDiagAl.mapRegion(RW_E, ctx, lrt),
        ctx,
        lrt
    );
    Array<floatType> origB(
        origBl.mapRegion(RW_E, ctx, lrt),
        ctx,
        lrt
    );
    //
    CopyMatrixDiagonal(A, origDiagA, ctx, lrt);
    CopyVector(origDiagA, exaggeratedDiagA, ctx, lrt);
    CopyVector(b, origB, ctx, lrt);

    // Modify the matrix diagonal to greatly exaggerate diagonal values.  CG
    // should converge in about 10 iterations for this problem, regardless of
    // problem size.
    const global_int_t *const AlocalToGlobalMap = A.localToGlobalMap->data();
    for (local_int_t i = 0; i < nrow; ++i) {
        global_int_t globalRowID = AlocalToGlobalMap[i];
        if (globalRowID < 9) {
            double scale = (globalRowID + 2) * 1.0e6;
            ScaleVectorValue(exaggeratedDiagA, i, scale, ctx, lrt);
            ScaleVectorValue(b, i, scale, ctx, lrt);
        }
        else {
            ScaleVectorValue(exaggeratedDiagA, i, 1.0e6, ctx, lrt);
            ScaleVectorValue(b, i, 1.0e6, ctx, lrt);
        }
    }
    ReplaceMatrixDiagonal(A, exaggeratedDiagA, ctx, lrt);

    int niters = 0;
    double normr = 0.0;
    double normr0 = 0.0;
    int maxIters = 50;
    int numberOfCgCalls = 2;
    // Set tolerance to reasonable value for grossly scaled diagonal terms.
    double tolerance = 1.0e-12;
    // For the unpreconditioned CG call, we should take about 10 iterations,
    // permit 12.
    testcg_data.expected_niters_no_prec = 12;
    // For the preconditioned case, we should take about 1 iteration, permit 2.
    testcg_data.expected_niters_prec = 2;
    testcg_data.niters_max_no_prec = 0;
    testcg_data.niters_max_prec = 0;
    //
    // This loop tests both unpreconditioned and preconditioned runs.
    for (int k = 0; k < 2; ++k) {
        int expected_niters = testcg_data.expected_niters_no_prec;
        if (k == 1) expected_niters = testcg_data.expected_niters_prec;
        for (int i = 0; i < numberOfCgCalls; ++i) {
            ZeroVector(x, ctx, lrt); // Zero out x
            int ierr = CG(A,
                          data,
                          b,
                          x,
                          maxIters,
                          tolerance,
                          niters,
                          normr,
                          normr0,
                          &times[0],
                          k == 1,
                          ctx,
                          lrt
                       );
            if (ierr) cerr << "Error in call to CG: " << ierr << ".\n" << endl;
            if (niters <= expected_niters) {
                ++testcg_data.count_pass;
            }
            else {
                ++testcg_data.count_fail;
            }
            if (k == 0 && niters > testcg_data.niters_max_no_prec) {
                // Keep track of largest iter count.
                testcg_data.niters_max_no_prec = niters;
            }
            if (k == 1 && niters > testcg_data.niters_max_prec) {
                // Same for preconditioned run.
                testcg_data.niters_max_prec = niters;
            }
            if (A.geom->data()->rank == 0) {
                cout << "Call [" << i << "] Number of Iterations [" << niters
                     << "] Scaled Residual [" << normr/normr0 << "]" << endl;
                if (niters > expected_niters) {
                    cout << " Expected " << expected_niters << " iterations.  Performed "
                         << niters << "." << endl;
                }
            }
        }
    }
    // Restore matrix diagonal and RHS.
    ReplaceMatrixDiagonal(A, origDiagA, ctx, lrt);
    //
    CopyVector(origB, b, ctx, lrt);
    // Delete vectors.
    // TODO make sure you do this elsewhere.
    origDiagAl.deallocate(ctx, lrt);
    origDiagAl.unmapRegion(ctx, lrt);
    exaggeratedDiagAl.deallocate(ctx, lrt);
    exaggeratedDiagAl.unmapRegion(ctx, lrt);
    origBl.deallocate(ctx, lrt);
    origBl.unmapRegion(ctx, lrt);
    //
    testcg_data.normr = normr;
    //
    return 0;
}
