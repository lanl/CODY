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
    @file TestSymmetry.hpp

    HPCG routine
 */

#pragma once

#include "hpcg.hpp"
#include "LegionStuff.hpp"
#include "LegionArrays.hpp"
#include "LegionMatrices.hpp"
#include "ComputeSPMV.hpp"
#include "ComputeMG.hpp"
#include "ComputeDotProduct.hpp"
#include "ComputeResidual.hpp"
#include "Geometry.hpp"
#include "VectorOps.hpp"

#include <fstream>
#include <iostream>
#include <cfloat>
#include <vector>
#include <cmath>

struct TestSymmetryData {
    //!< Departure from symmetry for the SPMV kernel.
    double depsym_spmv;
    //!< Departure from symmetry for the SYMGS kernel.
    double depsym_mg;
    //!< Number of failures in the symmetry tests.
    int count_fail;
};

/*!
    Tests symmetry-preserving properties of the sparse matrix vector multiply
    and symmetric Gauss-Siedel routines.

    @param[in]    geom   The description of the problem's geometry.

    @param[in]    A      The known system matrix.

    @param[in]    b      The known right hand side vector.

    @param[in]    xexact The exact solution vector.

    @param[inout] testSymmetryData The data structure with the results of the
                  CG symmetry test including pass/fail information.

    @return returns 0 upon success and non-zero otherwise.

    @see ComputeDotProduct
    @see ComputeSPMV
    @see ComputeMG
*/
inline int
TestSymmetry(
    SparseMatrix &A,
    Array<floatType> &b,
    Array<floatType> &xexact,
    TestSymmetryData &testSymmetryData,
    Context ctx,
    HighLevelRuntime *lrt
) {
    using namespace std;

    const auto *const Asclrs = A.sclrs->data();
    const Geometry *const Ageom = A.geom->data();
    //
    local_int_t nrow = Asclrs->localNumberOfRows;
    local_int_t ncol = Asclrs->localNumberOfColumns;

    LogicalArray<floatType> x_ncoll, y_ncoll, z_ncoll;
    //
    x_ncoll.allocate("x_ncol", ncol, ctx, lrt);
    y_ncoll.allocate("y_ncol", ncol, ctx, lrt);
    z_ncoll.allocate("z_ncol", ncol, ctx, lrt);
    //
    Array<floatType> x_ncol(
        x_ncoll.mapRegion(RW_E, ctx, lrt), ctx, lrt
    );
    Array<floatType> y_ncol(
        y_ncoll.mapRegion(RW_E, ctx, lrt), ctx, lrt
    );
    Array<floatType> z_ncol(
        z_ncoll.mapRegion(RW_E, ctx, lrt), ctx, lrt
    );
    //
    Item< DynColl<floatType> > &dcFT = *A.dcAllRedSumFT;
    // Needed for dot-product call, otherwise unused.
    double t4 = 0.0;
    testSymmetryData.count_fail = 0;

    ////////////////////////////////////////////////////////////////////////////
    // Test symmetry of matrix.
    ////////////////////////////////////////////////////////////////////////////
    // First load vectors with random values
    FillRandomVector(x_ncol, ctx, lrt);
    FillRandomVector(y_ncol, ctx, lrt);

    double xNorm2, yNorm2;
    double ANorm = 2 * 26.0;

    // Next, compute x'*A*y
    ComputeDotProduct(nrow, y_ncol, y_ncol, yNorm2, t4, dcFT, ctx, lrt);
    //
    // z_nrow = A*y_overlap
    int ierr = ComputeSPMV(A, y_ncol, z_ncol, ctx, lrt);
    if (ierr) cerr << "Error in call to SpMV: " << ierr << ".\n" << endl;
    // x'*A*y
    double xtAy = 0.0;
    ierr = ComputeDotProduct(nrow, x_ncol, z_ncol, xtAy, t4, dcFT, ctx, lrt);
    if (ierr) cerr << "Error in call to dot: " << ierr << ".\n" << endl;
    // Next, compute y'*A*x
    ComputeDotProduct(nrow, x_ncol, x_ncol, xNorm2, t4, dcFT, ctx, lrt); 
    // b_computed = A*x_overlap
    ierr = ComputeSPMV(A, x_ncol, z_ncol, ctx, lrt);
    if (ierr) cerr << "Error in call to SpMV: " << ierr << ".\n" << endl;
    double ytAx = 0.0;
    // y'*A*x
    ierr = ComputeDotProduct(nrow, y_ncol, z_ncol, ytAx, t4, dcFT, ctx, lrt); 
    if (ierr) cerr << "Error in call to dot: " << ierr << ".\n" << endl;
    testSymmetryData.depsym_spmv = std::fabs((long double)(xtAy - ytAx))
                                 / ((xNorm2 * ANorm * yNorm2
                                     + yNorm2 * ANorm * xNorm2)
                                 *  (DBL_EPSILON));
    if (testSymmetryData.depsym_spmv > 1.0) {
        // If the difference is > 1, count it wrong.
        ++testSymmetryData.count_fail;
    }
    if (Ageom->rank == 0) {
        cout << "Departure from symmetry (scaled) for SpMV abs(x'*A*y - y'*A*x) = "
             << testSymmetryData.depsym_spmv << endl;
    }

    ////////////////////////////////////////////////////////////////////////////
    // Test symmetry of symmetric Gauss-Seidel.
    ////////////////////////////////////////////////////////////////////////////
    // Compute x'*Minv*y
    ierr = ComputeMG(A, y_ncol, z_ncol); // z_ncol = Minv*y_ncol
    if (ierr) HPCG_fout << "Error in call to MG: " << ierr << ".\n" << endl;
    double xtMinvy = 0.0;
    ierr = ComputeDotProduct(nrow, x_ncol, z_ncol, xtMinvy, t4, A.isDotProductOptimized); // x'*Minv*y
    if (ierr) HPCG_fout << "Error in call to dot: " << ierr << ".\n" << endl;
#if 0
    // Next, compute z'*Minv*x
    ierr = ComputeMG(A, x_ncol, z_ncol); // z_ncol = Minv*x_ncol
    if (ierr) HPCG_fout << "Error in call to MG: " << ierr << ".\n" << endl;
    double ytMinvx = 0.0;
    ierr = ComputeDotProduct(nrow, y_ncol, z_ncol, ytMinvx, t4, A.isDotProductOptimized); // y'*Minv*x
    if (ierr) HPCG_fout << "Error in call to dot: " << ierr << ".\n" << endl;

    testSymmetryData.depsym_mg = std::fabs((long double) (xtMinvy - ytMinvx))/((xNorm2*ANorm*yNorm2 + yNorm2*ANorm*xNorm2) * (DBL_EPSILON));
    if (testSymmetryData.depsym_mg > 1.0) ++testSymmetryData.count_fail;  // If the difference is > 1, count it wrong
    if (A.geom->rank==0) HPCG_fout << "Departure from symmetry (scaled) for MG abs(x'*Minv*y - y'*Minv*x) = " << testSymmetryData.depsym_mg << endl;

    CopyVector(xexact, x_ncol); // Copy exact answer into overlap vector

    int numberOfCalls = 2;
    double residual = 0.0;
    for (int i=0; i< numberOfCalls; ++i) {
        ierr = ComputeSPMV(A, x_ncol, z_ncol); // b_computed = A*x_overlap
        if (ierr) HPCG_fout << "Error in call to SpMV: " << ierr << ".\n" << endl;
        if ((ierr = ComputeResidual(A.localNumberOfRows, b, z_ncol, residual)))
            HPCG_fout << "Error in call to compute_residual: " << ierr << ".\n" << endl;
        if (A.geom->rank==0) HPCG_fout << "SpMV call [" << i << "] Residual [" << residual << "]" << endl;
    }
    DeleteVector(x_ncol);
    DeleteVector(y_ncol);
    DeleteVector(z_ncol);
#endif
    return 0;
}

