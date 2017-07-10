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
 @file CG.hpp

 HPCG routine
 */

#include <fstream>
#include <cmath>
#include <unistd.h>

#include "hpcg.hpp"
#include "mytimer.hpp"

#include "LegionArrays.hpp"
#include "LegionMatrices.hpp"
#include "LegionCGData.hpp"
#include "VectorOps.hpp"

#include "ComputeSPMV.hpp"
#include "ComputeWAXPBY.hpp"
#include "ComputeDotProduct.hpp"

#if 0
#include "ComputeMG_ref.hpp"
#endif

// Use TICK and TOCK to time a code section in MATLAB-like fashion
#define TICK()  t0 = mytimer() //!< record current time in 't0'
#define TOCK(t) t += mytimer() - t0 //!< store time difference in 't' using time in 't0'

/*!
  Reference routine to compute an approximate solution to Ax = b

  @param[inout] A    The known system matrix
  @param[inout] data The data structure with all necessary CG vectors
                     preallocated
  @param[in]    b    The known right hand side vector
  @param[inout] x    On entry: the initial guess; on exit: the new approximate
                     solution
  @param[in]    max_iter  The maximum number of iterations to perform, even if
                tolerance is not met.
  @param[in]    tolerance The stopping criterion to assert convergence: if norm
                of residual is <= to tolerance.
  @param[out]   niters    The number of iterations actually performed.
  @param[out]   normr     The 2-norm of the residual vector after the last
                          iteration.
  @param[out]   normr0    The 2-norm of the residual vector before the first
                          iteration.
  @param[out]   times     The 7-element vector of the timing information
                          accumulated during all of the iterations.
  @param[in]    doPreconditioning The flag to indicate whether the
                preconditioner should be invoked at each iteration.

  @return Returns zero on success and a non-zero value otherwise.

  @see CG()
*/
inline int
CG(
    SparseMatrix &A,
    CGData &data,
    const Array<floatType> &b,
    Array<floatType> &x,
    const int max_iter,
    const double tolerance,
    int &niters,
    double &normr,
    double &normr0,
    double *times,
    bool doPreconditioning,
    LegionRuntime::HighLevel::Context ctx,
    LegionRuntime::HighLevel::Runtime *lrt
) {
    using namespace std;

    double t_begin = mytimer();  // Start timing right away
    //
    const int print_freq = 50;
    const int rank = A.geom->data()->rank;
    const local_int_t nrow = A.sclrs->data()->localNumberOfRows;
    //
    double rtz = 0.0, oldrtz = 0.0, alpha = 0.0, beta = 0.0, pAp = 0.0;
    double t0 = 0.0, t1 = 0.0, t2 = 0.0, t3 = 0.0, t4 = 0.0, t5 = 0.0, t6 = 0.0;

    normr = 0.0;

    Array<floatType> &r  = *(data.r); // Residual vector.
    Array<floatType> &z  = *(data.z); // Preconditioned residual vector.
    Array<floatType> &p  = *(data.p); // Direction vector (ncol >= nrow).
    Array<floatType> &Ap = *(data.Ap);// Holds result from A * p.

    Item< DynColl<floatType> > &dcFT = *A.dcAllRedSumFT;

    if (!doPreconditioning && rank == 0) {
        std::cout << "WARNING: PERFORMING UNPRECONDITIONED ITERATIONS" << std::endl;
    }
    // p is of length ncols, copy x to p for sparse MV operation
    CopyVector(x, p, ctx, lrt);
    //
    TICK(); ComputeSPMV(A, p, Ap, ctx, lrt);  TOCK(t3); // Ap = A*p

    // r = b - Ax (x stored in p)
    TICK(); ComputeWAXPBY(nrow, 1.0, b, -1.0, Ap, r); TOCK(t2);

    TICK();
    ComputeDotProduct(nrow, r, r, normr, t4, dcFT, ctx, lrt);
    TOCK(t1);

    normr = sqrt(normr);

    if (rank == 0) std::cout << "Initial Residual = "<< normr << std::endl;

    // Record initial residual for convergence testing
    normr0 = normr;

  // Start iterations

    for (int k = 1; k <= max_iter && normr / normr0 > tolerance; k++ ) {
        TICK();
        if (doPreconditioning) {
#if 0
            ComputeMG_ref(A, r, z); // Apply preconditioner
#endif
        }
        else {
            // copy r to z (no preconditioning)
            ComputeWAXPBY(nrow, 1.0, r, 0.0, r, z);
        }
        TOCK(t5); // Preconditioner apply time

        if (k == 1) {
            TICK(); CopyVector(z, p, ctx, lrt); TOCK(t2); // Copy Mr to p
            // rtz = r'*z
            TICK();
            ComputeDotProduct(nrow, r, z, rtz, t4, dcFT, ctx, lrt);
            TOCK(t1);
        }
        else {
            oldrtz = rtz;
            //
            TICK(); // rtz = r'*z
            ComputeDotProduct(nrow, r, z, rtz, t4, dcFT, ctx, lrt);
            TOCK(t1);
            beta = rtz / oldrtz;
            //
            TICK(); // p = beta*p + z
            ComputeWAXPBY(nrow, 1.0, z, beta, p, p);
            TOCK(t2);
        }

        TICK(); // Ap = A*p
        ComputeSPMV(A, p, Ap, ctx, lrt);
        TOCK(t3);
        //
        TICK(); // alpha = p'*Ap
        ComputeDotProduct(nrow, p, Ap, pAp, t4, dcFT, ctx, lrt);
        TOCK(t1);
        //
        alpha = rtz / pAp;
        //
        TICK(); // x = x + alpha*p
        ComputeWAXPBY(nrow, 1.0, x, alpha, p, x);
        // r = r - alpha*Ap
        ComputeWAXPBY(nrow, 1.0, r, -alpha, Ap, r);
        TOCK(t2);
        //
        TICK();
        ComputeDotProduct(nrow, r, r, normr, t4, dcFT, ctx, lrt);
        TOCK(t1);
        normr = sqrt(normr);
        //
        if (rank == 0 && ( k % print_freq == 0 || k == max_iter)) {
            cout << "Iteration = "<< k << "   Scaled Residual = "
                << normr / normr0 << std::endl;
        }
        niters = k;
  }
    // Store times
    times[1] += t1; // dot product time
    times[2] += t2; // WAXPBY time
    times[3] += t3; // SPMV time
    times[4] += t4; // AllReduce time
    times[5] += t5; // preconditioner apply time
    times[6] += t6; // exchange halo time
    times[0] += mytimer() - t_begin;  // Total time. All done...

    return 0;
}
