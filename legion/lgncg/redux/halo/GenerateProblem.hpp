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
 @file GenerateProblem_ref.hpp

 HPCG routine
 */

#pragma once

#include "hpcg.hpp"
#include "LegionMatrices.hpp"

#include <cassert>

// TODO RM
#define HPCG_NO_MPI

/*!
    Reference version of GenerateProblem to generate the sparse matrix, right
    hand side, initial guess, and exact solution.

    @param[in]  A        The known system matrix
    @param[inout] b      The newly allocated and generated right hand side
    vector (if b!=0 on entry)
    @param[inout] x      The newly allocated solution vector with entries set to
    0.0 (if x!=0 on entry)
    @param[inout] xexact The newly allocated solution vector with entries set to
    the exact solution (if the xexact!=0 non-zero on entry)

    @see GenerateGeometry
*/
void
GenerateProblem(
    SparseMatrix &A,
    void *b,
    void *x,
    void *xexact
) {
    (void)b;
    (void)x;
    (void)xexact;
    // Make local copies of geometry information.  Use global_int_t since the
    // RHS products in the calculations below may result in global range values.
    global_int_t nx  = A.geom->nx;
    global_int_t ny  = A.geom->ny;
    global_int_t nz  = A.geom->nz;
    global_int_t npx = A.geom->npx;
    global_int_t npy = A.geom->npy;
    global_int_t npz = A.geom->npz;
    global_int_t ipx = A.geom->ipx;
    global_int_t ipy = A.geom->ipy;
    global_int_t ipz = A.geom->ipz;
    global_int_t gnx = nx*npx;
    global_int_t gny = ny*npy;
    global_int_t gnz = nz*npz;
    // This is the size of our subblock
    local_int_t localNumberOfRows = nx*ny*nz;
    // If this assert fails, it most likely means that the local_int_t is set to
    // int and should be set to.  Throw an exception of the number of rows is
    // less than zero (can happen if int overflow)long long
    assert(localNumberOfRows>0);
    // We are approximating a 27-point finite element/volume/difference 3D
    // stencil
    local_int_t numberOfNonzerosPerRow = 27;

    // Total number of grid points in mesh
    global_int_t totalNumberOfRows =
        ((global_int_t)localNumberOfRows)*((global_int_t)A.geom->size);
    // If this assert fails, it most likely means that the global_int_t is set
    // to int and should be set to long long
    assert(totalNumberOfRows>0);
    // Allocate arrays that are of length localNumberOfRows
    //char *nonzerosInRow     = A.nonzerosInRow;
#if 0
    global_int_t **mtxIndG  = new global_int_t*[localNumberOfRows];
    local_int_t  **mtxIndL  = new local_int_t*[localNumberOfRows];
    double **matrixValues   = new double*[localNumberOfRows];
    double **matrixDiagonal = new double*[localNumberOfRows];
    //
    if (b!=0)      InitializeVector(*b,      localNumberOfRows);
    if (x!=0)      InitializeVector(*x,      localNumberOfRows);
    if (xexact!=0) InitializeVector(*xexact, localNumberOfRows);
    //
    double *bv = 0;
    double *xv = 0;
    double *xexactv = 0;
    if (b!=0) bv = b->values; // Only compute exact solution if requested
    if (x!=0) xv = x->values; // Only compute exact solution if requested
    // Only compute exact solution if requested
    if (xexact!=0) xexactv = xexact->values;
    A.localToGlobalMap.resize(localNumberOfRows);

    for (local_int_t i=0; i< localNumberOfRows; ++i) {
        matrixValues[i] = 0;
        matrixDiagonal[i] = 0;
        mtxIndG[i] = 0;
        mtxIndL[i] = 0;
    }
    // Now allocate the arrays pointed to
    for (local_int_t i=0; i< localNumberOfRows; ++i) {
        mtxIndL[i] = new local_int_t[numberOfNonzerosPerRow];
        matrixValues[i] = new double[numberOfNonzerosPerRow];
        mtxIndG[i] = new global_int_t[numberOfNonzerosPerRow];
    }
    //
    local_int_t localNumberOfNonzeros = 0;
    for (local_int_t iz=0; iz<nz; iz++) {
        global_int_t giz = ipz*nz+iz;
        for (local_int_t iy=0; iy<ny; iy++) {
            global_int_t giy = ipy*ny+iy;
            for (local_int_t ix=0; ix<nx; ix++) {
                global_int_t gix = ipx*nx+ix;
                local_int_t currentLocalRow = iz*nx*ny+iy*nx+ix;
                global_int_t currentGlobalRow = giz*gnx*gny+giy*gnx+gix;
                A.globalToLocalMap[currentGlobalRow] = currentLocalRow;
                A.localToGlobalMap[currentLocalRow] = currentGlobalRow;
                char numberOfNonzerosInRow = 0;
                // Pointer to current value in current row
                double *currentValuePointer = matrixValues[currentLocalRow];
                // Pointer to current index in current row
                global_int_t * currentIndexPointerG = mtxIndG[currentLocalRow];
                for (int sz=-1; sz<=1; sz++) {
                    if (giz+sz>-1 && giz+sz<gnz) {
                        for (int sy=-1; sy<=1; sy++) {
                            if (giy+sy>-1 && giy+sy<gny) {
                                for (int sx=-1; sx<=1; sx++) {
                                    if (gix+sx>-1 && gix+sx<gnx) {
                                        global_int_t curcol =
                                            currentGlobalRow+sz*gnx*gny+sy*gnx+sx;
                                        if (curcol==currentGlobalRow) {
                                            matrixDiagonal[currentLocalRow] =
                                                currentValuePointer;
                                            *currentValuePointer++ = 26.0;
                                        } else {
                                            *currentValuePointer++ = -1.0;
                                        }
                                        *currentIndexPointerG++ = curcol;
                                        numberOfNonzerosInRow++;
                                    } // end x bounds test
                                } // end sx loop
                            } // end y bounds test
                        } // end sy loop
                    } // end z bounds test
                } // end sz loop
                nonzerosInRow[currentLocalRow] = numberOfNonzerosInRow;
                localNumberOfNonzeros += numberOfNonzerosInRow;
                if (b!=0) {
                    bv[currentLocalRow] =
                        26.0 - ((double)(numberOfNonzerosInRow-1));
                }
                if (x!=0) {
                    xv[currentLocalRow] = 0.0;
                }
                if (xexact!=0) {
                    xexactv[currentLocalRow] = 1.0;
                }
            } // end ix loop
        } // end iy loop
    } // end iz loop
#ifdef HPCG_DETAILED_DEBUG
    HPCG_fout     << "Process " << A.geom->rank << " of " << A.geom->size <<"
        has " << localNumberOfRows    << " rows."     << endl << "Process " <<
        A.geom->rank << " of " << A.geom->size <<" has " <<
        localNumberOfNonzeros<< " nonzeros." <<endl;
#endif

    global_int_t totalNumberOfNonzeros = 0;
#ifndef HPCG_NO_MPI
    // Use MPI's reduce function to sum all nonzeros
    // convert to 64 bit for MPI call
    long long lnnz = localNumberOfNonzeros, gnnz = 0;
    MPI_Allreduce(&lnnz, &gnnz, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
    totalNumberOfNonzeros = gnnz; // Copy back
#else
  totalNumberOfNonzeros = localNumberOfNonzeros;
#endif
    // If this assert fails, it most likely means that the global_int_t is set
    // to int and should be set to long long This assert is usually the first to
    // fail as problem size increases beyond the 32-bit integer range.  Throw an
    // exception of the number of nonzeros is less than zero (can happen if int
    // overflow)
    assert(totalNumberOfNonzeros>0);

    A.title = 0;
    A.totalNumberOfRows = totalNumberOfRows;
    A.totalNumberOfNonzeros = totalNumberOfNonzeros;
    A.localNumberOfRows = localNumberOfRows;
    A.localNumberOfColumns = localNumberOfRows;
    A.localNumberOfNonzeros = localNumberOfNonzeros;
    A.nonzerosInRow = nonzerosInRow;
    A.mtxIndG = mtxIndG;
    A.mtxIndL = mtxIndL;
    A.matrixValues = matrixValues;
    A.matrixDiagonal = matrixDiagonal;

    return;
#endif
}
