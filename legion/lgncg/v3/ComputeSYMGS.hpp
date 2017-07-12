
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
      @file ComputeSYMGS.hpp

      HPCG routine
 */

#include "LegionStuff.hpp"
#include "LegionArrays.hpp"
#include "LegionMatrices.hpp"
#include "ExchangeHalo.hpp"

#include <cassert>

/*!
    Computes one step of symmetric Gauss-Seidel:

    Assumption about the structure of matrix A:
    - Each row 'i' of the matrix has nonzero diagonal value whose address is
      matrixDiagonal[i]
    - Entries in row 'i' are ordered such that:
         - Lower triangular terms are stored before the diagonal element.
         - Upper triangular terms are stored after the diagonal element.
         - No other assumptions are made about entry ordering.

    Symmetric Gauss-Seidel notes:
    - We use the input vector x as the RHS and start with an initial guess for y
    of all zeros.
    - We perform one forward sweep.  x should be initially zero on the first GS
      sweep, but we do not attempt to exploit this fact.
    - We then perform one back sweep.
    - For simplicity we include the diagonal contribution in the for-j loop,
      then correct the sum after

    @param[in] A the known system matrix.

    @param[in] r the input vector.

    @param[inout] x On entry, x should contain relevant values, on exit x
                  contains the result of one symmetric GS sweep with r as the
                  RHS.

    @warning Early versions of this kernel (Version 1.1 and earlier) had the r
    and x arguments in reverse order, and out of sync with other kernels.

    @return returns 0 upon success and non-zero otherwise.

    @see ComputeSYMGS
*/
inline int
ComputeSYMGS(
    SparseMatrix &A,
    Array<floatType> &r,
    Array<floatType> &x,
    Context ctx,
    Runtime *lrt
) {
    const SparseMatrixScalars *const Asclrs = A.sclrs->data();
    // Make sure x contain space for halo values.
    assert(x.length() == size_t(Asclrs->localNumberOfColumns));
    //
    ExchangeHalo(A, x, ctx, lrt);
    //
    const local_int_t nrow = Asclrs->localNumberOfRows;
    const local_int_t nnpr = A.geom->data()->stencilSize;
    const floatType *const matrixDiagonal = A.matrixDiagonal->data();
    const floatType *const rv = r.data();
    floatType *const xv = x.data();
    // Interpreted as 2D array
    Array2D<floatType> matrixValues(
        nrow, nnpr, A.matrixValues->data()
    );
    // Interpreted as 2D array
    Array2D<local_int_t> mtxIndL(
        nrow, nnpr, A.mtxIndL->data()
    );
    const char *const nonzerosInRow = A.nonzerosInRow->data();
    //
    for (local_int_t i = 0; i < nrow; i++) {
        const floatType *const currentValues = matrixValues(i);
        const local_int_t *const currentColIndices = mtxIndL(i);
        const int currentNumberOfNonzeros = nonzerosInRow[i];
        const floatType currentDiagonal = matrixDiagonal[i];
        double sum = rv[i]; // RHS value
        //
        for (int j = 0; j < currentNumberOfNonzeros; j++) {
            const local_int_t curCol = currentColIndices[j];
            sum -= currentValues[j] * xv[curCol];
        }
        // Remove diagonal contribution from previous loop.
        sum += xv[i] * currentDiagonal;
        //
        xv[i] = sum / currentDiagonal;
    }
    // Now the back sweep.
    for (local_int_t i = nrow - 1; i >= 0; i--) {
        const floatType *const currentValues = matrixValues(i);
        const local_int_t *const currentColIndices = mtxIndL(i);
        const int currentNumberOfNonzeros = nonzerosInRow[i];
        const floatType currentDiagonal = matrixDiagonal[i];
        double sum = rv[i]; // RHS value
        //
        for (int j = 0; j < currentNumberOfNonzeros; j++) {
            const local_int_t curCol = currentColIndices[j];
            sum -= currentValues[j] * xv[curCol];
        }
        // Remove diagonal contribution from previous loop.
        sum += xv[i] * currentDiagonal;
        xv[i] = sum / currentDiagonal;
    }
    //
    return 0;
}
