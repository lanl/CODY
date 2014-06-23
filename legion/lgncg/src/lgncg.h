/**
 * Copyright (c) 2014      Los Alamos National Security, LLC
 *                         All rights reserved.
 */

#ifndef LGNCG_H_INCLUDED
#define LGNCG_H_INCLUDED

#include "vector.h"
#include "sparsemat.h"

#include "legion.h"

#include <stdint.h>

struct SparseMatrix;

namespace lgncg {
void
init(void);

void
finalize(void);

void
cgSolv(SparseMatrix &A,
       Vector &b,
       double tolerance,
       int64_t maxIters,
       Vector &x,
       int64_t nParts,
       bool doPreconditioning,
       LegionRuntime::HighLevel::Context ctx,
       LegionRuntime::HighLevel::HighLevelRuntime *lrt);

} // end lgncg namespace

#endif
