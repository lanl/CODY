/**
 * Copyright (c) 2014      Los Alamos National Security, LLC
 *                         All rights reserved.
 */

#ifndef LGNCG_TIDS_H_INCLUDED
#define LGNCG_TIDS_H_INCLUDED

namespace lgncg {

enum {
    // FIXME - central registration thing again...
    LGNCG_VECCP_TID = 32,
    LGNCG_VEC_ZERO_TID,
    LGNCG_SPMV_TID,
    LGNCG_WAXPBY_TID,
    LGNCG_DOTPROD_TID,
    LGNCG_SYMGS_TID,
    LGNCG_RESTRICTION_TID,
    LGNCG_PROLONGATION_TID
};

} // end lgncg namespace

#endif
