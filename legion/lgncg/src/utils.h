/**
 * Copyright (c) 2014      Los Alamos National Security, LLC
 *                         All rights reserved.
 */

#ifndef LGNCG_UTILS_H_INCLUDED
#define LGNCG_UTILS_H_INCLUDED

#include "legion.h"

namespace lgncg {

/**
 * convenience routine to get a task's ID
 */
static inline int
getTaskID(const LegionRuntime::HighLevel::Task *task)
{
    return task->index_point.point_data[0];
}

/**
 * courtesy of some other legion code.
 */
template <unsigned DIM, typename T>
static inline bool
offsetsAreDense(const Rect<DIM> &bounds,
                const LegionRuntime::Accessor::ByteOffset *offset)
{
    off_t exp_offset = sizeof(T);
    for (unsigned i = 0; i < DIM; i++) {
        bool found = false;
        for (unsigned j = 0; j < DIM; j++)
            if (offset[j].offset == exp_offset) {
                found = true;
                exp_offset *= (bounds.hi[j] - bounds.lo[j] + 1);
                break;
            }
        if (!found) return false;
    }
    return true;
}

} // end lgncg namespace

#endif
