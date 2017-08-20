/**
 * Copyright (c) 2016-2017 Los Alamos National Security, LLC
 *                         All rights reserved.
 *
 * Copyright (c) 2016      Stanford University
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

#pragma once

////////////////////////////////////////////////////////////////////////////////
// Task IDs
////////////////////////////////////////////////////////////////////////////////
enum {
    MAIN_TID = 0,
    GEN_PROB_TID,
    START_BENCHMARK_TID,
    REGION_TO_REGION_COPY_TID,
    DYN_COLL_TASK_CONTRIB_GIT_TID,
    DYN_COLL_TASK_CONTRIB_FT_TID,
    FLOAT_REDUCE_SUM_TID,
    FLOAT_REDUCE_MIN_TID,
    FLOAT_REDUCE_MAX_TID,
    INT_REDUCE_SUM_TID,
    COPY_VECTOR_TID,
    ZERO_VECTOR_TID,
    FILLRAND_VECTOR_TID,
    WAXPBY_TID,
    SPMV_TID,
    DDOT_TID,
    SYMGS_TID,
    PROLONGATION_TID,
    RESTRICTION_TID,
    FUTURE_MATH_TID,
    COMPUTE_RESIDUAL_TID,
    EXCHANGE_HALO_TID
};
