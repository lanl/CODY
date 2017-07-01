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

#include "ReduceSum.hpp"

const floatType FloatReduceSumAccumulate::identity = 0.0;

template<>
void
FloatReduceSumAccumulate::apply<true>(LHS &lhs, RHS rhs) {
    lhs += rhs;
}

template<>
void
FloatReduceSumAccumulate::apply<false>(LHS &lhs, RHS rhs) {
    int64_t *target = (int64_t *)&lhs;
    union { int64_t as_int; floatType as_T; } oldval, newval;
    do {
        oldval.as_int = *target;
        newval.as_T = oldval.as_T + rhs;
    } while (!__sync_bool_compare_and_swap(target, oldval.as_int, newval.as_int));
}

template<>
void
FloatReduceSumAccumulate::fold<true>(RHS &rhs1, RHS rhs2) {
    rhs1 += rhs2;
}

template<>
void
FloatReduceSumAccumulate::fold<false>(RHS &rhs1, RHS rhs2) {
    int64_t *target = (int64_t *)&rhs1;
    union { int64_t as_int; floatType as_T; } oldval, newval;
    do {
        oldval.as_int = *target;
        newval.as_T = oldval.as_T + rhs2;
    } while (!__sync_bool_compare_and_swap(target, oldval.as_int, newval.as_int));
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
const global_int_t IntReduceSumAccumulate::identity = 0;

template<>
void
IntReduceSumAccumulate::apply<true>(LHS &lhs, RHS rhs) {
    lhs += rhs;
}

template<>
void
IntReduceSumAccumulate::apply<false>(LHS &lhs, RHS rhs) {
    int64_t *target = (int64_t *)&lhs;
    union { int64_t as_int; global_int_t as_T; } oldval, newval;
    do {
        oldval.as_int = *target;
        newval.as_T = oldval.as_T + rhs;
    } while (!__sync_bool_compare_and_swap(target, oldval.as_int, newval.as_int));
}

template<>
void
IntReduceSumAccumulate::fold<true>(RHS &rhs1, RHS rhs2) {
    rhs1 += rhs2;
}

template<>
void
IntReduceSumAccumulate::fold<false>(RHS &rhs1, RHS rhs2) {
    int64_t *target = (int64_t *)&rhs1;
    union { int64_t as_int; global_int_t as_T; } oldval, newval;
    do {
        oldval.as_int = *target;
        newval.as_T = oldval.as_T + rhs2;
    } while (!__sync_bool_compare_and_swap(target, oldval.as_int, newval.as_int));
}
