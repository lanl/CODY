#
# Copyright (c) 2014      Los Alamos National Security, LLC
#                         All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 
# 1. Redistributions of source code must retain the above copyright notice, this
# list of conditions and the following disclaimer.
# 
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# LA-CC 10-123
#

# cg implementation in octave (maybe matlab)
# in octave run with: run octcg.m

1;

function x = cg(A, b, tol, maxIts, x)
    r = b - A * x;
    p = r;
    rsold = r' * r;
    for i = 1:maxIts
        Ap = A * p;
        alpha = rsold / (p' * Ap)
        x = x + alpha * p;
        r = r - alpha * Ap;
        rsnew = r' * r;
        if sqrt(rsnew) < tol
            break
        end
        p = r + (rsnew / rsold * p);
        rsold = rsnew
    end
endfunction

A = [4 1; 1 3]
x = [1 ; 1]
b = [1 ; 2]

cganswer = cg(A, b, 0.1, 128, x)

A = [4 1; 1 3];
x = [1 ; 1];
b = [1 ; 2];

pcranswer = pcr(A, b, 0.001, 128)

# cganswer and pcranswer should be about the same
