#!/usr/bin/env rdmd

/**
 * Copyright (c) 2014-2015 Los Alamos National Security, LLC
 *                         All rights reserved.
 *
 * This software was produced under U.S. Government contract DE-AC52-06NA25396
 * for Los Alamos National Laboratory (LANL), which is operated by Los Alamos
 * National Security, LLC for the U.S. Department of Energy. The U.S. Government
 * has rights to use, reproduce, and distribute this software.  NEITHER THE
 * GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS
 * OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If
 * software is modified to produce derivative works, such modified software
 * should be clearly marked, so as not to confuse it with the version available
 * from LANL.
 *
 * Additionally, redistribution and use in source and binary forms, with or
 * without modification, are permitted provided that the following conditions
 * are met:
 *
 * . Redistributions of source code must retain the above copyright notice, this
 *   list of conditions and the following disclaimer.
 *
 * . Redistributions in binary form must reproduce the above copyright notice,
 *   this list of conditions and the following disclaimer in the documentation
 *   and/or other materials provided with the distribution.
 *
 * . Neither the name of Los Alamos National Security, LLC, Los Alamos National
 *   Laboratory, LANL, the U.S. Government, nor the names of its contributors
 *   may be used to endorse or promote products derived from this software
 *   without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT
 * NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL
 * SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 * OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
 * OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * LA-CC 10-123
 */

// build optimized version:
// dmd -O -release -inline -noboundscheck ./heattx.d
// run unoptimized: ./heattx.d

import std.stdio;
import std.string;

immutable string APP_NAME = "heattx-d";
immutable string APP_VER = "0.1";
immutable ulong T_MAX = 1024;
immutable ulong N = 512;
immutable double THERM_COND = 0.6;
immutable double K = 0.4;

class Mesh {
    ulong nx, ny;
    double[][] cells;

    this(ulong n) {
        ny = nx = n;
        cells = new double[][](nx, ny);
        foreach (i ; 0 .. ny) {
            foreach (j ; 0 .. nx) {
                cells[i][j] = 0.0;
            }
        }
    }
public:
    void setInitialConds() {
        auto x0 = nx / 2;
        auto y0 = ny / 2;
        auto x  = nx / 4, y = 0;
        long radius_err = 1 - x;
        while (x >= y) {
            cells[ x + x0][ y + y0] = K * .50;
            cells[ y + x0][ x + y0] = K * .60;
            cells[-x + x0][ y + y0] = K * .70;
            cells[-y + x0][ x + y0] = K * .80;
            cells[-x + x0][-y + y0] = K * .70;
            cells[-y + x0][-x + y0] = K * .60;
            cells[ x + x0][-y + y0] = K * .50;
            cells[ y + x0][-x + y0] = K;
            y++;
            if (radius_err < 0) radius_err += 2 * y + 1;
            else {
                --x;
                radius_err += 2 * (y - x + 1);
            }
        }
    }
};

class Simulation {
    Mesh oldMesh, newMesh;
    double c, deltaS, deltaT;
    ulong maxT;

    this(ulong n, double c, ulong maxT) {
        this.c = c;
        this.deltaS = 1.0 / cast(double)(n + 1);
        this.deltaT = (deltaS * deltaS) / (4.0 * c);
        this.maxT = maxT;
        oldMesh = new Mesh(n);
        newMesh = new Mesh(n);
    }
public:
    void run() {
        auto ny  = oldMesh.ny - 1;
        auto nx  = oldMesh.nx - 1;
        auto ds2 = deltaS * deltaS;
        auto cdtods2 = (c * deltaT) / ds2;
        writeln("o starting simulation...");
        foreach (t ; 0 .. maxT) {
            if (0 == t % 100) {
                writeln(". starting iteration ", t, " of ", maxT);
            }
            foreach (i ; 1 .. nx) {
                auto nci  = newMesh.cells[i];
                auto oci  = oldMesh.cells[i];
                auto ocip = oldMesh.cells[i - 1];
                auto ocin = oldMesh.cells[i + 1];
                foreach (j ; 1 .. ny) {
                    nci[j] = oci[j] + (cdtods2 * (ocin[j] + ocip[j] - 4.0 *
                                       oci[j] + oci[j + 1] + oci[j - 1]));
                }
            }
            // swap old and new
            auto tmp = newMesh; newMesh = oldMesh; oldMesh = tmp;
            // Constant heat source
            oldMesh.setInitialConds();
        }
    }
    void dump() {
        auto f = File("heat-img.dat", "w");
        foreach (i ; 0 .. newMesh.nx) {
            foreach (j ; 0 .. newMesh.ny) {
                f.writef("%f", newMesh.cells[i][j]);
                if (j != newMesh.ny - 1) {
                    f.write(" ");
                }
            }
            f.writeln();
        }
        // fis automatically closed
    }
}

void main()
{
    writeln(APP_NAME, " ", APP_VER);
    Simulation sim = new Simulation(N, THERM_COND, T_MAX);
    sim.oldMesh.setInitialConds();
    sim.run();
    sim.dump();
}
