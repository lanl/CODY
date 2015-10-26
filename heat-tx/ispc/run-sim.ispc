/**
 * Copyright (c) 2015, Los Alamos National Security, LLC All rights reserved.
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

/* A simple 2D heat transfer simulation in C by Samuel K. Gutierrez and
 * ported to ISPC by Scott Pakin */

#include "heat-tx.h"

extern "C" uniform int
set_initial_conds(uniform mesh_t *uniform sim);

export uniform int
run_simulation(uniform simulation_t *uniform sim)
{
    uniform int rc = FAILURE;
    uniform uint64_t t;
    const uniform uint64_t nx = sim->old_mesh->nx;
    const uniform uint64_t ny = sim->old_mesh->ny;
    const uniform double ds2 = sim->params->delta_s * sim->params->delta_s;
    const uniform double cdtods2 = (sim->params->c * sim->params->delta_t) / ds2;
    const uniform uint64_t t_max = sim->params->max_t;
    uniform mesh_t *uniform new_mesh = sim->new_mesh;
    uniform mesh_t *uniform old_mesh = sim->old_mesh;

    print("o starting simulation...\n");
    for (t = 0; t < t_max; ++t) {
        if (0 == t % 100) {
            print(". starting iteration % of %\n", t, t_max);
        }
        double **uniform new_cells = new_mesh->cells;
        double **uniform old_cells = old_mesh->cells;
        for (uniform int i = 1; i < nx - 1; ++i) {
            double *uniform nci =  new_cells[i];
            double *uniform oci =  old_cells[i];
            double *uniform ocip = old_cells[i - 1];
            double *uniform ocin = old_cells[i + 1];
            foreach (j = 1 ... ny - 1) {
                double ocij = oci[j];
                nci[j] = ocij + (cdtods2 * (ocin[j] + ocip[j] - DOUBLE_C(4.0) *
                                            ocij + oci[j + 1] + oci[j - 1]));
            }
        }
        /* swap the mesh pointers */
        uniform mesh_t *uniform tmp_meshp = old_mesh; old_mesh = new_mesh; new_mesh = tmp_meshp;
        /* constant heat source */
        if (SUCCESS != (rc = set_initial_conds(old_mesh))) return rc;
    }
    return SUCCESS;
}
