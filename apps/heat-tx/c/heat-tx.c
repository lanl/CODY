/**
 * Copyright (c) 2014, Los Alamos National Security, LLC All rights reserved.
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

/* A simple 2D heat transfer simulation in C by Samuel K. Gutierrez */

/* what we are solving
 *
 * u_t = c * (u_xx * u_yy), 0 <= x,y <= NX, t >= 0
 */

/* http://www.cosy.sbg.ac.at/events/parnum05/book/horak1.pdf */

/*
 * NOTES
 * see also: Crank-Nicolson method
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <inttypes.h>
#include <math.h>

/* return codes */
enum {
    SUCCESS = 0,
    FAILURE,
    FAILURE_OOR,
    FAILURE_IO,
    FAILURE_INVALID_ARG
};

static char *app_name = "c-heat-tx";
static char *app_ver = "0.2";

/* max simulation time */
#define T_MAX 1024
/* nx and ny */
#define N 512
/* thermal conductivity */
#define THERM_COND 0.6
/* some constant */
#define K 0.4

typedef struct mesh_t {
    /* mesh size in x and y */
    uint64_t nx, ny; 
    /* mesh cells */
    double **cells;
} mesh_t;

/* simulation parameters */
typedef struct simulation_params_t {
    /* thermal conductivity */
    double c;
    double delta_s;
    /* time interval */
    double delta_t;
    /* max simulation time */
    uint64_t max_t;
} simulation_params_t;

typedef struct simulation_t {
    /* the meshes */
    mesh_t *old_mesh, *new_mesh;
    /* simulation parameters */
    simulation_params_t *params;
} simulation_t;

/* static forward declarations */
static int
mesh_construct(mesh_t **new_mesh,
               int x /* number of rows */,
               int y /* number of columns */);

static int
mesh_destruct(mesh_t *mesh);

static int
params_construct(simulation_params_t **params);

static int
params_destruct(simulation_params_t *params);

static int
set_initial_conds(mesh_t *sim);

/* ////////////////////////////////////////////////////////////////////////// */
static int
sim_param_cp(const simulation_params_t *from,
             simulation_params_t *to)
{
    if (!from || !to) return FAILURE_INVALID_ARG;
    (void)memcpy(to, from, sizeof(*from));
    return SUCCESS;
}

/* ////////////////////////////////////////////////////////////////////////// */
static int
params_construct(simulation_params_t **params)
{
    simulation_params_t *tmp = NULL;
    if (!params) return FAILURE_INVALID_ARG;
    if (NULL == (tmp = calloc(1, sizeof(*tmp)))) {
        fprintf(stderr, "out of resources @ %s:%d\n", __FILE__, __LINE__);
        return FAILURE_OOR;
    }
    *params = tmp;
    return SUCCESS;
}

/* ////////////////////////////////////////////////////////////////////////// */
static int
params_destruct(simulation_params_t *params)
{
    if (!params) return FAILURE_INVALID_ARG;
    free(params);
    return SUCCESS;
}

/* ////////////////////////////////////////////////////////////////////////// */
static int
mesh_construct(mesh_t **new_mesh,
               int x /* number of rows */,
               int y /* number of columns */)
{
    mesh_t *tmp_mesh = NULL;
    int i;

    if (NULL == new_mesh) return FAILURE_INVALID_ARG;

    tmp_mesh = (mesh_t *)calloc(1, sizeof(*tmp_mesh));
    if (NULL == tmp_mesh) {
        fprintf(stderr, "out of resources @ %s:%d\n", __FILE__, __LINE__);
        /* just bail */
        return FAILURE_OOR;
    }
    tmp_mesh->cells = (double **)calloc(x, sizeof(double *));
    if (NULL == tmp_mesh->cells) {
        fprintf(stderr, "out of resources @ %s:%d\n", __FILE__, __LINE__);
        goto error;
    }
    for (i = 0; i < y; ++i) {
        tmp_mesh->cells[i] = (double *)calloc(y, sizeof(double));
        if (NULL == tmp_mesh->cells[i]) {
            fprintf(stderr, "out of resources @ %s:%d\n", __FILE__, __LINE__);
            goto error;
        }
    }
    tmp_mesh->nx = x;
    tmp_mesh->ny = y;

    *new_mesh = tmp_mesh;
    return SUCCESS;

error:
    if (NULL != tmp_mesh) free(tmp_mesh);
    mesh_destruct(tmp_mesh);
    return FAILURE_OOR;
}

/* ////////////////////////////////////////////////////////////////////////// */
static int
mesh_destruct(mesh_t *mesh)
{
    uint64_t i;

    if (!mesh) return FAILURE_INVALID_ARG;
    if (mesh->cells) {
        for (i = 0; i < mesh->ny; ++i) {
            if (NULL != mesh->cells[i]) {
                free(mesh->cells[i]);
                mesh->cells[i] = NULL;
            }
        }
        free(mesh->cells);
        mesh->cells = NULL;
    }
    free(mesh);
    return SUCCESS;
}

/* ////////////////////////////////////////////////////////////////////////// */
static int
gen_meshes(simulation_t *sim, uint64_t nx, uint64_t ny)
{
    int rc = FAILURE;

    if (SUCCESS != (rc = mesh_construct(&sim->old_mesh, nx, ny))) {
        fprintf(stderr, "\nmesh_construct failure @ %s:%d\n", __FILE__,
                __LINE__);
        goto out;
    }
    if (SUCCESS != (rc = mesh_construct(&sim->new_mesh, nx, ny))) {
        fprintf(stderr, "\nmesh_construct failure @ %s:%d\n", __FILE__,
                __LINE__);
        goto out;
    }
out:
    if (SUCCESS != rc) {
        mesh_destruct(sim->new_mesh);
        mesh_destruct(sim->old_mesh);
    }
    return rc;
}

/* ////////////////////////////////////////////////////////////////////////// */
static int
simulation_destruct(simulation_t *sim)
{
    if (!sim) return FAILURE_INVALID_ARG;
    (void)params_destruct(sim->params);
    (void)mesh_destruct(sim->new_mesh);
    (void)mesh_destruct(sim->old_mesh);
    free(sim);
    return SUCCESS;
}
/* ////////////////////////////////////////////////////////////////////////// */
static int
simulation_construct(simulation_t **new_sim,
               const simulation_params_t *params)
{
    simulation_t *sim = NULL;
    int rc = FAILURE;

    if (!new_sim || !params) return FAILURE_INVALID_ARG;

    if (NULL == (sim = (simulation_t *)calloc(1, sizeof(*sim)))) {
        fprintf(stderr, "out of resources @ %s:%d\n", __FILE__, __LINE__);
        return FAILURE_OOR;
    }
    if (SUCCESS != (rc = params_construct(&sim->params))) {
        fprintf(stderr, "params_construct failure @ %s:%d: rc = %d\n", __FILE__,
                __LINE__, rc);
        goto out;
    }
    if (SUCCESS != (rc = sim_param_cp(params, sim->params))) {
        fprintf(stderr, "sim_param_cp failure @ %s:%d: rc = %d\n", __FILE__,
                __LINE__, rc);
        goto out;
    }
    if (SUCCESS != (rc = gen_meshes(sim, N, N))) {
        fprintf(stderr, "gen_meshes failure @ %s:%d: rc = %d\n", __FILE__,
                __LINE__, rc);
        /* on failure, gen_meshes cleans up after itself */
        goto out;
    }
    *new_sim = sim;
out:
    return rc;
}

/* ////////////////////////////////////////////////////////////////////////// */
static int
init_params(simulation_params_t *params,
            int n,
            double c,
            int max_t)
{
    if (NULL == params) return FAILURE_INVALID_ARG;

    printf("o initializing simulation parameters...\n");

    params->c = c;
    params->max_t = max_t;
    params->delta_s = 1.0 / (double)(n + 1);
    /* we know from theory that we have to obey the restriction:
     * delta_t <= (delta_s)^2/2c. so just make them equal.
     */
    params->delta_t = pow(params->delta_s, 2.0) / (4.0 * params->c);

    printf(". max_t: %"PRIu64"\n", params->max_t);
    printf(". c: %lf\n", params->c);
    printf(". delta_s: %lf\n", params->delta_s);
    printf(". delta_t: %lf\n\n", params->delta_t);

    return SUCCESS;
}

/* ////////////////////////////////////////////////////////////////////////// */
static int
dump(const simulation_t *sim)
{
    FILE *imgfp = NULL;
    uint64_t i, j;

    if (NULL == (imgfp = fopen("heat-img.dat", "wb"))) {
        fprintf(stderr, "fopen failure @ %s:%d\n", __FILE__, __LINE__);
        return FAILURE_IO;
    }
    /* write the matrix */
    for (i = 0; i < sim->new_mesh->nx; ++i) {
        for (j = 0; j < sim->new_mesh->ny; ++j) {
            fprintf(imgfp, "%lf%s", sim->new_mesh->cells[i][j],
                    (j == sim->new_mesh->ny - 1) ? "" : " ");
        }
        fprintf(imgfp, "\n");
    }
    fflush(imgfp);
    fclose(imgfp);
    return SUCCESS;
}

/* ////////////////////////////////////////////////////////////////////////// */
static int
run_simulation(simulation_t *sim)
{
    int rc = FAILURE;
    uint64_t t, i, j;
    uint64_t nx = sim->old_mesh->nx;
    uint64_t ny = sim->old_mesh->ny;
    double ds2 = sim->params->delta_s * sim->params->delta_s;
    double cdtods2 = (sim->params->c * sim->params->delta_t) / ds2;
    uint64_t t_max = sim->params->max_t;
    mesh_t *new_mesh = sim->new_mesh;
    mesh_t *old_mesh = sim->old_mesh;

    printf("o starting simulation...\n");
    for (t = 0; t < t_max; ++t) {
        if (0 == t % 100) {
            printf(". starting iteration %"PRIu64" of %"PRIu64"\n", t, t_max); 
        }
        for (i = 1; i < nx - 1; ++i) {
            for (j = 1; j < ny - 1; ++j) {
                new_mesh->cells[i][j] =
                    old_mesh->cells[i][j] +
                    (cdtods2 *
                     (old_mesh->cells[i + 1][j] +
                      old_mesh->cells[i - 1][j] -
                      4.0 * old_mesh->cells[i][j] +
                      old_mesh->cells[i][j + 1] +
                      old_mesh->cells[i][j - 1]));
            }
        }
        /* swap the mesh pointers */
        mesh_t *tmp_meshp = old_mesh; old_mesh = new_mesh; new_mesh = tmp_meshp;
        /* constant heat source */
        if (SUCCESS != (rc = set_initial_conds(old_mesh))) return rc;
    }
    return SUCCESS;
}

/* ////////////////////////////////////////////////////////////////////////// */
static int
set_initial_conds(mesh_t *mesh)
{
    if (NULL == mesh) return FAILURE_INVALID_ARG;

    int x0 = mesh->nx / 2;
    int y0 = mesh->ny / 2;
    int x = mesh->nx / 4, y = 0;
    int radius_err = 1 - x;

    while (x >= y) {
        mesh->cells[ x + x0][ y + y0] = K;
        mesh->cells[ x + x0][ y + y0] = K * .50;
        mesh->cells[ y + x0][ x + y0] = K * .60;
        mesh->cells[-x + x0][ y + y0] = K * .70;
        mesh->cells[-y + x0][ x + y0] = K * .80;
        mesh->cells[-x + x0][-y + y0] = K * .70;
        mesh->cells[-y + x0][-x + y0] = K * .60;
        mesh->cells[ x + x0][-y + y0] = K * .50;
        mesh->cells[ y + x0][-x + y0] = K;
        y++;
        if (radius_err < 0) radius_err += 2 * y + 1;
        else {
            --x;
            radius_err += 2 * (y - x + 1);
        }
    }
    return SUCCESS;
}

/* ////////////////////////////////////////////////////////////////////////// */
int
main(void)
{
    int rc = FAILURE;
    int erc = EXIT_FAILURE;
    simulation_params_t *params = NULL;
    simulation_t *sim = NULL;

    /* print application banner */
    printf("o %s %s\n", app_name, app_ver);

    if (SUCCESS != (rc = params_construct(&params))) {
        fprintf(stderr, "params_construct failure @ %s:%d: rc = %d\n", __FILE__,
                __LINE__, rc);
        goto cleanup;
    }
    if (SUCCESS != (rc = init_params(params, N, THERM_COND, T_MAX))) {
        fprintf(stderr, "init_params failure @ %s:%d: rc = %d\n", __FILE__,
                __LINE__, rc);
        goto cleanup;
    }
    if (SUCCESS != (rc = simulation_construct(&sim, params))) {
        fprintf(stderr, "simulation_construct failure @ %s:%d: rc = %d\n",
                __FILE__, __LINE__, rc);
        goto cleanup;
    }
    if (SUCCESS != (rc = set_initial_conds(sim->old_mesh))) {
        fprintf(stderr, "set_initial_conds failure @ %s:%d: rc = %d\n",
                __FILE__, __LINE__, rc);
        goto cleanup;
    }
    if (SUCCESS != (rc = run_simulation(sim))) {
        fprintf(stderr, "run_simulation failure @ %s:%d: rc = %d\n",
                __FILE__, __LINE__, rc);
        goto cleanup;
    }
    if (SUCCESS != dump(sim)) {
        fprintf(stderr, "dump failure @ %s:%d: rc = %d\n",
                __FILE__, __LINE__, rc);
        goto cleanup;
    }
    /* all is well */
    erc = EXIT_SUCCESS;

cleanup:
    (void)simulation_destruct(sim);
    (void)params_destruct(params);
    return erc;
}
