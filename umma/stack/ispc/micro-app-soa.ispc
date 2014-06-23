#define NPOINTS 10
#define NEDGES  10

struct graph {
    int v0[NEDGES];
    int v1[NEDGES];
    float v0_data[NEDGES][3];
    float v1_data[NEDGES][3];
    float data[NEDGES];
};

export void edge_gather(uniform int nedges, 
        uniform struct graph * uniform g, 
        uniform float pt_data[][3],
        uniform float edge_data[]) {
    int v0;
    int v1;

    foreach (i = 0 ... nedges) {

        if (i < nedges) {
            v0 = g->v0[i];
            v1 = g->v1[i];

            g->v0_data[i][0] = pt_data[v0][0];
            g->v0_data[i][1] = pt_data[v0][1];
            g->v0_data[i][2] = pt_data[v0][2];

            g->v1_data[i][0] = pt_data[v1][0];
            g->v1_data[i][1] = pt_data[v1][1];
            g->v1_data[i][2] = pt_data[v1][2];

            g->data[i] = edge_data[i];
        }
    }
}

export void edge_compute(uniform int nedges, uniform struct graph * uniform g) {
    float v0_p0, v0_p1, v0_p2;
    float v1_p0, v1_p1, v1_p2;
    float x0, x1, x2;
    float e_data;

    foreach (i = 0 ... nedges) {
        v0_p0 = g->v0_data[i][0];
        v0_p1 = g->v0_data[i][1];
        v0_p2 = g->v0_data[i][2];

        v1_p0 = g->v1_data[i][0];
        v1_p1 = g->v1_data[i][1];
        v1_p2 = g->v1_data[i][2];

        e_data = g->data[i];

        x0 = (v0_p0 + v1_p0) * e_data;
        x1 = (v0_p1 + v1_p1) * e_data;
        x2 = (v0_p2 + v1_p2) * e_data;

        g->v0_data[i][0] = x0;
        g->v0_data[i][1] = x1;
        g->v0_data[i][2] = x2;

        g->v1_data[i][0] = x0;
        g->v1_data[i][1] = x1;
        g->v1_data[i][2] = x2;
    }
}

export void edge_scatter(uniform int nedges, uniform struct graph * uniform g,
        uniform float pt_data[][3]) {
    int v0;
    int v1;

    foreach (i = 0 ... nedges) {

        foreach_active(j) {
            v0 = g->v0[i];
            v1 = g->v1[i];

            pt_data[v0][0] += g->v0_data[i][0];
            pt_data[v0][1] += g->v0_data[i][1];
            pt_data[v0][2] += g->v0_data[i][2];

            pt_data[v1][0] += g->v1_data[i][0];
            pt_data[v1][1] += g->v1_data[i][1];
            pt_data[v1][2] += g->v1_data[i][2];
        }
    }
}

