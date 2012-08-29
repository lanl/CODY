struct edge {
    int v0;
    int v1;
    float data;
    float v0_pt_data[3];
    float v1_pt_data[3];
};

export void edge_gather(uniform int nedges, 
        uniform struct edge edges[],
        uniform float pt_data[][3], 
        uniform float edge_data[]) {
    int v0;
    int v1;

    foreach (i = 0 ... nedges) {
        v0 = edges[i].v0;
        v1 = edges[i].v1;

        edges[i].v0_pt_data[0] = pt_data[v0][0];
        edges[i].v0_pt_data[1] = pt_data[v0][1];
        edges[i].v0_pt_data[2] = pt_data[v0][2];

        edges[i].v1_pt_data[0] = pt_data[v1][0];
        edges[i].v1_pt_data[1] = pt_data[v1][1];
        edges[i].v1_pt_data[2] = pt_data[v1][2];

        edges[i].data = edge_data[i];
    }
}

export void edge_compute(uniform int nedges, 
        uniform struct edge edges[]) {
    float v0_p0, v0_p1, v0_p2;
    float v1_p0, v1_p1, v1_p2;
    float x0, x1, x2;
    float e_data;

    foreach (i = 0 ... nedges) {
        v0_p0 = edges[i].v0_pt_data[0];
        v0_p1 = edges[i].v0_pt_data[1];
        v0_p2 = edges[i].v0_pt_data[2];

        v1_p0 = edges[i].v1_pt_data[0];
        v1_p1 = edges[i].v1_pt_data[1];
        v1_p2 = edges[i].v1_pt_data[2];

        e_data = edges[i].data;

        x0 = (v0_p0 + v1_p0) * e_data;
        x1 = (v0_p1 + v1_p1) * e_data;
        x2 = (v0_p2 + v1_p2) * e_data;

        edges[i].v0_pt_data[0] = x0;
        edges[i].v0_pt_data[1] = x1;
        edges[i].v0_pt_data[2] = x2;

        edges[i].v1_pt_data[0] = x0;
        edges[i].v1_pt_data[1] = x1;
        edges[i].v1_pt_data[2] = x2;
    }
}

export void edge_scatter(uniform int nedges, 
        uniform struct edge edges[],
        uniform float pt_data[][3]) {
    int v0;
    int v1;

    foreach (i = 0 ... nedges) {
        v0 = edges[i].v0;
        v1 = edges[i].v1;

        pt_data[v0][0] += edges[i].v0_pt_data[0];
        pt_data[v0][1] += edges[i].v0_pt_data[1];
        pt_data[v0][2] += edges[i].v0_pt_data[2];

        pt_data[v1][0] += edges[i].v1_pt_data[0];
        pt_data[v1][1] += edges[i].v1_pt_data[1];
        pt_data[v1][2] += edges[i].v1_pt_data[2];
    }
}

