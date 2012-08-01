#ifndef _utilities_h_
#define _utilities_h_

#ifdef __cplusplus
extern "C" {
#endif 

#define NEDGES   10
#define NPOINTS  10
#define NTHREADS 128

struct edge {
    int v0;
    int v1;
    float data;
    float v0_pt_data[3];
    float v1_pt_data[3];
};

struct graph {
    int v0[NEDGES];
    int v1[NEDGES];
    float v0_data[NEDGES][3];
    float v1_data[NEDGES][3];
    float data[NEDGES];
};

int graph_init_aos(char* graph_type, int npoints, int nedges,
        struct edge* edges);

int graph_init_soa(char* graph_type, int npoints, int nedges,
        struct graph* gr);

#ifdef __cplusplus
}
#endif 

#endif
