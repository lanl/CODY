#ifndef _micro_app_cuda_h_
#define __micro_app_cuda_h_

#ifdef __cplusplus
extern "C" {
#endif 

#define NEDGES   10000
#define NPOINTS  10000
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
        struct edge* edges, char* fname);

int graph_init_soa(char* graph_type, int npoints, int nedges,
        struct graph* gr, char* fname);

#ifdef __cplusplus
}
#endif 

#endif
