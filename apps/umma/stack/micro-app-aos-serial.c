#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <getopt.h>
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>

#define NPOINTS 10000
#define NEDGES  10000

struct edge {
    int v0;
    int v1;
    float data;
    float v0_pt_data[3];
    float v1_pt_data[3];
};

struct edge edges[NEDGES];
float pt_data[NPOINTS][3];
float edge_data[NEDGES];

void print_help() {
    printf("Usage: \n");
    printf("\t --help print this message and exit \n");
    printf("\t --type Type of graph, must be one of:\n");
    printf("\t\t\t pure_random \n");
    printf("\t\t\t regular_random \n");
    printf("\t\t\t contiguous \n");
    printf("\t\t\t file \n");
    printf("\t --nloops Number of repetitions, must be \n");
    printf("\t          at least one. \n");
    printf("\t --file File from which to read graph \n");
}

double timer() {
    struct timeval tp;
    struct timezone tzp;
    long i;

    i = gettimeofday(&tp, &tzp);
    return ((double)tp.tv_sec) + ((double) tp.tv_usec) * 1e-6;
}

int graph_init(char* graph_type, char* fname) {
    lua_State *L;
    int i, j, k, v;

    L = luaL_newstate();
    luaL_openlibs(L);
    luaL_loadfile(L, "graph.lua");
    lua_pcall(L, 0, 0, 0);

    lua_getglobal(L, "create_graph");
    lua_pushstring(L, graph_type);
    lua_pushinteger(L, NPOINTS);
    lua_pushinteger(L, NEDGES);

    if (fname != NULL) {
        lua_pushstring(L, fname);
        lua_call(L, 4, 1);
    } else {
        lua_call(L, 3, 1);
    }

    /* Table is now sitting at the top of the stack */
    i = 0;
    lua_pushnil(L);  /* Make sure lua_next starts at beginning */
    while (lua_next(L, -2) != 0) {
        // fetch first key
        k = lua_tointeger(L, -2);
        lua_pushnil(L);
        while (lua_next(L, -2) != 0) { // loop over neighbors
            lua_pop(L,1);
            v = lua_tointeger(L, -1);

            // build edges array here
            edges[i].v0 = k - 1;
            edges[i].v1 = v - 1;

            for (j = 0; j < 3; j++) {
                edges[i].v0_pt_data[j] = 0;
                edges[i].v1_pt_data[j] = 0;
            }

            i++;
        }
        lua_pop(L,1);
    }

    lua_close(L);

    if (i == 0) {
        return -1;
    } else {
        return 0;
    }
}

int data_init() {
    int i;

    for (i = 0; i < NPOINTS; i++) {
        pt_data[i][0] = 1;
        pt_data[i][1] = 1;
        pt_data[i][2] = 1;
    }

    return 0;
}

int edge_data_init() {
    int i;

    for (i = 0; i < NEDGES; i++) {
        edge_data[i] = 1;
    }

    return 0;
}

int edge_gather() {
    int i,j;
    int v0;
    int v1;

    for (i = 0; i < NEDGES; i++) {
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

    return 0;
}

int edge_compute() {
    int i, j;
    float v0_p0, v0_p1, v0_p2;
    float v1_p0, v1_p1, v1_p2;
    float x0, x1, x2;
    float e_data;

    for (i = 0; i < NEDGES; i++) {
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

    return 0;
}

int edge_scatter() {
    int i;
    int v0;
    int v1;

    for (i = 0; i < NEDGES; i++) {
        v0 = edges[i].v0;
        v1 = edges[i].v1;

        pt_data[v0][0] += edges[i].v0_pt_data[0];
        pt_data[v0][1] += edges[i].v0_pt_data[1];
        pt_data[v0][2] += edges[i].v0_pt_data[2];

        pt_data[v1][0] += edges[i].v1_pt_data[0];
        pt_data[v1][1] += edges[i].v1_pt_data[1];
        pt_data[v1][2] += edges[i].v1_pt_data[2];
    }

    return 0;
}

int main(int argc, char** argv) {
    int i;
    int rv;
    double time0, time1;
    int c, opt_i;
    int nloops = 0;
    char* gt = "";
    char* fname = "";

    static struct option long_opts[] = {
        {"help",   no_argument,       0, 0},
        {"type",   required_argument, 0, 0},
        {"nloops", required_argument, 0, 0},
        {"file",   required_argument, 0, 0}
    };

    /* Parse command-line arguments */
    while (1) {
        c = getopt_long(argc, argv, "", 
                long_opts, &opt_i);

        if (c == -1) {
            break;
        }

        if (c == 0) {
            switch (opt_i) {
                case 0:
                    print_help();
                    exit(0);
                case 1:
                    gt = optarg;
                    break;
                case 2:
                    nloops = atoi(optarg);
                    break;
                case 3:
                    fname = optarg;
                    break;
            }
        } else {
            print_help();
            exit(0);
        }
    }

    /* check for errors */
    if (gt == NULL || nloops < 1) {
        print_help();
        exit(0);
    }

    
    // initialize data structures
    rv = graph_init(gt, fname);
    if (rv < 0) {
        printf("Error creating graph. \n");
        exit(0);
    }

    data_init();
    edge_data_init();

    // loop
    time0 = timer();
    for (i = 0; i < nloops; i++) {
        edge_gather();
        edge_compute();
        edge_scatter();
    }
    time1 = timer();


    // print results
    for (i = 0; i < 10; i++) {
        printf("%i : %f %f %f \n", i, pt_data[i][0], pt_data[i][1], pt_data[i][2]);
    }

    printf("Time: %f s \n", (time1 - time0) / ((float) nloops));

    return 0;
}
