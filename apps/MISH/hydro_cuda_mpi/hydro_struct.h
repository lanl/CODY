#ifndef HYDRO_STRUCT_H_
#define HYDRO_STRUCT_H_

typedef struct __hydroProb{
    // Time
    double t;

    // Dimensions
    int nx, ny;

    // Physics
    int nvar;
    double dx, dy;
    double gamma;

    //Boundaries
    int bndL, bndR;
    int bndU, bndD;
} hydro_prob;

#define PREFIX_LEN 20
#define INIT_FN_LEN 50

typedef struct __hydroArgs{
    // I/O arguments
    double tend;
    int nstepmax;
    int noutput, nprtLine;
    double dtoutput;
    char outPre[PREFIX_LEN];
    char initFile[INIT_FN_LEN];

   //Dimensions
   int nx, ny;

    // Physics
    double sigma;
    double smallc, smallr;

    // Numerical scheme
    int niter_riemann;
    int iorder;
    double slope_type;
    int scheme;
} hydro_args;

#endif
