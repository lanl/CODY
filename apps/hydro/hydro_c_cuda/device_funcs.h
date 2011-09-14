#ifndef DEVICE_UTILS_H_INCL
#define DEVICE_UTILS_H_INCL

#include "hydro_structs.h"

void handle_cuda_err(char *,const char *,int,cudaError_t);

#define HANDLE_CUDA_ERROR(x) if(((x)=cudaGetLastError())!=cudaSuccess)  handle_cuda_err(__FILE__,__func__,__LINE__,(x))

void device_init(hydro_args *Ha, hydro_prob *Hp);
void symbol_test(hydro_args *Ha, hydro_prob *Hp);

__global__ void calc_denom(double *u,double *den);
__global__ void redu_max(double *arrIn, double *arrOut, int nVals);

__global__ void gen_bndX(double *u, double *xBnd, int bndL, int bndR);
__global__ void gen_bndY(double *u, double *yBnd, int bndU, int bndD);

#endif
