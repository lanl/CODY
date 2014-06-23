#ifndef DEV_FUNCS_H_
#define DEV_FUNCS_H_

#include "hydro_struct.h"

void print_cuda_err(char *file, const char *func, int line, cudaError_t err);

void handle_cuda_err(char * file, const char *func, int line, cudaError_t err);

#define HANDLE_CUDA_ERROR(x) if(((x)=cudaGetLastError())!=cudaSuccess)  handle_cuda_err(__FILE__,__func__,__LINE__,(x))

void device_init(hydro_prob *Hp, hydro_args *Ha);

__global__ void calc_denom(double *u, double *den);
__global__ void redu_max(double *arrIn, double *arrOut, int nVals);

__global__ void gen_bndXL(double *u, int bnd);
__global__ void gen_bndXU(double *u, int bnd);
__global__ void gen_bndYL(double *u, int bnd);
__global__ void gen_bndYU(double *u, int bnd);

__global__ void toPrimX(double *q, double *u);
__global__ void toPrimY(double *q, double *u);

__global__ void trace(double *ql, double *qr, double *q, double dtdx, int np, int nt);
__global__ void riemann(double *flx, double *qxm, double *qxp, int np, int nt);

__global__ void addFluxX(double *u, double *flx, double dtdx);
__global__ void addFluxY(double *u, double *flx, double dtdx);

#endif
