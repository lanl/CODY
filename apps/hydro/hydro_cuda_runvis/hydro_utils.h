#ifndef HYDRO_UTILS_H_INCL
#define HYDRO_UTILS_H_INCL

#include "hydro_structs.h"

void hydro_init(double **u, hydro_args *Ha, hydro_prob *Hp);
void print_array(char *desc, double *arr, int nx, int ny, int nvar);

#endif
