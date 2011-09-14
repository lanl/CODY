#ifndef HYDRO_DMP_H_INCL
#define HYDRO_DMP_H_INCL

#include "hydro_structs.h"

int wrt_dmp(char *filename, hydro_prob *Hp, double *u);

int rd_dmp(char *filename, hydro_prob *Hp, double **u);

void test_mesh(hydro_prob *Hp, double *u);

#endif
