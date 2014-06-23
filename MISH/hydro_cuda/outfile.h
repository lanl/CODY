#ifndef OUTFILE_H_
#define OUTFILE_H_

#include "hydro_struct.h"

void writeVis(char *name, double *mesh, double dx, double dy, int nvar, int nx, int ny);

#endif //OUTFILE_H_
