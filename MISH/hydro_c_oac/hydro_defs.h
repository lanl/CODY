#ifndef HYDRO_DEFS_H_
#define HYDRO_DEFS_H_

#define MAX(x,y) ((x)<(y))?(y):(x)
#define MIN(x,y) ((x)>(y))?(y):(x)

#define VARRHO 0
#define VARVX  1
#define VARVY  2
#define VARPR  3

#define BND_REFL 0
#define BND_PERM 1

#endif //HYDRO_DEFS_H_
