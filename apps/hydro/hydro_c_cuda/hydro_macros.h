#ifndef HYDRO_MACROS_H_INCL
#define HYDRO_MACROS_H_INCL

#define HSCHEME_MUSCL 0
#define HSCHEME_PLMDE 1
#define HSCHEME_COLLELA 2
#define N_HSCHEMES 3

#define ID 0
#define IU 1
#define IV 2
#define IP 3

#define BND_REFL 1
#define BND_PERM 2
#define BND_WRAP 3
#define BND_INTR 4

#define SQ(x) ((x)*(x))
#define MIN(x,y) ((x)<(y)?(x):(y))
#define BL(totTh, maxTh) ((totTh)+(maxTh)-1)/(maxTh)
#define BL_TH(totTh, maxTh) ((totTh)+(maxTh)-1)/(maxTh),(maxTh)

#endif
