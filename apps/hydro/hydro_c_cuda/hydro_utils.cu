#include <stdio.h>
#include <stdlib.h>

#include "hydro_utils.h"
#include "hydro_macros.h"
#include "hydro_dmp.h"

void hydro_init(double **u, hydro_args *Ha, hydro_prob *Hp){
    int i,j;
    int nx, ny;

    if(Ha->initFile[0]=='\0'){
     //Lower Left wave
      //Problem variables
      Hp->t=0.0;
      nx=Hp->nx=Ha->nx;
      ny=Hp->ny=Ha->ny;
      Hp->nvar=4;
      Hp->dx=5.0/Hp->nx;
      Hp->dy=5.0/Hp->ny;
      Hp->gamma=1.4;

      //Boundaries
      Hp->bndL=BND_REFL;
      Hp->bndR=BND_REFL;
      Hp->bndU=BND_REFL;
      Hp->bndD=BND_REFL;

      //Allocate mesh
      *u=(double *)calloc(Hp->nx*Hp->ny*Hp->nvar,sizeof(double));

     //Initialize mesh
      //Ambient
      for(j=0;j<ny;++j){
        for(i=0;i<nx;++i){
          (*u)[(ID*Hp->ny+j)*Hp->nx+i]=1.0;
          (*u)[(IU*Hp->ny+j)*Hp->nx+i]=0.0;
          (*u)[(IV*Hp->ny+j)*Hp->nx+i]=0.0;
          (*u)[(IP*Hp->ny+j)*Hp->nx+i]=1.0e-5;
        }
      }
      (*u)[(IP*Hp->ny+0)*Hp->nx+0]=1.0/(Hp->dx*Hp->dy);
    }else{
      if(rd_dmp(Ha->initFile,Hp,u)){
        fprintf(stderr,"Could not use dump file %s as init\n",Ha->initFile);
        exit(1);
      }
    }
}

void print_array(char *desc, double *arr, int nx, int ny, int nvar){
  int i,j,k;

  fprintf(stderr,"%s: %d variable %dx%d mesh\n",desc,nvar,nx,ny);
  for(k=0;k<nvar;++k){
    fprintf(stderr,"%s: State variable %d\n",desc,k);
    for(j=0;j<ny;++j){
      fprintf(stderr,"%s: ",desc);
      for(i=0;i<nx;++i){
        fprintf(stderr,"%lf ",arr[(k*ny+j)*nx+i]);
      }
      fprintf(stderr,"\n");
    }
    fprintf(stderr,"----\n");
  }
}

