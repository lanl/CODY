#include <stdio.h>
#include <stdlib.h>
#include "hydro.h"

int main(int argc, char* argv[]){
  double *mesh;
  hydro_prob Hp;
  hydro_args Ha;
  int i,j;
  int init=0;
  int pSMul=1;
  int iDiv=0, jDiv=0;

  //select from hardcoded inits

  if(argc<2){
    printf("No init supplied\n");
  } else {
    if(!strcmp(argv[1],"sod")) init=1;
    else if(!strcmp(argv[1],"crn")) init=2;
    else if(!strcmp(argv[1],"wsc")) init=3;
    else if(!strcmp(argv[1],"scs")) init=4;
    else printf("Unknown init\n");
  }

  //Get size multiplier
  if(init>2&&(argc<3||sscanf(argv[2],"%d",&pSMul)!=1)){
    printf("No problem size multiplier supplied\n");
    return 1;
  }
  if(pSMul<=0){
    return 1;
  }

  Ha.sigma=0.9;
  Ha.nprtLine=100;

  printf("INIT:%s\n",argv[1]);
  switch(init){
    case 1:
      Hp.nx=100;
      Hp.ny=1000;
      Hp.dx=0.25/Hp.nx;
      Hp.dy=1.0/Hp.ny;
      iDiv=100;
      jDiv=500;
      sprintf(Ha.outPre,"outDir/sod");
      Ha.tend=-1.0;
      Ha.dtoutput=-0.01;
      Ha.nstepmax=1000;
      Ha.noutput=-1;
      break;
    case 2:
      Hp.nx=1000;
      Hp.ny=1000;
      Hp.dx=1.0/Hp.nx;
      Hp.dy=1.0/Hp.ny;
      iDiv=500;
      jDiv=500;
      sprintf(Ha.outPre,"outDir/crn");
      Ha.tend=-1.0;
      Ha.dtoutput=-0.01;
      Ha.nstepmax=1000;
      Ha.noutput=-1;
      break;
    case 3:
      Hp.nx=100;
      Hp.ny=1000*pSMul;
      Hp.dx=0.25/Hp.nx;
      Hp.dy=1.0/Hp.ny;
      iDiv=100;
      jDiv=0.5*Hp.ny;
      sprintf(Ha.outPre,"outDir/wsc");
      Ha.tend=-1.0;
      Ha.dtoutput=-0.01;
      Ha.nstepmax=1000;
      Ha.noutput=-1;
      break;
    case 4:
      Hp.nx=10*pSMul;
      Hp.ny=10*pSMul;
      Hp.dx=1.0/Hp.nx;
      Hp.dy=1.0/Hp.ny;
      iDiv=4;
      jDiv=500;
      sprintf(Ha.outPre,"outDir/scs");
      Ha.tend=-1.0;
      Ha.dtoutput=-0.01;
      Ha.nstepmax=1000;
      Ha.noutput=-1;
      break;
    default:
      Hp.nx=100;
      Hp.ny=100;
      Hp.dx=0.1;
      Hp.dy=0.1;
      iDiv=0;
      jDiv=0;
      sprintf(Ha.outPre,"outDir/out");
      Ha.tend=-1.0;
      Ha.dtoutput=-0.01;
      Ha.nstepmax=100;
      Ha.noutput=-1;
      break;
  }
  Hp.t=0.0;

  Hp.nvar=4;

  mesh=(double*)malloc(Hp.nvar*Hp.nx*Hp.ny*sizeof(double));

  Hp.gamma=1.4;

  Hp.bndL=BND_REFL;
  Hp.bndR=BND_REFL;
  Hp.bndU=BND_REFL;
  Hp.bndD=BND_REFL;

  Ha.smallr=1e-10;
  Ha.smallc=1e-10;
  Ha.niter_riemann=10;
  
  for(j=0;j<Hp.ny;j++){
    for(i=0;i<Hp.nx;i++){
      mesh[i+Hp.nx*(j+Hp.ny*VARRHO)]=0.125;
      mesh[i+Hp.nx*(j+Hp.ny*VARVX )]=0.0;
      mesh[i+Hp.nx*(j+Hp.ny*VARVY )]=0.0;
      mesh[i+Hp.nx*(j+Hp.ny*VARPR )]=0.25;
    }
  }

  for(j=0;j<jDiv;j++){
    for(i=0;i<iDiv;i++){
      mesh[i+Hp.nx*(j+Hp.ny*VARRHO)]=1.0;
      mesh[i+Hp.nx*(j+Hp.ny*VARVX )]=0.0;
      mesh[i+Hp.nx*(j+Hp.ny*VARVY )]=0.0;
      mesh[i+Hp.nx*(j+Hp.ny*VARPR )]=2.5;
    }
  }

  engine(&argc,&argv,mesh,&Hp,&Ha);
  free(mesh);
}
