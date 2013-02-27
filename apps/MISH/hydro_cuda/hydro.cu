#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "hydro.h"
#include "dev_funcs.h"
#include "outfile.h"

#define CDT_REGS 64
#define STEP_REGS 64

hydro_args *Ha;
hydro_prob *Hp;
size_t varSize, meshSize, primSize, qSize, flxSize;
double *bndLS, *bndLR, *bndHS, *bndHR;
double *d_u;
double *d_q;
double *d_qr, *d_ql;
double *d_flx;

//temp variable to store data for checking
double *h_ref;

//CUDA vars
int nTh;

//MPI Vars
int bndT, bndB;

void printArray(char* label,double *arr, int nvar, int nx, int ny, int nHx, int nHy){
  int nV,i,j;
  printf("Array %s\n",label);
  for(nV=0;nV<nvar;nV++){
    printf("Variable %d:\n",nV);
    printf("%s-<%3d>%3s",label,nV,"");
    for(i=0;i<nx;i++){
      printf("|i%8di",i);
    }
    printf("|\n");
    for(j=0;j<ny;j++){
      printf("%s-<%3d>%3d",label,nV,j);
      for(i=0;i<nx;i++){
	printf("|%10g",arr[i+nHx+(nx+2*nHx)*(j+nHy+(ny+2*nHy)*nV)]);
      }
      printf("|\n");
    }
  }
}

void setHHalo(int LBnd, int RBnd){
  gen_bndXL<<<BL_TH(2*Hp->ny,nTh)>>>(d_u,LBnd);
  gen_bndXU<<<BL_TH(2*Hp->ny,nTh)>>>(d_u,RBnd);
}

void setVHalo(int TBnd, int BBnd){
  gen_bndYL<<<BL_TH(2*Hp->nx,nTh)>>>(d_u,BBnd);
  gen_bndYU<<<BL_TH(2*Hp->nx,nTh)>>>(d_u,TBnd);
}

double sumArray(double *mesh, int var, int nx, int ny, int nHx, int nHy){
  int lI,i,j;
  double sum, corr;
  double nsum, c_nxt;

  sum=0.0;
  corr=0.0;


  for(lI=0;lI<nx*ny;lI++){
    i=lI%(nx);
    j=lI/(nx);
    c_nxt=mesh[i+nHx+(nx+2*nHx)*(j+nHy+(ny+2*nHy)*var)];
    nsum=sum+c_nxt;
    corr=(nsum-sum)-c_nxt;
    sum=nsum;
  }
  return sum+corr;
}

int nansIn(double *mesh, int var, int nx, int ny, int nHx, int nHy);

void nanScan(int *ret, double *mesh, int nvar, int nx, int ny, int nHx, int nHy){
  int i;

  for(i=0;i<nvar;i++){
    ret[i]=nansIn(mesh,i,nx,ny,nHx,nHy);
  }
}

int nansIn(double *mesh, int var, int nx, int ny, int nHx, int nHy){
  int lI,i,j;
  int cnt=0;
  double v;
  
  for(lI=0;lI<nx*ny;lI++){
    i=lI%(nx);
    j=lI/(nx);
    v=mesh[i+nHx+(nx+2*nHx)*(j+nHy+(ny+2*nHy)*var)];
    if(isnan(v)||isinf(v)){
      cnt++;
      printf("NAN in var %d @ %d,%d: %g\n",var,i,j,v);
    }
  }
  return cnt;
}

void runPass(double dt, int dir){
  int np,nt;
  double dx,dy;
  char dCh;
  char outLab[30];


  if(dir==0){
    np=Hp->nx;
    nt=Hp->ny;
    dx=Hp->dx;
    dy=Hp->dy;
    dCh='x';
    setHHalo(Hp->bndL,Hp->bndR);
    toPrimX<<<BL_TH((np+4)*nt,nTh)>>>(d_q,d_u);
  }else{
    np=Hp->ny;
    nt=Hp->nx;
    dx=Hp->dy;
    dy=Hp->dx;
    dCh='y';
    setVHalo(bndT,bndB);
    toPrimY<<<BL_TH((np+4)*nt,nTh)>>>(d_q,d_u);
  }
  //cudaMemcpy(h_ref,d_u,meshSize*sizeof(double),cudaMemcpyDeviceToHost);
  //sprintf(outLab,"MESH-%c",dCh);
  //if(dCh=='x')
    //printArray(outLab,h_ref,Hp->nvar,Hp->nx+4,Hp->ny,0,2);
  //else
    //printArray(outLab,h_ref,Hp->nvar,Hp->nx,Hp->ny+4,2,0);
  //cudaMemcpy(h_ref,d_q,primSize*sizeof(double),cudaMemcpyDeviceToHost);
  //sprintf(outLab,"Q   -%c",dCh);
  //printArray(outLab,h_ref,Hp->nvar,np+4,nt,0,0);
  trace<<<BL_TH((np+2)*nt,nTh)>>>(d_ql,d_qr,d_q,dt/dx,np,nt);
  //cudaMemcpy(h_ref,d_ql,qSize*sizeof(double),cudaMemcpyDeviceToHost);
  //sprintf(outLab,"QL  -%c",dCh);
  //printArray(outLab,h_ref,Hp->nvar,np+2,nt,0,0);
  //cudaMemcpy(h_ref,d_qr,qSize*sizeof(double),cudaMemcpyDeviceToHost);
  //sprintf(outLab,"QR  -%c",dCh);
  //printArray(outLab,h_ref,Hp->nvar,np+2,nt,0,0);
  riemann<<<BL_TH((np+1)*nt,nTh)>>>(d_flx,d_ql,d_qr,np,nt);
  //cudaMemcpy(h_ref,d_flx,flxSize*sizeof(double),cudaMemcpyDeviceToHost);
  //sprintf(outLab,"FLX -%c",dCh);
  //printArray(outLab,h_ref,Hp->nvar,np+1,nt,0,0);
  if(dir==0){
    addFluxX<<<BL_TH((np)*nt,nTh)>>>(d_u,d_flx,dt/dx);
  }else{
    addFluxY<<<BL_TH((np)*nt,nTh)>>>(d_u,d_flx,dt/dx);
  }
  //cudaMemcpy(h_ref,d_u,meshSize*sizeof(double),cudaMemcpyDeviceToHost);
  //sprintf(outLab,"POST-%c",dCh);
  //printArray(outLab,h_ref,Hp->nvar,Hp->nx,Hp->ny,2,2);
}

void engine(int *argc, char **argv[], double *gMesh, hydro_prob *Hyp, hydro_args *Hya){
  int n, nV, lI, i,j,k;
  int bndL;
  int bndH;
  double dt, dt_denom;
  double cTime, nxttout;

  double volCell;
  double oTM,oTE;
  double TM,TE;
  int M_exp, E_exp;
  double M_prec, E_prec;
  double *lMesh;
  double *recvMesh;
  double *tmp, *d_denA, *d_denB;

  char outfile[30];
  char outLab[30];

  float runT;

  //Cuda vars
  int dev;
  cudaDeviceProp prop;
  int nDen, redBlocks;
  cudaError_t cuErrVar;
  int rpBl;
  int mxTh, thWp;
  int nBlockM;
  int nThCDT, nThStep;
  size_t shMpBl;
  size_t mem_reqd, mem_avail;

  cudaEvent_t start, end;

  Hp=Hyp;
  Ha=Hya;

  device_init(Hp,Ha);

  cudaGetDevice(&dev);
  HANDLE_CUDA_ERROR(cuErrVar);
  cudaGetDeviceProperties(&prop,dev);
  HANDLE_CUDA_ERROR(cuErrVar);

  //Set dev comm vars
  printf("Setting vars for kernel calls\n");
  printf("Device is %s with compute capability %d.%d\n",prop.name,prop.major,prop.minor);
  mxTh=prop.maxThreadsPerBlock;
  shMpBl=prop.sharedMemPerBlock;
  rpBl=prop.regsPerBlock;
  thWp=prop.warpSize;
  nThCDT=mxTh;//thread upper bound
  nThCDT=MIN(nThCDT,thWp*(rpBl/(CDT_REGS*thWp)));//register upper bound
  nThStep=thWp;//thread upper bound
  nThStep=MIN(nThStep,thWp*(rpBl/(STEP_REGS*thWp)));//register upper bound
  nTh=nThCDT;
  nBlockM=((Hp->ny*Hp->nx)+nTh-1)/nTh;
  printf("Per block: Max threads %d, regs %d\n", mxTh,rpBl, shMpBl);
  printf("Block size lims: cdt %d step %d\n",nThCDT,nThStep);

  n=0;
  cTime=0;
  nxttout=-1.0;

  if(Ha->nstepmax<0&&Ha->tend<0.0)return;


  //Calculate sizes for dispersal
  Hp->ny=Hp->ny;

  //Calculate arraysizes
  varSize=(Hp->nx+4)*(Hp->ny+4);
  meshSize=Hp->nvar*(Hp->nx+4)*(Hp->ny+4);
  if(Hp->ny>=Hp->nx){
    primSize=Hp->nvar*(Hp->nx+4)*Hp->ny;
    qSize   =Hp->nvar*(Hp->nx+2)*Hp->ny;
    flxSize =Hp->nvar*(Hp->nx+1)*Hp->ny;
  }else{
    primSize=Hp->nvar*(Hp->ny+4)*Hp->nx;
    qSize   =Hp->nvar*(Hp->ny+2)*Hp->nx;
    flxSize =Hp->nvar*(Hp->ny+1)*Hp->nx;
  }
  //printf("Done loc size calcs\n");

  //Print relative sizes of memory requirements
  mem_reqd=(meshSize+primSize+2*qSize+flxSize)*sizeof(double);
  mem_avail=prop.totalGlobalMem;
  printf("%u/%u of %f%% memory required for a %d var %dx%d mesh\n",mem_reqd,mem_avail,
         ((double)mem_reqd/(double)mem_avail)*100.0,Hp->nvar,Hp->nx,Hp->ny);

  if(mem_reqd>mem_avail){
    fprintf(stderr,"Not enough memory on GPU\n");
    exit(1);
  }

  //Allocate arrays
  recvMesh=(double*)malloc(Hp->nvar*Hp->nx*Hp->ny*sizeof(double));
  lMesh =(double*)malloc(meshSize*sizeof(double));
  h_ref =(double*)malloc(meshSize*sizeof(double));
  cudaMalloc(&d_u,meshSize*sizeof(double));
  HANDLE_CUDA_ERROR(cuErrVar);
  cudaMalloc(&d_q,primSize*sizeof(double));
  HANDLE_CUDA_ERROR(cuErrVar);
  cudaMalloc(&d_qr,qSize*sizeof(double));
  HANDLE_CUDA_ERROR(cuErrVar);
  cudaMalloc(&d_ql,qSize*sizeof(double));
  HANDLE_CUDA_ERROR(cuErrVar);
  cudaMalloc(&d_flx,flxSize*sizeof(double));
  HANDLE_CUDA_ERROR(cuErrVar);
  cudaMalloc(&d_denA,nBlockM*sizeof(double));
  HANDLE_CUDA_ERROR(cuErrVar);
  cudaMalloc(&d_denB,nBlockM*sizeof(double));
  HANDLE_CUDA_ERROR(cuErrVar);

  //printf("Arrays allocated\n");

  //Zero lMesh
  for(k=0;k<Hp->nvar;k++){
    for(i=0;i<Hp->nx;i++){
      for(j=0;j<Hp->ny;j++){
        lMesh[i+2+(Hp->nx+4)*(j+2+(Hp->ny+4)*k)]=gMesh[i+(Hp->nx)*(j+Hp->ny*k)];
      }
    }
  }

  if(Ha->tend>0.0){
    nxttout=Ha->tend;
  }
  if(Ha->dtoutput>0.0&&nxttout>Ha->dtoutput){
    nxttout=Ha->dtoutput;
  }

  volCell=Hp->dx*Hp->dy;
  oTM=sumArray(lMesh,VARRHO,Hp->nx,Hp->ny,2,2);
  oTE=sumArray(lMesh,VARPR ,Hp->nx,Hp->ny,2,2);
#ifdef M_PREC_CMP
  frexp(oTM,&M_exp);
  frexp(oTE,&E_exp);

  M_prec=ldexp(DBL_EPSILON,M_exp-1);
  E_prec=ldexp(DBL_EPSILON,E_exp-1);

  printf("TM:%g+-%g TE:%g+-%g\n",volCell*oTM,volCell*M_prec,volCell*oTE,volCell*E_prec);
#endif

  //Print initial condition
  snprintf(outfile,29,"%s%05d",Ha->outPre,n);
  writeVis(outfile,gMesh,Hp->dx,Hp->dy,Hp->nvar,Hp->nx,Hp->ny);
  printf("INIT: TM: %g TE: %g\n",volCell*oTM,volCell*oTE);

  //Move mesh onto GPU
  cudaMemcpy(d_u,lMesh,meshSize*sizeof(double),cudaMemcpyHostToDevice);
  HANDLE_CUDA_ERROR(cuErrVar);

  //Get start time  
  cudaEventCreate(&start);
  cudaEventCreate(&end);
  cudaEventRecord(start,0);

  while((n<Ha->nstepmax||Ha->nstepmax<0)&&(cTime<Ha->tend||Ha->tend<0)){
    //Calculate timestep
    dt=0.0;
    nDen=nBlockM;
    redBlocks=nDen;
    nTh=nThCDT;
    //cudaMemcpy(lMesh,d_u,meshSize*sizeof(double),cudaMemcpyDeviceToHost);
    //printArray("Mesh",lMesh,4,Hp->nx,Hp->ny,2,2);
    calc_denom<<<nBlockM,nTh,nTh*sizeof(double)>>>(d_u,d_denA);
    //cudaMemcpy(recvMesh,d_denA,nDen*sizeof(double),cudaMemcpyDeviceToHost);
    //printArray("Dens",recvMesh,1,nDen,1,0,0);
    while(redBlocks>1){
      redBlocks=(nDen+2*nTh-1)/(2*nTh);
      redu_max<<<redBlocks,nTh,nTh*sizeof(double)>>>(d_denA,d_denB,nDen);
      nDen=redBlocks;
      //cudaMemcpy(recvMesh,d_denA,nDen*sizeof(double),cudaMemcpyDeviceToHost);
      //printArray("Dens",recvMesh,1,nDen,1,0,0);
      tmp=d_denA;
      d_denA=d_denB;
      d_denB=tmp;
    }
    cudaMemcpy(&dt_denom,d_denA,sizeof(double),cudaMemcpyDeviceToHost);
    //printf("ITER %d denom=%g\n",n,dt_denom);
    dt=0.5*Ha->sigma/dt_denom;
    if(nxttout>0.0&&dt>(nxttout-cTime)){
      printf("Adjusting timestep from %g to %g for iter %d\n",dt,nxttout-cTime,n);
      dt=(nxttout-cTime);
    }
    //break;
    if(n%2==0){
      //X Dir
      runPass(dt,0);
      //Y Dir
      runPass(dt,1);
    }else{
      //Y Dir
      runPass(dt,1);
      //X Dir
      runPass(dt,0);
    }
    n+=1;
    cTime+=dt;
    if(n%Ha->nprtLine==0||(cTime>=nxttout&&nxttout>0)||(Ha->noutput>0&&n%Ha->noutput==0)){
      cudaMemcpy(lMesh,d_u,meshSize*sizeof(double),cudaMemcpyDeviceToHost);
      HANDLE_CUDA_ERROR(cuErrVar);
      if(n%Ha->nprtLine==0){
        TM=sumArray(lMesh,VARRHO,Hp->nx,Hp->ny,2,2);
        TE=sumArray(lMesh,VARPR ,Hp->nx,Hp->ny,2,2);
        printf("Iter %05d time %f dt %g TM: %g TE: %g\n",n,cTime,dt,volCell*TM,volCell*TE);
        if(0&&(fabs(TM-oTM)>0.0||fabs(TE-oTE)>0.0)){
          printf("ERR(%5s): Mass %g Ene %g\n","RAW",volCell*(TM-oTM),volCell*(TE-oTE));
          printf("ERR(%5s): Mass %g Ene %g\n","%",100.0*(TM-oTM)/oTM,100.0*(TE-oTE)/oTE);
#ifdef M_PREC_CMP
          printf("ERR(%5s): Mass %g Ene %g\n","M PRE",(TM-oTM)/M_prec,(TE-oTE)/E_prec);
#endif
        }
      }
      if((cTime>=nxttout&&nxttout>0)||(Ha->noutput>0&&n%Ha->noutput==0)){
        printf("Vis file @ time %f iter %d\n",cTime,n);
        if(cTime>=nxttout&&nxttout>0.0){
          nxttout+=Ha->dtoutput;
          if(nxttout>Ha->tend&&Ha->tend>0)nxttout=Ha->tend;
          printf("Next Vis Time: %f\n",nxttout);
        }
        for(k=0;k<Hp->nvar;k++){
          for(i=0;i<Hp->nx;i++){
	    for(j=0;j<Hp->ny;j++){
              gMesh[(k*Hp->ny+j)*Hp->nx+i]=lMesh[(k*(Hp->ny+4)+j+2)*(Hp->nx+4)+i+2];
	    }
	  }
        }
        snprintf(outfile,29,"%s%05d",Ha->outPre,n);
        writeVis(outfile,gMesh,Hp->dx,Hp->dy,Hp->nvar,Hp->nx,Hp->ny);
      }
    }
  }
  printf("time: %f, %d iters run\n",cTime,n);

  //Get end time
  cudaEventRecord(end,0);
  cudaEventElapsedTime(&runT,start,end);

  printf("TFMT:%s,%s,%s,%s,%s,%s,%s,%s\n","cType","mType","init","nproc","nth","niters","ncells","wRunt");
  printf("TIME:%s,%s,%s,%d,%d,%d,%d,%g\n","\"CUDA\"","\"GPU:?\"","\"Init\"",1,1,n,Hp->nx*Hp->ny,runT*1.0e-3);

  //Print final condition
  cudaMemcpy(lMesh,d_u,meshSize*sizeof(double),cudaMemcpyDeviceToHost);
  HANDLE_CUDA_ERROR(cuErrVar);
  for(k=0;k<Hp->nvar;k++){
    for(i=0;i<Hp->nx;i++){
      for(j=0;j<Hp->ny;j++){
        gMesh[(k*Hp->ny+j)*Hp->nx+i]=lMesh[(k*(Hp->ny+4)+j+2)*(Hp->nx+4)+i+2];
      }
    }
  }
  snprintf(outfile,29,"%s%05d",Ha->outPre,n);
  writeVis(outfile,gMesh,Hp->dx,Hp->dy,Hp->nvar,Hp->nx,Hp->ny);
  Hp->t+=cTime;

  free(recvMesh);
  free(lMesh);
  cudaFree(d_u  );
  cudaFree(d_q  );
  cudaFree(d_qr );
  cudaFree(d_ql );
  cudaFree(d_flx);
  cudaFree(d_denA);
  cudaFree(d_denB);

  printf("Returning from engine\n");
}
