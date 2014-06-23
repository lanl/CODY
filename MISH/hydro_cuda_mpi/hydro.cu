#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <mpi.h>
#include <float.h>
#include "hydro.h"
#include "dev_funcs.h"
#include "outfile.h"

#define CDT_REGS 64
#define STEP_REGS 64

hydro_args *Ha;
hydro_prob *Hp;
double *lMesh;
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
int pProc, nProc;
int rank, size;
int myNy;
int varSize;

void printArray(char* label,double *arr, int nvar, int nx, int ny){
  int nV,i,j;
  printf("N[%2d]: Array %s\n",rank,label);
  for(nV=0;nV<nvar;nV++){
    printf("N[%2d]:Variable %d:\n",rank,nV);
    printf("N[%2d]:%3s",rank,"");
    for(i=0;i<nx;i++){
      printf("|i%8di",i);
    }
    printf("|\n");
    for(j=0;j<ny;j++){
      printf("N[%2d]:%s-<%3d>%3d",rank,label,nV,j);
      for(i=0;i<nx;i++){
	printf("|%10g",arr[i+nx*(j+ny*nV)]);
	//printf("N[%2d]:V(%d,%d,%d)=%g\n",rank,nV,i,j,arr[i+nx*(j+ny*nV)]);
      }
      printf("|\n");
    }
  }
}

void setHHalo(int LBnd, int RBnd){
  gen_bndXL<<<BL_TH(2*myNy,nTh)>>>(d_u,LBnd);
  gen_bndXU<<<BL_TH(2*myNy,nTh)>>>(d_u,RBnd);
}

void setVHalo(int TBnd, int BBnd){
  int row=Hp->nx+4;
  int nV;
  MPI_Request reqs[16];
  MPI_Status stat[16];
  for(nV=0;nV<Hp->nvar;nV++){
    cudaMemcpy(lMesh+(2   )*row+nV*varSize,d_u+(2   )*row+nV*varSize,2*row*sizeof(double),cudaMemcpyDeviceToHost);
    cudaMemcpy(lMesh+(myNy)*row+nV*varSize,d_u+(myNy)*row+nV*varSize,2*row*sizeof(double),cudaMemcpyDeviceToHost);
    MPI_Irecv(lMesh+(     0)*row+nV*varSize,2*row,MPI_DOUBLE,pProc,1+2*nV,MPI_COMM_WORLD,reqs+ 0+nV);
    MPI_Isend(lMesh+(     2)*row+nV*varSize,2*row,MPI_DOUBLE,pProc,2+2*nV,MPI_COMM_WORLD,reqs+ 4+nV);
    MPI_Irecv(lMesh+(myNy+2)*row+nV*varSize,2*row,MPI_DOUBLE,nProc,2+2*nV,MPI_COMM_WORLD,reqs+ 8+nV);
    MPI_Isend(lMesh+(myNy  )*row+nV*varSize,2*row,MPI_DOUBLE,nProc,1+2*nV,MPI_COMM_WORLD,reqs+12+nV);
  }
  MPI_Waitall(16,reqs,stat);
  for(nV=0;nV<Hp->nvar;nV++){
    cudaMemcpy(d_u+(     0)*row+nV*varSize,lMesh+(     0)*row+nV*varSize,2*row*sizeof(double),cudaMemcpyHostToDevice);
    cudaMemcpy(d_u+(myNy+2)*row+nV*varSize,lMesh+(myNy+2)*row+nV*varSize,2*row*sizeof(double),cudaMemcpyHostToDevice);
  }
  gen_bndYL<<<BL_TH(2*Hp->nx,nTh)>>>(d_u,TBnd);
  gen_bndYU<<<BL_TH(2*Hp->nx,nTh)>>>(d_u,BBnd);
}

double sumArray(double *mesh, int var, int nx, int ny, int nHx, int nHy){
  int lI,i,j;
  double sum, corr;
  double nsum, c_nxt;
  double gsum, gcorr;

  sum=0.0;
  corr=0.0;

  //printf("N[%2d]: summing var %d over range (%d,%d) w/ halos of size (%d,%d)\n",rank,var,nx,ny,nHx,nHy);

  for(lI=0;lI<nx*ny;lI++){
    i=lI%(nx);
    j=lI/(nx);
    c_nxt=mesh[i+nHx+(nx+2*nHx)*(j+nHy+(ny+2*nHy)*var)];
    nsum=sum+c_nxt;
    corr=(nsum-sum)-c_nxt;
    sum=nsum;
    //printf("N[%2d]: @(%d,%d) %g\n",rank,i,j,sum+corr);
  }
  //printf("N[%2d]: lsum %g+%g\n",rank,sum,corr);
  MPI_Allreduce(&sum,&gsum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&corr,&gcorr,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  //printf("N[%2d]: gsum %g+%g\n",rank,gsum,gcorr);
  return gsum+gcorr;
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
      printf("N[%2d]: NAN in var %d @ %d,%d: %g\n",rank,var,i,j,v);
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
    nt=myNy;
    dx=Hp->dx;
    dy=Hp->dy;
    dCh='x';
    //printf("N[%2d]:X-pass\n",rank);
    setHHalo(Hp->bndL,Hp->bndR);
    toPrimX<<<BL_TH((np+4)*nt,nTh)>>>(d_q,d_u);
  }else{
    np=myNy;
    nt=Hp->nx;
    dx=Hp->dy;
    dy=Hp->dx;
    dCh='y';
    //printf("N[%2d]:Y-pass\n",rank);
    setVHalo(bndT,bndB);
    toPrimY<<<BL_TH((np+4)*nt,nTh)>>>(d_q,d_u);
  }
  trace<<<BL_TH((np+2)*nt,nTh)>>>(d_ql,d_qr,d_q,dt/dx,np,nt);
  riemann<<<BL_TH((np+1)*nt,nTh)>>>(d_flx,d_ql,d_qr,np,nt);
  if(dir==0){
    addFluxX<<<BL_TH((np)*nt,nTh)>>>(d_u,d_flx,dt/dx);
  }else{
    addFluxY<<<BL_TH((np)*nt,nTh)>>>(d_u,d_flx,dt/dx);
  }
}

void engine(int *argc, char **argv[], double *gMesh, hydro_prob *Hyp, hydro_args *Hya){
  int n, nV, lI, i,j;
  int bndL;
  int bndH;
  double dt, dt_denom, gDenom;
  double cTime, nxttout;

  double volCell;
  double ogTM,ogTE;
  double lTM,lTE;
  double gTM,gTE;
  int M_exp, E_exp;
  double M_prec, E_prec;
  double *recvMesh;
  double *tmp, *d_denA, *d_denB;

  char outfile[30];
  char outLab[30];

  size_t meshSize, primSize, qSize, flxSize;

  double initT, endT;

  //MPI vars
  int mpi_err;
  int *counts, *dspls;

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

  int nanCnt[4];

  mpi_err=MPI_Init(argc,argv);

  if(mpi_err!=MPI_SUCCESS){
    printf("Error initializing MPI\n");
  }

  Hp=Hyp;
  Ha=Hya;

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
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  //if(rank==0)printf("Before loc size calcs\n");
  counts=(int *)malloc(size*sizeof(int));
  dspls=(int *)malloc(size*sizeof(int));
  dspls[0]=0;
  counts[0]=Hp->ny/size;
  if((Hp->ny%size)!=0)counts[0]++;
  counts[0]*=Hp->nx;
  for(i=1;i<size;i++){
    counts[i]=Hp->ny/size;
    if((Hp->ny%size)>i)counts[i]++;
    counts[i]*=Hp->nx;
    dspls[i]=dspls[i-1]+counts[i-1];
  }
  myNy=counts[rank]/Hp->nx;
  if(rank==0){
    pProc=MPI_PROC_NULL;
    bndT=Hp->bndU;
  }else{
    pProc=rank-1;
    bndT=BND_INT;
  }
  if(rank==size-1){
    nProc=MPI_PROC_NULL;
    bndB=Hp->bndD;
  }else{
    nProc=rank+1;
    bndB=BND_INT;
  }
  
  device_init(Hp,Ha,myNy);

  //Calculate arraysizes
  varSize=(Hp->nx+4)*(myNy+4);
  meshSize=Hp->nvar*(Hp->nx+4)*(myNy+4);
  if(myNy>=Hp->nx){
    primSize=Hp->nvar*(Hp->nx+4)*myNy;
    qSize   =Hp->nvar*(Hp->nx+2)*myNy;
    flxSize =Hp->nvar*(Hp->nx+1)*myNy;
  }else{
    primSize=Hp->nvar*(myNy+4)*Hp->nx;
    qSize   =Hp->nvar*(myNy+2)*Hp->nx;
    flxSize =Hp->nvar*(myNy+1)*Hp->nx;
  }
  //if(rank==0)printf("Done loc size calcs\n");

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
  recvMesh=(double*)malloc(Hp->nvar*Hp->nx*myNy*sizeof(double));
  lMesh =(double*)malloc(Hp->nvar*varSize*sizeof(double));
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

  //if(rank==0)printf("Arrays allocated\n");

  //Zero lMesh
  for(i=0;i<Hp->nvar*varSize;i++){
    lMesh[i]=0.0;
  }

  //distribute over processors
  for(nV=0;nV<Hp->nvar;nV++){
    mpi_err=MPI_Scatterv(gMesh+nV*Hp->nx*Hp->ny,counts,dspls,MPI_DOUBLE,recvMesh+nV*Hp->nx*myNy,Hp->nx*myNy,MPI_DOUBLE,0,MPI_COMM_WORLD);
    if(mpi_err!=MPI_SUCCESS){
      printf("Error scattering data to other processors\n");
    }
    for(lI=0;lI<Hp->nx*myNy;lI++){
      i=lI%(Hp->nx);
      j=lI/(Hp->nx);
      lMesh[i+2+(Hp->nx+4)*(j+2+(myNy+4)*nV)]=recvMesh[i+Hp->nx*(j+myNy*nV)];
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  //sprintf(outLab,"Loc mesh n %d",rank);
  //printArray(outLab,lMesh,4,Hp->nx+4,myNy+4);

  if(rank==0)printf("Initial conditions distributed\n");

  if(Ha->tend>0.0){
    nxttout=Ha->tend;
  }
  if(Ha->dtoutput>0.0&&nxttout>Ha->dtoutput){
    nxttout=Ha->dtoutput;
  }

  volCell=Hp->dx*Hp->dy;
  lTM=sumArray(lMesh,VARRHO,Hp->nx,myNy,2,2);
  lTE=sumArray(lMesh,VARPR ,Hp->nx,myNy,2,2);
  MPI_Allreduce(&lTM,&ogTM,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD); 
  MPI_Allreduce(&lTE,&ogTE,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#ifdef M_PREC_CMP
  frexp(ogTM,&M_exp);
  frexp(ogTE,&E_exp);

  M_prec=ldexp(DBL_EPSILON,M_exp-1);
  E_prec=ldexp(DBL_EPSILON,E_exp-1);

  if(rank==0)printf("TM:%g+-%g TE:%g+-%g\n",volCell*ogTM,volCell*M_prec,volCell*ogTE,volCell*E_prec);
#endif

  //Print initial condition
  for(nV=0;nV<Hp->nvar;nV++){
    for(lI=0;lI<Hp->nx*myNy;lI++){
      i=lI%(Hp->nx);
      j=lI/(Hp->nx);
      recvMesh[i+Hp->nx*(j+myNy*nV)]=lMesh[i+2+(Hp->nx+4)*(j+2+(myNy+4)*nV)];
    }
    MPI_Gatherv(recvMesh+nV*Hp->nx*myNy,Hp->nx*myNy,MPI_DOUBLE,gMesh+nV*Hp->nx*Hp->ny,counts,dspls,MPI_DOUBLE,0,MPI_COMM_WORLD);
  }
  if(rank==0){
    snprintf(outfile,29,"%s%05d",Ha->outPre,n);
    writeVis(outfile,gMesh,Hp->dx,Hp->dy,Hp->nvar,Hp->nx,Hp->ny);
    printf("INIT: TM: %g TE: %g\n",volCell*ogTM,volCell*ogTE);
  }
 
  cudaMemcpy(d_u,lMesh,meshSize*sizeof(double),cudaMemcpyHostToDevice);
 
  initT=MPI_Wtime();

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
    //printf("N[%2d] ITER %d l denom=%g\n",rank,n,dt_denom);
    MPI_Allreduce(&dt_denom,&gDenom,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    //printf("N[%2d] ITER %d g denom=%g\n",rank,n,gDenom);
    dt=0.5*Ha->sigma/gDenom;
    if(nxttout>0.0&&dt>(nxttout-cTime)){
      printf("Adjusting timestep from %g to %g for iter %d\n",dt,nxttout-cTime,n);
      dt=(nxttout-cTime);
    }
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
    if(n%Ha->nprtLine==0||((cTime>=nxttout&&nxttout>0)||(Ha->noutput>0&&n%Ha->noutput==0))){
      cudaMemcpy(lMesh,d_u,meshSize*sizeof(double),cudaMemcpyDeviceToHost);
      if(n%Ha->nprtLine==0){
        //nanScan(nanCnt,lMesh,Hp->nvar,Hp->nx,myNy,2,2);
        lTM=sumArray(lMesh,VARRHO,Hp->nx,myNy,2,2);
        lTE=sumArray(lMesh,VARPR ,Hp->nx,myNy,2,2);
        MPI_Allreduce(&lTM,&gTM,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        MPI_Allreduce(&lTE,&gTE,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        //printf("N[%2d] Nans: %d %d %d %d\n",rank,nanCnt[0],nanCnt[1],nanCnt[2],nanCnt[3]);
        if(rank==0){
          printf("Iter %05d time %f dt %g TM: %g TE: %g\n",n,cTime,dt,volCell*gTM,volCell*gTE); 
          if(0&&(fabs(gTM-ogTM)>0.0||fabs(gTE-ogTE)>0.0)){
            printf("ERR(%5s): Mass %g Ene %g\n","RAW",volCell*(gTM-ogTM),volCell*(gTE-ogTE));
            printf("ERR(%5s): Mass %g Ene %g\n","%",100.0*(gTM-ogTM)/ogTM,100.0*(gTE-ogTE)/ogTE);
#ifdef M_PREC_CMP
            printf("ERR(%5s): Mass %g Ene %g\n","M PRE",(TM-oTM)/M_prec,(TE-oTE)/E_prec);
#endif
          }
        }
      }
      if((cTime>=nxttout&&nxttout>0)||(Ha->noutput>0&&n%Ha->noutput==0)){
        if(rank==0)printf("Vis file @ time %f iter %d\n",cTime,n);
        if(cTime>=nxttout&&nxttout>0.0){
          nxttout+=Ha->dtoutput;
          if(nxttout>Ha->tend&&Ha->tend>0)nxttout=Ha->tend;
          //if(rank==0)printf("Next Vis Time: %f\n",nxttout);
        }
        for(nV=0;nV<Hp->nvar;nV++){
          for(lI=0;lI<Hp->nx*myNy;lI++){
            i=lI%(Hp->nx);
            j=lI/(Hp->nx);
            recvMesh[i+Hp->nx*(j+myNy*nV)]=lMesh[i+2+(Hp->nx+4)*(j+2+(myNy+4)*nV)];
          }
          MPI_Gatherv(recvMesh+nV*Hp->nx*myNy,Hp->nx*myNy,MPI_DOUBLE,gMesh+nV*Hp->nx*Hp->ny,counts,dspls,MPI_DOUBLE,0,MPI_COMM_WORLD);
        }
        if(rank==0){
          snprintf(outfile,29,"%s%05d",Ha->outPre,n);
          writeVis(outfile,gMesh,Hp->dx,Hp->dy,Hp->nvar,Hp->nx,Hp->ny);
        }
      }
    }
  }
  if(rank==0)printf("time: %f, %d iters run\n",cTime,n);

  endT=MPI_Wtime();

  if(rank==0){
    printf("TFMT:%s,%s,%s,%s.%s,%s,%s,%s\n","cType","mType","init","nproc","nth","niters","ncells","wRunt");
    printf("TIME:%s,%s,%s,%d,%d,%d,%d,%g\n","\"MPI/CUDA\"","\"GPU:?\"","\"Init\"",size,1,n,Hp->nx*Hp->ny,(endT-initT));
  }
  
  cudaMemcpy(lMesh,d_u,meshSize*sizeof(double),cudaMemcpyDeviceToHost);

  //Print final condition
  for(nV=0;nV<Hp->nvar;nV++){
    for(lI=0;lI<Hp->nx*myNy;lI++){
      i=lI%(Hp->nx);
      j=lI/(Hp->nx);
      recvMesh[i+Hp->nx*(j+myNy*nV)]=lMesh[i+2+(Hp->nx+4)*(j+2+(myNy+4)*nV)];
    }
    MPI_Gatherv(recvMesh+nV*Hp->nx*myNy,Hp->nx*myNy,MPI_DOUBLE,gMesh+nV*Hp->nx*Hp->ny,counts,dspls,MPI_DOUBLE,0,MPI_COMM_WORLD);
  }
  if(rank==0){
    snprintf(outfile,29,"%s%05d",Ha->outPre,n);
    writeVis(outfile,gMesh,Hp->dx,Hp->dy,Hp->nvar,Hp->nx,Hp->ny);
  }
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

  printf("NODE %d: Finalizing MPI\n",rank);
  MPI_Finalize();

  printf("N %d: Returning from engine\n",rank);
}
