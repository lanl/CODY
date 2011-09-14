#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>

#include "engine.h"
#include "hydro_utils.h"
#include "hydro_macros.h"
#include "device_funcs.h"
#include "hydro_dmp.h"

#define FN_LEN 50

#define CDT_REGS 11
#define STEP_REGS 16

__global__ void halfX(double *upx, double *u, double *bnd, double dt);
__global__ void halfY(double *upy, double *u, double *bnd, double dt);
__global__ void cmp_flxX(double *flx, double *u, double *bnd, double *upx, double dtdx);
__global__ void cmp_flxY(double *flx, double *u, double *bnd, double *upy, double dtdy);
__global__ void cmp_dStX(double *du, double *flx);
__global__ void cmp_dStY(double *du, double *flx);
__global__ void add_dSt(double *u, double *du);

void engine(hydro_args Ha){
    // Host problem state vars
    hydro_prob Hp;
    double *hst_uold;
    // Host vars
    int nstep;
    double dt, dt_denom;
    int nDen, redBlocks;
    double *tmp;
    // Device vars
    double *dev_uold, *dev_du;
    double *dev_up, *dev_flx;
    double *dev_denA, *dev_denB;
    double *dev_bnd;
    // Cuda vars
    int dev;
    cudaDeviceProp prop;
    cudaError_t cuErrVar;
    size_t meshVarSize, bndVarSize, hStepVarSize;
    int nxy;
    int nTh, rpBl;
    int mxTh, thWp;
    int nThCDT, nThStep;
    size_t shMpBl;
    int nBlockM;

    // I/O vars
    double lastTOut;
    char filename[FN_LEN];

    //debug printing vars
    //int i,j,k;
    size_t mem_avail,mem_used;

    //initialize mesh/device
    hydro_init(&hst_uold,&Ha,&Hp);
    device_init(&Ha,&Hp);

    //Device settings
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
    nThStep=mxTh;//thread upper bound
    nThStep=MIN(nThStep,thWp*(rpBl/(STEP_REGS*thWp)));//register upper bound
    nTh=nThCDT;
    nBlockM=((Hp.ny*Hp.nx)+nTh-1)/nTh;
    printf("Per block: Max threads %d, regs %d\n", mxTh,rpBl, shMpBl);
    printf("Block size lims: cdt %d step %d\n",nThCDT,nThStep);

    //calculate size of mesh vars
    nxy=(Hp.nx>Hp.ny)?Hp.nx:Hp.ny;
    bndVarSize=Hp.nvar*nxy*4*sizeof(double);
    meshVarSize=Hp.nvar*Hp.ny*Hp.nx*sizeof(double);
    hStepVarSize=Hp.nvar*(Hp.nx+1)*(Hp.ny+1)*sizeof(double);

    //print relative size of global memory
    mem_used=2*meshVarSize+bndVarSize+hStepVarSize+2*nBlockM*sizeof(double);
    mem_avail=prop.totalGlobalMem;
    printf("%u/%u or %f%% memory used for a %d var %dx%d mesh\n",mem_used, mem_avail, (double)mem_used/(double)mem_avail*100.0, Hp.nvar, Hp.nx, Hp.ny);

    //allocate device vars
    cudaMalloc(&dev_uold,meshVarSize);
    HANDLE_CUDA_ERROR(cuErrVar);
    cudaMalloc(&dev_bnd,bndVarSize);
    HANDLE_CUDA_ERROR(cuErrVar);
    cudaMalloc(&dev_up,hStepVarSize);
    HANDLE_CUDA_ERROR(cuErrVar);
    cudaMalloc(&dev_flx,hStepVarSize);
    HANDLE_CUDA_ERROR(cuErrVar);
    cudaMalloc(&dev_du,meshVarSize);
    HANDLE_CUDA_ERROR(cuErrVar);
    cudaMalloc(&dev_denA,nBlockM*sizeof(double));
    HANDLE_CUDA_ERROR(cuErrVar);
    cudaMalloc(&dev_denB,nBlockM*sizeof(double));
    HANDLE_CUDA_ERROR(cuErrVar);

    //Set initial state on device
    cudaMemcpy(dev_uold,hst_uold,meshVarSize,cudaMemcpyHostToDevice);
    HANDLE_CUDA_ERROR(cuErrVar);

    //initialize main loop/output vars 
    nstep=0;
    lastTOut=Hp.t;
    //Write init conditions
    sprintf(filename,"%s-lw-i%07d.dmp",Ha.dmpPre,nstep);
    wrt_dmp(filename,&Hp,hst_uold);
    //main loop
    while((nstep<Ha.nstepmax)&&((Hp.t<Ha.tend)||(Ha.tend<0))){
      //calculate timestep
      dt=0.0;
      nDen=nBlockM;
      redBlocks=nDen;
      nTh=nThCDT;
      calc_denom<<<BL_TH(Hp.nx*Hp.ny,nTh)>>>(dev_uold, dev_denA);
      HANDLE_CUDA_ERROR(cuErrVar);
      while(redBlocks>1){
        redBlocks=(nDen+2*nTh-1)/(2*nTh);
        redu_max<<<redBlocks,nTh>>>(dev_denA,dev_denB,nDen);
        HANDLE_CUDA_ERROR(cuErrVar);
        nDen=redBlocks;
        tmp=dev_denA;
        dev_denA=dev_denB;
        dev_denB=tmp;
      }
      cudaMemcpy(&dt_denom,dev_denA,sizeof(double),cudaMemcpyDeviceToHost);
      HANDLE_CUDA_ERROR(cuErrVar);//invalid argument error in above line
      dt=Ha.sigma/dt_denom;
      //Shift dt to account for output times or tend
      if((lastTOut+Ha.dtoutput<Hp.t+dt)&&Ha.dtoutput>0.0)
        dt=lastTOut+Ha.dtoutput-Hp.t;
      if((Ha.tend<Hp.t+dt)&&Ha.tend>0.0) dt=Ha.tend-Hp.t;
      //Step
      nTh=nThStep;
      cudaMemset(dev_up,0,hStepVarSize);
      cudaMemset(dev_flx,0,hStepVarSize);
      cudaMemset(dev_du,0,meshVarSize);//zero change in state var
      gen_bndX<<<BL_TH(Hp.ny,nTh)>>>(dev_uold,dev_bnd,Hp.bndL,Hp.bndR);
      halfX   <<<BL_TH((Hp.nx+1)*Hp.ny,nTh)>>>(dev_up,dev_uold,dev_bnd,dt);
      cmp_flxX<<<BL_TH((Hp.nx+1)*Hp.ny,nTh)>>>(dev_flx,dev_uold,dev_bnd,dev_up,dt/Hp.dx);
      cmp_dStX<<<BL_TH((Hp.nx  )*Hp.ny,nTh)>>>(dev_du,dev_flx);

      gen_bndY<<<BL_TH(Hp.nx,nTh)>>>(dev_uold,dev_bnd,Hp.bndL,Hp.bndR);
      halfY   <<<BL_TH((Hp.ny+1)*Hp.nx,nTh)>>>(dev_up,dev_uold,dev_bnd,dt);
      cmp_flxY<<<BL_TH((Hp.ny+1)*Hp.nx,nTh)>>>(dev_flx,dev_uold,dev_bnd,dev_up,dt/Hp.dy);
      cmp_dStY<<<BL_TH((Hp.nx  )*Hp.ny,nTh)>>>(dev_du,dev_flx);

      add_dSt <<<BL_TH(Hp.nx*Hp.ny,nTh)>>>(dev_uold,dev_du);
      HANDLE_CUDA_ERROR(cuErrVar);
      //Finish main loop
      ++nstep;
      Hp.t+=dt;
      if((nstep%Ha.nprtLine)==0){
        fprintf(stdout,"Run Iter=%6d t=%12.5g dt=%12g\n", nstep, Hp.t, dt);
      }
      //Output
      filename[0]='\0';
      if(Ha.dtoutput<=0){
        if(nstep%Ha.noutput==0){
          sprintf(filename,"%s-lw-i%07d.dmp",Ha.dmpPre,nstep);
        }
      }else{
        if(Hp.t>=lastTOut+Ha.dtoutput){
          sprintf(filename,"%s-lw-i%07d.dmp",Ha.dmpPre,nstep);
          lastTOut=Hp.t;
        }
      }
      if(filename[0]){
        //Fetch data from device
        cudaMemcpy(hst_uold,dev_uold,meshVarSize,cudaMemcpyDeviceToHost);
        HANDLE_CUDA_ERROR(cuErrVar);
        //Write dump file iter line
        fprintf(stdout,"Dmp Iter=%6d t=%12.5g dt=%12g\n", nstep, Hp.t, dt);
        //Write dump file
        wrt_dmp(filename,&Hp,hst_uold);
      }
    }
    //Fetch state from device
    cudaMemcpy(hst_uold,dev_uold,meshVarSize,cudaMemcpyDeviceToHost);
    HANDLE_CUDA_ERROR(cuErrVar);
    //Write final iter line
    fprintf(stdout,"End Iter=%6d t=%12.5g dt=%12g\n", nstep, Hp.t, dt);
    //Write dump file w/ final state
    sprintf(filename,"%s-lw-i%07d.dmp",Ha.dmpPre,nstep);
    wrt_dmp(filename,&Hp,hst_uold);

    //Free cuda vars
    cudaFree(dev_uold);
    cudaFree(dev_bnd);
    cudaFree(dev_flx);
    cudaFree(dev_du);
    cudaFree(dev_denA);
    cudaFree(dev_denB);
    //Finalize
    free(hst_uold);
}
