#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>

#include "engine.h"
#include "hydro_utils.h"
#include "hydro_macros.h"
#include "device_funcs.h"
#include "hydro_dmp.h"

#define FN_LEN 50

#define CDT_REGS 16
#define STEP_REGS 16

__global__ void primX(double *q, double *u, double *bnd);
__global__ void primY(double *q, double *u, double *bnd);
__global__ void trace(double *qxm, double *qxp, double *q, int np, int nt, double dtdx);
__global__ void riemann(double *flx, double *qxm, double *qxp, int np, int nt);
__global__ void cmp_dStX(double *du, double *qgdnv, double dtdx);
__global__ void cmp_dStY(double *du, double *qgdnv, double dtdx);
__global__ void add_dSt(double *u, double *du);

#ifdef FINE_RES_TIMER
#define TIMER_INIT(eI,eF) cudaEventCreate(&(eI));cudaEventCreate(&(eF));cudaEventRecord(eI, 0);HANDLE_CUDA_ERROR(cuErrVar);
#define TIMER_FINAL(rec,eI,eF) cudaEventRecord((eF), 0);cudaEventSynchronize(eF);cudaEventElapsedTime(&(rec),(eI),(eF)); cudaEventDestroy(eI);cudaEventDestroy(eF);HANDLE_CUDA_ERROR(cuErrVar);
#define TIMER_CALC(c) c
#else
#define TIMER_INIT(eI,eF)
#define TIMER_FINAL(rec,eI,eF) HANDLE_CUDA_ERROR(cuErrVar);
#define TIMER_CALC(c)
#endif //FINE_RES_TIMER

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
    double *dev_uold, *dev_q;
    double *dev_qxm, *dev_qxp;
    double *dev_flx;
    double *dev_denA, *dev_denB;
    double *dev_bnd;
    // Cuda vars
    int dev;
    cudaDeviceProp prop;
    cudaError_t cuErrVar;
    size_t meshVarSize, bndVarSize, wkVarSize;
    int nxy;
    int nTh, rpBl;
    int mxTh, thWp;
    int nThCDT, nThStep;
    size_t shMpBl;
    int nBlockM;

    // I/O vars
    double lastTOut;
    char filename[FN_LEN];

    //Timing vars
    cudaEvent_t start, end;
    float t_rec;
    double rTime;
#ifdef FINE_RES_TIMER
    cudaEvent_t timeA, timeB;
    double t_transfer;
    double t_CD, t_red;
    double t_bndX, t_priX, t_traX, t_rieX, t_dStX;
    double t_bndY, t_priY, t_traY, t_rieY, t_dStY;
#endif //FINE_RES_TIMER

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
    nThStep=thWp;//thread upper bound
    nThStep=MIN(nThStep,thWp*(rpBl/(STEP_REGS*thWp)));//register upper bound
    nTh=nThCDT;
    nBlockM=((Hp.ny*Hp.nx)+nTh-1)/nTh;
    printf("Per block: Max threads %d, regs %d\n", mxTh,rpBl, shMpBl);
    printf("Block size lims: cdt %d step %d\n",nThCDT,nThStep);
    dim3 prKX(BL(Hp.nx+4,nThStep),Hp.ny), prKY(BL(Hp.nx  ,nThStep),Hp.ny+4);
    dim3 trKX(BL(Hp.nx+2,nThStep),Hp.ny), trKY(BL(Hp.nx+2,nThStep),Hp.nx);
    dim3 riKX(BL(Hp.nx+1,nThStep),Hp.ny), riKY(BL(Hp.nx+1,nThStep),Hp.nx);
    dim3 mesh(BL(Hp.nx  ,nThStep),Hp.ny);

    //calculate size of mesh vars
    nxy=(Hp.nx>Hp.ny)?Hp.nx:Hp.ny;
    meshVarSize=Hp.nvar*Hp.ny*Hp.nx*sizeof(double);
    wkVarSize=Hp.nvar*(Hp.ny+4)*(Hp.nx+4)*sizeof(double);
    bndVarSize=Hp.nvar*nxy*4*sizeof(double);

    //print relative size of global memory
    mem_used=meshVarSize+bndVarSize+4*wkVarSize+2*nBlockM*sizeof(double);
    mem_avail=prop.totalGlobalMem;
    printf("%u/%u or %f%% memory used for a %d var %dx%d mesh\n",mem_used, mem_avail, (double)mem_used/(double)mem_avail*100.0, Hp.nvar, Hp.nx, Hp.ny);
    if(mem_used>mem_avail){
      fprintf(stderr,"Not enough memory on GPU\n");
      exit(0);
    }

    //allocate device vars
    cudaMalloc(&dev_uold,meshVarSize);
    HANDLE_CUDA_ERROR(cuErrVar);
    cudaMalloc(&dev_q,wkVarSize);
    HANDLE_CUDA_ERROR(cuErrVar);
    cudaMalloc(&dev_qxm,wkVarSize);
    HANDLE_CUDA_ERROR(cuErrVar);
    cudaMalloc(&dev_qxp,wkVarSize);
    HANDLE_CUDA_ERROR(cuErrVar);
    cudaMalloc(&dev_flx,wkVarSize);
    HANDLE_CUDA_ERROR(cuErrVar);
    cudaMalloc(&dev_bnd,bndVarSize);
    HANDLE_CUDA_ERROR(cuErrVar);
    cudaMalloc(&dev_denA,nBlockM*sizeof(double));
    HANDLE_CUDA_ERROR(cuErrVar);
    cudaMalloc(&dev_denB,nBlockM*sizeof(double));
    HANDLE_CUDA_ERROR(cuErrVar);


    //initialize main loop/output vars 
    nstep=0;
    lastTOut=Hp.t;

    //Write init conditions
    sprintf(filename,"%s-god-i%07d.dmp",Ha.dmpPre,nstep);
    wrt_dmp(filename,&Hp,hst_uold);

    //Init timer vars
#ifdef FINE_RES_TIMER
    t_transfer=0.0;
    t_CD=0.0;t_red=0.0;
    t_bndX=0.0;t_priX=0.0;t_traX=0.0;t_rieX=0.0; t_dStX=0.0;
    t_bndY=0.0;t_priY=0.0;t_traY=0.0;t_rieY=0.0; t_dStY=0.0;
#endif //FINE_RES_TIMER

    //begin timer
    cudaEventCreate(&start);
    cudaEventCreate(&end);
    cudaEventRecord(start,0);

    //Set initial state on device
    TIMER_INIT(timeA,timeB);
    cudaMemcpy(dev_uold,hst_uold,meshVarSize,cudaMemcpyHostToDevice);
    TIMER_FINAL(t_rec,timeA,timeB);
    TIMER_CALC(t_transfer+=t_rec);
    HANDLE_CUDA_ERROR(cuErrVar);

    //main loop
    while((nstep<Ha.nstepmax)&&((Hp.t<Ha.tend)||(Ha.tend<0))){
      //calculate timestep
      if((nstep%2)==0){
        dt=0.0;
        nDen=nBlockM;
        redBlocks=nDen;
        nTh=nThCDT;
        TIMER_INIT(timeA,timeB);
        calc_denom<<<nBlockM,nTh,nTh*sizeof(double)>>>(dev_uold, dev_denA);
        TIMER_FINAL(t_rec,timeA,timeB);
        TIMER_CALC(t_CD+=t_rec);
        while(redBlocks>1){
          redBlocks=(nDen+2*nTh-1)/(2*nTh);
          TIMER_INIT(timeA,timeB);
          redu_max<<<redBlocks,nTh,nTh*sizeof(double)>>>(dev_denA,dev_denB,nDen);
          TIMER_FINAL(t_rec,timeA,timeB);
          TIMER_CALC(t_red+=t_rec);
          HANDLE_CUDA_ERROR(cuErrVar);
          nDen=redBlocks;
          tmp=dev_denA;
          dev_denA=dev_denB;
          dev_denB=tmp;
        }
        cudaMemcpy(&dt_denom,dev_denA,sizeof(double),cudaMemcpyDeviceToHost);
        HANDLE_CUDA_ERROR(cuErrVar);
        dt=0.5*Ha.sigma/dt_denom;
        //Shift dt to account for output times or tend
        if((lastTOut+Ha.dtoutput<Hp.t+dt)&&Ha.dtoutput>0.0)
          dt=lastTOut+Ha.dtoutput-Hp.t;
        else if((lastTOut+Ha.dtoutput<Hp.t+2.0*dt)&&Ha.dtoutput>0.0)
          dt=0.5*(lastTOut+Ha.dtoutput-Hp.t);
        if((Ha.tend<Hp.t+dt)&&Ha.tend>0.0) dt=Ha.tend-Hp.t;
        else if((Ha.tend<Hp.t+2.0*dt)&&Ha.tend>0.0) dt=0.5*(Ha.tend-Hp.t);
      }
      //Step
      nTh=nThStep;
      if(nstep%2==0){
        //x pass
        TIMER_INIT(timeA,timeB);
        gen_bndX<<<BL_TH(Hp.ny,nTh)>>>(dev_uold,dev_bnd,Hp.bndL,Hp.bndR);
        TIMER_FINAL(t_rec,timeA,timeB);
        TIMER_CALC(t_bndX+=t_rec);
        TIMER_INIT(timeA,timeB);
        primX   <<<prKX,nTh>>>(dev_q,dev_uold,dev_bnd);
        TIMER_FINAL(t_rec,timeA,timeB);
        TIMER_CALC(t_priX+=t_rec);
        TIMER_INIT(timeA,timeB);
        trace   <<<trKX,nTh>>>(dev_qxm, dev_qxp, dev_q, Hp.nx+2, Hp.ny, dt/Hp.dx);
        TIMER_FINAL(t_rec,timeA,timeB);
        TIMER_CALC(t_traX+=t_rec);
        TIMER_INIT(timeA,timeB);
        riemann <<<riKX,nTh>>>(dev_flx, dev_qxm, dev_qxp, Hp.nx+1, Hp.ny);
        TIMER_FINAL(t_rec,timeA,timeB);
        TIMER_CALC(t_rieX+=t_rec);
        TIMER_INIT(timeA,timeB);
        cmp_dStX<<<mesh,nTh>>>(dev_uold, dev_flx, dt/Hp.dx);
        TIMER_FINAL(t_rec,timeA,timeB);
        TIMER_CALC(t_dStX+=t_rec);
        HANDLE_CUDA_ERROR(cuErrVar);
        //y pass
        TIMER_INIT(timeA,timeB);
        gen_bndY<<<BL_TH(Hp.nx,nTh)>>>(dev_uold,dev_bnd,Hp.bndL,Hp.bndR);
        TIMER_FINAL(t_rec,timeA,timeB);
        TIMER_CALC(t_bndY+=t_rec);
        TIMER_INIT(timeA,timeB);
        primY   <<<prKY,nTh>>>(dev_q,dev_uold,dev_bnd);
        TIMER_FINAL(t_rec,timeA,timeB);
        TIMER_CALC(t_priY+=t_rec);
        TIMER_INIT(timeA,timeB);
        trace   <<<trKY,nTh>>>(dev_qxm, dev_qxp, dev_q, Hp.ny+2, Hp.nx, dt/Hp.dy);
        TIMER_FINAL(t_rec,timeA,timeB);
        TIMER_CALC(t_traY+=t_rec);
        TIMER_INIT(timeA,timeB);
        riemann <<<riKY,nTh>>>(dev_flx, dev_qxm, dev_qxp, Hp.ny+1, Hp.nx);
        TIMER_FINAL(t_rec,timeA,timeB);
        TIMER_CALC(t_rieY+=t_rec);
        TIMER_INIT(timeA,timeB);
        cmp_dStY<<<mesh,nTh>>>(dev_uold, dev_flx, dt/Hp.dy);
        TIMER_FINAL(t_rec,timeA,timeB);
        TIMER_CALC(t_dStY+=t_rec);
        TIMER_INIT(timeA,timeB);
      }else{
        //y pass
        TIMER_INIT(timeA,timeB);
        gen_bndY<<<BL_TH(Hp.nx,nTh)>>>(dev_uold,dev_bnd,Hp.bndL,Hp.bndR);
        TIMER_FINAL(t_rec,timeA,timeB);
        TIMER_CALC(t_bndY+=t_rec);
        TIMER_INIT(timeA,timeB);
        primY   <<<prKY,nTh>>>(dev_q,dev_uold,dev_bnd);
        TIMER_FINAL(t_rec,timeA,timeB);
        TIMER_CALC(t_priY+=t_rec);
        TIMER_INIT(timeA,timeB);
        trace   <<<trKY,nTh>>>(dev_qxm, dev_qxp, dev_q, Hp.ny+2, Hp.nx, dt/Hp.dy);
        TIMER_FINAL(t_rec,timeA,timeB);
        TIMER_CALC(t_traY+=t_rec);
        TIMER_INIT(timeA,timeB);
        riemann <<<riKY,nTh>>>(dev_flx, dev_qxm, dev_qxp, Hp.ny+1, Hp.nx);
        TIMER_FINAL(t_rec,timeA,timeB);
        TIMER_CALC(t_rieY+=t_rec);
        TIMER_INIT(timeA,timeB);
        cmp_dStY<<<mesh,nTh>>>(dev_uold, dev_flx, dt/Hp.dy);
        TIMER_FINAL(t_rec,timeA,timeB);
        TIMER_CALC(t_dStY+=t_rec);
        TIMER_INIT(timeA,timeB);
        //x pass
        TIMER_INIT(timeA,timeB);
        gen_bndX<<<BL_TH(Hp.ny,nTh)>>>(dev_uold,dev_bnd,Hp.bndL,Hp.bndR);
        TIMER_FINAL(t_rec,timeA,timeB);
        TIMER_CALC(t_bndX+=t_rec);
        TIMER_INIT(timeA,timeB);
        primX   <<<prKX,nTh>>>(dev_q,dev_uold,dev_bnd);
        TIMER_FINAL(t_rec,timeA,timeB);
        TIMER_CALC(t_priX+=t_rec);
        TIMER_INIT(timeA,timeB);
        trace   <<<trKX,nTh>>>(dev_qxm, dev_qxp, dev_q, Hp.nx+2, Hp.ny, dt/Hp.dx);
        TIMER_FINAL(t_rec,timeA,timeB);
        TIMER_CALC(t_traX+=t_rec);
        TIMER_INIT(timeA,timeB);
        riemann <<<riKX,nTh>>>(dev_flx, dev_qxm, dev_qxp, Hp.nx+1, Hp.ny);
        TIMER_FINAL(t_rec,timeA,timeB);
        TIMER_CALC(t_rieX+=t_rec);
        TIMER_INIT(timeA,timeB);
        cmp_dStX<<<mesh,nTh>>>(dev_uold, dev_flx, dt/Hp.dx);
        TIMER_FINAL(t_rec,timeA,timeB);
        TIMER_CALC(t_dStX+=t_rec);
      }
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
          sprintf(filename,"%s-god-i%07d.dmp",Ha.dmpPre,nstep);
        }
      }else{
        if(Hp.t>=lastTOut+Ha.dtoutput){
          sprintf(filename,"%s-god-i%07d.dmp",Ha.dmpPre,nstep);
          lastTOut=Hp.t;
        }
      }
      if(filename[0]){
        //Fetch data from device
        TIMER_INIT(timeA,timeB);
        cudaMemcpy(hst_uold,dev_uold,meshVarSize,cudaMemcpyDeviceToHost);
        TIMER_FINAL(t_rec,timeA,timeB);
        TIMER_CALC(t_transfer+=t_rec);
        HANDLE_CUDA_ERROR(cuErrVar);
        //Write dump file iter line
        fprintf(stdout,"Dmp Iter=%6d t=%12.5g dt=%12g\n", nstep, Hp.t, dt);
        //Write dump file
        wrt_dmp(filename,&Hp,hst_uold);
      }
    }
    //Fetch state from device
    TIMER_INIT(timeA,timeB);
    cudaMemcpy(hst_uold,dev_uold,meshVarSize,cudaMemcpyDeviceToHost);
    TIMER_FINAL(t_rec,timeA,timeB);
    TIMER_CALC(t_transfer+=t_rec);
    HANDLE_CUDA_ERROR(cuErrVar);
    //Get end time and calc runtime
    cudaEventRecord(end,0);
    HANDLE_CUDA_ERROR(cuErrVar);
    cudaEventSynchronize(end);
    HANDLE_CUDA_ERROR(cuErrVar);
    cudaEventElapsedTime(&t_rec,start,end);
    rTime=t_rec/1000.0;
    cudaEventDestroy(start);
    cudaEventDestroy(end);
    //Write final iter line
    fprintf(stdout,"End Iter=%6d t=%12.5g dt=%12g\n", nstep, Hp.t, dt);
    //Write dump file w/ final state
    sprintf(filename,"%s-god-i%07d.dmp",Ha.dmpPre,nstep);
    wrt_dmp(filename,&Hp,hst_uold);

    //Print run end info
    printf("Run completed successfully\n");
    printf("Computation took %lf seconds\n", rTime);
#ifdef FINE_RES_TIMER
    printf("%15s\t%lf ms(%lf%%)\n","transfers", t_transfer,t_transfer/(10.0*rTime));
    printf("%15s\t%lf ms(%lf%%)\n","calc_denom",t_CD,t_CD/(10.0*rTime));
    printf("%15s\t%lf ms(%lf%%)\n","reduce_max",t_red,t_red/(10.0*rTime));
    printf("%15s\t%lf ms(%lf%%)\n","gen_bndX",  t_bndX,t_bndX/(10.0*rTime));
    printf("%15s\t%lf ms(%lf%%)\n","gen_bndY",  t_bndY,t_bndY/(10.0*rTime));
    printf("%15s\t%lf ms(%lf%%)\n","primX",     t_priX,t_priX/(10.0*rTime));
    printf("%15s\t%lf ms(%lf%%)\n","primY",     t_priY,t_priY/(10.0*rTime));
    printf("%15s\t%lf ms(%lf%%)\n","trace X",   t_traX,t_traX/(10.0*rTime));
    printf("%15s\t%lf ms(%lf%%)\n","trace Y",   t_traY,t_traY/(10.0*rTime));
    printf("%15s\t%lf ms(%lf%%)\n","riemann X", t_rieX,t_rieX/(10.0*rTime));
    printf("%15s\t%lf ms(%lf%%)\n","riemann Y", t_rieY,t_rieY/(10.0*rTime));
    printf("%15s\t%lf ms(%lf%%)\n","cmp_dStX",  t_dStX,t_dStX/(10.0*rTime));
    printf("%15s\t%lf ms(%lf%%)\n","cmp_dStY",  t_dStY,t_dStY/(10.0*rTime));
#endif //FINE_RES_TIMER
    //Free cuda vars
    cudaFree(dev_uold);
    cudaFree(dev_q);
    cudaFree(dev_qxm);
    cudaFree(dev_qxp);
    cudaFree(dev_flx);
    cudaFree(dev_bnd);
    cudaFree(dev_denA);
    cudaFree(dev_denB);
    //Finalize
    free(hst_uold);
}
