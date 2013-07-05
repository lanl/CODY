#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>

#include "dev_funcs.h"
#include "hydro_defs.h"

__device__ __constant__ int d_nx, d_ny, d_nvar;
__device__ __constant__ double d_gamma, d_dx, d_dy;
__device__ __constant__ double d_smallr, d_smallc, d_smallp;
__device__ __constant__ int d_niterR;

void print_cuda_err(char *file, const char *func, int line, cudaError_t err){
  fprintf(stderr,"CUDA error in %s: %s(%d):\n",file,func,line);
  fprintf(stderr,"Cuda error description:\n\t%s\n",cudaGetErrorString(err));
}

void handle_cuda_err(char *file, const char *func, int line, cudaError_t err){
  if(err!=cudaSuccess){
    print_cuda_err(file,func,line,err);
    exit(1);
  }
}

void device_init(hydro_prob *Hp, hydro_args *Ha){
  cudaError_t errVar;

  double smallp;

  smallp=Ha->smallc*Ha->smallc/Hp->gamma; 

  cudaMemcpyToSymbol(d_nx,&(Hp->nx),sizeof(int));
  HANDLE_CUDA_ERROR(errVar);
  cudaMemcpyToSymbol(d_ny,&(Hp->ny),sizeof(int));
  HANDLE_CUDA_ERROR(errVar);
  cudaMemcpyToSymbol(d_dx,&(Hp->dx),sizeof(double));
  HANDLE_CUDA_ERROR(errVar);
  cudaMemcpyToSymbol(d_dy,&(Hp->dy),sizeof(double));
  HANDLE_CUDA_ERROR(errVar);
  cudaMemcpyToSymbol(d_gamma,&(Hp->gamma),sizeof(double));
  HANDLE_CUDA_ERROR(errVar);
  cudaMemcpyToSymbol(d_smallr,&(Ha->smallr),sizeof(double));
  HANDLE_CUDA_ERROR(errVar);
  cudaMemcpyToSymbol(d_smallc,&(Ha->smallc),sizeof(double));
  HANDLE_CUDA_ERROR(errVar);
  cudaMemcpyToSymbol(d_smallp,&smallp,sizeof(double));
  HANDLE_CUDA_ERROR(errVar);
  cudaMemcpyToSymbol(d_niterR,&(Ha->niter_riemann),sizeof(int));
  HANDLE_CUDA_ERROR(errVar);
}

extern __shared__ double dynVar[];
__global__ void calc_denom(double *u, double *den){
double *denom=dynVar;
    uint stInd, i,j;
    int thInd;
    int mxInd, stride;
    double rho, vx, vy, eken, eint;
    double p, c;
    double courx, coury;

    //Calculate indicies
    thInd=threadIdx.x;
    stInd=threadIdx.x+blockDim.x*blockIdx.x;
    denom[thInd]=0.0;
    //Calculate denominator for my cell
    if(stInd<d_nx*d_ny){
      i=stInd%d_nx;
      j=stInd/d_nx;
      rho =u[(VARRHO*(d_ny+4)+j+2)*(d_nx+4)+i+2];
      vx  =u[(VARVX *(d_ny+4)+j+2)*(d_nx+4)+i+2];
      vy  =u[(VARVY *(d_ny+4)+j+2)*(d_nx+4)+i+2];
      eint=u[(VARPR *(d_ny+4)+j+2)*(d_nx+4)+i+2];
      rho=fmax(rho,d_smallr);
      vx=vx/rho;
      vy=vy/rho;
      eken=0.5*(vx*vx+vy*vy);
      eint=eint/rho-eken;

      p=fmax((d_gamma-(double)1.0)*rho*eint,rho*d_smallp);
      c=sqrt(d_gamma*p/rho);

      courx=(c+fabs(vx))/d_dx;
      coury=(c+fabs(vy))/d_dy;
      denom[thInd]=courx+coury;
    }
    //Reduce to maximum value of denomenator over whole mesh.
    __syncthreads();

    mxInd=blockDim.x;
    stride=mxInd;
    while(stride>0){
      stride>>=1;
      if(thInd+stride<blockDim.x){
        denom[thInd]=fmax(denom[thInd],denom[thInd+stride]);
      }
      __syncthreads();
    }
    //Have first thread write max denom to output array
    if(thInd==0)den[blockIdx.x]=denom[0];
}

__global__ void redu_max(double *arrIn, double *arrOut, int nVals){
    double *arr=dynVar;
    int arrInd, thInd, mxInd, stride;

    //Calculate indicies
    thInd=threadIdx.x;
    arrInd=thInd+2*blockDim.x*blockIdx.x;
    mxInd=blockDim.x;
    //Get max of two vals from input array
    stride=mxInd;
    arr[thInd]=0.0;
    if(arrInd<nVals){
      arr[thInd]=arrIn[arrInd];
    }
    if(arrInd+stride<nVals){
      arr[thInd]=fmax(arrIn[arrInd+stride],arr[thInd]);
    }
    //Find max over block
    while(stride>0){
      __syncthreads();
      stride>>=1;
      if(thInd+stride<mxInd){
        arr[thInd]=fmax(arr[thInd],arr[thInd+stride]);
      }
    }
    //Have first thread write max to out array
    if(thInd==0)arrOut[blockIdx.x]=arr[0];
}

__global__ void gen_bndXL(double *u, int bndT){
  int thInd=threadIdx.x+blockDim.x*blockIdx.x;
  int i=thInd%2;
  int j=thInd/2;
  int wInd, rInd;
  int vSize;

  wInd=  i+(d_nx+4)*(j+2);
  rInd=3-i+(d_nx+4)*(j+2);
  vSize=(d_nx+4)*(d_ny+4);
  if(j<d_ny){
    if(bndT==BND_REFL){
      u[wInd+VARRHO*vSize]= u[rInd+VARRHO*vSize];
      u[wInd+VARVX *vSize]=-u[rInd+VARVX *vSize];
      u[wInd+VARVY *vSize]= u[rInd+VARVY *vSize];
      u[wInd+VARPR *vSize]= u[rInd+VARPR *vSize];
    }else if(bndT==BND_PERM){
      u[wInd+VARRHO*vSize]= u[rInd+VARRHO*vSize];
      u[wInd+VARVX *vSize]= u[rInd+VARVX *vSize];
      u[wInd+VARVY *vSize]= u[rInd+VARVY *vSize];
      u[wInd+VARPR *vSize]= u[rInd+VARPR *vSize]; 
    }
  } 
}

__global__ void gen_bndXU(double *u, int bndT){
  int thInd=threadIdx.x+blockDim.x*blockIdx.x;
  int i=thInd%2;
  int j=thInd/2;
  int wInd, rInd;
  int vSize;

  wInd=d_nx+2+i+(d_nx+4)*(j+2);
  rInd=d_nx+1-i+(d_nx+4)*(j+2);
  vSize=(d_nx+4)*(d_ny+4);
  if(j<d_ny){
    if(bndT==BND_REFL){
      u[wInd+VARRHO*vSize]= u[rInd+VARRHO*vSize];
      u[wInd+VARVX *vSize]=-u[rInd+VARVX *vSize];
      u[wInd+VARVY *vSize]= u[rInd+VARVY *vSize];
      u[wInd+VARPR *vSize]= u[rInd+VARPR *vSize];
    }else if(bndT==BND_PERM){
      u[wInd+VARRHO*vSize]= u[rInd+VARRHO*vSize];
      u[wInd+VARVX *vSize]= u[rInd+VARVX *vSize];
      u[wInd+VARVY *vSize]= u[rInd+VARVY *vSize];
      u[wInd+VARPR *vSize]= u[rInd+VARPR *vSize]; 
    }
  } 
}

__global__ void gen_bndYL(double *u, int bndT){
  int thInd=threadIdx.x+blockDim.x*blockIdx.x;
  int i=thInd/2;
  int j=thInd%2;
  int wInd, rInd;
  int vSize;

  wInd=i+2+(d_nx+4)*(  j);
  rInd=i+2+(d_nx+4)*(3-j);
  vSize=(d_nx+4)*(d_ny+4);
  if(i<d_nx){
    if(bndT==BND_REFL){
      u[wInd+VARRHO*vSize]= u[rInd+VARRHO*vSize];
      u[wInd+VARVX *vSize]= u[rInd+VARVX *vSize];
      u[wInd+VARVY *vSize]=-u[rInd+VARVY *vSize];
      u[wInd+VARPR *vSize]= u[rInd+VARPR *vSize];
    }else if(bndT==BND_PERM){
      u[wInd+VARRHO*vSize]= u[rInd+VARRHO*vSize];
      u[wInd+VARVX *vSize]= u[rInd+VARVX *vSize];
      u[wInd+VARVY *vSize]= u[rInd+VARVY *vSize];
      u[wInd+VARPR *vSize]= u[rInd+VARPR *vSize]; 
    }
  } 
}

__global__ void gen_bndYU(double *u, int bndT){
  int thInd=threadIdx.x+blockDim.x*blockIdx.x;
  int i=thInd/2;
  int j=thInd%2;
  int wInd, rInd;
  int vSize;

  wInd=i+2+(d_nx+4)*(d_ny+2+j);
  rInd=i+2+(d_nx+4)*(d_ny+1-j);
  vSize=(d_nx+4)*(d_ny+4);
  if(i<d_nx){
    if(bndT==BND_REFL){
      u[wInd+VARRHO*vSize]= u[rInd+VARRHO*vSize];
      u[wInd+VARVX *vSize]= u[rInd+VARVX *vSize];
      u[wInd+VARVY *vSize]=-u[rInd+VARVY *vSize];
      u[wInd+VARPR *vSize]= u[rInd+VARPR *vSize];
    }else if(bndT==BND_PERM){
      u[wInd+VARRHO*vSize]= u[rInd+VARRHO*vSize];
      u[wInd+VARVX *vSize]= u[rInd+VARVX *vSize];
      u[wInd+VARVY *vSize]= u[rInd+VARVY *vSize];
      u[wInd+VARPR *vSize]= u[rInd+VARPR *vSize]; 
    }
  } 
}

__global__ void toPrimX(double *q, double *u){
  int thInd=threadIdx.x+blockDim.x*blockIdx.x;
  int i=thInd%(d_nx+4);
  int j=thInd/(d_nx+4);
  double r, vx, vy, eint, p;

  if(j<d_ny){
    r   =fmax(u[i+(d_nx+4)*(j+2+(d_ny+4)*VARRHO)],d_smallr);
    vx  =     u[i+(d_nx+4)*(j+2+(d_ny+4)*VARVX )]/r;
    vy  =     u[i+(d_nx+4)*(j+2+(d_ny+4)*VARVY )]/r;
    eint=     u[i+(d_nx+4)*(j+2+(d_ny+4)*VARPR )]-0.5*r*(vx*vx+vy*vy);
    p   =fmax((d_gamma-1)*r*eint,d_smallp);
    q[i+(d_nx+4)*(j+d_ny*VARRHO)]=r;
    q[i+(d_nx+4)*(j+d_ny*VARVX )]=vx;
    q[i+(d_nx+4)*(j+d_ny*VARVY )]=vy;
    q[i+(d_nx+4)*(j+d_ny*VARPR )]=p;
  }
}

__global__ void toPrimY(double *q, double *u){
  int thInd=threadIdx.x+blockDim.x*blockIdx.x;
  int i=thInd/(d_ny+4);
  int j=thInd%(d_ny+4);
  double r, vx, vy, eint, p;

  if(i<d_nx){
    r   =fmax(u[i+2+(d_nx+4)*(j+(d_ny+4)*VARRHO)],d_smallr);
    vx  =     u[i+2+(d_nx+4)*(j+(d_ny+4)*VARVX )]/r;
    vy  =     u[i+2+(d_nx+4)*(j+(d_ny+4)*VARVY )]/r;
    eint=     u[i+2+(d_nx+4)*(j+(d_ny+4)*VARPR )]-0.5*r*(vx*vx+vy*vy);
    p   =fmax((d_gamma-1)i*r*eint,d_smallp);
    q[j+(d_ny+4)*(i+d_nx*VARRHO)]=r;
    q[j+(d_ny+4)*(i+d_nx*VARVX )]=vy;
    q[j+(d_ny+4)*(i+d_nx*VARVY )]=vx;
    q[j+(d_ny+4)*(i+d_nx*VARPR )]=p;
  }
}

__device__ double slope(double *q, int ind){
  double dlft, drgt, dcen, dsgn, dlim;
  dlft=q[ind  ]-q[ind-1];
  drgt=q[ind+1]-q[ind  ];
  dcen=0.5*(dlft+drgt);
  dsgn=(dcen>=0)?1.0:-1.0;
  dlim=fmin(fabs(dlft),fabs(drgt));
  if(dlft*drgt<0)dlim=0.0;
  return dsgn*fmin(dlim,fabs(dcen));
}

__global__ void trace(double *ql, double *qr, double *q, double dtdx, int np, int nt){
  int thInd=threadIdx.x+blockDim.x*blockIdx.x;
  int i=thInd%(np+2);
  int j=thInd/(np+2);
  double  r,  u,  v1,  p;
  double dr, du, dv1, dp;
  double cc, csq;
  double alpham, alphap, alphazr;
  double spplus, spzero, spminus;
  double ap, am, azr, azv1;

  if(j<nt){
    r =q[i+1+(np+4)*(j+nt*VARRHO)];
    u =q[i+1+(np+4)*(j+nt*VARVX )];
    v1=q[i+1+(np+4)*(j+nt*VARVY )];
    p =q[i+1+(np+4)*(j+nt*VARPR )];

    csq=d_gamma*p/r;
    cc=sqrt(csq);

    dr =slope(q,i+1+(np+4)*(j+nt*VARRHO));
    du =slope(q,i+1+(np+4)*(j+nt*VARVX ));
    dv1=slope(q,i+1+(np+4)*(j+nt*VARVY ));
    dp =slope(q,i+1+(np+4)*(j+nt*VARPR ));
   

    alpham  = 0.5*(dp/(r*cc)-du)*r/cc;
    alphap  = 0.5*(dp/(r*cc)+du)*r/cc;
    alphazr = dr-dp/csq;

    //Right
    spminus=((u-cc)>=0.0)?0.0:(u-cc)*dtdx+1.0;
    spzero =((u   )>=0.0)?0.0:(u   )*dtdx+1.0;
    spplus =((u+cc)>=0.0)?0.0:(u+cc)*dtdx+1.0;
    ap  =-0.5*spplus *alphap;
    am  =-0.5*spminus*alpham;
    azr =-0.5*spzero *alphazr;
    azv1=-0.5*spzero *dv1;
    qr[i+(np+2)*(j+nt*VARRHO)]=r +(ap+am+azr);
    qr[i+(np+2)*(j+nt*VARVX )]=u +(am-am    )*cc/r;
    qr[i+(np+2)*(j+nt*VARVY )]=v1+(azv1     );
    qr[i+(np+2)*(j+nt*VARPR )]=p +(ap+am    )*csq;

    //Left
    spminus=((u-cc)<=0.0)?0.0:(u-cc)*dtdx-1.0;
    spzero =((u   )<=0.0)?0.0:(u   )*dtdx-1.0;
    spplus =((u+cc)<=0.0)?0.0:(u+cc)*dtdx-1.0;
    ap  =-0.5*spplus *alphap;
    am  =-0.5*spminus*alpham;
    azr =-0.5*spzero *alphazr;
    azv1=-0.5*spzero *dv1;
    ql[i+(np+2)*(j+nt*VARRHO)]=r +(ap+am+azr);
    ql[i+(np+2)*(j+nt*VARVX )]=u +(am-am    )*cc/r;
    ql[i+(np+2)*(j+nt*VARVY )]=v1+(azv1     );
    ql[i+(np+2)*(j+nt*VARPR )]=p +(ap+am    )*csq;
  }
}

__global__ void riemann(double *flx, double *qxm, double *qxp, int np, int nt){
  int thInd=threadIdx.x+blockDim.x*blockIdx.x;
  int i=thInd%(np+1);
  int j=thInd/(np+1);
  int n;
  double gmma6, entho, smallpp;
  double qgdnvR, qgdnvVX, qgdnvVY, qgdnvP;
  double rl,vxl,vyl,pl,cl,wl,ql,vsl;
  double rr,vxr,vyr,pr,cr,wr,qr,vsr;
  double ro,vxo,po,wo,co;
  double rx,vxx,px,wx,cx;
  double sgnm,scr,frac;
  double spout,spin,ushk;
  double ekin,etot,delp;

  smallpp=d_smallr*d_smallp;
  gmma6=(d_gamma+1)/(2.0*d_gamma);
  entho=1.0/(d_gamma-1.0);

  if(j<nt){
    rl =fmax(qxm[i  +(np+2)*(j+nt*VARRHO)],d_smallr);
    vxl=     qxm[i  +(np+2)*(j+nt*VARVX )];
    vyl=     qxm[i  +(np+2)*(j+nt*VARVY )];
    pl =fmax(qxm[i  +(np+2)*(j+nt*VARPR )],rl*d_smallp);

    rr =fmax(qxp[i+1+(np+2)*(j+nt*VARRHO)],d_smallr);
    vxr=     qxp[i+1+(np+2)*(j+nt*VARVX )];
    vyr=     qxp[i+1+(np+2)*(j+nt*VARVY )];
    pr =fmax(qxp[i+1+(np+2)*(j+nt*VARPR )],rl*d_smallp);

    cl=d_gamma*pl*rl;
    cr=d_gamma*pr*rr;

    wl=sqrt(cl);
    wr=sqrt(cr);

    px=fmax(0.0,((wr*pl+wl*pr)+wl*wr*(vxl-vxr))/(wl+wr));
    for(n=0;n<d_niterR;n++){
      wl=sqrt(cl*(1.0+gmma6*(px-pl)/pl));
      wr=sqrt(cr*(1.0+gmma6*(px-pr)/pr));
      ql=2.0*wl*wl*wl/(wl*wl+cl);
      qr=2.0*wr*wr*wr/(wr*wr+cr);
      vsl=vxl-(px-pl)/wl;
      vsr=vxr+(px-pr)/wr;
      delp=fmax(-px,qr*ql/(qr+ql)*(vsl-vsr));
      px+=delp;
      vxo=fabs(delp/(px+smallpp));
      if(vxo<1.0e-6)break;
    }
    wl=sqrt(cl*(1.0+gmma6*(px-pl)/pl));
    wr=sqrt(cr*(1.0+gmma6*(px-pr)/pr));
    vxx=0.5*(vxl+(pl-px)/wl+
             vxr-(pr-px)/wr);
    if(vxx>=0.0){
      sgnm=1.0;
      ro = rl;
      vxo=vxl;
      po = pl;
      wo = wl;
      qgdnvVY=vyl;
    }else{
      sgnm=-1.0;
      ro = rr;
      vxo=vxr;
      po = pr;
      wo = wr;
      qgdnvVY=vyr;
    }
    co=fmax(d_smallc,sqrt(fabs(d_gamma*po/ro)));
    rx=fmax(d_smallr,ro/(1.0+ro*(po-px)/(wo*wo)));
    cx=fmax(d_smallc,sqrt(fabs(d_gamma*px/rx)));

    spout=co   -sgnm*vxo;
    spin =cx   -sgnm*vxx;
    ushk =wo/ro-sgnm*vxo;

    if(px>=po){
      spin=ushk;
      spout=ushk;
    }

    scr=fmax(spout-spin,d_smallc+fabs(spout+spin));

    frac=0.5*(1.0+(spout+spin)/scr);
    frac=fmax(0.0,fmin(1.0,frac));
    qgdnvR =frac* rx+(1.0-frac)* ro;
    qgdnvVX=frac*vxx+(1.0-frac)*vxo;
    qgdnvP=frac* px+(1.0-frac)* po;
    if(spout<0.0){
      qgdnvR = ro;
      qgdnvVX=vxo;
      qgdnvP = po;
    }
    if(spin>0.0){
      qgdnvR = rx;
      qgdnvVX=vxx;
      qgdnvP = px;
    }

    flx[i+(np+1)*(j+nt*VARRHO)]=qgdnvR*qgdnvVX;
    flx[i+(np+1)*(j+nt*VARVX )]=qgdnvR*qgdnvVX*qgdnvVX+qgdnvP;
    flx[i+(np+1)*(j+nt*VARVY )]=qgdnvVY;//qgdnvR*qgdnvVX*qgdnvVY;
    ekin=0.5*qgdnvR*(qgdnvVX*qgdnvVX+qgdnvVY*qgdnvVY);
    etot=qgdnvP*entho+ekin;
    flx[i+(np+1)*(j+nt*VARPR )]=qgdnvVX*(etot+qgdnvP);
  }
}

__global__ void addFluxX(double *u, double *flx, double dtdx){
  int thInd=threadIdx.x+blockDim.x*blockIdx.x;
  int i=thInd%d_nx;
  int j=thInd/d_nx;

  if(j<d_ny){
    u[i+2+(d_nx+4)*(j+2+(d_ny+4)*VARRHO)]+=dtdx*(flx[i  +(d_nx+1)*(j+d_ny*VARRHO)]-
                                                 flx[i+1+(d_nx+1)*(j+d_ny*VARRHO)]);
    u[i+2+(d_nx+4)*(j+2+(d_ny+4)*VARVX )]+=dtdx*(flx[i  +(d_nx+1)*(j+d_ny*VARVX )]-
                                                 flx[i+1+(d_nx+1)*(j+d_ny*VARVX )]);
    u[i+2+(d_nx+4)*(j+2+(d_ny+4)*VARVY )]+=dtdx*(flx[i  +(d_nx+1)*(j+d_ny*VARVY )]-
                                                 flx[i+1+(d_nx+1)*(j+d_ny*VARVY )]);
    u[i+2+(d_nx+4)*(j+2+(d_ny+4)*VARPR )]+=dtdx*(flx[i  +(d_nx+1)*(j+d_ny*VARPR )]-
                                                 flx[i+1+(d_nx+1)*(j+d_ny*VARPR )]);
  }
}

__global__ void addFluxY(double *u, double *flx, double dtdx){
  int thInd=threadIdx.x+blockDim.x*blockIdx.x;
  int i=thInd/d_ny;
  int j=thInd%d_ny;

  if(i<d_nx){
    u[i+2+(d_nx+4)*(j+2+(d_ny+4)*VARRHO)]+=dtdx*(flx[j  +(d_ny+1)*(i+d_nx*VARRHO)]-
                                                 flx[j+1+(d_ny+1)*(i+d_nx*VARRHO)]);
    u[i+2+(d_nx+4)*(j+2+(d_ny+4)*VARVX )]+=dtdx*(flx[j  +(d_ny+1)*(i+d_nx*VARVY )]-
                                                 flx[j+1+(d_ny+1)*(i+d_nx*VARVY )]);
    u[i+2+(d_nx+4)*(j+2+(d_ny+4)*VARVY )]+=dtdx*(flx[j  +(d_ny+1)*(i+d_nx*VARVX )]-
                                                 flx[j+1+(d_ny+1)*(i+d_nx*VARVX )]);
    u[i+2+(d_nx+4)*(j+2+(d_ny+4)*VARPR )]+=dtdx*(flx[j  +(d_ny+1)*(i+d_nx*VARPR )]-
                                                 flx[j+1+(d_ny+1)*(i+d_nx*VARPR )]);
  }
}

