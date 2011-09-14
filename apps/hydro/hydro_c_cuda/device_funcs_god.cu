#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>

#include "device_funcs.h"
#include "hydro_macros.h"

#define MAX_THREADS 1024

//Constants
 //Problem constants
__device__ __constant__ int dev_nx, dev_ny, dev_nvar;
__device__ __constant__ int dev_nTh;
__device__ __constant__ double dev_dx, dev_dy;
__device__ __constant__ double dev_gamma;
 //Model constants
__device__ __constant__ double dev_smallr, dev_smallc, dev_smallp;
__device__ __constant__ int dev_iter_RS;
__device__ __constant__ double dev_zL, dev_zR, dev_proj;
__device__ __constant__ int dev_iorder;
__device__ __constant__ double dev_sl_type;

//error handlers

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

void device_init(hydro_args *Ha, hydro_prob *Hp){
    cudaError_t errVar;
    double smallp;
    //Scheme arrays, defines value for each scheme, removes need for 3 if
    //statements in trace
    double zerol;
    double zeror;
    double project;


    if(Ha->scheme==HSCHEME_MUSCL){
      zerol=-100.0;
      zeror=+100.0;
      project=1.0;
    }else if(Ha->scheme==HSCHEME_PLMDE){
      zerol=0.0;
      zeror=0.0;
      project=1.0;
    }else {
      zerol=0.0;
      zeror=0.0;
      project=0.0;
    }

    //Initialize constants
    smallp=Ha->smallc*Ha->smallc/Hp->gamma;
     //Problem constants
    cudaMemcpyToSymbol("dev_nx",&(Hp->nx),sizeof(int));
    HANDLE_CUDA_ERROR(errVar);
    cudaMemcpyToSymbol("dev_ny",&(Hp->ny),sizeof(int));
    HANDLE_CUDA_ERROR(errVar);
    cudaMemcpyToSymbol("dev_nvar",&(Hp->nvar),sizeof(int));
    HANDLE_CUDA_ERROR(errVar);
    cudaMemcpyToSymbol("dev_dx",&(Hp->dx),sizeof(double));
    HANDLE_CUDA_ERROR(errVar);
    cudaMemcpyToSymbol("dev_dy",&(Hp->dy),sizeof(double));
    HANDLE_CUDA_ERROR(errVar);
    cudaMemcpyToSymbol("dev_gamma",&(Hp->gamma),sizeof(double));
    HANDLE_CUDA_ERROR(errVar);
     //Model constants
    cudaMemcpyToSymbol("dev_smallr",&(Ha->smallr),sizeof(double));
    HANDLE_CUDA_ERROR(errVar);
    cudaMemcpyToSymbol("dev_smallc",&(Ha->smallc),sizeof(double));
    HANDLE_CUDA_ERROR(errVar);
    cudaMemcpyToSymbol("dev_smallp",&(smallp),sizeof(double));
    HANDLE_CUDA_ERROR(errVar);
    cudaMemcpyToSymbol("dev_sl_type",&(Ha->slope_type),sizeof(double));
    HANDLE_CUDA_ERROR(errVar);
    cudaMemcpyToSymbol("dev_iter_RS",&(Ha->niter_riemann),sizeof(int));
    HANDLE_CUDA_ERROR(errVar);
    cudaMemcpyToSymbol("dev_zL",&zerol,sizeof(double));
    HANDLE_CUDA_ERROR(errVar);
    cudaMemcpyToSymbol("dev_zR",&zeror,sizeof(double));
    HANDLE_CUDA_ERROR(errVar);
    cudaMemcpyToSymbol("dev_proj",&project,sizeof(double));
    HANDLE_CUDA_ERROR(errVar);
    cudaMemcpyToSymbol("dev_iorder",&(Ha->iorder),sizeof(int));
    HANDLE_CUDA_ERROR(errVar);
}

__global__ void cpy_consts(int *i_c, double *d_c){
  //Integer constants
  i_c[0]=dev_nx;
  i_c[1]=dev_ny;
  i_c[2]=dev_nvar;
  i_c[3]=dev_iter_RS;
  i_c[4]=dev_iorder;
  //Float constants
  d_c[0]=dev_dx;
  d_c[1]=dev_dy;
  d_c[2]=dev_gamma;
  d_c[3]=dev_smallr;
  d_c[4]=dev_smallc;
  d_c[5]=dev_smallp;
  d_c[6]=dev_sl_type;
}

void symbol_test(hydro_args *Ha, hydro_prob *Hp){
  cudaError_t errVar;
  int i_c[5], *dev_i_c;//variable to recieve integer constants
  double d_c[7], *dev_d_c; //variable to recieve double constants

  errVar=cudaMalloc(&dev_i_c,5*sizeof(int));
  HANDLE_CUDA_ERROR(errVar);
  errVar=cudaMalloc(&dev_d_c,7*sizeof(double));
  HANDLE_CUDA_ERROR(errVar);

  cpy_consts<<<1,1>>>(dev_i_c,dev_d_c);
  errVar=cudaGetLastError();
  HANDLE_CUDA_ERROR(errVar);

  errVar=cudaMemcpy(i_c,dev_i_c,5*sizeof(int),cudaMemcpyDeviceToHost);
  HANDLE_CUDA_ERROR(errVar);
  errVar=cudaMemcpy(d_c,dev_d_c,7*sizeof(double),cudaMemcpyDeviceToHost);
  HANDLE_CUDA_ERROR(errVar);

  errVar=cudaFree(dev_i_c);
  HANDLE_CUDA_ERROR(errVar);
  errVar=cudaFree(dev_d_c);
  HANDLE_CUDA_ERROR(errVar);

  fprintf(stderr,"Device Variables:\n"                 );
  fprintf(stderr,"----------------------------------\n");
  fprintf(stderr,"%-15s=%12d\n","nx",i_c[0]            );
  fprintf(stderr,"%-15s=%12d\n","ny",i_c[1]            );
  fprintf(stderr,"%-15s=%12d\n","nvar",i_c[2]          );
  fprintf(stderr,"%-15s=%12.5g\n","dx",d_c[0]          );
  fprintf(stderr,"%-15s=%12.5g\n","dy",d_c[1]          );
  fprintf(stderr,"%-15s=%12.5g\n","gamma",d_c[2]       );
  fprintf(stderr,"----------------------------------\n");
  fprintf(stderr,"%-15s=%12.5g\n","smallr",d_c[3]      );
  fprintf(stderr,"%-15s=%12.5g\n","smallc",d_c[4]      );
  fprintf(stderr,"%-15s=%12.5g\n","smallp",d_c[5]      );
  fprintf(stderr,"%-15s=%12d\n","niter_riemann",i_c[3] );
  fprintf(stderr,"%-15s=%12d\n","iorder",i_c[4]        );
  fprintf(stderr,"%-15s=%12.5g\n","slope_type",d_c[6]  );
  fprintf(stderr,"----------------------------------\n");

}

extern __shared__ double dynVar[];
__global__ void calc_denom(double *u, double *den){
    double *denom=dynVar;
    uint stInd;
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
    if(stInd<dev_nx*dev_ny){
      rho= u[(ID*dev_ny*dev_nx)+stInd];
      vx=  u[(IU*dev_ny*dev_nx)+stInd];
      vy=  u[(IV*dev_ny*dev_nx)+stInd];
      eint=u[(IP*dev_ny*dev_nx)+stInd];
      rho=fmax(rho,dev_smallr);
      vx=vx/rho;
      vy=vy/rho;
      eken=0.5*(vx*vx+vy*vy);
      eint=eint/rho-eken;

      p=fmax((dev_gamma-(double)1.0)*rho*eint,rho*dev_smallp);
      c=sqrt(dev_gamma*p/rho);

      courx=(c+fabs(vx))/dev_dx;
      coury=(c+fabs(vy))/dev_dy;
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

//Generate boundary conditions
__global__ void gen_bndX(double *u, double *bnd, int bndL, int bndR){
    int stIndY=threadIdx.x+blockDim.x*blockIdx.x;
    int i,k;
    double sgn;
    int ind;

    if(stIndY<dev_ny){
      //Left boundary
      for(k=0;k<dev_nvar;++k){
        for(i=0;i<2;++i){
          sgn=1.0;
          if(bndL==BND_REFL){
            ind=1-i;
            if(k==IU)sgn=-1.0;
          }else if(bndL==BND_PERM){
            ind=0;
          }else if(bndL==BND_WRAP){
            ind=dev_nx-2+i;
          }else if(bndL==BND_INTR){
            continue;
          }
          bnd[(k*dev_ny + stIndY)*4 +i]=sgn*u[(k*dev_ny + stIndY)*dev_nx + ind];
        }
      }
      //Right boundary
      for(k=0;k<dev_nvar;++k){
        for(i=2;i<4;++i){
          sgn=1.0;
          if(bndR==BND_REFL){
            ind=dev_nx+1-i;
            if(k==IU)sgn=-1.0;
          }else if(bndR==BND_PERM){
            ind=dev_nx-1;
          }else if(bndR==BND_WRAP){
            ind=i-2;
          }else if(bndR==BND_INTR){
            continue;
          }
          bnd[(k*dev_ny + stIndY)*4 +i]=sgn*u[(k*dev_ny + stIndY)*dev_nx + ind];
        }
      }
    }
}
__global__ void gen_bndY(double *u, double *bnd, int bndU, int bndD){
    int stIndX=threadIdx.x+blockDim.x*blockIdx.x;
    int j,k;
    double sgn;
    int ind;

    if(stIndX<dev_nx){
      //Down boundary
      for(k=0;k<dev_nvar;++k){
        for(j=0;j<2;++j){
          sgn=1.0;
          if(bndD==BND_REFL){
            ind=1-j;
            if(k==IV)sgn=-1.0;
          }else if(bndD==BND_PERM){
            ind=0;
          }else if(bndD==BND_WRAP){
            ind=dev_ny-2+j;
          }else if(bndD==BND_INTR){
            continue;
          }
          bnd[(k*dev_nx + stIndX)*4 +j]=sgn*u[(k*dev_ny +ind)*dev_nx +stIndX];
        }
      }
      //Up boundary
      for(k=0;k<dev_nvar;++k){
        for(j=2;j<4;++j){
          sgn=1.0;
          if(bndU==BND_REFL){
            ind=dev_ny+1-j;
            if(k==IV)sgn=-1.0;
          }else if(bndU==BND_PERM){
            ind=dev_ny-1;
          }else if(bndU==BND_WRAP){
            ind=j-2;
          }else if(bndU==BND_INTR){
            continue;
          }
          bnd[(k*dev_nx + stIndX)*4 +j]=sgn*u[(k*dev_ny +ind)*dev_nx +stIndX];
        }
      }
    }
}

__device__ double meshv(double *u, double *bnd, int x, int y, int nv){
  if(-2<=x&&x<=dev_nx+1&&-2<=y&&y<=dev_ny+1){
    if(x<0)
      return bnd[(nv*dev_ny+y)*4+x+2];//return from left bound array
    if(y<0)
      return bnd[(nv*dev_nx+x)*4+y+2];//return from down bound array
    if(x>=dev_nx)
      return bnd[(nv*dev_ny+y)*4+(x-dev_nx)+2];//return from right bound array
    if(y>=dev_ny)
      return bnd[(nv*dev_nx+x)*4+(y-dev_ny)+2];//return from up bound array
    return u[(nv*dev_ny+y)*dev_nx+x];
  }
  return 0.0;
}

//Prepare standard mesh of primitive state variables for
//computations X-dir
__global__ void primX(double *q, double *u, double *bnd){
    int i=threadIdx.x+blockDim.x*blockIdx.x;
    int j=blockIdx.y;
    int k;
    double rho, vx, vy, eken, eint, p;

    if(i<dev_nx+4&&j<dev_ny){
      rho= meshv(u,bnd,i-2,j,ID);
      vx=  meshv(u,bnd,i-2,j,IU);
      vy=  meshv(u,bnd,i-2,j,IV);
      eint=meshv(u,bnd,i-2,j,IP);
      rho=fmax(rho,dev_smallr);
      vx=vx/rho;
      vy=vy/rho;
      eken=0.5*(vx*vx+vy*vy);
      eint=eint/rho-eken;
      p=fmax((dev_gamma-1.0)*rho*eint,rho*dev_smallp); 

      q[(ID*dev_ny+j)*(dev_nx+4)+i]=rho;
      q[(IU*dev_ny+j)*(dev_nx+4)+i]=vx;
      q[(IV*dev_ny+j)*(dev_nx+4)+i]=vy;
      q[(IP*dev_ny+j)*(dev_nx+4)+i]=p;
      for(k=4;k<dev_nvar;++k)
        q[(k*dev_ny+j)*(dev_nx+4)+i]=meshv(u,bnd,i-2,j,k);
    }
}

//Prepare standard mesh of primitive state variables for
//computations X-dir
__global__ void primY(double *q, double *u, double *bnd){
    int i=threadIdx.x+blockDim.x*blockIdx.x;
    int j=blockIdx.y;
    int k;
    double rho, vx, vy, eken, eint, p;

    if(i<dev_nx&&j<dev_ny+4){
      rho= meshv(u,bnd,i,j-2,ID);
      vx=  meshv(u,bnd,i,j-2,IU);
      vy=  meshv(u,bnd,i,j-2,IV);
      eint=meshv(u,bnd,i,j-2,IP);
      rho=fmax(rho,dev_smallr);
      vx=vx/rho;
      vy=vy/rho;
      eken=0.5*(vx*vx+vy*vy);
      eint=eint/rho-eken;
      p=fmax((dev_gamma-1.0)*rho*eint,rho*dev_smallp); 

      q[(ID*dev_nx+i)*(dev_ny+4)+j]=rho;
      q[(IU*dev_nx+i)*(dev_ny+4)+j]=vy;
      q[(IV*dev_nx+i)*(dev_ny+4)+j]=vx;
      q[(IP*dev_nx+i)*(dev_ny+4)+j]=p;
      for(k=4;k<dev_nvar;++k)
        q[(k*dev_nx+i)*(dev_ny+4)+j]=meshv(u,bnd,i,j-2,k);
    }
}

//calculate slope of state var nv @ i,j
__device__ double slope(double *q, int i, int j, int nv, int np, int nt){
    double dlft, drgt, dcen, dsgn, dlim;
    dlft=dev_sl_type*(q[(nv*nt+j)*np+i  ]-q[(nv*nt+j)*np+i-1]);
    drgt=dev_sl_type*(q[(nv*nt+j)*np+i+1]-q[(nv*nt+j)*np+i  ]);
    dcen=0.5*(dlft+drgt)/dev_sl_type;
    if(dcen>=0)dsgn=1.0;
    else dsgn=-1.0;
    dlim=fmin(fabs(dlft),fabs(drgt));
    if((dlft*drgt)<=0.0)dlim=0.0;
    return dsgn*fmin(dlim,fabs(dcen));
}

//Produce interpolations
__global__ void trace(double *qxm, double *qxp, double *q, int np, int nt, double dtdx){
    int i=threadIdx.x+blockDim.x*blockIdx.x;
    int j=blockIdx.y;
    int k;

    double cc, csq;
    double  r, u, v, p, a;
    double dr,du,dv,dp,da;
    double alpham,alphap, alpha0r, alpha0v;
    double spminus,spzero,spplus;
    double ap, am, azr, azv1, acmp;
    double zerol, zeror, project;

    if(i<np&&j<nt){
      zerol=dev_zL/dtdx;
      zeror=dev_zR/dtdx;
      project=dev_proj;

      //collect primitive state vars
      r=q[(ID*nt+j)*(np+2)+i+1];
      u=q[(IU*nt+j)*(np+2)+i+1];
      v=q[(IV*nt+j)*(np+2)+i+1];
      p=q[(IP*nt+j)*(np+2)+i+1];

      //calculate sound speeds
      cc=sqrt(dev_gamma*p/r);
      csq=cc*cc;

      //calculate slopes
      dr=slope(q,i+1,j,ID,np+2,nt);
      du=slope(q,i+1,j,IU,np+2,nt);
      dv=slope(q,i+1,j,IV,np+2,nt);
      dp=slope(q,i+1,j,IP,np+2,nt);

      //Calculate ??
      alpham = 0.5*(dp/(r*cc)-du)*r/cc;
      alphap = 0.5*(dp/(r*cc)+du)*r/cc;
      alpha0r = dr - dp/csq;
      alpha0v = dv;

      //Right state
      spminus = (u-cc)*dtdx+1.0;
      spzero  = (u   )*dtdx+1.0;
      spplus  = (u+cc)*dtdx+1.0;
      if((u-cc)>=zeror)spminus=project;
      if((u   )>=zeror)spzero =project;
      if((u+cc)>=zeror)spplus =project;
      ap  =-0.5*spplus *alphap;
      am  =-0.5*spminus*alpham;
      azr =-0.5*spzero *alpha0r;
      azv1=-0.5*spzero *alpha0v;
      qxp[(ID*nt+j)*np+i]=r+(ap+am+azr);
      qxp[(IU*nt+j)*np+i]=u+(ap-am    )*cc/r;
      qxp[(IV*nt+j)*np+i]=v+(azv1     );
      qxp[(IP*nt+j)*np+i]=p+(ap+am    )*csq;

      //Left state
      spminus = (u-cc)*dtdx-1.0;
      spzero  = (u   )*dtdx-1.0;
      spplus  = (u+cc)*dtdx-1.0;
      if((u-cc)<=zerol)spminus=-project;
      if((u   )<=zerol)spzero =-project;
      if((u+cc)<=zerol)spplus =-project;
      ap  =-0.5*spplus *alphap;
      am  =-0.5*spminus*alpham;
      azr =-0.5*spzero *alpha0r;
      azv1=-0.5*spzero *alpha0v;
      qxm[(ID*nt+j)*np+i]=r+(ap+am+azr);
      qxm[(IU*nt+j)*np+i]=u+(ap-am    )*cc/r;
      qxm[(IV*nt+j)*np+i]=v+(azv1     );
      qxm[(IP*nt+j)*np+i]=p+(ap+am    )*csq;

      for(k=4;k<dev_nvar;++k){
        a=q[(k*nt+j)*(np+2)+i+1];
        da=slope(q,i+1,j,k,np+2,nt);

        //right state
        spzero=u*dtdx+1.0;
        if(u>=zeror)spzero=project;
        acmp=-0.5*spzero*da;
        qxp[(k*nt+j)*np+i]=a+acmp;

        //left state
        spzero=u*dtdx+1.0;
        if(u<=zeror)spzero=-project;
        acmp=-0.5*spzero*da;
        qxm[(k*nt+j)*np+i]=a+acmp;
      }
    }
}

__global__ void riemann(double *flx, double *qxm, double *qxp, int np, int nt){
    int i=threadIdx.x+blockDim.x*blockIdx.x;
    int j=blockIdx.y;
    int k;
    int n;

    double rl, ul, vl, pl, cl, wl;//left
    double rr, ur, vr, pr, cr, wr;//right
    double ro, uo, po, co, wo;//??
    double rx, ux, px, cx;//star
    double sgnm,spin,spout,ushk;
    double frac, scr,delp;
    double smallpp,gmma6,ql,qr,usr,usl;

    double qgdnv[4];

    double entho, etot, ekin;

    if(i<np&&j<nt){
      smallpp=dev_smallr*dev_smallp;
      gmma6=(dev_gamma+1.0)/(2.0*dev_gamma);
	
      //density, pressure, and velocity values
      rl=fmax(qxm[(ID*nt+j)*(np+1)+i  ],dev_smallr);
      ul=     qxm[(IU*nt+j)*(np+1)+i  ];
      vl=     qxm[(IV*nt+j)*(np+1)+i  ];
      pl=fmax(qxm[(IP*nt+j)*(np+1)+i  ],rl*dev_smallp);
      rr=fmax(qxp[(ID*nt+j)*(np+1)+i+1],dev_smallr);
      ur=     qxp[(IU*nt+j)*(np+1)+i+1];
      vr=     qxp[(IV*nt+j)*(np+1)+i+1];
      pr=fmax(qxp[(IP*nt+j)*(np+1)+i+1],rr*dev_smallp);

      //Lagrangian sound speed
      cl=dev_gamma*pl*rl;
      cr=dev_gamma*pr*rr;

      //init guess
      wl=sqrt(cl);
      wr=sqrt(cr);
      px=((wr*pl+wl*pr)+wl*wr*(ul-ur))/(wl+wr);
      px=fmax(px,0.0);
      //Run Newton-Raphson iterations
      for(n=0;n<dev_iter_RS;++n){
        wl=sqrt(cl*(1.0+gmma6*(px-pl)/pl));
        wr=sqrt(cr*(1.0+gmma6*(px-pr)/pr));
        ql=2.0*wl*wl*wl/(wl*wl+cl);
        qr=2.0*wr*wr*wr/(wr*wr+cr);
        usl=ul-(px-pl)/wl;
        usr=ur+(px-pr)/wr;
        delp=fmax(qr*ql/(qr+ql)*(usl-usr),-px);
        px+=delp;
        uo=fabs(delp/(px+smallpp));
        if(uo<1.0e-6)break;
      }
      //Do calculations for actual values
      wl=sqrt(cl*(1.0+gmma6*(px-pl)/pl));
      wr=sqrt(cr*(1.0+gmma6*(px-pr)/pr));
      ux=0.5*(ul+(pl-px)/wl
             +ur-(pr-px)/wr);
      //??
      if(ux>=0){
        sgnm=1.0;
        ro=rl;
        uo=ul;
        po=pl;
        wo=wl;
        qgdnv[IV]=vl;
      }else{
        sgnm=-1.0;
        ro=rr;
        uo=ur;
        po=pr;
        wo=wr;
        qgdnv[IV]=vr;
      }
      //??
      co=fmax(dev_smallc,sqrt(fabs(dev_gamma*po/ro)));
      //??
      rx=ro/(1.0+ro*(po-px)/(wo*wo));
      rx=fmax(rx,dev_smallr);
      //??
      cx=fmax(dev_smallc,sqrt(fabs(dev_gamma*px/rx)));
      //??
      spout=co   -sgnm*uo;
      spin =cx   -sgnm*ux;
      ushk =wo/ro-sgnm*uo;
      //??
      if(px>=po){
        spin=ushk;
        spout=ushk;
      }
      //??
      scr=fmax(spout-spin,dev_smallc+fabs(spout+spin));
      //Calculate interface values?
      frac=(1.0+(spout+spin)/scr)*0.5;
      frac=fmax((double)0.0,fmin((double)1.0,frac));
      qgdnv[ID]=frac*rx+(1.0-frac)*ro;
      qgdnv[IU]=frac*ux+(1.0-frac)*uo;
      qgdnv[IP]=frac*px+(1.0-frac)*po;
      if(spout<0.0){
        qgdnv[ID]=ro;
        qgdnv[IU]=uo;
        qgdnv[IP]=po;
      }
      if(spin>0.0){
        qgdnv[ID]=rx;
        qgdnv[IU]=ux;
        qgdnv[IP]=px;
      }
      //Calculate fluxes

      entho=1.0/(dev_gamma-1.0);

      //Mass density flux term
      flx[(ID*nt+j)*np+i]=qgdnv[ID]*qgdnv[IU];
      //Momentum flux term
      flx[(IU*nt+j)*np+i]=qgdnv[ID]*qgdnv[IU]*qgdnv[IU]+qgdnv[IP];
      //Transverse momentum flux term
      flx[(IV*nt+j)*np+i]=qgdnv[ID]*qgdnv[IU]*qgdnv[IV];
      //Total energy flux term
      ekin=0.5*qgdnv[ID]*(SQ(qgdnv[IU])+SQ(qgdnv[IV]));
      etot=qgdnv[IP]*entho+ekin;
      flx[(IP*nt+j)*np+i]=qgdnv[IU]*(etot+qgdnv[IP]);

      for(k=4;k<dev_nvar;++k){
        if(sgnm==1.0){
          flx[(k*nt+j)*np+i]=flx[(ID*nt+j)*np+i]*qxm[(k*nt+j)*(np+1)+i  ];
        }else{
          flx[(k*nt+j)*np+i]=flx[(ID*nt+j)*np+i]*qxp[(k*nt+j)*(np+1)+i+1];
        }
      }
    }
}

//Compute change in state variables for X-dir pass
__global__ void cmp_dStX(double *du, double *flx, double dtdx){
    int i=threadIdx.x+blockDim.x*blockIdx.x;
    int j=blockIdx.y;
    int k;

    if(i<dev_nx&&j<dev_ny){
      du[(ID*dev_ny+j)*dev_nx+i]+= dtdx*(flx[(ID*dev_ny+j)*(dev_nx+1)+i  ]
                                        -flx[(ID*dev_ny+j)*(dev_nx+1)+i+1]);
      du[(IU*dev_ny+j)*dev_nx+i]+= dtdx*(flx[(IU*dev_ny+j)*(dev_nx+1)+i  ]
                                        -flx[(IU*dev_ny+j)*(dev_nx+1)+i+1]);
      du[(IV*dev_ny+j)*dev_nx+i]+= dtdx*(flx[(IV*dev_ny+j)*(dev_nx+1)+i  ]
                                        -flx[(IV*dev_ny+j)*(dev_nx+1)+i+1]);
      du[(IP*dev_ny+j)*dev_nx+i]+= dtdx*(flx[(IP*dev_ny+j)*(dev_nx+1)+i  ]
                                        -flx[(IP*dev_ny+j)*(dev_nx+1)+i+1]);
      for(k=4;k<dev_nvar;++k)
        du[(k*dev_ny+j)*dev_nx+i]+=dtdx*(flx[( k*dev_ny+j)*(dev_nx+1)+i  ]
                                        -flx[( k*dev_ny+j)*(dev_nx+1)+i+1]);
    }
}

//Compute change in state variables for Y-dir pass
__global__ void cmp_dStY(double *du, double *flx, double dtdx){
    int i=threadIdx.x+blockDim.x*blockIdx.x;
    int j=blockIdx.y;
    int k;

    if(i<dev_nx&&j<dev_ny){
      du[(ID*dev_ny+j)*dev_nx+i]+= dtdx*(flx[(ID*dev_nx+i)*(dev_ny+1)+j  ]
                                        -flx[(ID*dev_nx+i)*(dev_ny+1)+j+1]);
      du[(IU*dev_ny+j)*dev_nx+i]+= dtdx*(flx[(IV*dev_nx+i)*(dev_ny+1)+j  ]
                                        -flx[(IV*dev_nx+i)*(dev_ny+1)+j+1]);
      du[(IV*dev_ny+j)*dev_nx+i]+= dtdx*(flx[(IU*dev_nx+i)*(dev_ny+1)+j  ]
                                        -flx[(IU*dev_nx+i)*(dev_ny+1)+j+1]);
      du[(IP*dev_ny+j)*dev_nx+i]+= dtdx*(flx[(IP*dev_nx+i)*(dev_ny+1)+j  ]
                                        -flx[(IP*dev_nx+i)*(dev_ny+1)+j+1]);
      for(k=4;k<dev_nvar;++k)
        du[(k*dev_ny+j)*dev_nx+i]+=dtdx*(flx[( k*dev_nx+i)*(dev_ny+1)+j  ]
                                        -flx[( k*dev_nx+i)*(dev_ny+1)+j+1]);
    }
}

//Add change in state variables
__global__ void add_dSt(double *u, double *du){
    int i=threadIdx.x+blockDim.x*blockIdx.x;
    int j=blockIdx.y;
    int k;

    if(i<dev_nx&&j<dev_ny)
      for(k=0;k<dev_nvar;++k)
        u[(k*dev_ny+j)*dev_nx+i]+=du[(k*dev_ny+j)*dev_nx+i];
}

