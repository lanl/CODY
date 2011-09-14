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
__device__ __constant__ int dev_scheme;
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
    cudaMemcpyToSymbol("dev_scheme",&(Ha->scheme),sizeof(int));
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
  i_c[4]=dev_scheme;
  i_c[5]=dev_iorder;
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
  int i_c[6], *dev_i_c;//variable to recieve integer constants
  double d_c[7], *dev_d_c; //variable to recieve double constants

  errVar=cudaMalloc(&dev_i_c,6*sizeof(int));
  HANDLE_CUDA_ERROR(errVar);
  errVar=cudaMalloc(&dev_d_c,7*sizeof(double));
  HANDLE_CUDA_ERROR(errVar);

  cpy_consts<<<1,1>>>(dev_i_c,dev_d_c);
  errVar=cudaGetLastError();
  HANDLE_CUDA_ERROR(errVar);

  errVar=cudaMemcpy(i_c,dev_i_c,6*sizeof(int),cudaMemcpyDeviceToHost);
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
  if(i_c[4]==HSCHEME_MUSCL)
    fprintf(stderr,"%-15s=%12s\n","scheme","muscl"     );
  else if(i_c[4]==HSCHEME_PLMDE)
    fprintf(stderr,"%-15s=%12s\n","scheme","plmde"     );
  else
    fprintf(stderr,"%-15s=%12s\n","scheme","collela"   );
  fprintf(stderr,"%-15s=%12d\n","iorder",i_c[5]        );
  fprintf(stderr,"%-15s=%12.5g\n","slope_type",d_c[6]  );
  fprintf(stderr,"----------------------------------\n");

}

__global__ void calc_denom(double *u, double *den){
    __shared__ double denom[MAX_THREADS];
    uint stInd;
    int thInd;
    int mxInd;
    double rho, vx, vy, eken, eint;
    double p, c;
    double courx, coury;

    thInd=threadIdx.x;
    stInd=threadIdx.x+blockDim.x*blockIdx.x;
    if(stInd<dev_nx*dev_ny){
      rho=fmax(u[(ID*dev_ny*dev_nx)+stInd],dev_smallr);
      vx=u[(IU*dev_ny*dev_nx)+stInd]/rho;
      vy=u[(IV*dev_ny*dev_nx)+stInd]/rho;
      eken=0.5*(vx*vx+vy*vy);
      eint=u[(IP*dev_ny*dev_nx)+stInd]/rho-eken;

      p=fmax((dev_gamma-(double)1.0)*rho*eint,rho*dev_smallp);
      c=sqrt(dev_gamma*p/rho);

      courx=(c+fabs(vx))/dev_dx;
      coury=(c+fabs(vy))/dev_dy;
      denom[thInd]=courx+coury;
    }else{
      denom[thInd]=0.0;
    }

    __syncthreads();

    mxInd=blockDim.x;
    while(mxInd>1){
      mxInd=(mxInd+1)/2;
      if(thInd+mxInd<blockDim.x){
        denom[thInd]=fmax(denom[thInd],denom[thInd+mxInd]);
      }
      __syncthreads();
    }

    if(thInd==0)den[blockIdx.x]=denom[0];
}

__global__ void redu_max(double *arrIn, double *arrOut, int nVals){
    __shared__ double arr[MAX_THREADS];
    int arrInd, thInd, mxInd, stride;

    thInd=threadIdx.x;
    arrInd=thInd+2*blockDim.x*blockIdx.x;
    mxInd=blockDim.x;
    stride=mxInd;
    if(arrInd+stride<nVals){
      arr[thInd]=fmax(arrIn[arrInd+stride],arrIn[arrInd]);
    }else if(arrInd<nVals){
      arr[thInd]=arrIn[arrInd];
    }else{
      arr[thInd]=0.0;
    }
    while(stride>0){
      __syncthreads();
      stride>>=1;
      if(thInd+stride<mxInd){
        arr[thInd]=fmax(arr[thInd],arr[thInd+stride]);
      }else{
        arr[thInd]=0.0;
      }
    }
    if(thInd==0)arrOut[blockIdx.x]=arr[0];
}

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
          }else{
            ind=dev_nx-2+i;
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
          }else{
            ind=i-2;
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
          }else{
            ind=dev_ny-2+j;
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
          }else{
            ind=j-2;
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

__device__ double p_ind(double *u, double *bnd, int x, int y){
    double rho, vx, vy, eken, eint;

    rho=fmax(meshv(u,bnd,x,y,ID),dev_smallr);
    vx=meshv(u,bnd,x,y,IU)/rho;
    vy=meshv(u,bnd,x,y,IV)/rho;
    eken=0.5*(vx*vx+vy*vy);
    eint=meshv(u,bnd,x,y,IP)/rho-eken;

    return fmax((dev_gamma-(double)1.0)*rho*eint,rho*dev_smallp);
}

__global__ void halfX(double *upx, double *u, double *bnd, double dt){
    uint thInd=threadIdx.x+blockDim.x*blockIdx.x;
    int i=thInd%(dev_nx+1);
    int j=thInd/(dev_nx+1);
    int k;

    

    if(i<(dev_nx+1)&&j<dev_ny){
      //run half step
      //density calculation
      upx[(ID*dev_ny+j)*(dev_nx+1)+i]=
        0.5*(meshv(u,bnd,i,j,ID)+meshv(u,bnd,i-1,j,ID))-dt/(2.0*dev_dx)*
        ((meshv(u,bnd,i  ,j,IU))
        -(meshv(u,bnd,i-1,j,IU)));
      //momentum x calculation
      upx[(IU*dev_ny+j)*(dev_nx+1)+i]=
        0.5*(meshv(u,bnd,i,j,IU)+meshv(u,bnd,i-1,j,IU))-dt/(2.0*dev_dx)*
        ((meshv(u,bnd,i  ,j,IU)*meshv(u,bnd,i  ,j,IU)/meshv(u,bnd,i  ,j,ID)+p_ind(u,bnd,i  ,j))
        -(meshv(u,bnd,i-1,j,IU)*meshv(u,bnd,i-1,j,IU)/meshv(u,bnd,i-1,j,ID)+p_ind(u,bnd,i-1,j)));
      //momentum y calculation
      upx[(IV*dev_ny+j)*(dev_nx+1)+i]=
        0.5*(meshv(u,bnd,i,j,IV)+meshv(u,bnd,i-1,j,IV))-dt/(2.0*dev_dx)*
        ((meshv(u,bnd,i  ,j,IV)*meshv(u,bnd,i  ,j,IU)/meshv(u,bnd,i  ,j,ID))
        -(meshv(u,bnd,i-1,j,IV)*meshv(u,bnd,i-1,j,IU)/meshv(u,bnd,i-1,j,ID)));
      //energy calculation
      upx[(IP*dev_ny+j)*(dev_nx+1)+i]=
        0.5*(meshv(u,bnd,i,j,IP)+meshv(u,bnd,i-1,j,IP))-dt/(2.0*dev_dx)*
        (((meshv(u,bnd,i  ,j,IU)/meshv(u,bnd,i  ,j,ID))*(meshv(u,bnd,i  ,j,IP)+p_ind(u,bnd,i  ,j)))
        -((meshv(u,bnd,i-1,j,IU)/meshv(u,bnd,i-1,j,ID))*(meshv(u,bnd,i-1,j,IP)+p_ind(u,bnd,i-1,j))));
      //advect other variables
      for(k=4;k<dev_nvar;++k)
        upx[(k*dev_ny+j)*(dev_nx+1)+i]=
          0.5*(meshv(u,bnd,i,j,k)+meshv(u,bnd,i-1,j,k))-dt/(2.0*dev_dx)*
          ((meshv(u,bnd,i  ,j,k)*meshv(u,bnd,i  ,j,IU)/meshv(u,bnd,i  ,j,ID))
          -(meshv(u,bnd,i-1,j,k)*meshv(u,bnd,i-1,j,IU)/meshv(u,bnd,i-1,j,ID)));
    }
}

__global__ void halfY(double *upy, double *u, double *bnd, double dt){
    uint thInd=threadIdx.x+blockDim.x*blockIdx.x;
    int i=thInd%dev_nx;
    int j=thInd/dev_nx;
    int k;

    if(i<dev_nx&&j<(dev_ny+1)){
      //run half step
      //density calculation
      upy[(ID*(dev_ny+1)+j)*dev_nx+i]=
        0.5*(meshv(u,bnd,i,j,ID)+meshv(u,bnd,i,j-1,ID))-dt/(2.0*dev_dy)*
        ((meshv(u,bnd,i,j  ,IV))
        -(meshv(u,bnd,i,j-1,IV)));
      //momentum x calculation
      upy[(IU*(dev_ny+1)+j)*dev_nx+i]=
        0.5*(meshv(u,bnd,i,j,IU)+meshv(u,bnd,i,j-1,IU))-dt/(2.0*dev_dy)*
        ((meshv(u,bnd,i,j  ,IU)*meshv(u,bnd,i,j  ,IV)/meshv(u,bnd,i,j  ,ID))
        -(meshv(u,bnd,i,j-1,IU)*meshv(u,bnd,i,j-1,IV)/meshv(u,bnd,i,j-1,ID)));
      //momentum y calculation
      upy[(IV*(dev_ny+1)+j)*dev_nx+i]=
        0.5*(meshv(u,bnd,i,j,IV)+meshv(u,bnd,i,j-1,IV))-dt/(2.0*dev_dy)*
        ((meshv(u,bnd,i,j  ,IV)*meshv(u,bnd,i,j  ,IV)/meshv(u,bnd,i,j  ,ID)+p_ind(u,bnd,i,j  ))
        -(meshv(u,bnd,i,j-1,IV)*meshv(u,bnd,i,j-1,IV)/meshv(u,bnd,i,j-1,ID)+p_ind(u,bnd,i,j-1)));
      //energy calculation
      upy[(IP*(dev_ny+1)+j)*dev_nx+i]=
        0.5*(meshv(u,bnd,i,j,IP)+meshv(u,bnd,i,j-1,IP))-dt/(2.0*dev_dy)*
        (((meshv(u,bnd,i,j  ,IV)/meshv(u,bnd,i,j  ,ID))*(meshv(u,bnd,i,j  ,IP)+p_ind(u,bnd,i,j  )))
        -((meshv(u,bnd,i,j-1,IV)/meshv(u,bnd,i,j-1,ID))*(meshv(u,bnd,i,j-1,IP)+p_ind(u,bnd,i,j-1))));
      //Advect other variables
      for(k=4;k<dev_nvar;++k)
        upy[(k*(dev_ny+1)+j)*dev_nx+i]=
          0.5*(meshv(u,bnd,i,j,k)+meshv(u,bnd,i,j-1,k))-dt/(2.0*dev_dy)*
          ((meshv(u,bnd,i,j  ,k)*meshv(u,bnd,i,j  ,IV)/meshv(u,bnd,i,j  ,ID))
          -(meshv(u,bnd,i,j-1,k)*meshv(u,bnd,i,j-1,IV)/meshv(u,bnd,i,j-1,ID)));
    }
}

//LIMITERS
#define LIMITER limb
__device__ double lima(double p, double m){
    return fmax(fmin((double) 1.0, p),(double) 1.0)
          +fmax(fmin((double) 1.0, m),(double) 1.0)-1.0;
}

__device__ double limb(double p, double m){
    return fmax(fmin(fmin((double)1.0,m),p),(double)0.0);
}

__device__ double libc(double p, double m){
    return fmax(fmin(fmin(2.0,2.0*m),fmin(2.0*p,0.5*(m+p))),0.0);
}


__device__ double qxcalc(double *u,double *bnd,int i, int j){
    double rplus, rminus, rdenom;
    double dumr, dumu, dumv, dume;
    double duor, duou, duov, duoe;
    double dupr, dupu, dupv, dupe;
    double eps=1.0e-30;

    if(i-2>=-2){
      dumr=meshv(u,bnd,i-1,j,ID)-meshv(u,bnd,i-2,j,ID);
      dumu=meshv(u,bnd,i-1,j,IU)-meshv(u,bnd,i-2,j,IU);
      dumv=meshv(u,bnd,i-1,j,IV)-meshv(u,bnd,i-2,j,IV);
      dume=meshv(u,bnd,i-1,j,IP)-meshv(u,bnd,i-2,j,IP);
    }else{
      dumr=0.0;
      dumu=0.0;
      dumv=0.0;
      dume=0.0;
    }

    duor=meshv(u,bnd,i  ,j,ID)-meshv(u,bnd,i-1,j,ID);
    duou=meshv(u,bnd,i  ,j,IU)-meshv(u,bnd,i-1,j,IU);
    duov=meshv(u,bnd,i  ,j,IV)-meshv(u,bnd,i-1,j,IV);
    duoe=meshv(u,bnd,i  ,j,IP)-meshv(u,bnd,i-1,j,IP);

    if(i+1<=dev_nx+1){
      dupr=meshv(u,bnd,i+1,j,ID)-meshv(u,bnd,i  ,j,ID);
      dupu=meshv(u,bnd,i+1,j,IU)-meshv(u,bnd,i  ,j,IU);
      dupv=meshv(u,bnd,i+1,j,IV)-meshv(u,bnd,i  ,j,IV);
      dupe=meshv(u,bnd,i+1,j,IP)-meshv(u,bnd,i  ,j,IP);
    }else{
      dupr=0.0;
      dupu=0.0;
      dupv=0.0;
      dupe=0.0;
    }

    rdenom=duor*duor+duou*duou+duov*duov+duoe*duoe;
    rplus =dupr*duor+dupu*duou+dupv*duov+dupe*duoe;
    rminus=dumr*duor+dumu*duou+dumv*duov+dume*duoe;

    if(rdenom!=0){
      rplus/=rdenom;
      rminus/=rdenom;
    }else{
      rplus=(rplus+eps)/(eps);
      rminus=(rminus+eps)/(eps);
    }

    return LIMITER(rplus,rminus);
}

__global__ void cmp_flxX(double *flx, double *u, double *bnd, double *upx, double dtdx){
    uint thInd=threadIdx.x+blockDim.x*blockIdx.x;
    int i=thInd%(dev_nx+1);
    int j=thInd/(dev_nx+1);
    int k;

    double rpx, vxpx, vypx, ekenpx, eintpx,ppx;
    double cs, nu, q, cv, w;

    if(i<(dev_nx+1)&&j<dev_ny){
      rpx=fmax(upx[(ID*dev_ny+j)*(dev_nx+1)+i],dev_smallr);
      vxpx=upx[(IU*dev_ny+j)*(dev_nx+1)+i]/rpx;
      vypx=upx[(IV*dev_ny+j)*(dev_nx+1)+i]/rpx;
      ekenpx=0.5*(vxpx*vxpx+vypx*vypx);
      eintpx=upx[(IP*dev_ny+j)*(dev_nx+1)+i]/rpx-ekenpx;
      ppx=fmax((dev_gamma-1)*rpx*(eintpx),rpx*dev_smallp);

      cs=sqrt(dev_gamma*ppx/rpx);
      nu=(fabs(vxpx)+cs)*dtdx;
      q=qxcalc(u,bnd,i,j);
      cv=nu*(1.0-nu);
      w=0.5*cv*(1.0-q);

      //Density
      flx[(ID*dev_ny+j)*(dev_nx+1)+i]=dtdx*(rpx*vxpx)
         -w*(meshv(u,bnd,i  ,j,ID)-meshv(u,bnd,i-1,j,ID));
      //Momentum x dir
      flx[(IU*dev_ny+j)*(dev_nx+1)+i]=dtdx*(rpx*vxpx*vxpx+ppx)
         -w*(meshv(u,bnd,i  ,j,IU)-meshv(u,bnd,i-1,j,IU));
      //Momentum y dir
      flx[(IV*dev_ny+j)*(dev_nx+1)+i]=dtdx*(rpx*vypx*vxpx)
         -w*(meshv(u,bnd,i  ,j,IV)-meshv(u,bnd,i-1,j,IV));
      //Energy
      flx[(IP*dev_ny+j)*(dev_nx+1)+i]=dtdx*(vxpx*(upx[(IP*dev_ny+j)*(dev_nx+1)+i]+ppx))
         -w*(meshv(u,bnd,i  ,j,IP)-meshv(u,bnd,i-1,j,IP));
      //Other advected vars
      for(k=4;k<dev_nvar;++k)
        flx[(k*dev_ny+j)*(dev_nx+1)+i]=dtdx*(upx[(k*dev_ny+j)*(dev_nx+1)+i]*vxpx/rpx)
           -w*(meshv(u,bnd,i  ,j,k)-meshv(u,bnd,i-1,j,k));
    }
}

__device__ double qycalc(double *u, double *bnd, int i, int j){
    double rplus, rminus, rdenom;
    double dumr, dumu, dumv, dume;
    double duor, duou, duov, duoe;
    double dupr, dupu, dupv, dupe;
    double eps=1.0e-30;

    if(i-2>=-2){
      dumr=meshv(u,bnd,i,j-1,ID)-meshv(u,bnd,i,j-2,ID);
      dumu=0.0;//meshv(u,bnd,i,j-1,IU)-meshv(u,bnd,i,j-2,IU);
      dumv=meshv(u,bnd,i,j-1,IV)-meshv(u,bnd,i,j-2,IV);
      dume=meshv(u,bnd,i,j-1,IP)-meshv(u,bnd,i,j-2,IP);
    }else{
      dumr=0.0;
      dumu=0.0;
      dumv=0.0;
      dume=0.0;
    }

    duor=meshv(u,bnd,i,j  ,ID)-meshv(u,bnd,i,j-1,ID);
    duou=0.0;//meshv(u,bnd,i,j  ,IU)-meshv(u,bnd,i,j-1,IU);
    duov=meshv(u,bnd,i,j  ,IV)-meshv(u,bnd,i,j-1,IV);
    duoe=meshv(u,bnd,i,j  ,IP)-meshv(u,bnd,i,j-1,IP);

    if(j+1<=dev_ny+1){
      dupr=meshv(u,bnd,i,j+1,ID)-meshv(u,bnd,i,j  ,ID);
      dupu=0.0;//meshv(u,bnd,i,j+1,IU)-meshv(u,bnd,i,j  ,IU);
      dupv=meshv(u,bnd,i,j+1,IV)-meshv(u,bnd,i,j  ,IV);
      dupe=meshv(u,bnd,i,j+1,IP)-meshv(u,bnd,i,j  ,IP);
    }else{
      dupr=0.0;
      dupu=0.0;
      dupv=0.0;
      dupe=0.0;
    }

    rdenom=duor*duor+duou*duou+duov*duov+duoe*duoe;
    rplus =dupr*duor+dupu*duou+dupv*duov+dupe*duoe;
    rminus=dumr*duor+dumu*duou+dumv*duov+dume*duoe;

    if(rdenom!=0){
      rplus/=rdenom;
      rminus/=rdenom;
    }else{
      rplus=(rplus+eps)/(eps);
      rminus=(rminus+eps)/(eps);
    }

    return LIMITER(rplus,rminus);
}

__global__ void cmp_flxY(double *flx, double *u, double *bnd, double *upy, double dtdy){
    uint thInd=threadIdx.x+blockDim.x*blockIdx.x;
    int i=thInd%dev_nx;
    int j=thInd/dev_nx;
    int k;

    double rpy, vxpy, vypy, ekenpy, eintpy, ppy;
    double cs, nu, q, cv, w;

    if(i<dev_nx&&j<dev_ny+1){
      rpy=fmax(upy[(ID*(dev_ny+1)+j)*dev_nx+i],dev_smallr);
      vxpy=upy[(IU*(dev_ny+1)+j)*dev_nx+i]/rpy;
      vypy=upy[(IV*(dev_ny+1)+j)*dev_nx+i]/rpy;
      ekenpy=0.5*(vxpy*vxpy+vypy*vypy);
      eintpy=upy[(IP*(dev_ny+1)+j)*dev_nx+i]/rpy-ekenpy;
      ppy=fmax((dev_gamma-1)*rpy*(eintpy),rpy*dev_smallp);

      cs=sqrt(dev_gamma*ppy/rpy);
      nu=(fabs(vypy)+cs)*dtdy;
      q=qycalc(u,bnd,i,j);
      cv=nu*(1.0-nu);
      w=0.5*cv*(1.0-q);

      //Density
      flx[(ID*(dev_ny+1)+j)*dev_nx+i]=dtdy*(rpy*vypy)
         -w*(meshv(u,bnd,i,j  ,ID)-meshv(u,bnd,i,j-1,ID));
      //Momentum x dir
      flx[(IU*(dev_ny+1)+j)*dev_nx+i]=dtdy*(rpy*vypy*vxpy)
         -w*(meshv(u,bnd,i,j  ,IU)-meshv(u,bnd,i,j-1,IU));
      //Momentum y dir
      flx[(IV*(dev_ny+1)+j)*dev_nx+i]=dtdy*(rpy*vypy*vypy+ppy)
         -w*(meshv(u,bnd,i,j  ,IV)-meshv(u,bnd,i,j-1,IV));
      //Energy
      flx[(IP*(dev_ny+1)+j)*dev_nx+i]=dtdy*(vypy*(upy[(IP*(dev_ny+1)+j)*dev_nx+i]+ppy))
         -w*(meshv(u,bnd,i,j  ,IP)-meshv(u,bnd,i,j-1,IP));
      //Other advected vars
      for(k=4;k<dev_nvar;++k)
        flx[(k*(dev_ny+1)+j)*dev_nx+i]=dtdy*(upy[(k*(dev_ny+1)+j)*dev_nx+i]*vypy/rpy)
           -w*(meshv(u,bnd,i,j  ,k)-meshv(u,bnd,i,j-1,k));
    }
}

__global__ void cmp_dStX(double *du, double *flx){
    uint thInd=threadIdx.x+blockDim.x*blockIdx.x;
    int i=thInd%dev_nx;
    int j=thInd/dev_nx;
    int k;

    if(i<dev_nx&&j<dev_ny){
      for(k=0;k<dev_nvar;++k)
        du[(k*dev_ny+j)*dev_nx+i]+=flx[(k*dev_ny+j)*(dev_nx+1)+i  ]
                                  -flx[(k*dev_ny+j)*(dev_nx+1)+i+1];
    }
}

__global__ void cmp_dStY(double *du, double *flx){
    uint thInd=threadIdx.x+blockDim.x*blockIdx.x;
    int i=thInd%dev_nx;
    int j=thInd/dev_nx;
    int k;

    if(i<dev_nx&&j<dev_ny){
      for(k=0;k<dev_nvar;++k)
        du[(k*dev_ny+j)*dev_nx+i]+=flx[(k*(dev_ny+1)+j  )*dev_nx+i]
                                  -flx[(k*(dev_ny+1)+j+1)*dev_nx+i];
    }
}

__global__ void add_dSt(double *u, double *du){
    uint thInd=threadIdx.x+blockDim.x*blockIdx.x;
    int i=thInd%dev_nx;
    int j=thInd/dev_nx;
    int k;

    if(i<dev_nx&&j<dev_ny){
      for(k=0;k<dev_nvar;k++)
        u[(k*dev_ny+j)*dev_nx+i]+=du[(k*dev_ny+j)*dev_nx+i];
    }
}

