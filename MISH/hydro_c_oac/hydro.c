#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <openacc.h>
#include <time.h>
#include <sys/time.h>
#include "hydro.h"
#include "float.h"

hydro_args *Ha;
hydro_prob *Hp;
int nx, ny;
double *q;
double *qr, *ql;
double *flx;

size_t meshSize, primSize, qSize, flxSize;


double slope(double *q,int ind);

//Utility function to get wall runtime
double getNow(){
  struct timeval tv;
  struct timezone tz;
  double t;

  gettimeofday(&tv, &tz);

  t = (double)tv.tv_sec;
  t += ((double)tv.tv_usec)/1000000.0;

  return t;
}


static inline void printArray(char* label,double *arr, int nvar, int nx, int ny, int nHx, int nHy){
#if 0
  int nV,i,j;
  printf("Array %s\n",label);
  for(nV=0;nV<nvar;nV++){
    printf("Variable %d:\n",nV);
    printf("%3s","");
    for(i=0;i<nx;i++){
      printf("|i%8di",i);
    }
    printf("|\n");
    for(j=0;j<ny;j++){
      printf("%3d",j);
      for(i=0;i<nx;i++){
	printf("|%10g",arr[i+nHx+(nx+2*nHx)*(j+nHy+(ny+2*nHy)*nV)]);
      }
      printf("|\n");
    }
  }
#endif
}

double calcDT(double *mesh){
  int i;
  double max_denom;
  double r,vx,vy,eint,p;
  double c,cx,cy;
  double smallp,smallr;
  double gamma;
  double dx,dy;
  double denom;

  max_denom=Ha->smallc;
  gamma=Hp->gamma;
  dx=Hp->dx;
  dy=Hp->dy;
  smallr=Ha->smallr;

  smallp=Ha->smallc*Ha->smallc/Hp->gamma;
#pragma acc kernels pcopyin(mesh[0:meshSize]) copyin(smallp,smallr,gamma,nx,ny,dx,dy)
#pragma acc loop independent reduction(max:max_denom) private(r,vx,vy,eint,p,c,cx,cy,denom) 
  for(i=0; i<nx*ny; i++){
    r   =fmax(mesh[i+nx*ny*VARRHO],smallr);
    vx  =     mesh[i+nx*ny*VARVX ]/r;
    vy  =     mesh[i+nx*ny*VARVY ]/r;
    eint=     mesh[i+nx*ny*VARPR ]-0.5*r*(vx*vx+vy*vy);
    p   =fmax((gamma-1.0)*eint,r*smallp);
    
    c=sqrt((gamma*p/r));
    cx=(c+fabs(vx))/dx;
    cy=(c+fabs(vy))/dy;
    denom=cx+cy;
    max_denom=fmax(denom,max_denom);
  }
  return 0.5/max_denom;
}

void toPrimX(double *q, double *mesh){
  int i;
  int xI, yI;
  double r,vx,vy,eint,p;
  double gamma, smallr, smallp;

  gamma=Hp->gamma;
  smallr=Ha->smallr;
  smallp=Ha->smallc*Ha->smallc/Hp->gamma;

#pragma acc kernels \
  pcopyin(mesh[0:meshSize])\
  pcopyout(q[0:primSize])\
  copyin(nx,ny,gamma,smallr,smallp)
#pragma acc loop independent private(i,xI,yI,r,vx,vy,eint,p)
  for(i=0;i<ny*nx;i++){
    xI=i%nx;
    yI=i/nx;
    r   =fmax(mesh[i+nx*ny*VARRHO],smallr);
    vx  =mesh[i+nx*ny*VARVX ]/r;
    vy  =mesh[i+nx*ny*VARVY ]/r;
    eint=mesh[i+nx*ny*VARPR ]-0.5*r*(vx*vx+vy*vy);
    p   =fmax((gamma-1)*r*eint,smallp);
    q[xI+2+(nx+4)*(yI+ny*VARRHO)]=1.0;//r;
    q[xI+2+(nx+4)*(yI+ny*VARVX )]=vx;
    q[xI+2+(nx+4)*(yI+ny*VARVY )]=vy;
    q[xI+2+(nx+4)*(yI+ny*VARPR )]=p;
  }
}

void toPrimY(double *q, double *mesh){
  int i;
  int xI, yI;
  double r,vx,vy,eint,p;
  double gamma, smallr, smallp;

  gamma=Hp->gamma;
  smallr=Ha->smallr;
  smallp=Ha->smallc*Ha->smallc/Hp->gamma;

#pragma acc kernels pcopyin(mesh[0:meshSize])\
  pcopyout(q[0:primSize])\
  copyin(gamma,smallr,smallp)
#pragma acc loop independent private(xI,yI,r,vx,vy,eint,p)
  for(i=0;i<ny*nx;i++){
    xI=i%nx;
    yI=i/nx;
    r   =fmax(mesh[i+nx*ny*VARRHO],smallr);
    vx  =mesh[i+nx*ny*VARVX ]/r;
    vy  =mesh[i+nx*ny*VARVY ]/r;
    eint=mesh[i+nx*ny*VARPR ]-0.5*r*(vx*vx+vy*vy);
    p   =fmax((gamma-1)*r*eint,smallp);
    q[yI+2+(ny+4)*(xI+nx*VARRHO)]=r;
    q[yI+2+(ny+4)*(xI+nx*VARVX )]=vy;
    q[yI+2+(ny+4)*(xI+nx*VARVY )]=vx;
    q[yI+2+(ny+4)*(xI+nx*VARPR )]=p;
  }
}

void setBndCnd(double* q, int cndL, int cndH, int np, int nt){
  int i;
  int pI,tI;
  int wInd, rInd;


  //printf("Running bnd cnds for %dx%d prims\n",np,nt);

#pragma acc kernels \
  pcopy(q[0:primSize])
#pragma acc loop independent private(pI,tI,wInd,rInd)
  for(i=0;i<2*nt;i++){
    pI=i%2;
    tI=i/2;
    wInd=pI+(np+4)*tI;
    rInd=3-pI+(np+4)*tI;
    //printf("BND: %d refs %d\n",wInd,rInd);
    if(cndL==BND_REFL){
      q[wInd+(np+4)*nt*VARRHO]= q[rInd+(np+4)*nt*VARRHO];
      q[wInd+(np+4)*nt*VARVX ]=-q[rInd+(np+4)*nt*VARVX ];
      q[wInd+(np+4)*nt*VARVY ]= q[rInd+(np+4)*nt*VARVY ];
      q[wInd+(np+4)*nt*VARPR ]= q[rInd+(np+4)*nt*VARPR ];
    }else if(cndL==BND_PERM){
      q[wInd+(np+4)*nt*VARRHO]= q[rInd+(np+4)*nt*VARRHO];
      q[wInd+(np+4)*nt*VARVX ]= q[rInd+(np+4)*nt*VARVX ];
      q[wInd+(np+4)*nt*VARVY ]= q[rInd+(np+4)*nt*VARVY ];
      q[wInd+(np+4)*nt*VARPR ]= q[rInd+(np+4)*nt*VARPR ];
    }
    wInd=np+2+pI+(np+4)*tI;
    rInd=np+1-pI+(np+4)*tI;
    //printf("BND: %d refs %d\n",wInd,rInd);
    if(cndH==BND_REFL){
      q[wInd+(np+4)*nt*VARRHO]= q[rInd+(np+4)*nt*VARRHO];
      q[wInd+(np+4)*nt*VARVX ]=-q[rInd+(np+4)*nt*VARVX ];
      q[wInd+(np+4)*nt*VARVY ]= q[rInd+(np+4)*nt*VARVY ];
      q[wInd+(np+4)*nt*VARPR ]= q[rInd+(np+4)*nt*VARPR ];
    }else if(cndH==BND_PERM){
      q[wInd+(np+4)*nt*VARRHO]= q[rInd+(np+4)*nt*VARRHO];
      q[wInd+(np+4)*nt*VARVX ]= q[rInd+(np+4)*nt*VARVX ];
      q[wInd+(np+4)*nt*VARVY ]= q[rInd+(np+4)*nt*VARVY ];
      q[wInd+(np+4)*nt*VARPR ]= q[rInd+(np+4)*nt*VARPR ];
    }
  }
}

void trace(double *ql, double *qr, double *q, double dtdx, int np, int nt){
  int lI;
  int i,j;
  double  r, u, v1, p, a;
  double dlft, drgt, dcen, dsgn, dlim;
  double dr,du,dv1,dp,da;
  double cc,csq;
  double alpham,alphap,alphazr;
  double spplus,spzerol,spzeror,spminus;
  double ap,am,azr,azv1,acmp;
  double gamma;

  gamma=Hp->gamma;

#pragma acc kernels pcopyin(q[0:primSize]) pcopyout(qr[0:qSize],ql[0:qSize]) copyin(gamma,dtdx,np,nt)
#pragma acc loop independent private(i,j,r,u,v1,p,csq,cc,dr,du,dv1,dp,ap,am,azr,azv1,acmp)
  for(lI=0;lI<(np+2)*nt;lI++){
    i=lI%(np+2);
    j=lI/(np+2);
    r =q[i+1+(np+4)*(j+nt*VARRHO)];
    u =q[i+1+(np+4)*(j+nt*VARVX )];
    v1=q[i+1+(np+4)*(j+nt*VARVY )];
    p =q[i+1+(np+4)*(j+nt*VARPR )];
    
    csq=gamma*p/r;
    cc=sqrt(csq);
    
    dlft=q[i+1+(np+4)*(j+nt*VARRHO)  ]-q[i+1+(np+4)*(j+nt*VARRHO)-1];
    drgt=q[i+1+(np+4)*(j+nt*VARRHO)+1]-q[i+1+(np+4)*(j+nt*VARRHO)  ];
    dcen=0.5*(dlft+drgt);
    dlim=(dlft*drgt<=0)?0.0:MIN(fabs(dlft),fabs(drgt));
    dr= (dcen>=0?1.0:-1.0)*MIN(dlim,fabs(dcen));

    dlft=q[i+1+(np+4)*(j+nt*VARVX )  ]-q[i+1+(np+4)*(j+nt*VARVX )-1];
    drgt=q[i+1+(np+4)*(j+nt*VARVX )+1]-q[i+1+(np+4)*(j+nt*VARVX )  ];
    dcen=0.5*(dlft+drgt);
    dlim=(dlft*drgt<=0)?0.0:MIN(fabs(dlft),fabs(drgt));
    du= (dcen>=0?1.0:-1.0)*MIN(dlim,fabs(dcen));

    dlft=q[i+1+(np+4)*(j+nt*VARVY )  ]-q[i+1+(np+4)*(j+nt*VARVY )-1];
    drgt=q[i+1+(np+4)*(j+nt*VARVY )+1]-q[i+1+(np+4)*(j+nt*VARVY )  ];
    dcen=0.5*(dlft+drgt);
    dlim=(dlft*drgt<=0)?0.0:MIN(fabs(dlft),fabs(drgt));
    dv1= (dcen>=0?1.0:-1.0)*MIN(dlim,fabs(dcen));

    dlft=q[i+1+(np+4)*(j+nt*VARPR )  ]-q[i+1+(np+4)*(j+nt*VARPR )-1];
    drgt=q[i+1+(np+4)*(j+nt*VARPR )+1]-q[i+1+(np+4)*(j+nt*VARPR )  ];
    dcen=0.5*(dlft+drgt);
    dlim=(dlft*drgt<=0)?0.0:MIN(fabs(dlft),fabs(drgt));
    dp= (dcen>=0?1.0:-1.0)*MIN(dlim,fabs(dcen));
    
    alpham = 0.5*(dp/(r*cc)-du)*r/cc;
    alphap = 0.5*(dp/(r*cc)+du)*r/cc;
    alphazr= dr-dp/csq;
    
    //right
    spminus=((u-cc)>=0.0)?0.0:(u-cc)*dtdx+1.0;
    spzeror=((u   )>=0.0)?0.0:(u   )*dtdx+1.0;
    spplus =((u+cc)>=0.0)?0.0:(u+cc)*dtdx+1.0;
    ap  =-0.5*spplus *alphap ;
    am  =-0.5*spminus*alpham ;
    azr =-0.5*spzeror*alphazr;
    azv1=-0.5*spzeror*dv1;
    qr[i+(np+2)*(j+nt*VARRHO)]=r +(ap+am+azr);
    qr[i+(np+2)*(j+nt*VARVX )]=u +(ap-am    )*cc/r;
    qr[i+(np+2)*(j+nt*VARVY )]=v1+(azv1     );
    qr[i+(np+2)*(j+nt*VARPR )]=p +(ap+am    )*csq;
    
    //left
    spminus=((u-cc)<=0.0)?0.0:(u-cc)*dtdx-1.0;
    spzerol=((u   )<=0.0)?0.0:(u   )*dtdx-1.0;
    spplus =((u+cc)<=0.0)?0.0:(u+cc)*dtdx-1.0;
    ap  =-0.5*spplus *alphap ;
    am  =-0.5*spminus*alpham ;
    azr =-0.5*spzerol*alphazr;
    azv1=-0.5*spzerol*dv1;
    ql[i+(np+2)*(j+nt*VARRHO)]=r +(ap+am+azr);
    ql[i+(np+2)*(j+nt*VARVX )]=u +(ap-am    )*cc/r;
    ql[i+(np+2)*(j+nt*VARVY )]=v1+(azv1     );
    ql[i+(np+2)*(j+nt*VARPR )]=p +(ap+am    )*csq;
  }
}

double slope(double *q,int ind){
  double dlft, drgt, dcen, dsgn, dlim;
  dlft=q[ind  ]-q[ind-1];
  drgt=q[ind+1]-q[ind  ];
  dcen=0.5*(dlft+drgt);
  if(dcen>=0)dsgn=1.0;
  else dsgn=-1.0;
  if(dlft*drgt<=0)dlim=0.0;
  else dlim=MIN(fabs(dlft),fabs(drgt));
  return dsgn*MIN(dlim,fabs(dcen));
}

void riemann(double *flx, double *qxm, double *qxp, int np, int nt){
  int lI, i,j,n;
  double smallr, smallc, smallp, smallpp;
  double gamma, gmma6, entho;
  double qgdnvR,qgdnvVX,qgdnvVY,qgdnvP;
  double rl,vxl,vyl,pl,cl,wl,ql,vsl;
  double rr,vxr,vyr,pr,cr,wr,qr,vsr;
  double ro,vxo,po,wo,co;
  double rx,vxx,px,wx,cx;
  double sgnm, scr, frac;
  double spout,spin,ushk;
  double ekin,etot,delp;
  int niter_riemann;

  niter_riemann=Ha->niter_riemann;
  gamma=Hp->gamma;
  smallr=Ha->smallr;
  smallc=Ha->smallc;
  smallp=Ha->smallc*Ha->smallc/Hp->gamma;
  smallpp=Ha->smallr*smallp;
  gmma6=(Hp->gamma+1.0)/(2.0*Hp->gamma);
  entho=1.0/(Hp->gamma-1.0);

#pragma acc kernels pcopyin(qxm[0:qSize],qxp[0:qSize]) pcopyout(flx[0:flxSize]) pcopyin(np,nt,gamma,smallr,smallc,smallp,smallpp,gmma6,entho)
#pragma acc loop independent
  for(lI=0;lI<(np+1)*nt;lI++){
    i=lI%(np+1);
    j=lI/(np+1);
    
    rl =MAX(qxm[i  +(np+2)*(j+nt*VARRHO)],smallr);
    vxl=    qxm[i  +(np+2)*(j+nt*VARVX )];
    vyl=    qxm[i  +(np+2)*(j+nt*VARVY )];
    pl =MAX(qxm[i  +(np+2)*(j+nt*VARPR )],rl*smallp);

    rr =MAX(qxp[i+1+(np+2)*(j+nt*VARRHO)],smallr);
    vxr=    qxp[i+1+(np+2)*(j+nt*VARVX )];
    vyr=    qxp[i+1+(np+2)*(j+nt*VARVY )];
    pr =MAX(qxp[i+1+(np+2)*(j+nt*VARPR )],rr*smallp);
    
    cl=gamma*pl*rl;
    cr=gamma*pr*rr;
    
    wl=sqrt(cl);
    wr=sqrt(cr);
    
    px=((wr*pl+wl*pr)+wl*wr*(vxl-vxr))/(wl+wr);
    px=MAX(px,0.0);
    for(n=0;n<niter_riemann;n++){
      wl=sqrt(cl*(1.0+gmma6*(px-pl)/pl));
      wr=sqrt(cr*(1.0+gmma6*(px-pr)/pr));
      ql=2.0*wl*wl*wl/(wl*wl+cl);
      qr=2.0*wr*wr*wr/(wr*wr+cr);
      vsl=vxl-(px-pl)/wl;
      vsr=vxr+(px-pr)/wr;
      delp=qr*ql/(qr+ql)*(vsl-vsr);
      delp=MAX(delp,-px);
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
    co=MAX(smallc,sqrt(fabs(gamma*po/ro)));
    
    rx=MAX(smallr,ro/(1.0+ro*(po-px)/(wo*wo)));
    
    cx=MAX(smallc,sqrt(fabs(gamma*px/rx)));
      
    spout=co   -sgnm*vxo;
    spin =cx   -sgnm*vxx;
    ushk =wo/ro-sgnm*vxo;
    
    if(px>=po){
      spin=ushk;
      spout=ushk;
    }

    scr=MAX(spout-spin,smallc+fabs(spout+spin));

    frac=0.5*(1.0+(spout+spin)/scr);
    frac=MAX(0.0,MIN(1.0,frac));
    qgdnvR =frac* rx+(1.0-frac)* ro;
    qgdnvVX=frac*vxx+(1.0-frac)*vxo;
    qgdnvP =frac* px+(1.0-frac)* po;
    if(spout<0.0){
      qgdnvR =ro;
      qgdnvVX=vxo;
      qgdnvP =po;
    }
    if(spin>0.0){
      qgdnvR =rx;
      qgdnvVX=vxx;
      qgdnvP =px;
    }
    flx[i+(np+1)*(j+nt*VARRHO)]=0.025;//qgdnvR*qgdnvVX;
    flx[i+(np+1)*(j+nt*VARVX )]=qgdnvR*qgdnvVX*qgdnvVX+qgdnvP;
    flx[i+(np+1)*(j+nt*VARVY )]=qgdnvR*qgdnvVX*qgdnvVY;
    ekin=0.5*qgdnvR*(qgdnvVX*qgdnvVX+qgdnvVY*qgdnvVY);
    etot=qgdnvP*entho+ekin;
    flx[i+(np+1)*(j+nt*VARPR )]=qgdnvVX*(etot+qgdnvP);
  }
}

void addFluxX(double *mesh, double *flx, double dtdx, int np, int nt){
  int lI, i, j;

#pragma acc kernels pcopy(mesh[0:meshSize]) pcopyin(flx[0:flxSize]) pcopyin(dtdx,np,nt)
#pragma acc loop independent
  for(lI=0;lI<np*nt;lI++){
    i=lI%np;
    j=lI/np;
    mesh[i+np*(j+nt*VARRHO)]+=dtdx*(flx[i  +(np+1)*(j+nt*VARRHO)]-
				    flx[i+1+(np+1)*(j+nt*VARRHO)]);
    mesh[i+np*(j+nt*VARVX )]+=dtdx*(flx[i  +(np+1)*(j+nt*VARVX )]-
				    flx[i+1+(np+1)*(j+nt*VARVX )]);
    mesh[i+np*(j+nt*VARVY )]+=dtdx*(flx[i  +(np+1)*(j+nt*VARVY )]-
				    flx[i+1+(np+1)*(j+nt*VARVY )]);
    mesh[i+np*(j+nt*VARPR )]+=dtdx*(flx[i  +(np+1)*(j+nt*VARPR )]-
				    flx[i+1+(np+1)*(j+nt*VARPR )]);
  }
}

void addFluxY(double *mesh, double *flx, double dtdx, int np, int nt){
  int lI, i, j;

#pragma acc kernels pcopy(mesh[0:meshSize]) pcopyin(flx[0:flxSize]) pcopyin(dtdx,np,nt)
#pragma acc loop independent
  for(lI=0;lI<np*nt;lI++){
    i=lI%np;
    j=lI/np;
    mesh[j+nt*(i+np*VARRHO)]+=dtdx*(flx[i  +(np+1)*(j+nt*VARRHO)]-
				    flx[i+1+(np+1)*(j+nt*VARRHO)]);
    mesh[j+nt*(i+np*VARVX )]+=dtdx*(flx[i  +(np+1)*(j+nt*VARVY )]-
				    flx[i+1+(np+1)*(j+nt*VARVY )]);
    mesh[j+nt*(i+np*VARVY )]+=dtdx*(flx[i  +(np+1)*(j+nt*VARVX )]-
				    flx[i+1+(np+1)*(j+nt*VARVX )]);
    mesh[j+nt*(i+np*VARPR )]+=dtdx*(flx[i  +(np+1)*(j+nt*VARPR )]-
				    flx[i+1+(np+1)*(j+nt*VARPR )]);
  }
}

int nansIn(double *mesh, int var, int nx, int ny, int nHx, int nHy);

void nanScan(int *ret, double *mesh, int nvar, int nx, int ny, int nHx, int nHy)
{
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

double sumArray(double *mesh, int var, int nx, int ny){
  int i;
  double sum, corr;
  double nsum, c_nxt;

  sum=0.0;
  corr=0.0;

#pragma acc kernels \
  pcopyin(mesh[0:meshSize])
#pragma acc loop independent reduction(+:sum,corr)
  for(i=0;i<nx*ny;i++){
    c_nxt=mesh[i+nx*ny*var];
    nsum=sum+c_nxt;
    corr=(nsum-sum)-c_nxt;
    sum=nsum;
  }
  return sum+corr;
}

void runPass(double *mesh, double dt, int n, int dir){
  int bndL,bndH;
  int np,nt;
  double dx,dy;
  char dirCh, outfile[30];
  int nanList[4];

  if(dir==0){
    np=nx;
    nt=ny;
    bndL=Hp->bndL;
    bndH=Hp->bndR;
    dx=Hp->dx;
    dy=Hp->dy;
    dirCh='x';
  }else{
    np=ny;
    nt=nx;
    bndL=Hp->bndU;
    bndH=Hp->bndD;
    dx=Hp->dy;
    dy=Hp->dx;
    dirCh='y';
  }
  //printf("Passdir: %c\n",dirCh);
  //printArray("PRE :",mesh,4,nx,ny,0,0);
  //nanScan(nanList,mesh,4,nx,ny,0,0);
  //printf("Pre-pass mesh: %d %d %d %d\n",nanList[0],nanList[1],nanList[2],nanList[3]);
  if(dir==0){
    toPrimX(q,mesh);
  }else{
    toPrimY(q,mesh);
  }
  //printf("Primitives done\n");
  //#pragma acc update host(q[0:primSize])
  //nanScan(nanList,q,4,np,nt,2,0);
  //printArray("Q   :",q,4,np,nt,2,0);
  //printf("Prim: %d %d %d %d\n",nanList[0],nanList[1],nanList[2],nanList[3]);
  setBndCnd(q,bndL,bndH,np,nt);
  //printf("Bnd cnds set\n");
  //#pragma acc update host(q[0:primSize])
  //nanScan(nanList,q,4,np+4,nt,0,0);
  //printArray("QBND:",q,4,np+4,nt,0,0);
  //printf("Prim(bnd): %d %d %d %d\n",nanList[0],nanList[1],nanList[2],nanList[3]);
  trace(ql,qr,q,dt/dx,np,nt);
  //#pragma acc update host(qr[0:qSize],ql[0:qSize])
  //printf("Trace complete\n");
  //nanScan(nanList,ql,4,np+2,nt,0,0);
  //printArray("QL  :",ql,4,np+2,nt,0,0);
  //printf("QL: %d %d %d %d\n",nanList[0],nanList[1],nanList[2],nanList[3]);
  //nanScan(nanList,qr,4,np+2,nt,0,0);
  //printArray("QR  :",qr,4,np+2,nt,0,0);
  //printf("QR: %d %d %d %d\n",nanList[0],nanList[1],nanList[2],nanList[3]);
  riemann(flx,ql,qr,np,nt);
  //printf("Riemann and flx computation complete\n");
  //#pragma acc update host(flx[0:flxSize])
  //nanScan(nanList,flx,4,np+1,nt,0,0);
  //printArray("FLX :",flx,4,np+1,nt,0,0);
  //printf("FLX: %d %d %d %d\n",nanList[0],nanList[1],nanList[2],nanList[3]);
  if(dir==0){
    addFluxX(mesh,flx,dt/dx,np,nt);
  }else{
    addFluxY(mesh,flx,dt/dx,np,nt);
  }
  //printf("Flx added\n");
  //#pragma acc update host(mesh[0:meshSize])
  //nanScan(nanList,mesh,4,nx,ny,0,0);
  //printArray("POST:",mesh,4,nx,ny,0,0);
  //printf("Post-pass mesh: %d %d %d %d\n",nanList[1],nanList[1],nanList[2],nanList[3]);
}

void engine(double *mesh, hydro_prob *Hyp, hydro_args *Hya){
  int n, i,j;
  int bndL;
  int bndH;
  double dt;
  double cTime, nxttout;

  double volCell;
  double oTM,oTE;
  double TM,TE;
  int M_exp, E_exp;
  double M_prec, E_prec;

  char outfile[30];

  double initT, endT;

  Hp=Hyp;
  Ha=Hya;
  nx=Hp->nx;
  ny=Hp->ny;

  n=0;
  cTime=0;
  nxttout=-1.0;

#ifdef __OPENACC
  acc_init(acc_device_default);
#endif  

  meshSize=Hp->nvar*nx*ny;
  if(Hp->ny>=Hp->nx){
    primSize=Hp->nvar*(nx+4)*ny;
    qSize   =Hp->nvar*(nx+2)*ny;
    flxSize =Hp->nvar*(nx+1)*ny;
  }else{
    primSize=Hp->nvar*(ny+4)*nx;
    qSize   =Hp->nvar*(ny+2)*nx;
    flxSize =Hp->nvar*(ny+1)*nx;
  }

  if(Ha->nstepmax<0&&Ha->tend<0.0)return;

  q  =(double*)malloc(primSize*sizeof(double));
  qr =(double*)malloc(qSize*sizeof(double));
  ql =(double*)malloc(qSize*sizeof(double));
  flx=(double*)malloc(flxSize*sizeof(double));

  if(Ha->tend>0.0){
    nxttout=Ha->tend;
  }
  if(Ha->dtoutput>0.0&&nxttout>Ha->dtoutput){
    nxttout=Ha->dtoutput;
  }

  volCell=Hp->dx*Hp->dy;
  oTM=sumArray(mesh,VARRHO,Hp->nx,Hp->ny);
  oTE=sumArray(mesh,VARPR ,Hp->nx,Hp->ny);
#ifdef M_PREC_CMP
  frexp(oTM,&M_exp);
  frexp(oTE,&E_exp);

  M_prec=ldexp(DBL_EPSILON,M_exp-1);
  E_prec=ldexp(DBL_EPSILON,E_exp-1);

  printf("TM:%g+-%g TE:%g+-%g\n",volCell*oTM,volCell*M_prec,volCell*oTE,volCell*E_prec);
#endif

  //Print initial condition
  snprintf(outfile,29,"%s%05d",Ha->outPre,n);
  writeVis(outfile,mesh,Hp->dx,Hp->dy,Hp->nvar,Hp->nx,Hp->ny);
  printf("INIT: TM: %g TE: %g\n",volCell*oTM,volCell*oTE);
  
  initT=getNow();

#pragma acc data copy(mesh[0:meshSize]) create(q[0:primSize],qr[0:qSize],ql[0:qSize],flx[0:flxSize])
  {
    while((n<Ha->nstepmax||Ha->nstepmax<0)&&(cTime<Ha->tend||Ha->tend<0)){
      //Calculate timestep
      //printf("Timestep calculation\n");
      dt=Ha->sigma*calcDT(mesh);
      //printf("DT=%g\n",dt);
      if(nxttout>0.0&&dt>(nxttout-cTime)){
	printf("Adjusting timestep from %g to %g for iter %d\n",dt,nxttout-cTime,n);
	dt=(nxttout-cTime);
      }
      if(n%2==0){
	//X Dir
	runPass(mesh,dt,n,0);
	//Y Dir
	runPass(mesh,dt,n,1);
      }else{
	//Y Dir
	runPass(mesh,dt,n,1);
	//X Dir
	runPass(mesh,dt,n,0);
      }
      n+=1;
      cTime+=dt;
      if(n%Ha->nprtLine==0){
	TM=sumArray(mesh,VARRHO,Hp->nx,Hp->ny);
	TE=sumArray(mesh,VARPR ,Hp->nx,Hp->ny);
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
	  //printf("Next Vis Time: %f\n",nxttout);
	}
        //#pragma acc update host(mesh[0:meshSize])
	snprintf(outfile,29,"%s%05d",Ha->outPre,n);
	writeVis(outfile,mesh,Hp->dx,Hp->dy,Hp->nvar,Hp->nx,Hp->ny);
        printf("Vis. file \"%s\" written.\n",outfile);
      }
    }
  }
  printf("time: %f, %d iters run\n",cTime,n);

  endT=getNow();

  printf("TFMT:%s,%s,%s,%s,%s,%s,%s,%s\n","cType","mType","init","nproc","nth","niters","ncells","wRunt");
  printf("TIME:%s,%s,%s,%d,%d,%d,%d,%g\n","\"OAC\"","\"GPU:?\"","\"Init\"",1,1,n,Hp->nx*Hp->ny,endT-initT);

  //Print final condition
  snprintf(outfile,29,"%s%05d",Ha->outPre,n);
  writeVis(outfile,mesh,Hp->dx,Hp->dy,Hp->nvar,Hp->nx,Hp->ny);

  Hp->t+=cTime;
  free(q  );
  free(qr );
  free(ql );
  free(flx);
}
