#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <mpi.h>
#include "hydro.h"
#include "float.h"

hydro_args *Ha;
hydro_prob *Hp;
double *q, *bndLS, *bndLR, *bndHS, *bndHR;
double *qr, *ql;
double *flx;

//MPI Vars
int bndT, bndB;
int pProc, nProc;
int rank, size;
int myNy;
int varSize;

double slope(double *q,int ind);

void printArray(char* label,double *arr, int nvar, int nx, int ny){
  int nV,i,j;
  //printf("N[%2d]: Array %s\n",rank,label);
  for(nV=0;nV<nvar;nV++){
    //printf("N[%2d]:Variable %d:\n",rank,nV);
    //printf("N[%2d]:%3s",rank,"");
    for(i=0;i<nx;i++){
      //printf("|i%8di",i);
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

double calcDT(double *mesh){
  int lI, i, j;
  double denom, max_denom;
  double gmax;
  double r,vx,vy,eint,p;
  double c,cx,cy;
  double smallp;

  c=0.0;
  max_denom=Ha->smallc;

  smallp=Ha->smallc*Ha->smallc/Hp->gamma;
  for (lI=0; lI<Hp->nx*myNy; lI++){
    i=lI%Hp->nx;
    j=lI/Hp->nx;
    r   =MAX(mesh[i+2+(j+2)*(Hp->nx+4)+varSize*VARRHO],Ha->smallr);
    vx  =    mesh[i+2+(j+2)*(Hp->nx+4)+varSize*VARVX ]/r;
    vy  =    mesh[i+2+(j+2)*(Hp->nx+4)+varSize*VARVY ]/r;
    eint=    mesh[i+2+(j+2)*(Hp->nx+4)+varSize*VARPR ]-0.5*r*(vx*vx+vy*vy);
    p   =MAX((Hp->gamma-1.0)*eint,r*smallp);
    
    c=sqrt((Hp->gamma*p/r));
    cx=(c+fabs(vx))/Hp->dx;
    cy=(c+fabs(vy))/Hp->dy;
    denom=cx+cy;
    if(max_denom<denom)max_denom=denom;
  }
  //printf("NODE %d: denom %g\n",rank,max_denom);
  MPI_Allreduce(&max_denom,&gmax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  //printf("NODE %d: gdenom %g\n",rank,gmax);
  return 0.5/gmax;
}

void setHHalo(double *mesh, int LBnd, int RBnd){
  int lI, i,j;
  for(lI=0;lI<2*myNy;lI++){
    i=lI%2;
    j=lI/2;
    //Left Boundary
    if(LBnd==BND_REFL){
      mesh[i+(Hp->nx+4)*(j+2+(myNy+4)*VARRHO)]= mesh[3-i+(Hp->nx+4)*(j+2+(myNy+4)*VARRHO)];
      mesh[i+(Hp->nx+4)*(j+2+(myNy+4)*VARVX )]=-mesh[3-i+(Hp->nx+4)*(j+2+(myNy+4)*VARVX )];
      mesh[i+(Hp->nx+4)*(j+2+(myNy+4)*VARVY )]= mesh[3-i+(Hp->nx+4)*(j+2+(myNy+4)*VARVY )];
      mesh[i+(Hp->nx+4)*(j+2+(myNy+4)*VARPR )]= mesh[3-i+(Hp->nx+4)*(j+2+(myNy+4)*VARPR )];
    }else if(LBnd==BND_PERM){
      mesh[i+(Hp->nx+4)*(j+2+(myNy+4)*VARRHO)]= mesh[3-i+(Hp->nx+4)*(j+2+(myNy+4)*VARRHO)];
      mesh[i+(Hp->nx+4)*(j+2+(myNy+4)*VARVX )]= mesh[3-i+(Hp->nx+4)*(j+2+(myNy+4)*VARVX )];
      mesh[i+(Hp->nx+4)*(j+2+(myNy+4)*VARVY )]= mesh[3-i+(Hp->nx+4)*(j+2+(myNy+4)*VARVY )];
      mesh[i+(Hp->nx+4)*(j+2+(myNy+4)*VARPR )]= mesh[3-i+(Hp->nx+4)*(j+2+(myNy+4)*VARPR )];
    }
    //Right Boundary
    if(RBnd==BND_REFL){
      mesh[Hp->nx+2+i+(Hp->nx+4)*(j+2+(myNy+4)*VARRHO)]= mesh[Hp->nx+i+(Hp->nx+4)*(j+2+(myNy+4)*VARRHO)];
      mesh[Hp->nx+2+i+(Hp->nx+4)*(j+2+(myNy+4)*VARVX )]=-mesh[Hp->nx+i+(Hp->nx+4)*(j+2+(myNy+4)*VARVX )];
      mesh[Hp->nx+2+i+(Hp->nx+4)*(j+2+(myNy+4)*VARVY )]= mesh[Hp->nx+i+(Hp->nx+4)*(j+2+(myNy+4)*VARVY )];
      mesh[Hp->nx+2+i+(Hp->nx+4)*(j+2+(myNy+4)*VARPR )]= mesh[Hp->nx+i+(Hp->nx+4)*(j+2+(myNy+4)*VARPR )];
    }else if(RBnd==BND_PERM){
      mesh[Hp->nx+2+i+(Hp->nx+4)*(j+2+(myNy+4)*VARRHO)]= mesh[Hp->nx+i+(Hp->nx+4)*(j+2+(myNy+4)*VARRHO)];
      mesh[Hp->nx+2+i+(Hp->nx+4)*(j+2+(myNy+4)*VARVX )]= mesh[Hp->nx+i+(Hp->nx+4)*(j+2+(myNy+4)*VARVX )];
      mesh[Hp->nx+2+i+(Hp->nx+4)*(j+2+(myNy+4)*VARVY )]= mesh[Hp->nx+i+(Hp->nx+4)*(j+2+(myNy+4)*VARVY )];
      mesh[Hp->nx+2+i+(Hp->nx+4)*(j+2+(myNy+4)*VARPR )]= mesh[Hp->nx+i+(Hp->nx+4)*(j+2+(myNy+4)*VARPR )];
    }
  }
}

void setVHalo(double *mesh, int TBnd, int BBnd){
  int lI, i,j;
  int row=Hp->nx+4;
  MPI_Request reqs[16];
  MPI_Status stat[16];
  MPI_Irecv(mesh+(     0)*row+VARRHO*varSize,2*row,MPI_DOUBLE,pProc,1,MPI_COMM_WORLD,reqs+ 0);
  MPI_Irecv(mesh+(     0)*row+VARVX *varSize,2*row,MPI_DOUBLE,pProc,3,MPI_COMM_WORLD,reqs+ 1);
  MPI_Irecv(mesh+(     0)*row+VARVY *varSize,2*row,MPI_DOUBLE,pProc,5,MPI_COMM_WORLD,reqs+ 2);
  MPI_Irecv(mesh+(     0)*row+VARPR *varSize,2*row,MPI_DOUBLE,pProc,7,MPI_COMM_WORLD,reqs+ 3);
  MPI_Isend(mesh+(     2)*row+VARRHO*varSize,2*row,MPI_DOUBLE,pProc,2,MPI_COMM_WORLD,reqs+ 4);
  MPI_Isend(mesh+(     2)*row+VARVX *varSize,2*row,MPI_DOUBLE,pProc,4,MPI_COMM_WORLD,reqs+ 5);
  MPI_Isend(mesh+(     2)*row+VARVY *varSize,2*row,MPI_DOUBLE,pProc,6,MPI_COMM_WORLD,reqs+ 6);
  MPI_Isend(mesh+(     2)*row+VARPR *varSize,2*row,MPI_DOUBLE,pProc,8,MPI_COMM_WORLD,reqs+ 7);
  MPI_Isend(mesh+(myNy  )*row+VARRHO*varSize,2*row,MPI_DOUBLE,nProc,1,MPI_COMM_WORLD,reqs+ 8);
  MPI_Isend(mesh+(myNy  )*row+VARVX *varSize,2*row,MPI_DOUBLE,nProc,3,MPI_COMM_WORLD,reqs+ 9);
  MPI_Isend(mesh+(myNy  )*row+VARVY *varSize,2*row,MPI_DOUBLE,nProc,5,MPI_COMM_WORLD,reqs+10);
  MPI_Isend(mesh+(myNy  )*row+VARPR *varSize,2*row,MPI_DOUBLE,nProc,7,MPI_COMM_WORLD,reqs+11);
  MPI_Irecv(mesh+(myNy+2)*row+VARRHO*varSize,2*row,MPI_DOUBLE,nProc,2,MPI_COMM_WORLD,reqs+12);
  MPI_Irecv(mesh+(myNy+2)*row+VARVX *varSize,2*row,MPI_DOUBLE,nProc,4,MPI_COMM_WORLD,reqs+13);
  MPI_Irecv(mesh+(myNy+2)*row+VARVY *varSize,2*row,MPI_DOUBLE,nProc,6,MPI_COMM_WORLD,reqs+14);
  MPI_Irecv(mesh+(myNy+2)*row+VARPR *varSize,2*row,MPI_DOUBLE,nProc,8,MPI_COMM_WORLD,reqs+15);
  MPI_Waitall(16,reqs,stat);
  for(lI=0;lI<2*Hp->nx;lI++){
    i=lI/2;
    j=lI%2;
    //Top boundary
    if(TBnd==BND_REFL){
      mesh[i+2+(Hp->nx+4)*(j+(myNy+4)*VARRHO)]= mesh[i+2+(Hp->nx+4)*(3-j+(myNy+4)*VARRHO)];
      mesh[i+2+(Hp->nx+4)*(j+(myNy+4)*VARVX )]= mesh[i+2+(Hp->nx+4)*(3-j+(myNy+4)*VARVX )];
      mesh[i+2+(Hp->nx+4)*(j+(myNy+4)*VARVY )]=-mesh[i+2+(Hp->nx+4)*(3-j+(myNy+4)*VARVY )];
      mesh[i+2+(Hp->nx+4)*(j+(myNy+4)*VARPR )]= mesh[i+2+(Hp->nx+4)*(3-j+(myNy+4)*VARPR )];
    }else if(TBnd==BND_PERM){
      mesh[i+2+(Hp->nx+4)*(j+(myNy+4)*VARRHO)]= mesh[i+2+(Hp->nx+4)*(3-j+(myNy+4)*VARRHO)];
      mesh[i+2+(Hp->nx+4)*(j+(myNy+4)*VARVX )]= mesh[i+2+(Hp->nx+4)*(3-j+(myNy+4)*VARVX )];
      mesh[i+2+(Hp->nx+4)*(j+(myNy+4)*VARVY )]= mesh[i+2+(Hp->nx+4)*(3-j+(myNy+4)*VARVY )];
      mesh[i+2+(Hp->nx+4)*(j+(myNy+4)*VARPR )]= mesh[i+2+(Hp->nx+4)*(3-j+(myNy+4)*VARPR )];
    }
    //Bottom boundary
    if(BBnd==BND_REFL){
      mesh[i+2+(Hp->nx+4)*(myNy+2+j+(myNy+4)*VARRHO)]= mesh[i+2+(Hp->nx+4)*(myNy+1-j+(myNy+4)*VARRHO)];
      mesh[i+2+(Hp->nx+4)*(myNy+2+j+(myNy+4)*VARVX )]= mesh[i+2+(Hp->nx+4)*(myNy+1-j+(myNy+4)*VARVX )];
      mesh[i+2+(Hp->nx+4)*(myNy+2+j+(myNy+4)*VARVY )]=-mesh[i+2+(Hp->nx+4)*(myNy+1-j+(myNy+4)*VARVY )];
      mesh[i+2+(Hp->nx+4)*(myNy+2+j+(myNy+4)*VARPR )]= mesh[i+2+(Hp->nx+4)*(myNy+1-j+(myNy+4)*VARPR )];
    }else if(BBnd==BND_PERM){
      mesh[i+2+(Hp->nx+4)*(myNy+2+j+(myNy+4)*VARRHO)]= mesh[i+2+(Hp->nx+4)*(myNy+1-j+(myNy+4)*VARRHO)];
      mesh[i+2+(Hp->nx+4)*(myNy+2+j+(myNy+4)*VARVX )]= mesh[i+2+(Hp->nx+4)*(myNy+1-j+(myNy+4)*VARVX )];
      mesh[i+2+(Hp->nx+4)*(myNy+2+j+(myNy+4)*VARVY )]= mesh[i+2+(Hp->nx+4)*(myNy+1-j+(myNy+4)*VARVY )];
      mesh[i+2+(Hp->nx+4)*(myNy+2+j+(myNy+4)*VARPR )]= mesh[i+2+(Hp->nx+4)*(myNy+1-j+(myNy+4)*VARPR )];
    }
  }
}

void toPrimX(double *q, double *mesh){
  int i;
  int xI, yI;
  double r,vx,vy,eint,p;
  double smallp;

  smallp=Ha->smallc*Ha->smallc/Hp->gamma;
  for(i=0;i<myNy*(Hp->nx+4);i++){
    xI=i%(Hp->nx+4);
    yI=i/(Hp->nx+4);
    r   =MAX(mesh[xI+(yI+2)*(Hp->nx+4)+varSize*VARRHO],Ha->smallr);
    vx  =mesh[xI+(yI+2)*(Hp->nx+4)+varSize*VARVX ]/r;
    vy  =mesh[xI+(yI+2)*(Hp->nx+4)+varSize*VARVY ]/r;
    eint=mesh[xI+(yI+2)*(Hp->nx+4)+varSize*VARPR ]-0.5*r*(vx*vx+vy*vy);
    p   =MAX((Hp->gamma-1)*eint,r*smallp);
    q[xI+(Hp->nx+4)*(yI+myNy*VARRHO)]=r;
    q[xI+(Hp->nx+4)*(yI+myNy*VARVX )]=vx;
    q[xI+(Hp->nx+4)*(yI+myNy*VARVY )]=vy;
    q[xI+(Hp->nx+4)*(yI+myNy*VARPR )]=p;
  }
}

void toPrimY(double *q, double *mesh){
  int i;
  int xI, yI;
  double r,vx,vy,eint,p;
  double smallp;

  smallp=Ha->smallc*Ha->smallc/Hp->gamma;
  for(i=0;i<(myNy+4)*Hp->nx;i++){
    xI=i%Hp->nx;
    yI=i/Hp->nx;
    r   =MAX(mesh[xI+2+yI*(Hp->nx+4)+varSize*VARRHO],Ha->smallr);
    vx  =mesh[xI+2+yI*(Hp->nx+4)+varSize*VARVX ]/r;
    vy  =mesh[xI+2+yI*(Hp->nx+4)+varSize*VARVY ]/r;
    eint=mesh[xI+2+yI*(Hp->nx+4)+varSize*VARPR ]-0.5*r*(vx*vx+vy*vy);
    p   =MAX((Hp->gamma-1)*eint,r*smallp);
    q[yI+(myNy+4)*(xI+Hp->nx*VARRHO)]=r;
    q[yI+(myNy+4)*(xI+Hp->nx*VARVX )]=vy;
    q[yI+(myNy+4)*(xI+Hp->nx*VARVY )]=vx;
    q[yI+(myNy+4)*(xI+Hp->nx*VARPR )]=p;
  }
}

void trace(double *ql, double *qr, double *q, double dtdx, int np, int nt){
  int lI;
  int i,j;
  double  r, u, v1, p, a;
  double dr,du,dv1,dp,da;
  double cc,csq;
  double alpham,alphap,alphazr;
  double spplus,spzerol,spzeror,spminus;
  double ap,am,azr,azv1,acmp;

  //if(isnan(dtdx))printf("N[%2d]: dtdx isnan\n",rank);

  for(lI=0;lI<(np+2)*nt;lI++){
    i=lI%(np+2);
    j=lI/(np+2);
    r =q[i+1+(np+4)*(j+nt*VARRHO)];
    u =q[i+1+(np+4)*(j+nt*VARVX )];
    v1=q[i+1+(np+4)*(j+nt*VARVY )];
    p =q[i+1+(np+4)*(j+nt*VARPR )];
    
    csq=Hp->gamma*p/r;
    cc=sqrt(csq);

    if(isnan(cc)||cc==0.0){
      printf("N[%2d]: CC(%d,%d) %g\n\tgmma=%g,p=%g,r=%g\n",rank,i,j,cc,Hp->gamma,p,r);
    }
    dr =slope(q,i+1+(np+4)*(j+nt*VARRHO));
    du =slope(q,i+1+(np+4)*(j+nt*VARVX ));
    dv1=slope(q,i+1+(np+4)*(j+nt*VARVY ));
    dp =slope(q,i+1+(np+4)*(j+nt*VARPR ));
    
    if(isnan(dr)||isnan(du)||isnan(dv1)||isnan(dp))
      printf("N[%d]: dq nans@(%d,%d) %d %d %d %d\n",
	     rank,i,j,isnan(dr),isnan(du),isnan(dv1),isnan(dp));

    alpham = 0.5*(dp/(r*cc)-du)*r/cc;
    alphap = 0.5*(dp/(r*cc)+du)*r/cc;
    alphazr= dr-dp/csq;
    
    if(isnan(alpham)||isnan(alphap)||isnan(alphazr))
      printf("N[%2d]: alpha(%d,%d) nans %d %d %d\n",rank,i,j,isnan(alpham),isnan(alphap),isnan(alphazr));

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
  //  printf("Calc slope for %d refs: [%d,%d,%d]\n",ind,ind-1,ind,ind+1);
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
  double smallp, smallpp;
  double gmma6, entho;
  double qgdnvR,qgdnvVX,qgdnvVY,qgdnvP;
  double rl,vxl,vyl,pl,cl,wl,ql,vsl;
  double rr,vxr,vyr,pr,cr,wr,qr,vsr;
  double ro,vxo,po,wo,co;
  double rx,vxx,px,wx,cx;
  double sgnm, scr, frac;
  double spout,spin,ushk;
  double ekin,etot,delp;

  smallp=Ha->smallc*Ha->smallc/Hp->gamma;
  smallpp=Ha->smallr*smallp;
  gmma6=(Hp->gamma+1.0)/(2.0*Hp->gamma);
  entho=1.0/(Hp->gamma-1.0);
  for(lI=0;lI<(np+1)*nt;lI++){
    i=lI%(np+1);
    j=lI/(np+1);
    
    rl =MAX(qxm[i  +(np+2)*(j+nt*VARRHO)],Ha->smallr);
    vxl=    qxm[i  +(np+2)*(j+nt*VARVX )];
    vyl=    qxm[i  +(np+2)*(j+nt*VARVY )];
    pl =MAX(qxm[i  +(np+2)*(j+nt*VARPR )],rl*smallp);

    rr =MAX(qxp[i+1+(np+2)*(j+nt*VARRHO)],Ha->smallr);
    vxr=    qxp[i+1+(np+2)*(j+nt*VARVX )];
    vyr=    qxp[i+1+(np+2)*(j+nt*VARVY )];
    pr =MAX(qxp[i+1+(np+2)*(j+nt*VARPR )],rr*smallp);
    
    cl=Hp->gamma*pl*rl;
    cr=Hp->gamma*pr*rr;
    
    wl=sqrt(cl);
    wr=sqrt(cr);
    
    px=((wr*pl+wl*pr)+wl*wr*(vxl-vxr))/(wl+wr);
    px=MAX(px,0.0);
    for(n=0;n<Ha->niter_riemann;n++){
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
    co=MAX(Ha->smallc,sqrt(fabs(Hp->gamma*po/ro)));
    
    rx=MAX(Ha->smallr,ro/(1.0+ro*(po-px)/(wo*wo)));
    
    cx=MAX(Ha->smallc,sqrt(fabs(Hp->gamma*px/rx)));
      
    spout=co   -sgnm*vxo;
    spin =cx   -sgnm*vxx;
    ushk =wo/ro-sgnm*vxo;
    
    if(px>=po){
      spin=ushk;
      spout=ushk;
    }

    scr=MAX(spout-spin,Ha->smallc+fabs(spout+spin));

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
    flx[i+(np+1)*(j+nt*VARRHO)]=qgdnvR*qgdnvVX;
    flx[i+(np+1)*(j+nt*VARVX )]=qgdnvR*qgdnvVX*qgdnvVX+qgdnvP;
    flx[i+(np+1)*(j+nt*VARVY )]=qgdnvR*qgdnvVX*qgdnvVY;
    ekin=0.5*qgdnvR*(qgdnvVX*qgdnvVX+qgdnvVY*qgdnvVY);
    etot=qgdnvP*entho+ekin;
    flx[i+(np+1)*(j+nt*VARPR )]=qgdnvVX*(etot+qgdnvP);
  }
}

void addFluxX(double *mesh, double *flx, double dtdx, int np, int nt){
  int lI, i, j;

  for(lI=0;lI<np*nt;lI++){
    i=lI%np;
    j=lI/np;
    mesh[i+2+(np+4)*(j+2+(nt+4)*VARRHO)]+=dtdx*(flx[i  +(np+1)*(j+nt*VARRHO)]-
						flx[i+1+(np+1)*(j+nt*VARRHO)]);
    mesh[i+2+(np+4)*(j+2+(nt+4)*VARVX )]+=dtdx*(flx[i  +(np+1)*(j+nt*VARVX )]-
						flx[i+1+(np+1)*(j+nt*VARVX )]);
    mesh[i+2+(np+4)*(j+2+(nt+4)*VARVY )]+=dtdx*(flx[i  +(np+1)*(j+nt*VARVY )]-
						flx[i+1+(np+1)*(j+nt*VARVY )]);
    mesh[i+2+(np+4)*(j+2+(nt+4)*VARPR )]+=dtdx*(flx[i  +(np+1)*(j+nt*VARPR )]-
				    flx[i+1+(np+1)*(j+nt*VARPR )]);
  }
}

void addFluxY(double *mesh, double *flx, double dtdx, int np, int nt){
  int lI, i, j;

  for(lI=0;lI<np*nt;lI++){
    i=lI%np;
    j=lI/np;
    mesh[j+2+(nt+4)*(i+2+(np+4)*VARRHO)]+=dtdx*(flx[i  +(np+1)*(j+nt*VARRHO)]-
						flx[i+1+(np+1)*(j+nt*VARRHO)]);
    mesh[j+2+(nt+4)*(i+2+(np+4)*VARVX )]+=dtdx*(flx[i  +(np+1)*(j+nt*VARVY )]-
						flx[i+1+(np+1)*(j+nt*VARVY )]);
    mesh[j+2+(nt+4)*(i+2+(np+4)*VARVY )]+=dtdx*(flx[i  +(np+1)*(j+nt*VARVX )]-
						flx[i+1+(np+1)*(j+nt*VARVX )]);
    mesh[j+2+(nt+4)*(i+2+(np+4)*VARPR )]+=dtdx*(flx[i  +(np+1)*(j+nt*VARPR )]-
						flx[i+1+(np+1)*(j+nt*VARPR )]);
  }
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
    i=lI%(nx+2*nHx);
    j=lI/(nx+2*nHx);
    v=mesh[i+nHx+(nx+2*nHx)*(j+nHy+(ny+2*nHy)*var)];
    if(isnan(v)||isinf(v)){
      cnt++;
      printf("N[%2d]: NAN in var %d @ %d,%d: %g\n",rank,var,i,j,v);
    }
  }
  return cnt;
}

void runPass(double *mesh, double dt, int n, int dir){
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
    setHHalo(mesh,Hp->bndL,Hp->bndR);
    toPrimX(q,mesh);
  }else{
    np=myNy;
    nt=Hp->nx;
    dx=Hp->dy;
    dy=Hp->dx;
    dCh='y';
    setVHalo(mesh,bndT,bndB);
    toPrimY(q,mesh);
  }
  trace(ql,qr,q,dt/dx,np,nt);
  riemann(flx,ql,qr,np,nt);
  if(dir==0){
    addFluxX(mesh,flx,dt/dx,np,nt);
  }else{
    addFluxY(mesh,flx,dt/dx,np,nt);
  }
}

void engine(int *argc, char **argv[], double *gMesh, hydro_prob *Hyp, hydro_args *Hya){
  int n, nV, lI, i,j;
  int bndL;
  int bndH;
  double dt;
  double cTime, nxttout;

  double volCell;
  double oTM,oTE;
  double TM,TE;
  int M_exp, E_exp;
  double M_prec, E_prec;
  double *lMesh;
  double *recvMesh;

  char outfile[30];
  char outLab[30];

  size_t primSize, qSize, flxSize;

  double initT, endT;

  int mpi_err;

  int *counts, *dspls;

  mpi_err=MPI_Init(argc,argv);

  if(mpi_err!=MPI_SUCCESS){
    printf("Error initializing MPI\n");
  }

  Hp=Hyp;
  Ha=Hya;

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
  

  //Calculate arraysizes
  varSize=(Hp->nx+4)*(myNy+4);
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

  //Allocate arrays
  recvMesh=(double*)malloc(Hp->nvar*Hp->nx*myNy*sizeof(double));
  lMesh =(double*)malloc(Hp->nvar*varSize*sizeof(double));
  q  =(double*)malloc(primSize*sizeof(double));
  qr =(double*)malloc(qSize*sizeof(double));
  ql =(double*)malloc(qSize*sizeof(double));
  flx=(double*)malloc(flxSize*sizeof(double));

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
  oTM=sumArray(lMesh,VARRHO,Hp->nx,myNy,2,2);
  oTE=sumArray(lMesh,VARPR ,Hp->nx,myNy,2,2);
#ifdef M_PREC_CMP
  frexp(oTM,&M_exp);
  frexp(oTE,&E_exp);

  M_prec=ldexp(DBL_EPSILON,M_exp-1);
  E_prec=ldexp(DBL_EPSILON,E_exp-1);

  if(rank==0)printf("TM:%g+-%g TE:%g+-%g\n",volCell*oTM,volCell*M_prec,volCell*oTE,volCell*E_prec);
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
    printf("INIT: TM: %g TE: %g\n",volCell*oTM,volCell*oTE);
  }
  
  initT=MPI_Wtime();

  while((n<Ha->nstepmax||Ha->nstepmax<0)&&(cTime<Ha->tend||Ha->tend<0)){
    //Calculate timestep
    dt=Ha->sigma*calcDT(lMesh);
    if(nxttout>0.0&&dt>(nxttout-cTime)){
      if(rank==0)printf("Adjusting timestep from %g to %g for iter %d\n",dt,nxttout-cTime,n);
      dt=(nxttout-cTime);
    }
    if(n%2==0){
      //X Dir
      runPass(lMesh,dt,n,0);
      //Y Dir
      runPass(lMesh,dt,n,1);
    }else{
      //Y Dir
      runPass(lMesh,dt,n,1);
      //X Dir
      runPass(lMesh,dt,n,0);
    }
    n+=1;
    cTime+=dt;
    if(n%Ha->nprtLine==0){
      TM=sumArray(lMesh,VARRHO,Hp->nx,myNy,2,2);
      TE=sumArray(lMesh,VARPR ,Hp->nx,myNy,2,2);
      if(rank==0){
	printf("Iter %05d time %f dt %g TM: %g TE: %g\n",n,cTime,dt,volCell*TM,volCell*TE);
	if(0&&(fabs(TM-oTM)>0.0||fabs(TE-oTE)>0.0)){
	  printf("ERR(%5s): Mass %g Ene %g\n","RAW",volCell*(TM-oTM),volCell*(TE-oTE));
	  printf("ERR(%5s): Mass %g Ene %g\n","%",100.0*(TM-oTM)/oTM,100.0*(TE-oTE)/oTE);
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
  if(rank==0)printf("time: %f, %d iters run\n",cTime,n);

  endT=MPI_Wtime();

  if(rank==0){
    printf("TFMT:%s,%s,%s,%s,%s,%s,%s,%s\n","cType","mType","init","nproc","nth","niters","ncells","wRunt");
    printf("TIME:%s,%s,%s,%d,%d,%d,%d,%g\n","\"MPI\"","\"CPU:?\"","\"Init\"",size,1,n,Hp->nx*Hp->ny,(endT-initT));
  }

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
  free(q  );
  free(qr );
  free(ql );
  free(flx);

  printf("NODE %d: Finalizing MPI\n",rank);
  MPI_Finalize();

  printf("N %d: Returning from engine\n",rank);
}
