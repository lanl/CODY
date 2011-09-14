#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>

#define GL_GLEXT_PROTOTYPES
#include <GL/gl.h>
#include <GL/glext.h>
#include <GL/glut.h>
#include <GL/glu.h>
#include <cuda_gl_interop.h>

#include "engine.h"
#include "hydro_utils.h"
#include "hydro_macros.h"
#include "device_funcs.h"
#include "hydro_dmp.h"

#define FN_LEN 50

#define CDT_REGS 16
#define STEP_REGS 16

//Computation kernels
__global__ void primX(double *q, double *u, double *bnd);
__global__ void primY(double *q, double *u, double *bnd);
__global__ void trace(double *qxm, double *qxp, double *q, int np, int nt, double dtdx);
__global__ void riemann(double *flx, double *qxm, double *qxp, int np, int nt);
__global__ void cmp_dStX(double *du, double *qgdnv, double dtdx);
__global__ void cmp_dStY(double *du, double *qgdnv, double dtdx);
__global__ void add_dSt(double *u, double *du);
__global__ void crDispBuff(uchar4 *buff, double *u, int dispVar, float datN, float datX);

//OpenGL callbacks
void initGL(int *argc, char **argv);
void display();
void idleFunc();
void selDispVar(int);
void mouse(int button, int state, int x, int y);
void motion(int x, int y);
void reshape(int w, int h);
void keyboard(unsigned char, int, int);
void endRun();
void crLegTex();
void initDatTex();

// Host problem state vars
hydro_args Ha;
hydro_prob Hp;
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
cudaError_t cuErrVar;
size_t meshVarSize, bndVarSize, wkVarSize;
int nxy;
int mxTh, thWp;
int nTh, rpBl;
int nThCDT, nThStep;
int nBlockM;
dim3 prKX, prKY;
dim3 trKX, trKY;
dim3 riKX, riKY;
dim3 mesh;

//display vars
void *font = GLUT_BITMAP_9_BY_15;
float textW=0.02, textH=0.04;
int dispVar=0;
int winW=1000, winH=1000;
float dispRatio=1.0;
float cenX=0.0, cenY=0.0;
float disW=2.0, disH=2.0;
int panning=0;
int clX, clY;
float datN=0.0, datX=2.5;

int runStep=0;
int showVarN=1, showCen=0, showLeg=0, showTime=0;

char *varNames[4]={"Density","Momentum X","Momentum Y","Energy"};
#define nCol 256

cudaGraphicsResource_t datBuff_cu[1];
GLuint legTex, meshTex;
GLuint datBuff;

void engine(int argc, char *argv[], hydro_args HaIn){
    // Host problem state vars
    double *hst_uold;

    //Cuda vars
    int dev, devCount;
    cudaDeviceProp prop;
    size_t shMpBl;

    size_t mem_avail,mem_used;

    initGL(&argc, argv);

    //Device settings
    cudaGLSetGLDevice(0);
    cudaDeviceSynchronize();
    HANDLE_CUDA_ERROR(cuErrVar);
    cudaGetDevice(&dev);
    cudaGetDeviceCount(&devCount);
    HANDLE_CUDA_ERROR(cuErrVar);
    printf("Using device %d/%d\n", dev, devCount);
    cudaGetDeviceProperties(&prop,dev);
    HANDLE_CUDA_ERROR(cuErrVar);

    //initialize mesh/device
    Ha=HaIn;
    hydro_init(&hst_uold,&Ha,&Hp);
    device_init(&Ha,&Hp);

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

    //Set initial state on device
    cudaMemcpy(dev_uold,hst_uold,meshVarSize,cudaMemcpyHostToDevice);
    HANDLE_CUDA_ERROR(cuErrVar);
    //cudaMemcpy(hst_uold,dev_uold,meshVarSize,cudaMemcpyDeviceToHost);
    //print_array("Data",hst_uold,Hp.nx,Hp.ny,Hp.nvar);
    free(hst_uold);

    crLegTex();
    initDatTex();

    //begin mainloop
    glutMainLoop();

    endRun();
}

//Calculation
void iterCalc(){
    dim3 prKX(BL(Hp.nx+4,nThStep),Hp.ny), prKY(BL(Hp.nx  ,nThStep),Hp.ny+4);
    dim3 trKX(BL(Hp.nx+2,nThStep),Hp.ny), trKY(BL(Hp.nx+2,nThStep),Hp.nx);
    dim3 riKX(BL(Hp.nx+1,nThStep),Hp.ny), riKY(BL(Hp.nx+1,nThStep),Hp.nx);
    dim3 mesh(BL(Hp.nx  ,nThStep),Hp.ny);
    //calculate timestep
    if((nstep%2)==0){
      dt=0.0;
      nDen=nBlockM;
      redBlocks=nDen;
      nTh=nThCDT;
      calc_denom<<<nBlockM,nTh,nTh*sizeof(double)>>>(dev_uold, dev_denA);
      HANDLE_CUDA_ERROR(cuErrVar);
      while(redBlocks>1){
        redBlocks=(nDen+2*nTh-1)/(2*nTh);
        redu_max<<<redBlocks,nTh,nTh*sizeof(double)>>>(dev_denA,dev_denB,nDen);
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
      if((Ha.tend<Hp.t+dt)&&Ha.tend>0.0) dt=Ha.tend-Hp.t;
      else if((Ha.tend<Hp.t+2.0*dt)&&Ha.tend>0.0) dt=0.5*(Ha.tend-Hp.t);
    }
    //Step
    nTh=nThStep;
    if(nstep%2==0){
      //x pass
      gen_bndX<<<BL_TH(Hp.ny,nTh)>>>(dev_uold,dev_bnd,Hp.bndL,Hp.bndR);
      primX   <<<prKX,nTh>>>(dev_q,dev_uold,dev_bnd);
      trace   <<<trKX,nTh,Hp.nvar*(nTh+2)*sizeof(double)>>>(dev_qxm, dev_qxp, dev_q, Hp.nx+2, Hp.ny, dt/Hp.dx);
      riemann <<<riKX,nTh>>>(dev_flx, dev_qxm, dev_qxp, Hp.nx+1, Hp.ny);
      cmp_dStX<<<mesh,nTh>>>(dev_uold, dev_flx, dt/Hp.dx);
      //y pass
      gen_bndY<<<BL_TH(Hp.nx,nTh)>>>(dev_uold,dev_bnd,Hp.bndL,Hp.bndR);
      primY   <<<prKY,nTh>>>(dev_q,dev_uold,dev_bnd);
      trace   <<<trKY,nTh,Hp.nvar*(nTh+2)*sizeof(double)>>>(dev_qxm, dev_qxp, dev_q, Hp.ny+2, Hp.nx, dt/Hp.dy);
      riemann <<<riKY,nTh>>>(dev_flx, dev_qxm, dev_qxp, Hp.ny+1, Hp.nx);
      cmp_dStY<<<mesh,nTh>>>(dev_uold, dev_flx, dt/Hp.dy);
    }else{
      //y pass
      gen_bndY<<<BL_TH(Hp.nx,nTh)>>>(dev_uold,dev_bnd,Hp.bndL,Hp.bndR);
      primY   <<<prKY,nTh>>>(dev_q,dev_uold,dev_bnd);
      trace   <<<trKY,nTh,Hp.nvar*(nTh+2)*sizeof(double)>>>(dev_qxm, dev_qxp, dev_q, Hp.ny+2, Hp.nx, dt/Hp.dy);
      riemann <<<riKY,nTh>>>(dev_flx, dev_qxm, dev_qxp, Hp.ny+1, Hp.nx);
      cmp_dStY<<<mesh,nTh>>>(dev_uold, dev_flx, dt/Hp.dy);
      //x pass
      gen_bndX<<<BL_TH(Hp.ny,nTh)>>>(dev_uold,dev_bnd,Hp.bndL,Hp.bndR);
      primX   <<<prKX,nTh>>>(dev_q,dev_uold,dev_bnd);
      trace   <<<trKX,nTh,Hp.nvar*(nTh+2)*sizeof(double)>>>(dev_qxm, dev_qxp, dev_q, Hp.nx+2, Hp.ny, dt/Hp.dx);
      riemann <<<riKX,nTh>>>(dev_flx, dev_qxm, dev_qxp, Hp.nx+1, Hp.ny);
      cmp_dStX<<<mesh,nTh>>>(dev_uold, dev_flx, dt/Hp.dx);
    }
    HANDLE_CUDA_ERROR(cuErrVar);
    //Finish main loop
    ++nstep;
    Hp.t+=dt;
    if((nstep%Ha.nprtLine)==0){
      fprintf(stdout,"Run Iter=%6d t=%12.5g dt=%12g\n", nstep, Hp.t, dt);
    }
}

void datToTex(float minVal, float maxVal){
    uchar4 *dptr=NULL;
    int nThD=thWp;
    int nBlockD=BL(Hp.nx*Hp.ny,nThD);
    size_t datSize=Hp.nx*Hp.ny*4;

    cudaGraphicsMapResources(1, datBuff_cu, 0);
    cudaGraphicsResourceGetMappedPointer((void **) &dptr, &datSize, datBuff_cu[0]);
    cudaThreadSynchronize();
    crDispBuff<<<nBlockD,nThD>>>(dptr, dev_uold, dispVar, minVal, maxVal);
    cudaThreadSynchronize();
    cudaGraphicsUnmapResources(1, datBuff_cu, 0);
    HANDLE_CUDA_ERROR(cuErrVar);

    glBindBuffer(GL_PIXEL_UNPACK_BUFFER, datBuff);
    glBindTexture(GL_TEXTURE_2D, meshTex);
    glPixelStorei(GL_UNPACK_ALIGNMENT, GL_UNSIGNED_BYTE);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, Hp.nx, Hp.ny, 
		    GL_RGBA, GL_UNSIGNED_BYTE, NULL);
    glBindTexture(GL_TEXTURE_2D, 0);
    glBindBuffer(GL_PIXEL_UNPACK_BUFFER, 0);
}

void fitDisp(){
  //Fit to colormap to current display;
  double *u;
  double mx, mn;
  int i;

  u=(double *)malloc(meshVarSize);
  cudaMemcpy(u,dev_uold,meshVarSize,cudaMemcpyDeviceToHost);
  mx=mn=u[(dispVar*Hp.nx*Hp.ny)];
  for(i=0;i<Hp.nx*Hp.ny;++i){
    double cVal=u[(dispVar*Hp.nx*Hp.ny)+i];
    if(cVal>mx)mx=cVal;
    if(cVal<mn)mn=cVal;
  }
  datX=mx;
  datN=mn;
}

//Display functions
void initGL(int *argc, char **argv){
    int i;

    //Run GL program
    glutInit(argc, argv);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowSize(winW, winH);
    glutCreateWindow("CUDA Godunov Hydrocode");

    //set callbacks
    glutDisplayFunc(display);
    glutKeyboardFunc(keyboard);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    glutReshapeFunc(reshape);
    glutIdleFunc(idleFunc);

    glutCreateMenu(selDispVar);
    for(i=0;i<4;++i){
      glutAddMenuEntry(varNames[i],i);
    }
    glutAttachMenu(GLUT_RIGHT_BUTTON);
}

void crLegTex(){
    GLubyte texDat[1][nCol][4];
    int i;
    float dispVal;

    for(i=0;i<nCol;++i){
      dispVal=(float)i/(float)nCol;
      texDat[0][i][0]=(GLubyte) (0xff*dispVal);
      texDat[0][i][1]=(GLubyte) (0x88-0xff*abs(0.5-dispVal));
      texDat[0][i][2]=(GLubyte) (0xff-0xff*dispVal);
      texDat[0][i][3]=(GLubyte) 0xff;
    }
    glPixelStorei(GL_UNPACK_ALIGNMENT, GL_FLOAT);
    glGenTextures(1, &legTex);
    glBindTexture(GL_TEXTURE_2D, legTex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, 
                   GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, 
                   GL_NEAREST);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 1, nCol, 0, GL_RGBA, GL_UNSIGNED_BYTE, texDat);
}

void initDatTex(){
    int size_tex_data=sizeof(GLubyte)*Hp.nx*Hp.ny*4;
    GLubyte *texDat;

    texDat=(GLubyte *)malloc(size_tex_data);

    for(int i=0;i<size_tex_data;++i){
      texDat[i]=0;
    }

    glPixelStorei(GL_UNPACK_ALIGNMENT, GL_UNSIGNED_BYTE);
    glGenTextures(1, &meshTex);
    glBindTexture(GL_TEXTURE_2D, meshTex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, 
                   GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, 
                   GL_NEAREST);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, Hp.nx, Hp.ny, 0, GL_RGBA, GL_UNSIGNED_BYTE, texDat);
    glGenBuffers(1, &datBuff);
    glBindBuffer(GL_PIXEL_UNPACK_BUFFER,datBuff);
    glBufferData(GL_PIXEL_UNPACK_BUFFER, size_tex_data, NULL, GL_DYNAMIC_COPY);
    glBindBuffer(GL_PIXEL_UNPACK_BUFFER, 0);
    cudaGraphicsGLRegisterBuffer(datBuff_cu, datBuff,cudaGraphicsRegisterFlagsNone);
    HANDLE_CUDA_ERROR(cuErrVar);
}

void idleFunc(){
    glutPostRedisplay();
}

void drawMesh(){
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, meshTex);
    glColor4f(1.0,1.0,1.0,1.0);
    glBegin(GL_QUADS);
      glTexCoord2f(0.0,0.0); glVertex3f(-0.5*Hp.nx*Hp.dx*dispRatio,-0.5*Hp.ny*Hp.dy*dispRatio,0.0);
      glTexCoord2f(1.0,0.0); glVertex3f( 0.5*Hp.nx*Hp.dx*dispRatio,-0.5*Hp.ny*Hp.dy*dispRatio,0.0);
      glTexCoord2f(1.0,1.0); glVertex3f( 0.5*Hp.nx*Hp.dx*dispRatio, 0.5*Hp.ny*Hp.dy*dispRatio,0.0);
      glTexCoord2f(0.0,1.0); glVertex3f(-0.5*Hp.nx*Hp.dx*dispRatio, 0.5*Hp.ny*Hp.dy*dispRatio,0.0);
    glEnd();
    glDisable(GL_TEXTURE_2D);
    glColor4f(0.0,0.0,0.0,0.0);
    glBegin(GL_LINE_LOOP);
      glVertex3f(-0.5*Hp.nx*Hp.dx*dispRatio,-0.5*Hp.ny*Hp.dy*dispRatio,0.0);
      glVertex3f( 0.5*Hp.nx*Hp.dx*dispRatio,-0.5*Hp.ny*Hp.dy*dispRatio,0.0);
      glVertex3f( 0.5*Hp.nx*Hp.dx*dispRatio, 0.5*Hp.ny*Hp.dy*dispRatio,0.0);
      glVertex3f(-0.5*Hp.nx*Hp.dx*dispRatio, 0.5*Hp.ny*Hp.dy*dispRatio,0.0);
    glEnd();
}

void dispText(float x, float y, char *text){
    int len, i;
    float w, h;

    len = (int) strlen(text);
    w=textW*(len+1);
    h=textH*1.2;
    glColor3f(1.0,1.0,1.0);
    glBegin(GL_QUADS);
      glVertex3f(x  , y  , 0.0);
      glVertex3f(x+w, y  , 0.0);
      glVertex3f(x+w, y+h, 0.0);
      glVertex3f(x  , y+h, 0.0);
    glEnd();
    glColor3f(0.0,0.0,0.0);
    glBegin(GL_LINE_LOOP);
      glVertex3f(x  , y  , 0.0);
      glVertex3f(x+w, y  , 0.0);
      glVertex3f(x+w, y+h, 0.0);
      glVertex3f(x  , y+h, 0.0);
    glEnd();
    glRasterPos3f(x+0.5*textW,y+0.4*textH,0.0);
    for(i=0;i<len;++i){
      glutBitmapCharacter(font,text[i]);
    }
}

void drawLegend(float datN, float datX){
    float nw=textW*11;
    float nh=1.2*textH;
    char bnd[11];

    sprintf(bnd, "%10f", datX);
    dispText(cenX+0.5*disW-nw, cenY+0.5*disH-nh, bnd);
    sprintf(bnd, "%10f", datN);
    dispText(cenX+0.5*disW-nw, cenY-0.5*disH   , bnd);
    glColor4f(0.0,0.0,0.0,1.0);
    glBegin(GL_LINES);
      glVertex3f(cenX+0.5*disW-nw, cenY-0.5*disH, 0.0);
      glVertex3f(cenX+0.5*disW-nw, cenY+0.5*disH, 0.0);
    glEnd();
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, legTex);
    glColor4f(1.0,1.0,1.0,1.0);
    glBegin(GL_QUADS);
      glTexCoord2f(0.0,0.0); glVertex3f(cenX+0.5*disW-nw, cenY-0.5*disH+nh, 0.0);
      glTexCoord2f(1.0,0.0); glVertex3f(cenX+0.5*disW   , cenY-0.5*disH+nh, 0.0);
      glTexCoord2f(1.0,1.0); glVertex3f(cenX+0.5*disW   , cenY+0.5*disH-nh, 0.0);
      glTexCoord2f(0.0,1.0); glVertex3f(cenX+0.5*disW-nw, cenY+0.5*disH-nh, 0.0);
    glEnd();
    glDisable(GL_TEXTURE_2D);
}

void display(){
    glClearColor(0.0,0.0,0.0,1.0);
    if(Hp.t<Ha.tend&&runStep){
      iterCalc();
    }
    datToTex(datN, datX);

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(cenX-0.5*disW,cenX+0.5*disW,cenY-0.5*disH,cenY+0.5*disH);
    glMatrixMode(GL_MODELVIEW);
    drawMesh();
    //Run finished Signal
    if(Hp.t>=Ha.tend){
      dispText(cenX-1.75*textW,cenY-0.5*disH,"END");
    }
    //Variable name
    if(showVarN){
      dispText(cenX-0.5*disW,cenY+0.5*disH-1.2*textH,varNames[dispVar]);
    }
    //Center of display
    if(showCen){
      char center[40];
      sprintf(center,"(%f,%f) %fX%f",0.5*Hp.nx*Hp.dx+(cenX/dispRatio),0.5*Hp.ny*Hp.dy+(cenY/dispRatio),disW/dispRatio,disH/dispRatio);
      dispText(cenX-0.5*disW,cenY-0.5*disH,center);
    }
    //Current time
    if(showTime){
      char tCurr[25];
      float ofst;

      sprintf(tCurr,"t=%lf",Hp.t);
      ofst=((strlen(tCurr)+1.0)/2.0)*textW;
      dispText(cenX-ofst,cenY+0.5*disH-1.2*textH,tCurr);
    }
    //Color legend
    if(showLeg){
      drawLegend(datN,datX);
    }
    glutSwapBuffers();
}

//Callback functions
void keyboard(unsigned char key, int x, int y){
  switch(key){
    case 27:
      endRun();
      exit(0);
      break;
    case '+':
      dispRatio*=2;
      if(dispRatio>pow(2.0,100)) dispRatio=pow(2.0,100);
      break;
    case '-':
      dispRatio/=2;
      if(dispRatio<pow(2.0,-100)) dispRatio=pow(2.0,-100);
      break;
    case 'c':
      cenX=0.0; cenY=0.0;
      break;
    case 'f':
      fitDisp();
      break;
    case 'd':
      showVarN ^=1;
      break;
    case 'v':
      showCen ^=1;
      break;
    case 'l':
      showLeg ^=1;
      break;
    case 't':
      showTime ^=1;
      break;
    case ' ':
      runStep ^=1;
      break;
    default: break;
  }
}

void reshape(int w, int h){
  glViewport(0,0,w,h);
  disW=(2.0*w)/1000.0;
  disH=(2.0*h)/1000.0;
  winW=w;
  winH=h;
}

void selDispVar(int var){
  if(0<=var&&var<4)dispVar=var;
}

void mouse(int button, int state, int x, int y){
  if(state == GLUT_DOWN){
    panning=1;
    clX=x;
    clY=y;
  }
  if(state == GLUT_UP){
    panning=0;
  }
}

void motion(int x, int y){
  if(panning){
    cenX-=(disW*(x-clX))/winW;
    cenY+=(disH*(y-clY))/winH;
    clX=x;
    clY=y;
  }
}

void endRun(){
    cudaGraphicsUnregisterResource(*datBuff_cu);
    //delete textures
    glDeleteTextures(1, &legTex);
    glDeleteTextures(1, &meshTex);
    glDeleteBuffers(1, &datBuff);
    //Write final iter line
    fprintf(stdout,"End Iter=%6d t=%12.5g\n", nstep, Hp.t);

    //Print run end info
    printf("Run completed successfully\n");
    //Free cuda vars
    cudaFree(dev_uold);
    cudaFree(dev_q);
    cudaFree(dev_qxm);
    cudaFree(dev_qxp);
    cudaFree(dev_flx);
    cudaFree(dev_bnd);
    cudaFree(dev_denA);
    cudaFree(dev_denB);
}
