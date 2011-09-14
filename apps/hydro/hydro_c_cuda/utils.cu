#include <stdio.h>
#include <stdlib.h>

#include "hydro_structs.h"
#include "hydro_dmp.h"
#include "hydro_macros.h"

void get_data(int argc, char *argv[]);
void gen_init(int argc, char *argv[]);
void gen_vis(int argc, char *argv[]);
void get_usage(int argc, char *argv[]);
void test_dmp(int argc, char *argv[]);
void test_sym(int argc, char *argv[]);

int main(int argc, char *argv[]){
    if(argc<2){
      fprintf(stderr,"Insufficient arguments to %s\n",argv[0]);
      return 1;
    }
    if(!strcmp("info",argv[1])){
      get_data(argc,argv);
      return 0;
    }
    if(!strcmp("init",argv[1])){
      gen_init(argc,argv);
      return 0;
    }
    if(!strcmp("sym",argv[1])){
      test_sym(argc,argv);
      return 0;
    }
    if(!strcmp("test",argv[1])){
      test_dmp(argc,argv);
      return 0;
    }
    if(!strcmp("usage",argv[1])){
      get_usage(argc,argv);
      return 0;
    }
    if(!strcmp("vis",argv[1])){
      gen_vis(argc,argv);
      return 0;
    }
    printf("Unknown command %s\n",argv[1]);
}

void get_data(int argc, char *argv[]){
    hydro_prob Hp;
    double *u;
    int i,j,n;
    double totM, totE;

    for(n=2;n<argc;++n){
      if(rd_dmp(argv[n],&Hp,&u))continue;
      totM=0.0;
      totE=0.0;
      for(j=0;j<Hp.ny;++j){
        for(i=0;i<Hp.nx;++i){
          totM+=Hp.dx*Hp.dy*u[(ID*Hp.ny+j)*Hp.nx+i];
          totE+=Hp.dx*Hp.dy*u[(IP*Hp.ny+j)*Hp.nx+i];
        }
      }
      printf("----------------------\n");
      printf("%s:\n",argv[n]);
      printf("%d variable %dx%d mesh ",Hp.nvar,Hp.nx,Hp.ny);
      printf("covering a %4fx%4f area\n",Hp.dx*Hp.nx,Hp.dy*Hp.ny);
      printf("Current time: %lf\n",Hp.t);
      printf("Total Mass %lf\n",totM);
      printf("Total Energy %lf\n",totE);
      printf("----------------------\n\n");
      free(u);
    }
}

void test_dmp(int argc, char *argv[]){
    hydro_prob Hp;
    double *u;
    int n;

    for(n=2;n<argc;++n){
      if(rd_dmp(argv[n],&Hp,&u))continue;
      printf("%s: ",argv[n]);
      test_mesh(&Hp,u);
    }
}

void test_sym(int argc, char *argv[]){
    hydro_prob Hp;
    double *u;
    int i,j,k,n;
    int clean;
    double v1, v2, denom, err;
    double percen=0.001, eps=1.0e-10;

    for(n=2;n<argc;++n){
      if(rd_dmp(argv[n],&Hp,&u))continue;
      printf("%s: ",argv[n]);
      if(Hp.nx!=Hp.ny){
        printf("Mesh not square\n");
        continue;
      }
      v1=Hp.dx;
      v2=Hp.dy;
      denom=(fabs(v1)+fabs(v2));
      if(denom==0.0)err=0.0;
      else err=2*(v1-v2)/(fabs(v1)+fabs(v2));
      if(fabs(err)>percen){
        printf("Cells not square\n");
        continue;
      }
      printf("Running check\n");
      clean=1;
      for(k=0;k<Hp.nvar;++k){
        for(j=0;j<Hp.ny;++j){
          for(i=0;i<=j;++i){
            printf("Testing state var %d [%d][%d]\r",k,j,i);
            v1=u[(k*Hp.ny+j)*Hp.nx+i];
            if(k==IU)
              v2=u[(IV*Hp.ny+i)*Hp.nx+j];
            else if(k==IV)
              v2=u[(IU*Hp.ny+i)*Hp.nx+j];
            else
              v2=u[(k*Hp.ny+i)*Hp.nx+j];
            denom=(fabs(v1)+fabs(v2));
            if(denom==0.0)err=0.0;
            else err=2*(v1-v2)/(fabs(v1)+fabs(v2));
            if(fabs(err)>percen&&fabs(v1-v2)>eps){
              printf("Symmetry error in state var %d[%d,%d] and [%d,%d] of %f%%:\n", k, i, j, j, i,err*100.0);
              printf("\tv1:%lg v2: %lg diff: %lg\n",v1, v2, v1-v2);
              clean=0;
            }
          }
        }
      }
      if(clean)printf("                             \rSYMMETRIC\n");
    }
}

#define INIT_SOD 1
#define INIT_IBW 2
#define INIT_SIB 3
#define INIT_EXPL 4
#define INIT_CORN 5
#define INIT_CSHK 6

void gen_init(int argc, char *argv[]){
    hydro_prob Hp;
    double *u;
    int init;
    int i, j, n;
    char *outFile;
    char defOut[10];

    //Init args to reasonable values
    Hp.nx=10;
    Hp.ny=10;
    init=INIT_SOD;
    strcpy(defOut,"init.dmp");
    outFile=defOut;
    //parse command line args
    i=2;
    if(argc<2)fprintf(stderr,"Insufficient arguments to %s\n",argv[0]);
    while(i<argc){
      if(!strcmp("-o",argv[i])){//output file
        ++i;
        if(i<argc){
          outFile=argv[i];
          ++i;
        }else{
          fprintf(stderr,"Error: no filename after -o\n");
          exit(1);
        }
        continue;
      }
      if(!strcmp("--dim",argv[i])){
        ++i;
        if(i<argc){
          n=sscanf(argv[i],"%d,%d",&Hp.nx,&Hp.ny);
          if(n!=2){
            fprintf(stderr,"Error: did not find values for nx and ny after --dim\n");
          }
        }else{
          fprintf(stderr,"Error: no values for nx and ny after --dim\n");
          exit(1);
        }
      }
      if(!strcmp("--init",argv[i])){//output file
        ++i;
        if(i<argc){
          if(!strcmp(argv[i],"sod")){
            init=INIT_SOD;
            ++i;
            continue;
          }
          if(!strcmp(argv[i],"ibw")){
            init=INIT_IBW;
            ++i;
            continue;
          }
          if(!strcmp(argv[i],"sib")){
            init=INIT_SIB;
            ++i;
            continue;
          }
          if(!strcmp(argv[i],"expl")){
            init=INIT_EXPL;
            ++i;
            continue;
          }
          if(!strcmp(argv[i],"simp")){
            init=INIT_CORN;
            ++i;
            continue;
          }
          if(!strcmp(argv[i],"cshk")){
            init=INIT_CSHK;
            ++i;
            continue;
          }
          fprintf(stderr,"Unknown init %s\n",argv[i]);
          ++i;
          exit(1);
        }else{
          fprintf(stderr,"Error: no initialization specified\n");
          exit(1);
        }
        continue;
      }
      ++i;
    }
    //Print init info
    printf("Init file characteristics:\n");
    printf("-----------------------------\n");
    printf("%-10s=%10d\n","nx",Hp.nx);
    printf("%-10s=%10d\n","ny",Hp.ny);
    printf("%-10s=%5s%s\n","filename","",outFile);
    if(init==INIT_SOD)
      printf("%-10s=%5s%s\n","init","","Sod problem");
    else if(init==INIT_IBW)
      printf("%-10s=%5s%s\n","init","","Interacting blast waves");
    else if(init==INIT_SIB)
      printf("%-10s=%5s%s\n","init","","Shock in a box");
    else if(init==INIT_EXPL)
      printf("%-10s=%5s%s\n","init","","Point explosion");
    else if(init==INIT_CORN)
      printf("%-10s=%5s%s\n","init","","Corner explosion");
    else if(init==INIT_CSHK)
      printf("%-10s=%5s%s\n","init","","Corner explosion");
    printf("-----------------------------\n\n");
   //generate problem init
    if(Hp.nx<0||Hp.ny<0){
      fprintf(stderr,"Bad value for nx or ny\n");
      exit(1);
    }
    //Set header
    Hp.t=0.0;
    if(init==INIT_SOD){
      Hp.nvar=4;
      Hp.dx=1.0/Hp.nx;
      Hp.dy=0.25/Hp.ny;
      Hp.gamma=1.4;
      Hp.bndL=BND_REFL;
      Hp.bndR=BND_REFL;
      Hp.bndU=BND_REFL;
      Hp.bndD=BND_REFL;
    }
    if(init==INIT_IBW){
      Hp.nvar=4;
      Hp.dx=1.0/Hp.nx;
      Hp.dy=0.25/Hp.ny;
      Hp.gamma=1.4;
      Hp.bndL=BND_REFL;
      Hp.bndR=BND_REFL;
      Hp.bndU=BND_REFL;
      Hp.bndD=BND_REFL;
    }
    if(init==INIT_SIB){
      Hp.nvar=4;
      Hp.dx=1.0/Hp.nx;
      Hp.dy=1.0/Hp.ny;
      Hp.gamma=1.4;
      Hp.bndL=BND_REFL;
      Hp.bndR=BND_REFL;
      Hp.bndU=BND_REFL;
      Hp.bndD=BND_REFL;
    }
    if(init==INIT_EXPL){
      Hp.nvar=4;
      Hp.dx=0.8/Hp.nx;
      Hp.dy=1.0/Hp.ny;
      Hp.gamma=1.4;
      Hp.bndL=BND_REFL;
      Hp.bndR=BND_REFL;
      Hp.bndU=BND_REFL;
      Hp.bndD=BND_REFL;
    }
    if(init==INIT_CORN){
      Hp.nvar=4;
      Hp.dx=5.0/Hp.nx;
      Hp.dy=5.0/Hp.ny;
      Hp.gamma=1.4;
      Hp.bndL=BND_REFL;
      Hp.bndR=BND_REFL;
      Hp.bndU=BND_REFL;
      Hp.bndD=BND_REFL;
    }
    if(init==INIT_CSHK){
      Hp.nvar=4;
      Hp.dx=5.0/Hp.nx;
      Hp.dy=5.0/Hp.ny;
      Hp.gamma=1.4;
      Hp.bndL=BND_REFL;
      Hp.bndR=BND_REFL;
      Hp.bndU=BND_REFL;
      Hp.bndD=BND_REFL;
    }
    //Fill mesh
    u=(double *)calloc(Hp.nvar*Hp.ny*Hp.nx,sizeof(double));
    if(init==INIT_SOD){
      for(j=0;j<Hp.ny;++j){
        for(i=0;i<0.5*Hp.nx;++i){
          u[(ID*Hp.ny+j)*Hp.nx+i]=1.0;
          u[(IU*Hp.ny+j)*Hp.nx+i]=0.0;
          u[(IV*Hp.ny+j)*Hp.nx+i]=0.0;
          u[(IP*Hp.ny+j)*Hp.nx+i]=2.5;
        }
        for(i=0.5*Hp.nx;i<Hp.nx;++i){
          u[(ID*Hp.ny+j)*Hp.nx+i]=0.125;
          u[(IU*Hp.ny+j)*Hp.nx+i]=0.0;
          u[(IV*Hp.ny+j)*Hp.nx+i]=0.0;
          u[(IP*Hp.ny+j)*Hp.nx+i]=0.25;
        }
      }
    }
    if(init==INIT_IBW){
      for(j=0;j<Hp.ny;++j){
        for(i=0;i<0.1*Hp.nx;++i){
          u[(ID*Hp.ny+j)*Hp.nx+i]=1.0;
          u[(IU*Hp.ny+j)*Hp.nx+i]=0.0;
          u[(IV*Hp.ny+j)*Hp.nx+i]=0.0;
          u[(IP*Hp.ny+j)*Hp.nx+i]=1.0e3;
        }
        for(i=0.1*Hp.nx;i<0.9*Hp.nx;++i){
          u[(ID*Hp.ny+j)*Hp.nx+i]=1.0;
          u[(IU*Hp.ny+j)*Hp.nx+i]=0.0;
          u[(IV*Hp.ny+j)*Hp.nx+i]=0.0;
          u[(IP*Hp.ny+j)*Hp.nx+i]=1.0e-2;
        }
        for(i=0.9*Hp.nx;i<Hp.nx;++i){
          u[(ID*Hp.ny+j)*Hp.nx+i]=1.0;
          u[(IU*Hp.ny+j)*Hp.nx+i]=0.0;
          u[(IV*Hp.ny+j)*Hp.nx+i]=0.0;
          u[(IP*Hp.ny+j)*Hp.nx+i]=1.0e2;
        }
      }
    }
    if(init==INIT_SIB){
      for(j=0;j<Hp.ny;++j){
        for(i=0;i<Hp.nx;++i){
          u[(ID*Hp.ny+j)*Hp.nx+i]=1.0;
          u[(IU*Hp.ny+j)*Hp.nx+i]=0.0;
          u[(IV*Hp.ny+j)*Hp.nx+i]=0.0;
          u[(IP*Hp.ny+j)*Hp.nx+i]=1.0e-5;
        }
      }
      for(j=(Hp.ny/2)-1+(Hp.ny%2);j<=(Hp.ny/2);++j){
        for(i=(Hp.nx/2)-1+(Hp.nx%2);i<=(Hp.nx/2);++i){
          u[(IP*Hp.ny+j)*Hp.nx+i]=1.0/((Hp.nx%2+1)*Hp.dx*(Hp.ny%2+1)*Hp.dy);
        }
      }
    }
    if(init==INIT_EXPL){
      for(j=0;j<Hp.ny;++j){
        for(i=0;i<Hp.nx;++i){
          u[(ID*Hp.ny+j)*Hp.nx+i]=1.0;
          u[(IU*Hp.ny+j)*Hp.nx+i]=0.0;
          u[(IV*Hp.ny+j)*Hp.nx+i]=0.0;
          u[(IP*Hp.ny+j)*Hp.nx+i]=1.0e-5;
        }
      }
      for(j=(Hp.ny/2)-1+(Hp.ny%2);j<=(Hp.ny/2);++j){
        for(i=(Hp.nx/2)-1+(Hp.nx%2);i<=(Hp.nx/2);++i){
          u[(IP*Hp.ny+j)*Hp.nx+i]=1.0/((Hp.nx%2+1)*Hp.dx*(Hp.ny%2+1)*Hp.dy);
        }
      }
    }
    if(init==INIT_CORN){
      for(j=0;j<Hp.ny;++j){
        for(i=0;i<Hp.nx;++i){
          u[(ID*Hp.ny+j)*Hp.nx+i]=1.0;
          u[(IU*Hp.ny+j)*Hp.nx+i]=0.0;
          u[(IV*Hp.ny+j)*Hp.nx+i]=0.0;
          u[(IP*Hp.ny+j)*Hp.nx+i]=1.0e-5;
        }
      }
      u[(IP*Hp.ny)*Hp.nx]=1.0/(Hp.dx*Hp.dy);
    }
    if(init==INIT_CSHK){
      for(j=0;j<0.5*Hp.ny;++j){
        for(i=0;i<0.5*Hp.nx;++i){
          u[(ID*Hp.ny+j)*Hp.nx+i]=1.0;
          u[(IU*Hp.ny+j)*Hp.nx+i]=0.0;
          u[(IV*Hp.ny+j)*Hp.nx+i]=0.0;
          u[(IP*Hp.ny+j)*Hp.nx+i]=2.5;
        }
        for(i=0.5*Hp.nx;i<Hp.nx;++i){
          u[(ID*Hp.ny+j)*Hp.nx+i]=1.0;
          u[(IU*Hp.ny+j)*Hp.nx+i]=0.0;
          u[(IV*Hp.ny+j)*Hp.nx+i]=0.0;
          u[(IP*Hp.ny+j)*Hp.nx+i]=0.25;
        }
      }
      for(j=0.5*Hp.ny;j<Hp.ny;++j){
        for(i=0;i<Hp.nx;++i){
          u[(ID*Hp.ny+j)*Hp.nx+i]=1.0;
          u[(IU*Hp.ny+j)*Hp.nx+i]=0.0;
          u[(IV*Hp.ny+j)*Hp.nx+i]=0.0;
          u[(IP*Hp.ny+j)*Hp.nx+i]=0.25;
        }
      }
    }
    //Write init file
    wrt_dmp(outFile,&Hp,u);
    free(u);
}

void gen_vis(int argc, char *argv[]){
    hydro_prob Hp;
    double *u;
    int i,j,nv,argInd;
    FILE *vis;
    char outName[53], *ext;
    char name[30];

    for(argInd=2;argInd<argc;++argInd){
      if(rd_dmp(argv[argInd],&Hp,&u))continue;
      strncpy(outName,argv[argInd],50);
      outName[49]='\0';
      if(!(ext=strrchr(outName,'.'))) ext=strchr(outName,'\0');
      sprintf(ext,".vts");
      vis=fopen(outName,"w");
      if(vis==NULL){
        fprintf(stderr,"Could not open file %s\n",outName);
        free(u);
        continue;
      }
      fprintf(vis, "<?xml version=\"1.0\"?>\n");
      fprintf(vis, "<VTKFile type=\"StructuredGrid\">\n");
      fprintf(vis, "<StructuredGrid WholeExtent=\"%ld %ld %ld %ld %ld %ld\">\n", (long)0, (long)Hp.nx, (long)0, (long) Hp.ny, (long)0, (long)0);
      fprintf(vis, "<Piece Extent=\"%ld %ld %ld %ld %ld %ld\">\n", (long)0, Hp.nx, (long)0, Hp.ny, (long)0, (long)0);
      fprintf(vis, "<PointData></PointData>\n");
      name[0] = 0;

      for (nv = 0; nv <= IP; nv++)
      {
        if (nv == ID)
          sprintf(name, "%s varID", name);
        if (nv == IU)
          sprintf(name, "%s varIU", name);
        if (nv == IV)
          sprintf(name, "%s varIV", name);
        if (nv == IP)
          sprintf(name, "%s varIP", name);
      }

      // declaration of the variable list
      fprintf(vis, "<CellData Scalars=\"%s\">\n", name);

      name[0] = 0;

      for (nv = 0; nv <= IP; nv++){
        if (nv == ID)
          sprintf(name, "varID");
        if (nv == IU)
          sprintf(name, "varIU");
        if (nv == IV)
          sprintf(name, "varIV");
        if (nv == IP)
          sprintf(name, "varIP");

        //Definition of the cell values
        fprintf(vis, "<DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n", name);

        // the image is the interior of the computed domain
        for (j = 0; j < Hp.ny; j++){
          for (i = 0; i < Hp.nx; i++){
                fprintf(vis, "%lf ", u[(nv*Hp.ny+j)*Hp.nx+i]);
          }
          fprintf(vis, "\n");
        }
        fprintf(vis, "</DataArray>\n");
      }
      fprintf(vis, "</CellData>\n");
      fprintf(vis, "<Points>\n");
      fprintf(vis,"<DataArray type=\"Float32\" format=\"ascii\" NumberOfComponents=\"3\">\n");

      for (j = 0; j < Hp.ny + 1; j++)
      {
        for (i = 0; i < Hp.nx + 1; i++)
        {
            fprintf(vis, "%f %f %f\n", i * Hp.dx, j * Hp.dx, 0.0);
        }
      }

      fprintf(vis, "</DataArray>\n");
      fprintf(vis, "</Points>\n");
      fprintf(vis, "</Piece>\n");
      fprintf(vis, "</StructuredGrid>\n");
      fprintf(vis, "</VTKFile>\n");
      fclose(vis);
      free(u);
    }
}

void get_usage(int argc, char *argv[]){
    int i;

    if(argc<2){
      printf("Tools available are:\n");
      printf("\t%7s-%s\n","info","Provides basic info from dump file");
      printf("\t%7s-%s\n","init","Creates intialization dump file");
      printf("\t%7s-%s\n","sym","Runs symmetry check");
      printf("\t%7s-%s\n","test","Checks for basic sanity of dump file values");
      printf("\t%7s-%s\n","usage","Print information on arguments");
      printf("\t%7s-%s\n","vis","Creates vis file from dump file");
    }

    for(i=2;i<argc;++i){
      printf("Usage for %s %s:\n",argv[0],argv[1]);
      if(!strcmp("info",argv[i])){
        printf("\t%s %s [<dump_file_a>] [<dump_file_b> ...\n", argv[0], argv[1]);
        continue;
      }
      if(!strcmp("init",argv[i])){
        printf("\t%s %s [-o <initFile>] [--dim <nx>,<ny>] [--init <init>]\n", argv[0], argv[1]);
        continue;
      }
      if(!strcmp("sym",argv[i])){
        printf("\t%s %s [<dump_file_a>] [<dump_file_b> ...\n", argv[0], argv[1]);
        continue;
      }
      if(!strcmp("test",argv[i])){
        printf("\t%s %s [<dump_file_a>] [<dump_file_b> ...\n", argv[0], argv[1]);
        continue;
      }
      if(!strcmp("usage",argv[i])){
        printf("\t%s %s [<util_a>] [<util_b> ...\n", argv[0], argv[1]);
        continue;
      }
      if(!strcmp("vis",argv[i])){
        printf("\t%s %s [<dump_file_a>] [<dump_file_b> ...\n", argv[0], argv[1]);
        continue;
      }
      printf("\tUnknown util\n");
    }
}

