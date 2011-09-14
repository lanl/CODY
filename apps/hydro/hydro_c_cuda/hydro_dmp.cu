#include <stdio.h>
#include <stdlib.h>

#include "hydro_dmp.h"
#include "hydro_macros.h"

int wrt_dmp(char *filename, hydro_prob *Hp, double *u){
    FILE *outFile;
    size_t count;
    size_t nCellMesh=Hp->nvar*Hp->ny*Hp->nx;

    outFile=fopen(filename,"wb");//open file
    if(outFile==NULL){
      fprintf(stderr,"Could not open dump file %s for writing\n",filename);
      return 1;
    }
    count=fwrite(Hp,sizeof(hydro_prob),1,outFile);//Write hydro_prob
    if(count!=1){
      fprintf(stderr,"Error while writing dump file header %s\n",filename);
      fclose(outFile);
      return 1;
    }
    count=fwrite(u,sizeof(double),nCellMesh,outFile);//Write mesh
    if(count!=nCellMesh){
      fprintf(stderr,"Error while writing mesh to file %s\n",filename);
      fclose(outFile);
      return 1;
    }
    fclose(outFile);
    return 0;
}

  int rd_dmp(char *filename, hydro_prob *Hp, double **u){
    FILE *inFile;
    size_t count;
    size_t nCellMesh;

    inFile=fopen(filename,"rb");//Open file
    if(inFile==NULL){
      fprintf(stderr,"Could not open dump file %s for reading\n",filename);
      return 1;
    }
    count=fread(Hp,sizeof(hydro_prob),1,inFile);
    if(count!=1){
      fprintf(stderr,"Error reading header data from dump file %s\n", filename);
      fclose(inFile);
      return 1;
    }
    //check for sane values
    if(Hp->nvar<4||Hp->ny<1||Hp->nx<1){
      fprintf(stderr,"Bad values for dimensions in file %s\n",filename);
      fclose(inFile);
      return 1;
    }
    if(Hp->dx<=0.0||Hp->dy<=0.0||Hp->gamma<=1.0){
      fprintf(stderr,"Bad physics vals in file %s\n",filename);
      fclose(inFile);
      return 1;
    }
    //read mesh
    nCellMesh=Hp->nvar*Hp->ny*Hp->nx;
    (*u)=(double *)calloc(nCellMesh,sizeof(double));
    count=fread((*u),sizeof(double),nCellMesh,inFile);
    if(count!=nCellMesh){
      fprintf(stderr,"Error reading mesh data from dump file %s\n",filename);
      fclose(inFile);
      return 1;
    }
    fclose(inFile);
    return 0;
}


void test_mesh(hydro_prob *Hp, double *u){
    int i, j, k;
    int clean=1;

    if(Hp->nx<=0||Hp->ny<=0){printf("\nBad mesh size");clean=0;}
    if(Hp->nvar<4){printf("\nNot enough state vars");clean=0;}
    if(Hp->dx<=0||Hp->dy<=0){printf("\nBad cell size");clean=0;}
    if(Hp->gamma<=1.0){printf("\nBad gamma val");clean=0;}
    if(!clean)return;
    for(j=0;j<Hp->ny;++j){
      for(i=0;i<Hp->nx;++i){
        if(u[(ID*Hp->ny+j)*Hp->nx+i]<=0){printf("\nBad density val[%d,%d]", i, j);clean=0;}
      }
    }
    for(j=0;j<Hp->ny;++j){
      for(i=0;i<Hp->nx;++i){
        if(u[(IP*Hp->ny+j)*Hp->nx+i]<=0){printf("\nBad energy val[%d,%d]", i, j);clean=0;}
      }
    }
    for(k=0;k<Hp->nvar;++k){
      for(j=0;j<Hp->ny;++j){
        for(i=0;i<Hp->nx;++i){
          if(isnan(u[(k*Hp->ny+j)*Hp->nx+i])){printf("\nNAN @ state var %d[%d,%d]",k,i,j);clean=0;}
        }
      }
    }
    if(clean)
      printf("OK\n");
    else
      printf("\n");
}
