#include <stdio.h>
#include <stdlib.h>

#include "hydro_args.h"
#include "hydro_macros.h"

void usage(){
    fprintf(stderr,"Possible options:\n");
    fprintf(stderr,"\t--help\t- Prints this message\n");
    fprintf(stderr,"\t-i [file]\t- input file\n"); 
}

void set_default(hydro_args *H){
    // I/O arguments
    H->tend=-1.0;
    H->nstepmax=100000;
    H->noutput=1000000;
    H->dtoutput=-1.0;
    H->nprtLine=1;
    strcpy(H->dmpPre,"out");
    H->initFile[0]='\0';

   //Dimensions
   H->nx=5;
   H->ny=5;

    // Physics
    H->sigma=0.5;
    H->smallc=1.0e-10;
    H->smallr=1.0e-10;

    // Numerical scheme
    H->niter_riemann=10;
    H->iorder=2;
    H->slope_type=1.0;
    H->scheme=HSCHEME_MUSCL;
}

void keyval(char *buf, char **pkey,char **pval){
  char *ptr;
  *pkey=buf;
  *pval=buf;

  //get rid of newline
  *pval=strchr(buf,'\n');
  if(*pval) **pval='\0';

  //Skip leading whitespace
  while((**pkey==' ')||(**pkey=='\t')) (*pkey)++;
  //Set val
  *pval=strchr(*pkey,'=');
  if(*pval){
    **pval='\0';
    (*pval)++;
  }
  //Strip key
  while((ptr=strchr(*pkey,' '))!=NULL)*ptr='\0';
  while((ptr=strchr(*pkey,'\t'))!=NULL)*ptr='\0';
}

void process_argFile(char *filename, hydro_args *H){
    FILE *fp=NULL;
    char buffer[1024];
    char *pkey, *pval;
    fp=fopen(filename,"r");
    if(fp==NULL){
      printf("Could not open input file %s\n",filename);
      exit(1);
    }
    while(fgets(buffer,1024,fp)==buffer){
      keyval(buffer,&pkey,&pval);

      //int params
      if(!strcmp(pkey,"nstepmax")){
        sscanf(pval,"%d",&H->nstepmax);
        continue;
      }
      if(!strcmp(pkey,"nx")){
        sscanf(pval,"%d",&H->nx);
        continue;
      }
      if(!strcmp(pkey,"ny")){
        sscanf(pval,"%d",&H->ny);
        continue;
      }
      if(!strcmp(pkey,"noutput")){
        sscanf(pval,"%d",&H->noutput);
        continue;
      }
      if(!strcmp(pkey,"nprtLine")){
        sscanf(pval,"%d",&H->nprtLine);
        continue;
      }
      if(!strcmp(pkey,"niter_riemann")){
        sscanf(pval,"%d",&H->niter_riemann);
        continue;
      }
      if(!strcmp(pkey,"iorder")){
        sscanf(pval,"%d",&H->iorder);
        continue;
      }
      //double params
      if(!strcmp(pkey,"tend")){
        sscanf(pval,"%lf",&H->tend);
        continue;
      }
      if(!strcmp(pkey,"dtoutput")){
        sscanf(pval,"%lf",&H->dtoutput);
        continue;
      }
      if(!strcmp(pkey,"courant_factor")){
        sscanf(pval,"%lf",&H->sigma);
        continue;
      }
      if(!strcmp(pkey,"smallr")){
        sscanf(pval,"%lf",&H->smallr);
        continue;
      }
      if(!strcmp(pkey,"smallc")){
        sscanf(pval,"%lf",&H->smallc);
        continue;
      }
      if(!strcmp(pkey,"slope_type")){
        sscanf(pval,"%lf",&H->slope_type);
        continue;
      }
      //string arg
      if(!strcmp(pkey,"scheme")){
        if(!strcmp(pval,"muscl"))H->scheme=HSCHEME_MUSCL;
        else if(!strcmp(pval,"plmde"))H->scheme=HSCHEME_PLMDE;
        else if(!strcmp(pval,"collela"))H->scheme=HSCHEME_COLLELA;
	else {
          printf("Unknown scheme %s, options are [muscl,plmde,collela]\n",pval);
          fclose(fp);
          exit(1);
        }
      }
      if(!strcmp(pkey,"dump_prefix")){
        strncpy(H->dmpPre,pval,PREFIX_LEN);
        H->dmpPre[PREFIX_LEN-1]='\0';
      }
      if(!strcmp(pkey,"init_file")){
        strncpy(H->initFile,pval,INIT_FN_LEN);
        H->initFile[INIT_FN_LEN-1]='\0';
      }
    }
    fclose(fp);
}

void print_args(hydro_args *H){
    printf( "-----------------------------------\n"         );
    printf( "Input variables.\n"                            );
    printf( "-----------------------------------\n"         );
    printf( "nstepmax       = %12lu\n",   H->nstepmax       );
    printf( "tend           = %12.5f\n",  H->tend           );
    if(H->dtoutput>0)
      printf( "dtoutput       = %12.5f\n",H->dtoutput       );
    else
      printf( "noutput        = %12lu\n", H->noutput        );
    printf( "nprtLine       = %12lu\n",   H->nprtLine       );
    printf( "dump_prefix    = %12s\n",    H->dmpPre         );
    printf( "init_file      = %12s\n",    H->initFile       );
    printf( "-----------------------------------\n"         );
    printf( "nx             = %12lu\n",   H->nx             );
    printf( "ny             = %12lu\n",   H->ny             );
    printf( "-----------------------------------\n"         );
    printf( "courant_factor = %12.5f\n",  H->sigma          );
    printf( "smallr         = %12.5g\n",  H->smallr         );
    printf( "smallc         = %12.5g\n",  H->smallc         );
    printf( "niter_riemann  = %12lu\n",   H->niter_riemann  );
    printf( "iorder         = %12lu\n",   H->iorder         );
    if ( H->scheme == HSCHEME_MUSCL   )
        printf( "scheme         =        muscl\n" );
    if ( H->scheme == HSCHEME_PLMDE   )
        printf( "scheme         =        plmde\n" );
    if ( H->scheme == HSCHEME_COLLELA )
        printf( "scheme         =      collela\n" );
    printf( "slope_type     = %12.5f\n",  H->slope_type     );
    printf( "-----------------------------------\n\n"       );
}

void process_args(int argc, char *argv[], hydro_args *H){
    int i;
    char inFile[512];

    i=1;
    inFile[0]='\0';
    set_default(H);
    while(i<argc){
      if(!strcmp(argv[i],"--help")){
        usage();
        i++;
        continue;
      }
      if(!strcmp(argv[i],"-i")){
        i++;
        if(i<argc){
          strncpy(inFile,argv[i],512);
          inFile[511]='\0';
          i++;
        }
        continue;
      }
      fprintf(stderr,"Unknown argument %s\n", argv[i]);
      i++;
    }
    if(inFile[0]!='\0'){
      process_argFile(inFile,H);
    }else{
      fprintf(stderr,"Missing argument -i is mandatory\n");
      exit(1);
    }
    print_args(H);
}
