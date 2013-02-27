#include <stdio.h>
#include <stdlib.h>
#include "hydro_struct.h"
#include "hydro_defs.h"

void writeVis(char* fname, double *u, double dx, double dy, int nvar, int nx, int ny){
  int i,j,nv,argInd;
  FILE *vis;
  char outName[53], *ext;
  char name[30];
  
  strncpy(outName,fname,50);
  outName[49]='\0';
  if(!(ext=strrchr(outName,'.'))) ext=strchr(outName,'\0');
  sprintf(ext,".vts");
  vis=fopen(outName,"w");
  if(vis==NULL){
    fprintf(stderr,"Could not open file %s\n",outName);
    return;
  }
  fprintf(vis, "<?xml version=\"1.0\"?>\n");
  fprintf(vis, "<VTKFile type=\"StructuredGrid\">\n");
  fprintf(vis, "<StructuredGrid WholeExtent=\"%ld %ld %ld %ld %ld %ld\">\n", (long)0, (long)nx, (long)0, (long) ny, (long)0, (long)0);
  fprintf(vis, "<Piece Extent=\"%ld %ld %ld %ld %ld %ld\">\n", (long)0, nx, (long)0, ny, (long)0, (long)0);
  fprintf(vis, "<PointData></PointData>\n");

  // declaration of the variable list

  fprintf(vis, "<CellData Scalars=\"");

  for (nv = 0; nv < nvar; nv++){
    switch(nv){
    case VARRHO:
	fprintf(vis, " Density");
	break;
    case VARVX:
	fprintf(vis, " MomX");
	break;
    case VARVY:
	fprintf(vis, " MomY");
	break;
    case VARPR:
	fprintf(vis, " ENE");
	break;
    default:
      fprintf(vis, " var%d", nv);
      break;
    }
  }
  fprintf(vis, "\">\n");
  
  name[0] = 0;

  for (nv = 0; nv < nvar; nv++){
    switch(nv){
    case VARRHO:
	sprintf(name, " Density");
	break;
    case VARVX:
	sprintf(name, " MomX");
	break;
    case VARVY:
	sprintf(name, " MomY");
	break;
    case VARPR:
	sprintf(name, " ENE");
	break;
    default:
      sprintf(name, " var%d", nv);
      break;
    }

    //Definition of the cell values
    fprintf(vis, "<DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n", name);

    // the image is the interior of the computed domain
    for (j = 0; j < ny; j++){
      for (i = 0; i < nx; i++){
	fprintf(vis, "%lf ", u[(nv*ny+j)*nx+i]);
      }
      fprintf(vis, "\n");
    }
    fprintf(vis, "</DataArray>\n");
  }
  fprintf(vis, "</CellData>\n");
  fprintf(vis, "<Points>\n");
  fprintf(vis,"<DataArray type=\"Float32\" format=\"ascii\" NumberOfComponents=\"3\">\n");

  for (j = 0; j <= ny; j++){
    for (i = 0; i <= nx; i++){
      fprintf(vis, "%f %f %f\n", i * dx, j * dy, 0.0);
    }
  }

  fprintf(vis, "</DataArray>\n");
  fprintf(vis, "</Points>\n");
  fprintf(vis, "</Piece>\n");
  fprintf(vis, "</StructuredGrid>\n");
  fprintf(vis, "</VTKFile>\n");
  fclose(vis);
}

