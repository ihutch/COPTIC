#include "visit_writer.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/* Fortran-callable wrapper for visit-writer routine to write out 
   a rectilinear array in vtk format that can be read into visit. 
   On entry:
      filename is the name of the file to be written, to which is appended vtk
      ifull[3] is the full dimensions of the fortran array
      ufortran(ifull(1),ifull(2),ifull(3))=u[ifull[2]][ifull[1]][ifull[0]]
      iuds[3] is the used dimensions of this array. 
      xn, yn, zn are the start addresses of the arrays of positions,
      ibinary=0 indicates write ascii; otherwise binary.
*/

void vtkwritescalar_(int *ifull, int *iuds, float *u, 
		     float *xn, float *yn, float *zn, 
		   int *ibinary, char *filename, char *varnames )
{ 
  int binary=1;
  int i,j,k;
  float *udata;
  udata=calloc(iuds[0]*iuds[1]*iuds[2],sizeof(float));
  /*  Testing diagnostics.
      fprintf(stderr,
	  "Allocated udata. ifull %d %d %d; iuds %d %d %d. Pointer = %p\n"
	  ,ifull[0],ifull[1],ifull[2],
	  iuds[0],iuds[1],iuds[2],udata);
  fprintf(stderr,"xn %f %f %f %f %f\n",*xn, *(xn+1), *(xn+2), *(xn+3), *(xn+4));
  fprintf(stderr,"Calling write. %f  %f  %f \n",*xn,*(xn+iuds[0]),*(xn+iuds[1])); */

  for (k=0;k<iuds[2];k++){
    for (j=0;j<iuds[1];j++){
      for (i=0;i<iuds[0];i++){
	/* Just compactify the data */
	*(udata+i+j*iuds[0]+k*iuds[1]*iuds[0])
	=*(u+i+j*ifull[0]+k*ifull[1]*ifull[0]);
      }
    }
  }

  int nvars=1;
  int vardims[] = { 1 };   // One scalar
  int centering[] = { 1 }; // node centered
  if(*ibinary==0) binary=0;
  write_rectilinear_mesh( (const char * const)filename , 
			  binary, iuds, xn, yn, zn,
			  nvars , vardims, centering,
			  (const char * const *)&varnames, &udata);
  free(udata);
}

/* Version to write a vector variable input as separate variables. */

void vtkwritevector_(int *ifull, int *iuds, float *u, 
		     float *xn, float *yn, float *zn,
	    int *ibinary,  int *ndims, char *filename, char *varnames )
{ // On entry u(ifull[1,2,3], ndims) is the vector (dimension ndims) on the grid
  int binary=1;
  int i,j,k,m;
  int nvars=1;
  int vardims[nvars];   // nvars vector of length ndims
  int centering[nvars]; // node centering

  float *udata;

  *vardims=*ndims;
  *centering=1;   // Centered.

  udata=calloc((*ndims)*iuds[0]*iuds[1]*iuds[2],sizeof(float));
  /*  Testing diagnostics.
      fprintf(stderr,
	  "Allocated udata. ifull %d %d %d; iuds %d %d %d. Pointer = %p\n"
	  ,ifull[0],ifull[1],ifull[2],
	  iuds[0],iuds[1],iuds[2],udata);
  fprintf(stderr,"xn %f %f %f %f %f\n",*xn, *(xn+1), *(xn+2), *(xn+3), *(xn+4));
  fprintf(stderr,"Calling write. %f  %f  %f \n",*xn,*(xn+iuds[0]),*(xn+iuds[1])); */

  for (k=0;k<iuds[2];k++){
    for (j=0;j<iuds[1];j++){
      for (i=0;i<iuds[0];i++){
	/* Compactify the data, and make the ndims triples the fast arg. */
	for (m=0;m<ndims[0];m++){
	  *(udata+m+ndims[0]*(i+iuds[0]*(j+iuds[1]*k)))
	    =*(u+i+ifull[0]*(j+ifull[1]*(k+ifull[2]*m)));
	}
      }
    }
  }

  if(*ibinary==0) binary=0;
  write_rectilinear_mesh( (const char * const)filename , 
			  binary, iuds, xn, yn, zn,
			  nvars , vardims, centering,
			  (const char * const *)&varnames, &udata);
  free(udata);
}
