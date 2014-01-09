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
  udata=(float*)calloc(iuds[0]*iuds[1]*iuds[2],sizeof(float));
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

void vtkwritevector(int *ifull, int *iuds, float *u,
		     float *xn, float *yn, float *zn,
		     int *ibinary,  int *ndims, char *filename, char *ivarnames,
		     int *stringlength, int *invars,
		     int *ivardims, int *ivardims_alloc,
		     int *icentering)
{ // On entry u(ifull[1,2,3], ndims) is the vector (dimension ndims) on the grid
  // Alternately u(ifull[1,2,3], vardims_alloc, nvars)
  // varnames(stringlength,nvars) assumed to contain null-terminated strings
  int binary=1;
  int i,j,k,m,n;
  int nvars=1;
  int *vardims;
  int vardims_storage[nvars];
  int *centering;
  int centering_storage[nvars]; // node centering
  char **varnames;
  int vardims_alloc;

  float *udata;
  float **vars;

  if (invars)
    nvars = *invars;

  if (ivardims) {
    vardims = ivardims;
  } else {
    *vardims_storage = *ndims;
    vardims = vardims_storage;
  }
  if (ivardims_alloc) {
    vardims_alloc = *ivardims_alloc;
  } else {
    vardims_alloc = *vardims;
  }

  if (icentering) {
    centering = icentering;
  } else {
    *centering_storage = 1;   // Node-centered.
    centering = centering_storage;
  }

  varnames = (char**)malloc(nvars*sizeof(char*));
  if (stringlength) {
    for (i=0; i<nvars; i++)
      varnames[i] = ivarnames + i*(*stringlength);
  } else {
    varnames[0] = ivarnames;
  }

  udata = (float*)calloc(nvars*vardims_alloc*iuds[0]*iuds[1]*iuds[2],sizeof(float));
  vars = (float**)malloc(nvars*sizeof(float*));
  /*  Testing diagnostics.
      fprintf(stderr,
	  "Allocated udata. ifull %d %d %d; iuds %d %d %d. Pointer = %p\n"
	  ,ifull[0],ifull[1],ifull[2],
	  iuds[0],iuds[1],iuds[2],udata);
  fprintf(stderr,"xn %f %f %f %f %f\n",*xn, *(xn+1), *(xn+2), *(xn+3), *(xn+4));
  fprintf(stderr,"Calling write. %f  %f  %f \n",*xn,*(xn+iuds[0]),*(xn+iuds[1])); */

  for (n=0; n<nvars; n++) {
    vars[n] = udata + n*vardims_alloc*iuds[0]*iuds[1]*iuds[2];
    for (k=0;k<iuds[2];k++){
      for (j=0;j<iuds[1];j++){
	for (i=0;i<iuds[0];i++){
	  /* Compactify the data, and make components the fast arg. */
	  for (m=0;m<vardims[n];m++){
	    *(vars[n]+m+vardims[n]*(i+iuds[0]*(j+iuds[1]*k)))
	      =*(u+n*vardims_alloc*ifull[0]*ifull[1]*ifull[2]
		 +i+ifull[0]*(j+ifull[1]*(k+ifull[2]*m)));
	  }
	}
      }
    }
  }

  if(*ibinary==0) binary=0;
  write_rectilinear_mesh( (const char * const)filename , 
			  binary, iuds, xn, yn, zn,
			  nvars , vardims, centering,
			  (const char * const *)varnames, vars);
  free(varnames);
  free(udata);
  free(vars);
}

void vtkwritevector_(int *ifull, int *iuds, float *u,
		     float *xn, float *yn, float *zn,
				int *ibinary,  int *ndims, char *filename, char *ivarnames) {
  vtkwritevector(ifull, iuds, u, xn, yn, zn, ibinary, ndims, filename,
		 ivarnames, NULL, NULL, NULL, NULL, NULL);
}

void vtkwritevars_(int *ifull, int *iuds, float *u,
		     float *xn, float *yn, float *zn,
		     int *ibinary,  int *ndims, char *filename, char *ivarnames,
		     int *stringlength, int *invars,
		     int *ivardims, int *ivardims_alloc,
		     int *icentering) {
  vtkwritevector(ifull, iuds, u, xn, yn, zn, ibinary, ndims, filename,
		 ivarnames, stringlength, invars, ivardims, ivardims_alloc,
		 icentering);
}


// Wrapper for writing vtk point mesh; see visit_writer.c for input spec.
void vtkwritescalarpoints_(int *npts, float *u, float *pts,
			   int *ibinary, char *filename, char *varnames)
{
  int binary=1;
  int nvars=1;
  int vardims[] = { 1 };   // One scalar
  if(*ibinary==0) binary=0;
  write_point_mesh( (const char * const)filename,
			  binary, *npts, pts,
			  nvars , vardims,
			  (const char * const *)&varnames, &u);
}

// Wrapper for writing unstruct. vtk mesh; see visit_writer.c for input spec.
void vtkwritescalarfacets_(int *npts, float *u, float *pts,
			   int *ncells, int *celltypes, int *conn,
			   int *centering,
			   int *ibinary, char *filename, char *varnames)
{
  int binary=1;
  int nvars=1;
  int vardims[] = { 1 };   // One scalar
  if(*ibinary==0) binary=0;
  write_unstructured_mesh((const char * const)filename, binary, *npts,
			  pts, *ncells, celltypes, conn,
			  nvars, vardims, centering,
			  (const char * const *)&varnames, &u);
}
