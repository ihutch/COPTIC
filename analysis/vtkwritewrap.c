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

// Version to write multicomp. variables; vtkwritevector_ for single vector,
// and vtkwritevars_ for multiple vectors

void vtkwritevector(int *ifull, int *iuds, float *u,
		    float *xn, float *yn, float *zn,
		    int *ibinary,  int *vardims, char *filename,
		    char *ivarnames, int *str_alloc, int *idims_alloc,
		    int *ivars, int *icentering)
{ // u(ifull[1,2,3], idims_alloc[, nvars]) contains the ivars vectors.
  // iuds(1:3) gives the number of used points in each spatial dimension.
  // vardims(ivars) gives the number of used components of each vector.
  // varnames(str_alloc,ivars) contains null-terminated strings.
  // icentering(ivars) specifies whether each vector is node-centered.
  int binary=1;
  int i,j,k,m,n;
  int nvars=1;
  int *centering;
  int centering_storage[nvars];
  char **varnames;
  int dims_alloc;

  float *udata;
  float **vars;

  // If called through vtkwritevars_ pointers are non-NULL
  if (ivars)
    nvars = *ivars;

  if (idims_alloc) {
    dims_alloc = *idims_alloc;
  } else {
    dims_alloc = *vardims;
  }

  if (icentering) {
    centering = icentering;
  } else {
    *centering_storage = 1;   // Node-centered.
    centering = centering_storage;
  }

  varnames = (char**)malloc(nvars*sizeof(char*));
  if (str_alloc) {
    for (i=0; i<nvars; i++)
      varnames[i] = ivarnames + i*(*str_alloc);
  } else {
    varnames[0] = ivarnames;
  }

  udata = (float*)calloc(nvars*dims_alloc*iuds[0]*iuds[1]*iuds[2],
			 sizeof(float));
  vars = (float**)malloc(nvars*sizeof(float*));
  /*  Testing diagnostics.
  fprintf(stderr,
	  "Allocated udata. ifull %d %d %d; iuds %d %d %d. Pointer = %p\n"
	  ,ifull[0],ifull[1],ifull[2],
	  iuds[0],iuds[1],iuds[2],udata);
  fprintf(stderr,"xn %f %f %f %f %f\n",*xn, *(xn+1), *(xn+2), *(xn+3), *(xn+4));
  for (i=0; i<nvars; i++) {
    fprintf(stderr,"u(0,0,0,0,%d)=%f, vardims(%d)=%d, varnames[%d]=%s\n", i,
	    *(u+i*dims_alloc*ifull[0]*ifull[1]*ifull[2]),
	    i, vardims[i], i, varnames[i]);
  }
  for (i=0; i<100; i++) {
    fprintf(stderr,"u(0,0,0,0,%d)=%f\n", i,
	    *(u+i+0*dims_alloc*ifull[0]*ifull[1]*ifull[2]));
  }
  fprintf(stderr,"Calling write. %f  %f  %f \n",*xn,*(xn+iuds[0]),*(xn+iuds[1]));
  */

  for (n=0; n<nvars; n++) {
    vars[n] = udata + n*dims_alloc*iuds[0]*iuds[1]*iuds[2];
    for (i=0; i<3; i++) {
      // NB: this doesn't work properly unless all or none centered
      if (centering[n]==0)
	iuds[i] -= 1;
    }
    for (k=0;k<iuds[2];k++){
      for (j=0;j<iuds[1];j++){
	for (i=0;i<iuds[0];i++){
	  /* Compactify the data, and make components the fast arg. */
	  for (m=0;m<vardims[n];m++){
	    *(vars[n]+m+vardims[n]*(i+iuds[0]*(j+iuds[1]*k)))
	      =*(u+n*dims_alloc*ifull[0]*ifull[1]*ifull[2]
		 +i+ifull[0]*(j+ifull[1]*(k+ifull[2]*m)));
	  }
	}
      }
    }
    for (i=0; i<3; i++) {
      // NB: this doesn't work properly unless all or none centered
      if (centering[n]==0)
	iuds[i] += 1;
    }
  }

  if(*ibinary==0) binary=0;
  write_rectilinear_mesh( (const char * const)filename,
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
		 ivarnames, NULL, NULL, NULL, NULL);
}

void vtkwritevars_(int *ifull, int *iuds, float *u,
		   float *xn, float *yn, float *zn,
		   int *ibinary,  int *vardims, char *filename,
		   char *ivarnames, int *str_alloc, int *ivars,
		   int *idims_alloc, int *icentering) {
  vtkwritevector(ifull, iuds, u, xn, yn, zn, ibinary, vardims, filename,
		 ivarnames, str_alloc, idims_alloc, ivars, icentering);
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
