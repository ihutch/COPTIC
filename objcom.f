
c Object-data storage, for keeping track of boundary conditions 
c at object surfaces for the potential sor solution. 
c Guess at needed total data size Lobjmax for each object.
c The object data consists of data enumerated as
c (2=forward/backward)*(ndims)*(ndata=fraction,B/A,C/A)
c + diagonal + potential contributors + region-code + inverse pointer
c + intersection-code

      integer ndims_sor,nobj_sor
c 3D choice needed for sharing common, I think. 
      parameter (ndims_sor=3,ndata_sor=ndims_sor)
c The total size of the structure member. Last= count of following
      parameter (nobj_sor=2*ndims_sor*ndata_sor+6)
c Pointer to diagonal contributions
      parameter (idgs_sor=2*ndims_sor*ndata_sor+1)
c Pointer to boundary contributions
      parameter (ibdy_sor=2*ndims_sor*ndata_sor+2)
c Pointer to flags
      parameter (iflag_sor=2*ndims_sor*ndata_sor+3)
c Pointer to region code of this node
      parameter (iregion_sor=2*ndims_sor*ndata_sor+4)
c Pointer to inverse pointer within u,cij (but relative to 2,2,2...).
      parameter (ipoint_sor=2*ndims_sor*ndata_sor+5)
c Pointer to intersection code of this node
      parameter (iinter_sor=2*ndims_sor*ndata_sor+6)
      parameter (Lobjmax=1000000)
      real dob_sor(nobj_sor,lobjmax)
      integer oi_sor
      common /objcom/oi_sor,dob_sor
c The following makes it possible to refer to dob_sor in integer form.
      integer idob_sor(nobj_sor,lobjmax)
      equivalence (idob_sor,dob_sor)
