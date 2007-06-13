c Object-data storage. Guess at needed size Lobjmax
c The object data consists of data enumerated as
c Inverse pointer.
c (2=forward/backward)*(ndims)*(ndata=fraction,B/A,C/A)
c + diagonal + potential contributors + flags.

      integer ndims_sor,nobj_sor
c 3D choice needed for sharing common, I think. 
      parameter (ndims_sor=3,ndata_sor=3)
c The total size of the structure member. Last= count of following
      parameter (nobj_sor=2*ndims_sor*ndata_sor+4)
c Pointer to diagonal contributions
      parameter (idgs_sor=2*ndims_sor*ndata_sor+1)
c Pointer to boundary contributions
      parameter (ibdy_sor=2*ndims_sor*ndata_sor+2)
c Pointer to flags
      parameter (iflag_sor=2*ndims_sor*ndata_sor+3)
c Pointer to region code
      parameter (iregion_sor=2*ndims_sor*ndata_sor+4)
      parameter (Lobjmax=100000)
      real dob_sor(nobj_sor,lobjmax)
      integer oi_sor
      common /objcom/oi_sor,dob_sor
c The following makes it possible to refer to dob_sor in integer form.
      integer idob_sor(nobj_sor,lobjmax)
      equivalence (idob_sor,dob_sor)
