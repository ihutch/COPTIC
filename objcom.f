
c Object-data storage, for keeping track of boundary conditions 
c at object surfaces for the potential sor solution. 
c Guess at needed total data size Lobjmax for each object.
c The object data consists of data enumerated as
c    (2=forward/backward)*(ndims)*(ndata=fraction,B/A,C/A)
c    + diagonal + potential contributors + region-code + inverse pointer
c    + intersection-object-code + extra-pointer.

      integer ndims_sor,nobj_sor,ndata_sor,idgs_sor,ibdy_sor,iflag_sor
      integer iregion_sor,ipoint_sor,iinter_sor,iextra_sor,Lobjmax
c 3D choice needed for sharing common, I think. 
      parameter (ndims_sor=3,ndata_sor=ndims_sor)
c The total size of the structure member. Last= count of following
      parameter (nobj_sor=2*ndims_sor*ndata_sor+7)
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
c Pointer to intersection object code of this node (used for color)
      parameter (iinter_sor=2*ndims_sor*ndata_sor+6)
c Pointer to additional intersection data for variable boundaries.
      parameter (iextra_sor=2*ndims_sor*ndata_sor+7)
      parameter (Lobjmax=1000000)
      real dob_sor(nobj_sor,Lobjmax)
      integer oi_sor
      common /objcom/oi_sor,dob_sor
c The following makes it possible to refer to dob_sor in integer form.
      integer idob_sor(nobj_sor,Lobjmax)
      equivalence (idob_sor,dob_sor)


c Extra chained object data when necessary is in a second frame the
c same size immediately following the first, but pointed to by iextra.
c Its data is
c    (2=forward/backward)*(ndims)*(ndata=3=iobj, ijbin, coefoa)
c    [succeeding data not yet used].
c    iobj is the geometry object number, ijbin the facet address
c    and coefoa is the coefficient by which to multiply c when
c    adding it to the ibdy contribution.
