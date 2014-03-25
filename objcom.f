c Object-data storage, for keeping track of boundary conditions 
c at object surfaces for the potential sor solution. 
c Requires ndimsdecl.f 
c Guess at needed total data size Lobjmax for each object.
c The object data consists of data enumerated as
c    (2=forward/backward)*(ndims)*(ndata=fraction,B/A,C/A)
c    + diagonal + potential contributors + region-code + inverse pointer
c    + intersection-object-code + extra-pointer.

      integer nobj_cij,ndata_cij,idgs_cij,ibdy_cij,iflag_cij
      integer iregion_cij,ipoint_cij,iinter_cij,iextra_cij,Lobjmax
c 3D choice needed for sharing common, I think. 
      parameter (ndata_cij=ndims)
c The total size of the structure member. Last= count of following
      parameter (nobj_cij=2*ndims*ndata_cij+7)
c Pointer to diagonal contributions
      parameter (idgs_cij=2*ndims*ndata_cij+1)
c Pointer to boundary contributions
      parameter (ibdy_cij=2*ndims*ndata_cij+2)
c Pointer to flags
      parameter (iflag_cij=2*ndims*ndata_cij+3)
c Pointer to region code of this node
      parameter (iregion_cij=2*ndims*ndata_cij+4)
c Pointer to inverse pointer within u,cij (but relative to 2,2,2...).
      parameter (ipoint_cij=2*ndims*ndata_cij+5)
c Pointer to intersection object code of this node (used for color)
      parameter (iinter_cij=2*ndims*ndata_cij+6)
c Pointer to additional intersection data for variable boundaries.
      parameter (iextra_cij=2*ndims*ndata_cij+7)
      parameter (Lobjmax=1000000)
      real dob_cij(nobj_cij,Lobjmax)
      integer oi_cij
      common /objcom/oi_cij,dob_cij
c The following makes it possible to refer to dob_cij in integer form.
      integer idob_cij(nobj_cij,Lobjmax)
      equivalence (idob_cij,dob_cij)


c Extra chained object data when necessary is in a second frame the
c same size immediately following the first, but pointed to by iextra.
c Its data is
c    (2=forward/backward)*(ndims)*(ndata=3=iobj, ijbin, coefoa)
c    [succeeding data not yet used].
c    iobj is the geometry object number, ijbin the facet address
c    and coefoa is the coefficient by which to multiply c when
c    adding it to the ibdy contribution.
