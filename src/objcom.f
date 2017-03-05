! Object-data storage, for keeping track of boundary conditions 
! at object surfaces for the potential sor solution. 
! Requires ndimsdecl.f 
! Guess at needed total data size Lobjmax for each object.
! The object data consists of data enumerated as
!    (2=forward/backward)*(ndims)*(ndata=fraction,B/A,C/A)
!    + diagonal + potential contributors + region-code + inverse pointer
!    + intersection-object-code + extra-pointer.

      integer nobj_cij,ndata_cij,idgs_cij,ibdy_cij,iflag_cij
      integer iregion_cij,ipoint_cij,iinter_cij,iextra_cij,Lobjmax
! 3D choice needed for sharing common, I think. 
      parameter (ndata_cij=ndims)
! The total size of the structure member. Last= count of following
      parameter (nobj_cij=2*ndims*ndata_cij+7)
! Pointer to diagonal contributions
      parameter (idgs_cij=2*ndims*ndata_cij+1)
! Pointer to boundary contributions
      parameter (ibdy_cij=2*ndims*ndata_cij+2)
! Pointer to flags
      parameter (iflag_cij=2*ndims*ndata_cij+3)
! Pointer to region code of this node
      parameter (iregion_cij=2*ndims*ndata_cij+4)
! Pointer to inverse pointer within u,cij (but relative to 2,2,2...).
      parameter (ipoint_cij=2*ndims*ndata_cij+5)
! Pointer to intersection object code of this node (used for color)
      parameter (iinter_cij=2*ndims*ndata_cij+6)
! Pointer to additional intersection data for variable boundaries.
      parameter (iextra_cij=2*ndims*ndata_cij+7)
      parameter (Lobjmax=10000000)
      real dob_cij(nobj_cij,Lobjmax)
      integer oi_cij
      common /objcom/oi_cij,dob_cij
! The following makes it possible to refer to dob_cij in integer form.
      integer idob_cij(nobj_cij,Lobjmax)
      equivalence (idob_cij,dob_cij)


! Extra chained object data when necessary is in a second frame the
! same size immediately following the first, but pointed to by iextra.
! Its data is
!    (2=forward/backward)*(ndims)*(ndata=3=iobj, ijbin, coefoa)
!    [succeeding data not yet used].
!    iobj is the geometry object number, ijbin the facet address
!    and coefoa is the coefficient by which to multiply c when
!    adding it to the ibdy contribution.
