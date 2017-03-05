! Declarations of the bbdycomm parameters for the call
!         call bbdy(iLs,ifull,iuds,u,nk,ndims,idims,
!     $        icoords,iLcoords,myside,myorig,
!     $        icommcart,mycartid,myid,lperiod)
! Number of dimensions and the number of blocks per dimension
! should be declared like this, before include
!       parameter (ndimsbbdy=3,idim1=3,idim2=2,idim3=2)
! The number of dimensions of the cartesian topology. (2 for 2d)
!      integer ndimsbbdy
! ifull full dimensions of u array
      integer ifull(ndimsbbdy)
! iuds used dimensions of u
      integer iuds(ndimsbbdy)
! nk      is iteration count.
! The length of each topology dimension
      integer idims(ndimsbbdy)
! For each topology dimension whether it is periodic or not.
      logical lperiod(ndimsbbdy)
! Cartesian topology coords of this process block (OUT)
      integer icoords(ndimsbbdy)
! structure of icoords (1,(idims(1)+1),(idims(1)+1)*(idims(2)+1),...)
      integer iLcoords(ndimsbbdy+1)
! My side length data, accounting for my position in the cluster.
      integer myside(ndimsbbdy)
! The origin of myblock such that u(myorig) is the first element.
      integer myorig
! icommcart is the id of the cartesian communicator (OUT).
      integer icommcart
! mycartid returns process id in cartesian topology communicator (OUT).
      integer mycartid,myid
! lperiod values must be set in the calling program after
! the above declarations. For lperiod, e.g.
!      data lperiod/ndimsbbdy*.false./
!      data ifull/ifd1,ifd2,ifd3/
! in bbdy initialization with the call to bbdydefine.
!--------------------------------------------------------------------
! Therefore the total example template for using bbdy is this
! 
!       parameter (ndimsbbdy=3,idim1=3,idim2=2,idim3=2)
! ndimsbbdy must be a parameter or all definitions including ndimsbbdy
! must be arguments of the subroutine. 
!       parameter (ifd1=100,ifd2=50,ifd3=200)
!       dimension u(ifd1,ifd2,ifd3)
!       include 'bbdydecl.f'
! And the values of ifull and iuds must be set either as data or 
! by assignment in the routine. For example:
!       integer ifull(ndimsbbdy),iLs(ndimsbbdy)
!       data ifull/ifd1,ifd2,ifd3/
!       data iuds/20,40,50/
!       data lperiod/ndimsbbdy*.false./
!       data idims/idim1,idim2,idim3/
!  ...
!       do nk=1,3
!         call bbdy(iLs,ifull,iuds,u,nk,ndimsbbdy,idims,
!     $        icoords,iLcoords,myside,myorig,
!     $        icommcart,mycartid,myid,lperiod)
!  ...
!       enddo
!--------------------------------------------------------------------
