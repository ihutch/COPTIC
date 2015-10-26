c Declarations of the bbdycomm parameters for the call
c         call bbdy(iLs,ifull,iuds,u,nk,ndims,idims,
c     $        icoords,iLcoords,myside,myorig,
c     $        icommcart,mycartid,myid,lperiod)
c Number of dimensions and the number of blocks per dimension
c should be declared like this, before include
c       parameter (ndimsbbdy=3,idim1=3,idim2=2,idim3=2)
c The number of dimensions of the cartesian topology. (2 for 2d)
c      integer ndimsbbdy
c ifull full dimensions of u array
      integer ifull(ndimsbbdy)
c iuds used dimensions of u
      integer iuds(ndimsbbdy)
c nk      is iteration count.
c The length of each topology dimension
      integer idims(ndimsbbdy)
c For each topology dimension whether it is periodic or not.
      logical lperiod(ndimsbbdy)
c Cartesian topology coords of this process block (OUT)
      integer icoords(ndimsbbdy)
c structure of icoords (1,(idims(1)+1),(idims(1)+1)*(idims(2)+1),...)
      integer iLcoords(ndimsbbdy+1)
c My side length data, accounting for my position in the cluster.
      integer myside(ndimsbbdy)
c The origin of myblock such that u(myorig) is the first element.
      integer myorig
c icommcart is the id of the cartesian communicator (OUT).
      integer icommcart
c mycartid returns process id in cartesian topology communicator (OUT).
      integer mycartid,myid
c lperiod values must be set in the calling program after
c the above declarations. For lperiod, e.g.
c      data lperiod/ndimsbbdy*.false./
c      data ifull/ifd1,ifd2,ifd3/
c in bbdy initialization with the call to bbdydefine.
c--------------------------------------------------------------------
c Therefore the total example template for using bbdy is this
c 
c       parameter (ndimsbbdy=3,idim1=3,idim2=2,idim3=2)
c ndimsbbdy must be a parameter or all definitions including ndimsbbdy
c must be arguments of the subroutine. 
c       parameter (ifd1=100,ifd2=50,ifd3=200)
c       dimension u(ifd1,ifd2,ifd3)
c       include 'bbdydecl.f'
c And the values of ifull and iuds must be set either as data or 
c by assignment in the routine. For example:
c       integer ifull(ndimsbbdy),iLs(ndimsbbdy)
c       data ifull/ifd1,ifd2,ifd3/
c       data iuds/20,40,50/
c       data lperiod/ndimsbbdy*.false./
c       data idims/idim1,idim2,idim3/
c  ...
c       do nk=1,3
c         call bbdy(iLs,ifull,iuds,u,nk,ndimsbbdy,idims,
c     $        icoords,iLcoords,myside,myorig,
c     $        icommcart,mycartid,myid,lperiod)
c  ...
c       enddo
c--------------------------------------------------------------------
