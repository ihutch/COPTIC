c Declarations of the bbdycomm parameters for the call
c         call bbdy(iLs,ifull,iuds,u,nk,ndims,idims,
c     $        icoords,iLcoords,myside,myorig,
c     $        icommcart,mycartid,myid,lperiod)
c Number of dimensions and the number of blocks per dimension
c should be declared like this, before include
c       parameter (ndimsdecl=3,idim1=3,idim2=2,idim3=2)
c The number of dimensions of the cartesian topology. (2 for 2d)
c      integer ndimsdecl
c Declared Dimensional structure of u must be (Li,Lj,Lk,..) passed using
c iLs, which embodies pointer steps. For 2d iLs=(1,Li,Li*Lj), 
c 3d (1,Li,Li*Lj,Li*Lj*Lk) etc 
      integer iLs(ndimsdecl+1)
c ifull full dimensions of u array
      integer ifull(ndimsdecl)
c iuds used dimensions of u
      integer iuds(ndimsdecl)
c nk      is iteration count.
c The length of each topology dimension
      integer idims(ndimsdecl)
c For each topology dimension whether it is periodic or not.
      logical lperiod(ndimsdecl)
c Cartesian topology coords of this process block (OUT)
      integer icoords(ndimsdecl)
c structure of icoords (1,(idims(1)+1),(idims(1)+1)*(idims(2)+1),...)
      integer iLcoords(ndimsdecl+1)
c My side length data, accounting for my position in the cluster.
      integer myside(ndimsdecl)
c The origin of myblock such that u(myorig) is the first element.
      integer myorig
c icommcart is the id of the cartesian communicator (OUT).
      integer icommcart
c mycartid returns process id in cartesian topology communicator (OUT).
      integer mycartid,myid
c lperiod values must be set in the calling program after
c the above declarations. For lperiod, e.g.
c      data lperiod/ndimsdecl*.false./
c iLs are normally set transparently using declared dimensions
c      data ifull/ifd1,ifd2,ifd3/
c in bbdy initialization with the call to bbdydefine.
c--------------------------------------------------------------------
c Therefore the total example template for using bbdy is this
c 
c       parameter (ndimsdecl=3,idim1=3,idim2=2,idim3=2)
c ndimsdecl must be a parameter.
c       parameter (ifd1=100,ifd2=50,ifd3=200)
c       dimension u(ifd1,ifd2,ifd3)
c       include 'bbdydecl.f'
c And the values of ifull and iuds must be set either as data or 
c by assignment in the routine. For example:
c       integer ifull(ndimsdecl)
c       data ifull/ifd1,ifd2,ifd3/
c       data iuds/20,40,50/
c       data lperiod/ndimsdecl*.false./
c       data idims/idim1,idim2,idim3/
c  ...
c       do nk=1,3
c         call bbdy(iLs,ifull,iuds,u,nk,ndimsdecl,idims,
c     $        icoords,iLcoords,myside,myorig,
c     $        icommcart,mycartid,myid,lperiod)
c  ...
c       enddo
c--------------------------------------------------------------------
