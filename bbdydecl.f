c Declarations of the bbdycomm parameters for the call
c         call bbdy(iLs,iuds,u,nk,iorig,ndims,idims,lperiod,
c     $        icoords,iLcoords,myside,myorig,
c     $        icommcart,mycartid,myid)
c Number of dimensions and the number of blocks per dimension
c should be declared like this, before include
c       parameter (ndims=3,idim1=3,idim2=2,idim3=2)
c Declared Dimensional structure of u must be (Li,Lj,Lk,..) passed using
c iLs, which embodies pointer steps. For 2d iLs=(1,Li,Li*Lj), 
c 3d (1,Li,Li*Lj,Li*Lj*Lk) etc 
      integer iLs(ndims+1)
c iuds used dimensions of u
      integer iuds(ndims)
c nk      is iteration count.
c iorig(idims(1)+1,idims(2)+1,...) provides origin of block(i,j,..) 
c within u. Blocks must be of equal size except for the uppermost
c The top of uppermost, with iblock=idims(n)) is indicated
c by a value pointing to 1 minus the length of u in that dimension.
c example only     integer iorig(idim1+1,idim2+1,idim3+1)
c The number of dimensions of the cartesian topology. (2 for 2d)
      integer ndims
c The length of each topology dimension
      integer idims(ndims)
c For each topology dimension whether it is periodic or not.
      logical lperiod(ndims)
c Cartesian topology coords of this process block (OUT)
      integer icoords(ndims)
c structure of icoords (1,(idims(1)+1),(idims(1)+1)*(idims(2)+1),...)
      integer iLcoords(ndims+1)
c My side length data, accounting for my position in the cluster.
      integer myside(ndims)
c The origin of myblock such that u(myorig) is the first element.
      integer myorig
c icommcart is the id of the cartesian communicator (OUT).
      integer icommcart
c mycartid returns process id in cartesian topology communicator (OUT).
      integer mycartid,myid
c iLs and lperiod values must be set in the calling program after
c the above declarations. For lperiod, e.g.
c      data lperiod/ndims*.false./
c iLs are normally set using already known declared dimensions
c      data ifull/ifd1,ifd2,ifd3/
c with the call which also defines iorig:
c      call bbdydefine(ndims,idims,ifull,iuds,iorig,iLs)
c--------------------------------------------------------------------
c Therefore the total example template for using bbdy is this
c 
c       parameter (ndims=3,idim1=3,idim2=2,idim3=2)
c ndims must be a parameter.
c idim1,idim2 etc must be parameters unless iorig is passed to caller.
c       integer iorig(idim1+1,idim2+1,idim3+1)
c       parameter (ifd1=100,ifd2=50,ifd3=200)
c       dimension u(ifd1,ifd2,ifd3)
c       include 'bbdydecl.f'
c And the values of ifull and iuds must be set either as data or 
c by assignment in the routine. For example:
c       integer ifull(ndims)
c       data ifull/ifd1,ifd2,ifd3/
c       data iuds/20,40,50/
c       data lperiod/ndims*.false./
c       data idims/idim1,idim2,idim3/
c  ...
c       call bbdydefine(ndims,idims,ifull,iuds,iorig,iLs)
c       do nk=1,3
c         call bbdy(iLs,iuds,u,nk,iorig,ndims,idims,lperiod,
c     $        icoords,iLcoords,myside,myorig,
c     $        icommcart,mycartid,myid)
c  ...
c       enddo
c--------------------------------------------------------------------
