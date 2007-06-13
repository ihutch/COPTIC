c Declarations of the blockbdy parameters for the call
c         call blockbdy(iLs,iuds,u,nk,iorig,nrd,idims,lperiod,
c     $        icoords,iLcoords,myside,
c     $        icommcart,mycartid,myid)
c Number of dimensions should be declared like this in the outer prog
c      parameter (ndims=2)
c Dimensional structure of u, for 2d should be (1,Li,Li*Lj), 
c 3d (1,Li,Li*Lj,Li*Lj*Lk) etc 
      integer iLs(ndims+1)
c ifull full dimensions of u
      integer ifull(ndims)
c iuds used dimensions of u
      integer iuds(ndims)
c nk      is iteration count.
c iorig(idims(1)+1,idims(2)+1,...) provides origin of block(i,j,..) 
c within u. Blocks must be of equal size except for the uppermost
c The top of uppermost, with iblock=idims(n)) is indicated
c by a value pointing to 1 minus the length of u in that dimension.
      integer iorig(idim1+1,idim2+1,idim3+1)
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
c icommcart is the id of the cartesian communicator (OUT).
      integer icommcart
c mycartid returns the process id in cartesian topology communicator (OUT).
      integer mycartid,myid

