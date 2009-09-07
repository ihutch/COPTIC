c***********************************************************************
c Block boundary communication.
      subroutine bbdy(iLs,ifull,iuds,u,kc,
     $     ndims,idims,lperiod,icoords,iLcoords,myside,myorig,
     $     icommcart,mycartid,myid)
c Dimensional structure of u, for 2d should be (1,Li,Lj), 
c 3d (1,Li,Li*Lj,Li*Lj*Lk) etc (last element may not be used)
      integer iLs(ndims+1)
c iuds used dimensions of u
      integer iuds(ndims)
c Inside this routine, u and iorig are referenced linearly.
      real u(*)
c kc      is iteration count which determines odd or even, and also 
c      if kc=-1 this is the final call: gather only.
c      if kc=-2 only (re)initialize the topology, cartesian communicator
c
      integer kc
c The number of dimensions of the cartesian topology. (2 for 2d) (IN)
      integer ndims
c The length of each topology dimension (number of blocks) (IN)
      integer idims(ndims)
c For each topology dimension whether it is periodic or not (IN)
      logical lperiod(ndims)
c Cartesian topology coords of this process block (OUT)
      integer icoords(ndims)
c structure of icoords (1,(idims(1)+1),(idims(1)+1)*(idims(2)+1),...)
      integer iLcoords(ndims+1)
c The origin of my block (OUT)
      integer myorig 
c My side length data, accounting for my position in the topology (OUT).
      integer myside(ndims)
c icommcart is the id of the cartesian topology communicator (OUT).
      integer icommcart
c mycartid is the process id in cartesian topology communicator (OUT).
      integer mycartid
c myid returns the process id in the MPI_COMM_WORLD (OUT)
      integer myid
c End of arguments.
c---------------------------------------------------------------
c The first time, set up the cartesian topology and return info 
c to calling process.
c Subsequently, do the boundary communication.
c---------------------------------------------------------------
c Start of local variables.
c      logical lalltoall
c      parameter (idebug=1,lalltoall=.false.)
      parameter (idebug=0)
c      parameter (lalltoall=.true.)
c Local storage:
      logical lreorder
c vector type ids for each dimension (maximum 10 dimensions)
      parameter(imds=10)
c facevector ids for each dimension; odd,even; bulk, top. 
c These are the handles to datatype that picks out data in the correct
c pattern, based on the u-address provided to the MPI call.
      integer iface(imds,2,2)
c Integer indication of whether we are bulk (1), or top (2)
      integer ibt(imds)
c stack pointers and lengths iodd(imds) points to the place in the 
c stack is where the odd vector starts. lodd is the vector length.
      integer iodd(imds),ieven(imds),lodd(imds),leven(imds)
c shift source and destination ids for each dimension left and right
      integer isdl(imds),iddl(imds),isdr(imds),iddr(imds)
c Right and left u-origins of each dimension for this block.
c      integer iobr(imds),iobl(imds)
c iside(2,imds) holds the side length of blocks in u for each dimension.
c iside(1,*) is the general value, iside(2,*) the uppermost value.
      integer iside(2,imds)
c irdims is the block count per dimension array with the order of the
c dimensions reversed to maintian fortran compatibility. 
c Ditto ircoords, lrperiod
      integer irdims(imds),ircoords(imds)
      logical lrperiod(imds)
c Scratch stack, whose length must be at least Prod_1^nd(iside(2,i)+1)
c Possibly this should be passed to the routine. 
      parameter (istacksize=100000)
      integer is(istacksize)
c      integer ktype(2**imds)
c Arrays for constructing ALLtoALL calls.
      parameter (maxprocs=1000)
      integer isdispls(maxprocs),irdispls(maxprocs)
      integer istypes(maxprocs),irtypes(maxprocs)
      integer iscounts(maxprocs),ircounts(maxprocs)
c      integer iconp(imds)
c      character*10 string
c Debugging arrays
c      parameter (ndebug=1000)
c      integer iaints(ndebug),iaadds(ndebug),iadats(ndebug)
c      integer isc(ndebug),isd(ndebug),ist(ndebug)
c      integer irc(ndebug),ird(ndebug),irt(ndebug)
c
c The iorig(idims(1)+1,idims(2)+1,...) provides origin of block(i,j,..) 
c within u. Blocks must be of equal size except for the uppermost
c The top of uppermost, with iblock=idims(n)) is indicated
c by a value pointing to 1 minus the length of u in that dimension.
c We declare it as a 1-d array for generality. This is the first time here:
c It must be of dimension greater than the number of processes (blocks)
      parameter (norigmax=1000)
      integer iorig(norigmax)
      common /iorigcom/iorig
      
      include 'mpif.h'
      integer status(MPI_STATUS_SIZE)
      integer iobindex
c Flag that we have called this routine once. Can't use multiple
c instances in one program. No way to reset this. Might include as
c an argument to allow us to reset. Not done yet.
      logical lflag
      data lflag/.false./
c Initialize iobindex to quite initialization warnings.
      data iobindex/1/
c This general save statement is necessary. But gives gfortran warnings.
      save

      if(kc.eq.-1) then 
c Just do the gather
         goto 100
      elseif(kc.eq.-2 .or. kc.eq.-3)then
c (Re)Initialize but don't communicate (-2) or do communicate (-3).
         lflag=.false.
      endif
      if(.not.lflag)then
c -----------------------------------------------------------------
c First time. Set up topology and calculate comm-types
c------------------------------------------------------------------
c------
c         write(*,*)'Setting up topology'
         if(ndims.gt.imds)then
            write(*,*)'MPI too many dimensions error',ndims
            goto 999
         endif
c Test if we fit in the number of processors in block storage.
         if((idims(1)+1)*(idims(2)+1).gt.norigmax)then
         write(*,*)'bbdy: Too many processes',idims(1),'x',idims(2),
     $        ' for norigmax=',norigmax
            goto 999
         endif
         nproc=1
         do n=1,ndims
c Count the processes needed for this topology
            nproc=nproc*idims(n)
         enddo
c Check the asked-for nproc
         if(nproc.gt.maxprocs) then
            write(*,*)'Too many processes',nproc,' Increase maxprocs.'
            goto 999
         endif
c End of safety checks
c------
c Initialize 
         call MPI_INITIALIZED(lflag,ierr)
         if(.not.lflag) call MPI_INIT( ierr )
         lflag=.true.
         call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
         call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
         if(nproc.ne.numprocs)then
            if(myid.eq.0)write(*,201)numprocs,
     $           nproc,(idims(n),n=1,ndims)
 201        format('WARNING: MPI processes',i4,
     $           ': don''t match this topology ',i4,':',6i3)
c Use MPI function to redimension block structure
            do ii=1,ndims
               idims(ii)=0
            enddo
            call MPI_DIMS_CREATE(numprocs,ndims,idims,ierr)
            if(ierr.eq.0)then
               nproc=numprocs
            else
               stop 'MPI_DIMS_CREATE error'
            endif
            if(myid.eq.0)write(*,'(''Reset to'',i4,'':'',6i3)')
     $           nproc,idims
         else
            if(myid.eq.0)write(*,'(''numprocs,idims()'',i4,10i4)')
     $           numprocs,(idims(n),n=1,ndims)
         endif
c End of topology idims resetting.
c-----
c Define mpi block structure.
c         write(*,*)'Calling bbdydefine'
         call bbdydefine(ndims,idims,ifull,iuds,iorig,iLs)
c Check roughly if the istacksize is enough.
c iLs is set in bbdydefine so it can't be earlier.
         if(2.*iLs(ndims).ge.istacksize) then
            write(*,*) 'Stack size',istacksize,
     $           ' possibly too small for arrays',iLs(ndims)
            write(*,*) 'Danger of overrun.'
            goto 999
         endif
c----
         mycartid=0
         iLcoords(1)=1
         do n=1,ndims
c Populate the block structure vector 
            if(n.gt.1) iLcoords(n)=iLcoords(n-1)*(idims(n-1)+1)
         enddo
c         write(*,*)'iLcoords',iLcoords

c Output some diagnostic data, perhaps.
         if(idebug.gt.0) call bbdyorigprint(ndims,idims,iorig)
         do n=1,ndims
c Calculate the block side lengths from the origins data.
            iside(1,n)=(iorig(1+iLcoords(n))-iorig(1))/iLs(n)+2
            kt=(1+ idims(n)*iLcoords(n))
            kn=(1+(idims(n)-1)*iLcoords(n))
            iside(2,n)=(iorig(kt)-iorig(kn))/iLs(n)+2
         enddo
c         write(*,*)'((iside(jj,ii),jj=1,2),ii=1,ndims)='
c         write(*,'(2i10)')((iside(jj,ii),jj=1,2),ii=1,ndims)
c-----

c Reverse the order of the dimensions in the calls to MPI_CART
c to compensate for the C-ordering in those.
         do n=1,ndims
            irdims(n)=idims(ndims-n+1)
            lrperiod(n)=lperiod(ndims-n+1)
c            irdims(n)=idims(n)
c            lrperiod(n)=lperiod(n)
         enddo
c Create topology
         lreorder=.true.
         call MPI_CART_CREATE(MPI_COMM_WORLD,ndims,irdims,lrperiod,
     $        lreorder,icommcart,ierr)
c Check if this process is unused
         if(icommcart.eq.MPI_COMM_NULL) goto 998
         call MPI_COMM_RANK(icommcart,mycartid, ierr )
c         write(*,*)'returned from create',ndims,idims,lperiod,ierr,
c     $        icommcart
c Determine my block cartesian coords, and hence start of my data.
c icoords are zero-origin not one-origin.
c Although we know idims,lperiod already, we want icoords.
c         call MPI_CART_GET(icommcart,ndims,irdims,lrperiod,
c     $        ircoords,ierr)
         call MPI_CART_COORDS(icommcart,mycartid,ndims,ircoords,ierr)
c Reverse the order of the dimensions in the calls to MPI_CART
         do n=1,ndims
c            icoords(n)=ircoords(n)
            icoords(n)=ircoords(ndims-n+1)
         enddo
        do n=1,ndims
            ibt(n)=1
            if(icoords(n).eq.idims(n)-1) ibt(n)=2
c Get my block side lengths, now knowing my cluster position.
            myside(n)=iside(ibt(n),n)
         enddo
         if(idebug.gt.0)then
            write(*,*)'u used dimensions',(iuds(k),k=1,ndims)
            write(*,*)'Blocks in each dim',(idims(k),k=1,ndims)
            write(*,*)'Rank=',myid,' icoords=',(icoords(k),k=1,ndims)
            write(*,*)'Block myside lengths',(myside(k),k=1,ndims)
            if(idebug.ge.2)then
               write(*,*)'u dim-structure',(iLs(k),k=1,3)
               write(*,*)'iLcoords=',iLcoords
            endif
         endif
         
c For all dimensions create vector types for communications.
c nn is the normal direction. 
         do nn=1,ndims
c Reuse the stack (drop previous data).
            ibeg=1
c Create the face vector indexes iodd(id),ieven(id), id=ndims.
c In iodd(1,...,ndims-1) are the subface indices 0...ndims-2,
c which however are not used externally.
            call bbdyfacecreate(nn,ndims,ibeg,is,iLs,myside,
     $           iodd,ieven,lodd,leven,id)
            if(idebug.gt.1)then
            write(*,*)'Face indices, direction nn=',nn,', block-dims'
     $           ,(myside(mod(k+nn-1,ndims)+1),k=1,2)
     $           ,' (group(1)= block-dim(1)/2)'
            write(*,*)'nn,id,iodd(id),lodd(id),ieven(id),leven(id)',
     $           nn,id,iodd(id),lodd(id),ieven(id),leven(id)
            write(*,*)' is(iodd (',nn,'))'
     $           ,(is(iodd(id)+ii),ii=0,lodd(id)-1)
            write(*,*)' is(ieven(',nn,'))'
     $           ,(is(ieven(id)+ii),ii=0,leven(id)-1)
            endif
c     Make a buffer of the correct length telling data lengths (1 each)
c     in scratch stack.
            iblens=ibeg
            do i=1,lodd(id)
               is(ibeg)=1
               ibeg=ibeg+1
            enddo
c Create the new data types, odd and even.
            call MPI_TYPE_INDEXED(lodd(id),is(iblens),is(iodd(id)),
     $           MPI_REAL,iface(nn,1,ibt(nn)),ierr)
            call MPI_TYPE_COMMIT(iface(nn,1,ibt(nn)),ierr)
            call MPI_TYPE_INDEXED(leven(id),is(iblens),is(ieven(id)),
     $           MPI_REAL,iface(nn,2,ibt(nn)),ierr)
            call MPI_TYPE_COMMIT(iface(nn,2,ibt(nn)),ierr)
         enddo

         do n=1,ndims
c Determine Shift ids.
c Shifts compensating for the dimension order reversal.
            call MPI_CART_SHIFT(icommcart,ndims-n,1,
     $           isdr(n),iddr(n),ierr)
            call MPI_CART_SHIFT(icommcart,ndims-n,-1,
     $           isdl(n),iddl(n),ierr)
         enddo
         iobindex=1
         ioffset=0
         do n=1,ndims
c Calculate my block origin index which is 
c   (1+icoords(1)*iLcoords(1)+icoords(2)*iLcoords(2),...)
            iobindex=iobindex+icoords(n)*iLcoords(n)
            ioffset=ioffset+iLs(n)
         enddo
         myorig=iorig(iobindex)
c--------------------------------------------------------------------
c Create the types, pointers and counts for block gathering.
c         write(*,*)'calling bbdyblockcreate'
         call bbdygatherparam(ndims,iLs,iside,
     $     mycartid,nproc,idims,iLcoords,ioffset,idebug,
     $     isdispls,istypes,iscounts,irdispls,irtypes,ircounts)
c--------------------------------------------------------------------
         if(idebug.ge.2) write(*,*)'End of initialization'
c         return
      endif
c If this is an unused process, do nothing.
      if(mycartid.eq.-1) return
c If this is a pure initialization, return
      if(kc.eq.-2) return
c--------------------------------------------------------------------
c (First and) Subsequent calls. Do the actual communication
c--------------------------------------------------------------------
c Even-odd pacing. We use all blocks each step but we transmit either
c the odd or the even data to the left, even or odd to right.

c If kc is odd, then odd parity (1,1,1, ...) is active this step.
c We need to read in the previous even data to be ready for it.
c Send odd/even data to right, receive odd/even data from left
      ko=mod(kc,2)+1
c Send even/odd data to left, receive even/odd data from right
      ke=mod(kc+1,2)+1
      itag=100
c      write(*,*)'ko,ke,iobindex',ko,ke,iobindex
c Origin of left face (the same for all dimensions)
      iolp=iorig(iobindex)
      do n=1,ndims
c         write(*,*)'iolp,n,ndims,iLcoords',iolp,n,ndims,iLcoords
c Origin of right face (i.e. of block to right)
         iorp=iorig(iobindex+iLcoords(n))
c      write(*,*)'iorp=',iorp
c iolp is for receiving +1 shift, iorp is for sending +1 shift.
         iolm=iolp+iLs(n)
         iorm=iorp+iLs(n)
c iorm is for receiving -1 shift, iolm is for sending -1 shift.
         if(idebug.ge.2)then
         if(n.eq.1) write(*,*)' n,ko,ke,    iolp,iorp,iolm,iorm,',
     $        'iddr,isdr,iddl(n),isdl(n)'
         write(*,'(3i3,8i6)') n,ko,ke,iolp,iorp,iolm,iorm,
     $        iddr(n),isdr(n),iddl(n),isdl(n)
         endif
c Send odd/even data to right, receive odd/even data from left
         if(iddr(n).ne.-1) call MPI_SEND(u(iorp),1,iface(n,ko,ibt(n)),
     $        iddr(n),itag,icommcart,ierr)
         if(isdr(n).ne.-1) call MPI_RECV(u(iolp),1,iface(n,ko,ibt(n)),
     $        isdr(n),itag,icommcart,status,ierr) 
c         call MPI_SENDRECV(
c     $        u(ior),1,iface(n,ko),iddr(n),itag,
c     $        u(iol),1,iface(n,ko),isdr(n),itag,
c     $        icommcart,status,ierr)
c         write(*,*)'Second n,ko,ke,iol,ior,isdl,isdr',
c     $        n,ko,ke,iol,ior,iddl(n),isdl(n)
c Send even/odd data to left, receive even/odd data from right
c shift in the direction left (-1).
         if(iddl(n).ne.-1) call MPI_SEND(u(iolm),1,iface(n,ke,ibt(n)),
     $        iddl(n),itag,icommcart,ierr)
         if(isdl(n).ne.-1) call MPI_RECV(u(iorm),1,iface(n,ke,ibt(n)),
     $        isdl(n),itag,icommcart,status,ierr)

      enddo
c Synchronize the processes. Did not help with segfaults.
c      call MPI_BARRIER(icommcart,ierr)
      return
c------------------------------------------------------------------
c Special cases determined by the value of kc.
 100  continue
c kc=-1. Do the block exchanging.
      if(idebug.ge.1)then
         write(*,*)'     np,sendcts,senddp,sendtp,recvcts,recvdp,recvtp'
         write(*,'(7i8)')(np,iscounts(np),isdispls(np)/4,
     $        mod(istypes(np),10000),
     $        ircounts(np),irdispls(np)/4,
     $        mod(irtypes(np),10000),
     $        np=1,nproc)

      endif

      call MPI_ALLTOALLW(u,iscounts,isdispls,istypes,
     $     u,ircounts,irdispls,irtypes,
     $     icommcart,ierr)

      return
c Unused process.
 998  continue
      mycartid=-1
      write(*,*)'Process',myid,' not used in cartesian communicator',
     $     (idims(i),i=1,ndims)
      return
c Exception stop:
 999  continue
      call MPI_FINALIZE(ierr)
      stop

      end
c******************************************************************
      subroutine bbdyalternate(is,ibeg,iodd,ieven,lodd,leven,id,
     $     idlen,iainc)
c Construct a new index array starting at ibeg in stack: is
c The index array is laminated from odd and even prior
c types which start at is(iodd(id-1)), is(ieven(id-1)) , 
c ibeg on return points to next unused stack position.

c storage index space
      integer is(*)
c current address of beginning of unused storage
      integer ibeg
c addresses and lengths of the odd and even types for dimensions nd.
c The values for i=id-1 must exist on input. Values for id are returned.
      integer iodd(*),ieven(*),lodd(*),leven(*)
c dimension at present issue
      integer id
c block length of present dimension.
      integer idlen
c full array index increment in this dimension
      integer iainc

      iodd(id)=ibeg
      iinc=0
      do i=1,idlen/2
c         write(*,*)'iinc=',iinc,iainc
         call bbdycatstep(is,ibeg,iodd(id-1),lodd(id-1),iinc)
         iinc=iinc+iainc
         call bbdycatstep(is,ibeg,ieven(id-1),leven(id-1),iinc)
         iinc=iinc+iainc
      enddo
      lodd(id)=ibeg-iodd(id)

      iinc=0
      ieven(id)=ibeg
      do i=1,idlen/2
         call bbdycatstep(is,ibeg,ieven(id-1),leven(id-1),iinc)
         iinc=iinc+iainc
         call bbdycatstep(is,ibeg,iodd(id-1),lodd(id-1),iinc)
         iinc=iinc+iainc
      enddo
      leven(id)=ibeg-ieven(id)

      end
c*****************************************************************
c Concatenate data from istart with length ilen into is at ibeg,
c with iinc added to it.
      subroutine bbdycatstep(is,ibeg,istart,ilen,iinc)
      integer is(*)
      do j=1,ilen
         is(ibeg)=is(istart-1+j)+iinc
         ibeg=ibeg+1
      enddo
      end
c****************************************************************
      subroutine bbdyfacecreate(nn,nd,ibeg,is,iLs,myside,
     $     iodd,ieven,lodd,leven,id)
c nn: normal dimension, nd: total dimensions.
      integer nn,nd,id
c ibeg: stack counter, is: stack array, iLs: u-structure,
c myside: block length in dimension n.
      integer ibeg,is(*),iLs(nd+1),myside(nd)
c iodd,ieven: returns the vector representing the face indices of is,
c in the elements iodd(nd), ieven(nd) (subfaces earlier).
      integer iodd(nd),ieven(nd),lodd(nd),leven(nd)

c      write(*,*)'Bbdyfacecreate: nn,nd,ibeg,myside',nn,nd,ibeg,myside
c      write(*,*) iLs
c Create odd/even vectors that get elements from face normal to nn.
c Iterate over all the dimensions

c Zeroth dimension addresses and lengths: 
      iodd(1)=ibeg
      is(ibeg)=0
      lodd(1)=1
      ibeg=ibeg+1
c even is a null length
      ieven(1)=ibeg-1
      leven(1)=0
      id=1
c If ndims>1 iterate to the higher face dimension.
      do  nc=nn,nn+nd-2
c     The actual dimension number: n
         n=mod(nc,nd)+1
c     Count of the dimension id starting at two.
         id=nc-nn+2
         call bbdyalternate(is,ibeg,iodd,ieven,lodd,leven,id,
     $           myside(n),iLs(n))
c         write(*,*)'ibeg,iodd(id),lodd(id),ieven(id),leven(id)',
c     $             ibeg,iodd(id),lodd(id),ieven(id),leven(id)
c     write(*,*)'iodd',(is(i),i=iodd(id),lodd(id)+iodd(id)-1)
c     write(*,*)'ieven',(is(i),i=ieven(id),leven(id)+ieven(id)-1)
c         else
      enddo
c     Now iodd(id),lodd(id),ieven(id),leven(id) 
c     are the even and odd vectors face normal to dimension nn.

c     write(*,*)'Bbdyfacecreate end: ibeg,iodd,lodd,ieven,leven',
c     $       nn,iodd,lodd,ieven,leven
c     write(*,*)'iodd',(is(i),i=iodd(id),lodd(id)+iodd(id)-1)
c     write(*,*)'ieven',(is(i),i=ieven(id),leven(id)+ieven(id)-1)
      end

c****************************************************************
      subroutine bbdyorigprint(ndims,idims,iorig)
      integer ndims
      integer idims(ndims),iorig(*)
      integer kk(3)
      character*40 string
c This diagnostic works only for up to 3-d arrays.
         if(ndims.le.3) then
            do j=3,1,-1
               if(ndims.lt.j)then 
                  kk(j)=1
               else
                  kk(j)=idims(j)+1
               endif
            enddo
            write(*,*)'j,k,Block origins(i)='
            write(string,'(''('',i4''i8)'')')kk(1)+2
            write(*,string)
     $           ((j,k,(iorig(i+kk(1)*((j-1)+kk(2)*(k-1))),
     $           i=1,kk(1)),j=1,kk(2)),k=1,kk(3))
         endif

         end
c*******************************************************************
c Return the block origins for a multidimensional block arrangement.
      subroutine bbdydefine(ndims,idims,ifull,iuds,iorig,iLs)
c ndims: number of dimensions, 
c idims(ndims) number of blocks in each dimension, 
c ifull(ndims) full lengths of declared array in each dimension. 
c iuds(ndims) used lengths in each dimension
c iorig: return block origin displacements. (OUT)
c iLs(ndims+1) full-length step structure of declared array (OUT)
      integer ndims,idims(ndims),ifull(ndims),iuds(ndims)
      integer iLs(ndims+1)
c Presumed effective dimensions of iorig are 
c    (idims(1)+1,idims(2)+1,...)
      integer iorig(*)
      character*20 string
c Define the block origin addresses
      iorig(1)=1
      ibeg=2
      ifn=1
      iLs(1)=1
      do n=1,ndims
         isz=2+(iuds(n)-2)/idims(n)
         if(mod(isz,2).eq.1)isz=isz-1
         istep=(isz-2)*ifn
         ilen=ibeg-1
         do i=1,idims(n)
            if(i.eq.idims(n))then
               iinc=(iuds(n)-2)*ifn
            else
               iinc=i*istep
            endif
c            write(*,*)'iinc,ifn,istep,isz=',iinc,ifn,istep,isz
            call bbdycatstep(iorig,ibeg,1,ilen,iinc)
         enddo
         ifn=ifull(n)*ifn
         iLs(n+1)=ifn
      enddo
c      write(*,*)'iorig='
      write(string,'(''('',i5,''i8)'')')(idims(1)+1)
c      write(*,string)(iorig(jj),jj=1,ibeg-1)
      end
c********************************************************************
      subroutine bbdyblockcreate(ndims,ktype,iLs,iside,iSIZEOFREAL)
c For a system of ndims dimensions, (IN)
c Create MPI type block handles in the array ktype (OUT)
c when the array structures are iLs (IN)
c the lengths of the block sides in dimension id are iside(2,id) (IN)
c with iside(1,id) referring to the bulk and iside(2,id) to the top.
      
c      parameter(imds=10)
      integer iLs(ndims+1)
c iside(1,*) is the general value, iside(2,*) the uppermost value.
      integer iside(2,ndims)
      
c ktype stores the type handles in the order
c 0-d; 1d bulk, 1d top; 2d bulk (1b,1t), 2d top (1b,1t), ...
c So for the nn-th level the start of the types is at 2**nn.
c ktype must have length (at least) 2**(ndims+1)-1
      integer ktype(*)
      include 'mpif.h'

      call MPI_TYPE_SIZE(MPI_REAL,iSIZEOFREAL,ierr)
c      write(*,*)'iSIZEOFREAL=',iSIZEOFREAL
c zeroth dimension type is MPI_REAL
      ktype(1)=MPI_REAL
      inew=2
      do nn=1,ndims
         iprior=2**(nn-1)
c         write(*,*)'bbdyblockcreate: nn,iprior,inew,ktype(iprior)'
c     $        ,nn,iprior,inew,ktype(iprior)
c laminate to higher dimension
         call bbdyblam(ndims,ktype,inew,iprior,nn,iside,
     $        iLs,iSIZEOFREAL)
c         write(*,*)'nn,ithis,inew,ktype(ithis...)',
c     $        nn,2**nn,inew,(ktype(ith),ith=(2**nn),2**(nn+1)-1)
      enddo
c The types we use now start at ktype(2**ndims), 2**ndims of them.
      end

c*********************************************************************
      subroutine bbdyblam(ndims,ktype,inew,iprior,nn,iside,
     $     iLs,iSIZEOFREAL)
c      parameter(imds=10)
c iside(1,*) is the general value, iside(2,*) the uppermost value.
      integer iside(2,ndims)
      integer iLs(ndims+1)
      integer ktype(*)
      include 'mpif.h'
c This was the cause of a lot of heartache. 
c Necessary for MPI_TYPE_CREATE_HVECTOR, wrong for TYPE_HVECTOR
c      integer(KIND=MPI_ADDRESS_KIND) istride
c The alternative for non F90 capable is integer*8 on 64 bit
c or integer*4 on 32 bit machines.
      ilen=1
      istride=iLs(nn)*iSIZEOFREAL
c For the bulk and top cases of this dimension
      do ibt=1,2
c Because we do not count the boundaries, 
c the count is 2 less than iside.
         icount=iside(ibt,nn)-2
c For all types in prior level, laminate enough together
c At each level there are 2**(nn-1) types in prior level
         do iold=iprior,iprior+2**(nn-1)-1
c            write(*,*)'iold,icount,istride,ktype(iold),inew',
c     $           iold,icount,istride,ktype(iold),inew
c This call needs istride to be MPI_ADDRESS_KIND
c            call MPI_TYPE_CREATE_HVECTOR(icount,ilen,istride,
c This call assumes istride is integer*4 (I think):
            call MPI_TYPE_HVECTOR(icount,ilen,istride,
     $           ktype(iold),ktype(inew),ierr)
c We only commit the top level that we are going to use.
            if(nn.eq.ndims)call MPI_TYPE_COMMIT(ktype(inew),ierr)
            inew=inew+1
         enddo
      enddo

      end
c*********************************************************************
      subroutine bbdycoords(nn,ndims,idims,icoords,ith,iLcoords,ionp)
c Obtain the cartesian coordinates of block for process nn,
c in ndims dimensions, whose lengths are idims(ndims).
c Return it in icoords. 
c Return the type index (zero-based), i.e.
c whether we are in bulk or on boundary in ith.
c Also the iorig index, ionp, for this block.
      integer nn,ndims
      integer icoords(ndims),idims(ndims)
      integer ith
      integer iLcoords(ndims+1)
      in=nn-1
      ith=0
      ionp=1
      do nd=1,ndims
         iquot=in/idims(nd)
         icoords(nd)=in-iquot*idims(nd)
         if(icoords(nd)+1.eq.idims(nd)) ith=ith+2**(nd-1)
c         write(*,*)'nd,icoords(nd),in,iquot',nd,icoords(nd),in,iquot
         in=iquot
c Since iorig has dimensions 1+idims this is needed:
         ionp=ionp+icoords(nd)*iLcoords(nd)
      enddo
c      write(*,*)'nn,ith,ndims,idims',nn,ith,ndims,idims
      end
c****************************************************************
      subroutine bbdygatherparam(ndims,iLs,iside,
     $     mycartid,nproc,idims,iLcoords,ioffset,idebug,
     $     isdispls,istypes,iscounts,irdispls,irtypes,ircounts)
c For ndims dimensions, with array structure iLs, and block lengths
c given by iside, set the send and receive alltoall(w) parameters:
c isdispls, istypes,iscounts, irdispls, irtypes, ircounts,
c for process mycartid, within total of nproc, cartesian communicator
c dimensions idims, iLcoords, ioffset, idebug.

      integer ndims,mycartid,ioffset
      integer iLs(ndims+1),iLcoords(ndims),idims(ndims)
      integer iside(2,ndims)
      integer isdispls(nproc),istypes(nproc),iscounts(nproc)
      integer irdispls(nproc),irtypes(nproc),ircounts(nproc)

      parameter (norigmax=1000)
      integer iorig(norigmax)
      common /iorigcom/iorig

      parameter (imds=10)
      integer ktype(2**(imds+1)),iconp(imds)
      logical lalltoall
      data lalltoall/.true./

c      include 'mpif.h'

      call bbdyblockcreate( ndims,ktype,iLs,iside,isizeofreal)
      ith0=2**ndims
      if(idebug.ge.1) write(*,*)'Block types:',(ktype(ith),
     $     ith=ith0,2*ith0-1)

c Create the required arrays for the ALLTOALL
c First find out where we are and return the type index in ithi.
      call bbdycoords(mycartid+1,ndims,idims,iconp,ithi,
     $     iLcoords,ionpmine)
      do np=1,nproc
c Here we need to send only the active length of data, so the origin
c is offset.
         isdispls(np)=(iorig(ionpmine)-1+ioffset)*isizeofreal
         istypes(np)=ktype(ith0+ithi)
c Get positions and type indexes for other processes
         call bbdycoords(np,ndims,idims,iconp,ithj,
     $        iLcoords,ionp)
         irdispls(np)=(iorig(ionp)-1+ioffset)*isizeofreal
         irtypes(np)=ktype(ith0+ithj)
         if(lalltoall) then
c All to all
            if(np.eq.mycartid+1)then
c Don't send to or receive from myself.
               iscounts(np)=0
               ircounts(np)=0
            else
               iscounts(np)=1
               ircounts(np)=1
            endif
         else
c Gather to process 0 (fortran index 1)
            if(np.eq.1.and.mycartid+1.ne.np)then
               iscounts(np)=1
            else
               iscounts(np)=0
            endif
            if(mycartid.eq.0.and.mycartid+1.ne.np)then
               ircounts(np)=1
            else
               ircounts(np)=0
            endif
         endif
      enddo
c Don't send to or receive from anyone test.
c Process 1 sends to process 0 test
c         if(mycartid.eq.0) ircounts(2)=1
c         if(mycartid.eq.1) iscounts(1)=1
c Process 2 sends to process 0 test
c         if(mycartid.eq.0) ircounts(3)=1
c         if(mycartid.eq.2) iscounts(1)=1
      end
c********************************************************************
c Abstraction to isolate mpi calls.
      subroutine mpifinalize(ierr)
      call MPI_FINALIZE(ierr)
      end
c*******************************************************************
      subroutine mpicommsize(numprocs,ierr)
      include 'mpif.h'
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
      end
