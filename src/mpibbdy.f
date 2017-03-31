!***********************************************************************
! Block boundary communication.
      subroutine bbdy(iLs,ifull,iuds,u,kc,
     $     ndims,idims,icoords,iLcoords,myside,myorig,
     $     icommcart,mycartid,myid,lperiod)
      implicit none
! The number of dimensions of the cartesian topology. (2 for 2d) (IN)
      integer ndims
! Dimensional structure of u, for 2d should be (1,Li,Lj), 
! 3d (1,Li,Li*Lj,Li*Lj*Lk) etc (last element may not be used)
      integer iLs(ndims+1),ifull(ndims)
! iuds used dimensions of u
      integer iuds(ndims)
! Inside this routine, u and iorig are referenced linearly.
      real u(*)
! kc      is iteration count which determines odd or even, and also 
!      if kc=-1 this is the final call: gather only.
!      if kc=-2 only (re)initialize the topology, cartesian communicator
!
      integer kc
! The length of each topology dimension (number of blocks) (IN/OUT)
      integer idims(ndims)
! For each topology dimension whether it is periodic or not (IN)
      logical lperiod(ndims)
! Cartesian topology coords of this process block (OUT)
      integer icoords(ndims)
! structure of icoords (1,(idims(1)+1),(idims(1)+1)*(idims(2)+1),...)
      integer iLcoords(ndims+1)
! The origin of my block (OUT)
      integer myorig 
! My side length data, accounting for my position in the topology (OUT).
      integer myside(ndims)
! icommcart is the id of the cartesian topology communicator (OUT).
      integer icommcart
! mycartid is the process id in cartesian topology communicator (OUT).
      integer mycartid
! myid returns the process id in the MPI_COMM_WORLD (OUT)
      integer myid
! End of arguments.
!---------------------------------------------------------------
! The first time, set up the cartesian topology and return info 
! to calling process.
! Subsequently, do the boundary communication.
!---------------------------------------------------------------
! Start of local variables.
      integer idebug
!      parameter (lalltoall=.true.)
      logical lreorder
! vector type ids for each dimension (maximum 10 dimensions)
      integer imds
      parameter(imds=10)
! facevector ids for each dimension; odd,even; bulk, top. 
! These are the handles to datatype that picks out data in the correct
! pattern, based on the u-address provided to the MPI call.
      integer iface(imds,2,2)
! Integer indication of whether we are bulk (1), or top (2)
      integer ibt(imds)
! stack pointers and lengths iodd(imds) points to the place in the 
! stack is where the odd vector starts. lodd is the vector length.
      integer iodd(imds),ieven(imds),lodd(imds),leven(imds)
! shift source and destination ids for each dimension left and right
      integer isdl(imds),iddl(imds),isdr(imds),iddr(imds)
! Right and left u-origins of each dimension for this block.
!      integer iobr(imds),iobl(imds)
! iside(2,imds) holds the side length of blocks in u for each dimension.
! iside(1,*) is the general value, iside(2,*) the uppermost value.
      integer iside(2,imds)
! irdims is the block count per dimension array with the order of the
! dimensions reversed to maintian fortran compatibility. 
! Ditto ircoords, lrperiod
      integer irdims(imds),ircoords(imds)
      logical lrperiod(imds)
! Scratch stack, whose length must be at least Prod_1^nd(iside(2,i)+1)
! Possibly this should be passed to the routine. 
      integer istacksize
      parameter (istacksize=3000000)
      integer is(istacksize)
      integer nstackneed
!      integer ktype(2**imds)
! Arrays for constructing ALLtoALL calls.
      integer maxprocs
      parameter (maxprocs=4096)
      integer isdispls(maxprocs),irdispls(maxprocs)
      integer istypes(maxprocs),irtypes(maxprocs)
      integer iscounts(maxprocs),ircounts(maxprocs)
!      integer iconp(imds)
!      character*10 string
! Debugging arrays
!      parameter (ndebug=1000)
!      integer iaints(ndebug),iaadds(ndebug),iadats(ndebug)
!      integer isc(ndebug),isd(ndebug),ist(ndebug)
!      integer irc(ndebug),ird(ndebug),irt(ndebug)
!
! The iorig(idims(1)+1,idims(2)+1,...) provides origin of block(i,j,..) 
! within u. That is, the bottom of the boundary cells in each dimension.
! Three blocks in the x-direction, one in y-direction on 16x4 grid:
! 1  ^----^
! 2      ^----^
! 3          ^------^
!
!  4 ................   ^
!  3 X...X...X.....X.   |
!  2 ................   |
!  1 X...X...X.....X.   v  1
!    1234567890123456
! X-Side-lengths 6[,6],8
! Blocks must be of equal size except for the uppermost
! The top of uppermost, with iblock=idims(n)) is indicated
! by a value pointing to the length of u in that dimension minus 1.
! Each side runs from e.g. iorig(i,j,...)+1 to iorig(i+1,j,...)
! We declare it as a 1-d array for generality. This is the first time here:
! It must be of dimension greater than the number of processes (blocks)
      integer norigmax
      parameter (norigmax=maxprocs)
      integer iorig(norigmax)
      common /iorigcom/iorig
      
      include 'mpif.h'
      integer status(MPI_STATUS_SIZE)
      integer iobindex
      data idebug/0/
! Flag that we have called this routine once.
      logical lflag
      data lflag/.false./
! Initialize iobindex to quiet initialization warnings.
      data iobindex/1/
! This general save statement gives gfortran warnings from the
! individual save statements in mpif.h. So instead save individually all
! the non-parameter or data variables defined above.
!      save
      save iface,ibt,iodd,ieven,lodd,leven,isdl,iddl,isdr,iddr
      save iside,irdims,ircoords,lrperiod,is,status
      save isdispls,irdispls,istypes,irtypes,iscounts,ircounts
! Various variables used as counters etc later. Not all saved now.
      integer i,ibeg,iblens,id,ierr,ii,ioffset,iolm,iolp,iorm,iorp
      integer itag,k,ke,kn,ko,kt,n,nn,np,nprcsses,nproc
      save nproc

      if(kc.eq.-1) then 
! Just do the gather
         goto 100
      elseif(kc.eq.-2 .or. kc.eq.-3)then
! (Re)Initialize but don't communicate (-2) or do communicate (-3).
         lflag=.false.
      endif
      if(.not.lflag)then
! -----------------------------------------------------------------
! First time. Set up topology and calculate comm-types
!------------------------------------------------------------------
!------
!         write(*,*)'Setting up topology'
         if(ndims.gt.imds)then
            write(*,*)'MPI too many dimensions error',ndims
            goto 999
         endif
         nproc=1
         do n=1,ndims
! Count the processes needed for this topology
            nproc=nproc*idims(n)
         enddo
! Initialize 
         call MPI_INITIALIZED(lflag,ierr)
         if(.not.lflag) call MPI_INIT( ierr )
         lflag=.true.
         call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
         call MPI_COMM_SIZE( MPI_COMM_WORLD, nprcsses, ierr )
         if(nprcsses.gt.maxprocs) then
            write(*,*)'Too many processes',nprcsses
     $           ,' Increase maxprocs parameter in mpibbdy.f'
            goto 999
         endif
! Check the asked-for nproc and if not equal to nprcsses, reapportion.
!         if(nproc.ne.nprcsses)then
! Only reapportion if too many processes were asked for.
! This fails in the coptic context because then some don't know phi!
         if(nproc.gt.nprcsses)then
            if(myid.eq.0 .and. idims(1).ne.99)write(*,201)nprcsses,
     $           nproc,(idims(n),n=1,ndims)
 201        format(' MPI processes',i4,
     $           ': don''t match this topology ',i7,':',6i4)
! Use MPI function to redimension block structure
            do ii=1,ndims
               if(iuds(ii).lt.8)then
                  idims(ii)=1
               else
                  idims(ii)=0
               endif
            enddo
            call MPI_DIMS_CREATE(nprcsses,ndims,idims,ierr)
            if(ierr.eq.0)then
               nproc=nprcsses
            else
               stop 'MPI_DIMS_CREATE error'
            endif
            if(myid.eq.0)write(*,'('' (Re)set MPI processes to'',i4
     $,'' :'',6i3)')nproc,idims
         else
            if(myid.eq.0 .and. nprcsses.gt.1)
     $           write(*,'('' MPI processes,idims()'',i4,'':'',10i4)')
     $           nprcsses,(idims(n),n=1,ndims)
         endif
! End of topology idims resetting.
!-----
! Define mpi block structure.
!         write(*,*)'Calling bbdydefine'
         call bbdydefine(ndims,idims,ifull,iuds,iorig,iLs)
!----
         mycartid=0
         iLcoords(1)=1
         do n=1,ndims
! Populate the block structure vector 
            if(n.gt.1) iLcoords(n)=iLcoords(n-1)*(idims(n-1)+1)
         enddo
!         write(*,*)'iLcoords',iLcoords

! Output some diagnostic data, perhaps.
         if(idebug.gt.0) call bbdyorigprint(ndims,idims,iorig)
         do n=1,ndims
! Calculate the block side lengths from the origins data.
            iside(1,n)=(iorig(1+iLcoords(n))-iorig(1))/iLs(n)+2
            kt=(1+ idims(n)*iLcoords(n))
            kn=(1+(idims(n)-1)*iLcoords(n))
            iside(2,n)=(iorig(kt)-iorig(kn))/iLs(n)+2
         enddo
!         write(*,*)'((iside(jj,ii),jj=1,2),ii=1,ndims)='
!         write(*,'(2i10)')((iside(jj,ii),jj=1,2),ii=1,ndims)
!-----

! Reverse the order of the dimensions in the calls to MPI_CART
! to compensate for the C-ordering in those.
         do n=1,ndims
            irdims(n)=idims(ndims-n+1)
            lrperiod(n)=lperiod(ndims-n+1)
!            irdims(n)=idims(n)
!            lrperiod(n)=lperiod(n)
         enddo
! Create topology
         lreorder=.true.
         call MPI_CART_CREATE(MPI_COMM_WORLD,ndims,irdims,lrperiod,
     $        lreorder,icommcart,ierr)
! Check if this process is unused
         if(icommcart.eq.MPI_COMM_NULL) goto 998
         call MPI_COMM_RANK(icommcart,mycartid, ierr )
!         write(*,*)'returned from create',ndims,idims,lperiod,ierr,
!     $        icommcart
! Determine my block cartesian coords, and hence start of my data.
! icoords are zero-origin not one-origin.
! Although we know idims,lperiod already, we want icoords.
!         call MPI_CART_GET(icommcart,ndims,irdims,lrperiod,
!     $        ircoords,ierr)
         call MPI_CART_COORDS(icommcart,mycartid,ndims,ircoords,ierr)
! Reverse the order of the dimensions in the calls to MPI_CART
         do n=1,ndims
!            icoords(n)=ircoords(n)
            icoords(n)=ircoords(ndims-n+1)
         enddo
         nstackneed=1
         do n=1,ndims
            ibt(n)=1
            if(icoords(n).eq.idims(n)-1) ibt(n)=2
! Get my block side lengths, now knowing my cluster position.
            myside(n)=iside(ibt(n),n)
            nstackneed=nstackneed*(myside(n)+1)
         enddo
! Check my stack size needs. If too great, stop.
         if(nstackneed.gt.istacksize)idebug=1
         if(idebug.gt.0)then
            write(*,*)'u used dimensions',(iuds(k),k=1,ndims)
            write(*,*)'Blocks in each dim',(idims(k),k=1,ndims)
            write(*,*)'Rank=',myid,' icoords=',(icoords(k),k=1,ndims)
            write(*,*)'Block myside lengths',(myside(k),k=1,ndims)
            if(idebug.ge.2)then
               write(*,*)'u dim-structure',(iLs(k),k=1,3)
               write(*,*)'iLcoords=',iLcoords
            endif
            if(nstackneed.gt.istacksize)then
               write(*,*)'nstackneed too large',nstackneed
     $              ,' larger than istacksize',istacksize
               write(*,*)'****Error****. Increase istacksize in mpibbdy'
               goto 999
            endif
         endif
         
! For all dimensions create vector types for communications.
! nn is the normal direction. 
         do nn=1,ndims
! Reuse the stack (drop previous data).
            ibeg=1
! Create the face vector indexes iodd(id),ieven(id), id=ndims.
! In iodd(1,...,ndims-1) are the subface indices 0...ndims-2,
! which however are not used externally.
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
!     Make a buffer of the correct length telling data lengths (1 each)
!     in scratch stack.
            iblens=ibeg
            do i=1,lodd(id)
               is(ibeg)=1
               ibeg=ibeg+1
            enddo
! Create the new data types, odd and even.
            call MPI_TYPE_INDEXED(lodd(id),is(iblens),is(iodd(id)),
     $           MPI_REAL,iface(nn,1,ibt(nn)),ierr)
            call MPI_TYPE_COMMIT(iface(nn,1,ibt(nn)),ierr)
            call MPI_TYPE_INDEXED(leven(id),is(iblens),is(ieven(id)),
     $           MPI_REAL,iface(nn,2,ibt(nn)),ierr)
            call MPI_TYPE_COMMIT(iface(nn,2,ibt(nn)),ierr)
         enddo

         do n=1,ndims
! Determine Shift ids.
! Shifts compensating for the dimension order reversal.
            call MPI_CART_SHIFT(icommcart,ndims-n,1,
     $           isdr(n),iddr(n),ierr)
            call MPI_CART_SHIFT(icommcart,ndims-n,-1,
     $           isdl(n),iddl(n),ierr)
         enddo
         iobindex=1
         ioffset=0
         do n=1,ndims
! Calculate my block origin index which is 
!   (1+icoords(1)*iLcoords(1)+icoords(2)*iLcoords(2),...)
            iobindex=iobindex+icoords(n)*iLcoords(n)
            ioffset=ioffset+iLs(n)
         enddo
         myorig=iorig(iobindex)
!--------------------------------------------------------------------
! Create the types, pointers and counts for block gathering.
         call bbdygatherparam(ndims,iLs,iside,
     $     mycartid,nproc,idims,iLcoords,ioffset,idebug,
     $     isdispls,istypes,iscounts,irdispls,irtypes,ircounts)
!--------------------------------------------------------------------
         if(idebug.ge.1) write(*,*)'End of initialization'
!         return
      endif
! If this is an unused process, do nothing.
      if(mycartid.eq.-1) return
! If this is a pure initialization, return
      if(kc.eq.-2) return
!--------------------------------------------------------------------
! (First and) Subsequent calls. Do the actual communication
!--------------------------------------------------------------------
! Even-odd pacing. We use all blocks each step but we transmit either
! the odd or the even data to the left, even or odd to right.

! If kc is odd, then odd parity (1,1,1, ...) is active this step.
! We need to read in the previous even data to be ready for it.
! Send odd/even data to right, receive odd/even data from left
      ko=mod(kc,2)+1
! Send even/odd data to left, receive even/odd data from right
      ke=mod(kc+1,2)+1
      itag=100
!      write(*,*)'ko,ke,iobindex',ko,ke,iobindex
! Origin of left face (the same for all dimensions)
      iolp=iorig(iobindex)
      do n=1,ndims
!         write(*,*)'iolp,n,ndims,iLcoords',iolp,n,ndims,iLcoords
! Origin of right face (i.e. of block to right)
         iorp=iorig(iobindex+iLcoords(n))
!      write(*,*)'iorp=',iorp
! iolp is for receiving +1 shift, iorp is for sending +1 shift.
         iolm=iolp+iLs(n)
         iorm=iorp+iLs(n)
! iorm is for receiving -1 shift, iolm is for sending -1 shift.
         if(idebug.ge.2)then
            if(n.eq.1) write(*,*)' n,ko,ke,    iolp,iorp,iolm,iorm,',
     $           'iddr,isdr,iddl(n),isdl(n)'
            write(*,'(''Communicate'',3i3,8i6)')n,ko,ke,iolp,iorp,iolm
     $           ,iorm,iddr(n),isdr(n),iddl(n),isdl(n)
         endif
! Send odd/even data to right, receive odd/even data from left
!         if(iddr(n).ne.-1) call MPI_SEND(u(iorp),1,iface(n,ko,ibt(n)),
!     $        iddr(n),itag,icommcart,ierr)
!         if(isdr(n).ne.-1) call MPI_RECV(u(iolp),1,iface(n,ko,ibt(n)),
!     $        isdr(n),itag,icommcart,status,ierr) 
! This type of call is necessary for periodic passes, e.g. to oneself.
         call MPI_SENDRECV(
     $        u(iorp),1,iface(n,ko,ibt(n)),iddr(n),itag,
     $        u(iolp),1,iface(n,ko,ibt(n)),isdr(n),itag,
     $        icommcart,status,ierr)
!         if(idebug.ge.2)write(*,*)'Second n,ko,ke,iol,ior,isdl,isdr',
!             n,ko,ke,iol,ior,iddl(n),isdl(n)
! Send even/odd data to left, receive even/odd data from right
! shift in the direction left (-1).
!         if(iddl(n).ne.-1) call MPI_SEND(u(iolm),1,iface(n,ke,ibt(n)),
!     $        iddl(n),itag,icommcart,ierr)
!         if(isdl(n).ne.-1) call MPI_RECV(u(iorm),1,iface(n,ke,ibt(n)),
!     $        isdl(n),itag,icommcart,status,ierr)
         call MPI_SENDRECV(
     $        u(iolm),1,iface(n,ke,ibt(n)),iddl(n),itag,
     $        u(iorm),1,iface(n,ke,ibt(n)),isdl(n),itag,
     $        icommcart,status,ierr)
      enddo
      return
!------------------------------------------------------------------
! Special cases determined by the value of kc.
 100  continue
! kc=-1. Do the block exchanging.
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
! Unused process.
 998  continue
      mycartid=-1
      write(*,'(a,i5,a,3i4)')' Process',myid
     $     ,' not in cartesian communicator',(idims(i),i=1,ndims)
      if(myid.eq.0)stop
      return
! Exception stop:
 999  continue
      call MPI_FINALIZE(ierr)
      stop

      end
!******************************************************************
      subroutine bbdyalternate(is,ibeg,iodd,ieven,lodd,leven,id,
     $     idlen,iainc)
! Construct a new index array starting at ibeg in stack: is
! The index array is laminated from odd and even prior
! types which start at is(iodd(id-1)), is(ieven(id-1)) , 
! ibeg on return points to next unused stack position.

! storage index space
      integer is(*)
! current address of beginning of unused storage
      integer ibeg
! addresses and lengths of the odd and even types for dimensions nd.
! The values for i=id-1 must exist on input. Values for id are returned.
      integer iodd(*),ieven(*),lodd(*),leven(*)
! dimension at present issue
      integer id
! block length of present dimension.
      integer idlen
! full array index increment in this dimension
      integer iainc

      iodd(id)=ibeg
      iinc=0
      do i=1,idlen/2
!         write(*,*)'iinc=',iinc,iainc
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
!*****************************************************************
! Concatenate data from istart with length ilen into is at ibeg,
! with iinc added to it.
      subroutine bbdycatstep(is,ibeg,istart,ilen,iinc)
      integer is(*)
      do j=1,ilen
         is(ibeg)=is(istart-1+j)+iinc
         ibeg=ibeg+1
      enddo
      end
!****************************************************************
      subroutine bbdyfacecreate(nn,nd,ibeg,is,iLs,myside,
     $     iodd,ieven,lodd,leven,id)
! nn: normal dimension, nd: total dimensions.
      integer nn,nd,id
! ibeg: stack counter, is: stack array, iLs: u-structure,
! myside: block length in dimension n.
      integer ibeg,is(*),iLs(nd+1),myside(nd)
! iodd,ieven: returns the vector representing the face indices of is,
! in the elements iodd(nd), ieven(nd) (subfaces earlier).
      integer iodd(nd),ieven(nd),lodd(nd),leven(nd)

!      write(*,*)'Bbdyfacecreate: nn,nd,ibeg,myside',nn,nd,ibeg,myside
!      write(*,*) iLs
! Create odd/even vectors that get elements from face normal to nn.
! Iterate over all the dimensions

! Zeroth dimension addresses and lengths: 
      iodd(1)=ibeg
      is(ibeg)=0
      lodd(1)=1
      ibeg=ibeg+1
! even is a null length
      ieven(1)=ibeg-1
      leven(1)=0
      id=1
! If ndims>1 iterate to the higher face dimension.
      do  nc=nn,nn+nd-2
!     The actual dimension number: n
         n=mod(nc,nd)+1
!     Count of the dimension id starting at two.
         id=nc-nn+2
         call bbdyalternate(is,ibeg,iodd,ieven,lodd,leven,id,
     $           myside(n),iLs(n))
!         write(*,*)'ibeg,iodd(id),lodd(id),ieven(id),leven(id)',
!     $             ibeg,iodd(id),lodd(id),ieven(id),leven(id)
!     write(*,*)'iodd',(is(i),i=iodd(id),lodd(id)+iodd(id)-1)
!     write(*,*)'ieven',(is(i),i=ieven(id),leven(id)+ieven(id)-1)
!         else
      enddo
!     Now iodd(id),lodd(id),ieven(id),leven(id) 
!     are the even and odd vectors face normal to dimension nn.

!     write(*,*)'Bbdyfacecreate end: ibeg,iodd,lodd,ieven,leven',
!     $       nn,iodd,lodd,ieven,leven
!     write(*,*)'iodd',(is(i),i=iodd(id),lodd(id)+iodd(id)-1)
!     write(*,*)'ieven',(is(i),i=ieven(id),leven(id)+ieven(id)-1)
      end

!****************************************************************
      subroutine bbdyorigprint(ndims,idims,iorig)
      integer ndims
      integer idims(ndims),iorig(*)
      integer kk(3)
      character*40 string
! This diagnostic works only for up to 3-d arrays.
         if(ndims.le.3) then
            do j=3,1,-1
               if(ndims.lt.j)then 
                  kk(j)=1
               else
                  kk(j)=idims(j)+1
               endif
            enddo
            write(*,*)'j,k,Block origins(i)='
            write(string,'(''('',i4''i10)'')')kk(1)+2
            write(*,string)
     $           ((j,k,(iorig(i+kk(1)*((j-1)+kk(2)*(k-1))),
     $           i=1,kk(1)),j=1,kk(2)),k=1,kk(3))
         endif

         end
!*******************************************************************
! Return the block origins for a multidimensional block arrangement.
      subroutine bbdydefine(ndims,idims,ifull,iuds,iorig,iLs)
! ndims: number of dimensions, 
! idims(ndims) number of blocks in each dimension, 
! ifull(ndims) full lengths of declared array in each dimension. 
! iuds(ndims) used lengths in each dimension
! iorig: return block origin displacements. (OUT)
! iLs(ndims+1) full-length step structure of declared array (OUT)
      integer ndims,idims(ndims),ifull(ndims),iuds(ndims)
      integer iLs(ndims+1)
! Presumed effective dimensions of iorig are 
!    (idims(1)+1,idims(2)+1,...)
      integer iorig(*)
      character*20 string
! Define the block origin addresses
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
!            write(*,*)'iinc,ifn,istep,isz=',iinc,ifn,istep,isz
            call bbdycatstep(iorig,ibeg,1,ilen,iinc)
         enddo
         ifn=ifull(n)*ifn
         iLs(n+1)=ifn
      enddo
!      write(*,*)'iorig='
      write(string,'(''('',i5,''i8)'')')(idims(1)+1)
!      write(*,string)(iorig(jj),jj=1,ibeg-1)
      end
!********************************************************************
      subroutine bbdyblockcreate(ndims,ktype,iLs,iside,iSIZEOFREAL)
! For a system of ndims dimensions, (IN)
! Create MPI type block handles in the array ktype (OUT)
! when the array structures are iLs (IN)
! the lengths of the block sides in dimension id are iside(2,id) (IN)
! with iside(1,id) referring to the bulk and iside(2,id) to the top.
      
!      parameter(imds=10)
      integer iLs(ndims+1)
! iside(1,*) is the general value, iside(2,*) the uppermost value.
      integer iside(2,ndims)
      
! ktype stores the type handles in the order
! 0-d; 1d bulk, 1d top; 2d bulk (1b,1t), 2d top (1b,1t), ...
! So for the nn-th level the start of the types is at 2**nn.
! ktype must have length (at least) 2**(ndims+1)-1
      integer ktype(*)
      include 'mpif.h'

      call MPI_TYPE_SIZE(MPI_REAL,iSIZEOFREAL,ierr)
!      write(*,*)'iSIZEOFREAL=',iSIZEOFREAL
! zeroth dimension type is MPI_REAL
      ktype(1)=MPI_REAL
      inew=2
      do nn=1,ndims
         iprior=2**(nn-1)
!         write(*,*)'bbdyblockcreate: nn,iprior,inew,ktype(iprior)'
!     $        ,nn,iprior,inew,ktype(iprior)
! laminate to higher dimension
         call bbdyblam(ndims,ktype,inew,iprior,nn,iside,
     $        iLs,iSIZEOFREAL)
!         write(*,*)'nn,ithis,inew,ktype(ithis...)',
!     $        nn,2**nn,inew,(ktype(ith),ith=(2**nn),2**(nn+1)-1)
      enddo
! The types we use now start at ktype(2**ndims), 2**ndims of them.
      end

!*********************************************************************
      subroutine bbdyblam(ndims,ktype,inew,iprior,nn,iside,
     $     iLs,iSIZEOFREAL)
!      parameter(imds=10)
! iside(1,*) is the general value, iside(2,*) the uppermost value.
      integer iside(2,ndims)
      integer iLs(ndims+1)
      integer ktype(*)
      include 'mpif.h'
! This was the cause of a lot of heartache. 
! Necessary for MPI_TYPE_CREATE_HVECTOR, wrong for TYPE_HVECTOR
!      integer(KIND=MPI_ADDRESS_KIND) istride
! The alternative for non F90 capable is integer*8 on 64 bit
! or integer*4 on 32 bit machines.
      ilen=1
      istride=iLs(nn)*iSIZEOFREAL
! For the bulk and top cases of this dimension
      do ibt=1,2
! Because we do not count the boundaries, 
! the count is 2 less than iside.
         icount=iside(ibt,nn)-2
! For all types in prior level, laminate enough together
! At each level there are 2**(nn-1) types in prior level
         do iold=iprior,iprior+2**(nn-1)-1
!            write(*,*)'iold,icount,istride,ktype(iold),inew',
!     $           iold,icount,istride,ktype(iold),inew
! This call needs istride to be MPI_ADDRESS_KIND
!            call MPI_TYPE_CREATE_HVECTOR(icount,ilen,istride,
! This call assumes istride is integer*4 (I think):
            call MPI_TYPE_HVECTOR(icount,ilen,istride,
     $           ktype(iold),ktype(inew),ierr)
! We only commit the top level that we are going to use.
            if(nn.eq.ndims)call MPI_TYPE_COMMIT(ktype(inew),ierr)
            inew=inew+1
         enddo
      enddo

      end
!*********************************************************************
      subroutine bbdycoords(nn,ndims,idims,icoords,ith,iLcoords,ionp)
! Obtain the cartesian coordinates of block for process nn,
! in ndims dimensions, whose lengths are idims(ndims).
! Return it in icoords. 
! Return the type index (zero-based), i.e.
! whether we are in bulk or on boundary in ith.
! Also the iorig index, ionp, for this block.
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
!         write(*,*)'nd,icoords(nd),in,iquot',nd,icoords(nd),in,iquot
         in=iquot
! Since iorig has dimensions 1+idims this is needed:
         ionp=ionp+icoords(nd)*iLcoords(nd)
      enddo
!      write(*,*)'nn,ith,ndims,idims',nn,ith,ndims,idims
      end
!****************************************************************
      subroutine bbdygatherparam(ndims,iLs,iside,
     $     mycartid,nproc,idims,iLcoords,ioffset,idebug,
     $     isdispls,istypes,iscounts,irdispls,irtypes,ircounts)
! For ndims dimensions, with array structure iLs, and block lengths
! given by iside, set the send and receive alltoall(w) parameters:
! isdispls, istypes,iscounts, irdispls, irtypes, ircounts,
! for process mycartid, within total of nproc, cartesian communicator
! dimensions idims, iLcoords, ioffset, idebug.

      integer ndims,mycartid,ioffset
      integer iLs(ndims+1),iLcoords(ndims),idims(ndims)
      integer iside(2,ndims)
      integer isdispls(nproc),istypes(nproc),iscounts(nproc)
      integer irdispls(nproc),irtypes(nproc),ircounts(nproc)

      integer norigmax
      parameter (norigmax=4096)
      integer iorig(norigmax)
      common /iorigcom/iorig

      parameter (imds=10)
      integer ktype(2**(imds+1)),iconp(imds)
      logical lalltoall
      data lalltoall/.true./

!      include 'mpif.h'

      call bbdyblockcreate( ndims,ktype,iLs,iside,iSIZEOFREAL)
      ith0=2**ndims
      if(idebug.ge.1) write(*,*)'Block types:',(ktype(ith),
     $     ith=ith0,2*ith0-1)

! Create the required arrays for the ALLTOALL
! First find out where we are and return the type index in ithi.
      call bbdycoords(mycartid+1,ndims,idims,iconp,ithi,
     $     iLcoords,ionpmine)
      do np=1,nproc
! Here we need to send only the active length of data, so the origin
! is offset.
         isdispls(np)=(iorig(ionpmine)-1+ioffset)*isizeofreal
         istypes(np)=ktype(ith0+ithi)
! Get positions and type indexes for other processes
         call bbdycoords(np,ndims,idims,iconp,ithj,
     $        iLcoords,ionp)
         irdispls(np)=(iorig(ionp)-1+ioffset)*isizeofreal
         irtypes(np)=ktype(ith0+ithj)
         if(lalltoall) then
! All to all
            if(np.eq.mycartid+1)then
! Don't send to or receive from myself.
               iscounts(np)=0
               ircounts(np)=0
            else
               iscounts(np)=1
               ircounts(np)=1
            endif
         else
! Gather to process 0 (fortran index 1)
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
! Don't send to or receive from anyone test.
! Process 1 sends to process 0 test
!         if(mycartid.eq.0) ircounts(2)=1
!         if(mycartid.eq.1) iscounts(1)=1
! Process 2 sends to process 0 test
!         if(mycartid.eq.0) ircounts(3)=1
!         if(mycartid.eq.2) iscounts(1)=1
      end
!********************************************************************
! Abstraction to isolate mpi calls.
      subroutine mpifinalize(ierr)
      call MPI_FINALIZE(ierr)
      end
!*******************************************************************
      subroutine mpicommsize(nprcsses,ierr)
      include 'mpif.h'
      call MPI_COMM_SIZE( MPI_COMM_WORLD, nprcsses, ierr )
      end
!********************************************************************
      subroutine mpigetmyid(myid,nprcsses,ierr)
! If necessary initialize the MPI system.
! Get my MPI id, and the number of processors.
      include 'mpif.h'
      logical lflag
      call MPI_INITIALIZED(lflag,ierr)
      if(.not.lflag) call MPI_INIT(ierr)
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, nprcsses, ierr )
      end
!********************************************************************
      subroutine mpiconvgreduce(convgd,icommcart,ierr)
      implicit none
      include 'mpif.h'
      real convgd(3),convgr(3)
      integer i,ierr,icommcart
      
! doing it in place:
!      write(*,*)'convgd,icommcart',convgd,icommcart
!........
! This is the MPI-2 version that is most convenient.
!      call MPI_ALLREDUCE(MPI_IN_PLACE,convgd,3,MPI_REAL,MPI_MAX,
!     $     icommcart,ierr)
! Since we use implicit none, if you try to use this version, then
! one ought to get an error if MPI_IN_PLACE is not in the headers.
! This should protect against MPI-1 problems.
! This is the version that one must use when MPI_IN_PLACE is absent.
      call MPI_ALLREDUCE(convgd,convgr,3,MPI_REAL,MPI_MAX,
     $     icommcart,ierr)
      do i=1,3
         convgd(i)=convgr(i)
      enddo
      end
!*********************************************************************
      subroutine mpibarrier(ierr)
      include 'mpif.h'
      integer ierr
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      end
!********************************************************************
! Broadcast the mesh quantity from the root.
      subroutine meshbroadcast(quant,ndims,ifull,iuds,iLs,iroot
     $     ,icommcart)
      implicit none
      integer ndims,ifull(ndims),iuds(ndims),iLs(ndims+1),iroot,ierr
     $     ,icommcart
      real quant(*)
      include 'mpif.h'
! Operator not really needed:
      external addarray_MPI
      logical lfirst,lneeded
      integer iaddtype,iaddop,iporig,nprcsses,ncartsize
      data lfirst/.true./lneeded/.true./
      save lfirst,iaddtype,iaddop,iporig,lneeded

      if(.not.lneeded)return
      
      if(lfirst)then
! Determine if this is all needed
         call MPI_COMM_SIZE(MPI_COMM_WORLD,nprcsses,ierr)
         if(icommcart.eq.MPI_COMM_NULL)then
            ncartsize=0
         else
            call MPI_COMM_SIZE(icommcart,ncartsize,ierr)
         endif
         if(nprcsses.gt.ncartsize)then
!         write(*,*)'Calling mpiopcreate',ndiags
! Create addtype and operator for reduce sum; full array.
            call mpiopcreate(ndims,ifull,iuds,addarray_MPI,
     $        iaddtype,iaddop)
!         write(*,*)'Returned from mpiopcreate',iaddtype
            lfirst=.false.
         else
            lneeded=.false.
            return
         endif
      endif
      iporig=1
!      write(*,*)'Bcast',icommcart
      call MPI_Bcast(quant(iporig),1,iaddtype,iroot,MPI_COMM_WORLD,ierr)
      end
!********************************************************************
