c***********************************************************************
c      subroutine mditerate(ndims,ifull,iused,routine,u,ipin)
c New argument order the same as mditerarg.
      subroutine mditerate(routine,ndims,ifull,iused,ipin,u)
c multidimensional iteration. For dimensions ndims, iterate over the
c array whose full dimensions are ifull(ndims) used iused(ndims).
c At each iteration, 
c         routine(inc,ipoint,indi,ndims,iused,u) is called,
c which accepts a pointer to the address in the array structure: ipoint
c referred to the full dimensions of u, so u(1+ipoint) is processed.  It
c returns the next increment in units of the lowest dimension: inc. The
c offsets for each dimension are passed in: indi(ndims), which have
c values in the range (0,iused(ndims)-1). They can (but need not) be
c adjusted in routine. If they are, ipoint must also be adjusted
c equivalently in routine. Again: in routine, ipoint addresses u(ifull).
c After the return of routine, ipoint and indi(1) are incremented by inc
c and wrapping is handled. Wrapping occurs at the _iused_ dimension
c length relative to the starting indices.  In other words, if the
c incremented indi(id)-indinp(id)+1 exceeds the used length iused(id),
c then the next higher dimension is incremented by 1, and indi(id)
c decremented by iused(id). Thus inc is the increment in the lowest
c dimension index, referenced to the _used_ dimensions, but ipoint is
c wrapped accounting for the knowledge of ifull so it points to the
c correct next position within the array. u is the data that is passed
c to the routine.
c
c On entry, ipin is the starting offset storage position from the start
c of u relative ifull array structure. The indinp(id) are the
c corresponding offset-indices of that position and are calculated
c internally.  Iteration does not extend past the used bounds of u
c (iused) plus the starting offset (indinp).  A shifted subarray can be
c treated by specifying ipoint=0, giving the address of the subarray as
c u and specifying in iused the lengths of the subarray used
c dimensions. Alternatively, it can be treated by passing the whole
c array and the appropriate pointer offset in ipin, and again the used
c subarray dimensions.
c
c Normal usage is to pass ipoint=0 to mditerate. Then the routine is called
c for the whole _used_ array at nodes spaced by inc in the 1st dimension.

      integer ndims
      integer ifull(ndims)
      integer iused(ndims)
      external routine
      real u(*)

      integer mdims
      parameter (mdims=10)
c Effective index in dimension, c-style (zero based)
      integer indi(mdims),indinp(mdims)
c Structure vector
      integer iLs(mdims+1)
      common /iLscom/iLs
      
c      write(*,*)'ndims',ndims,' ifull',ifull,' iused',iused

c Offset.
      call offsetexpand(ndims,ifull,ipin,indi)
      iLs(1)=1
      do n=1,ndims
         indinp(n)=indi(n)
         iLs(n+1)=iLs(n)*ifull(n)
         if(indinp(n)+iused(n).gt.ifull(n))then
            write(*,*) 'MDITERATE Error. Dimension',n,' ioff + iused',
     $           indinp(n),iused(n),' .gt. ifull',ifull(n)
         endif
      enddo
      ipoint=ipin
c      if(ipoint.ne.0)then
c         write(*,*)'MDITERATE', ipoint,
c     $        (indi(k),iused(k),ifull(k),k=1,ndims)
c      endif
      inc=1
c      icount=0
      n=1
c Iteration over the multidimensional array
 101  continue
c      write(*,'(''('',i1,i4,'') '',$)')n,indi(n)
      if(indi(n)-indinp(n).gt.iused(n)-1)then
c     Overflow. Subtract off enough (inm) of next dimension
c     and move ipoint to appropriate position in full array.
         inm=0
 102     inm=inm+1
         ipoint=ipoint+iLs(n+1)-iused(n)*iLs(n)
         indi(n)=indi(n)-iused(n)
         if(indi(n)-indinp(n).gt.iused(n)-1)goto 102
c Increment the next level.
         n=n+1
         if(n.gt.ndims)goto 201
         indi(n)=indi(n)+inm
         goto 101
      elseif(n.gt.1)then
c We've carried and cleared an increment.
c Return stepwise to base level
         n=n-1
         goto 101
      else
c We're at the base level and have succeeded in incrementing.
c Do whatever we need to and increment indi(1) and ipoint
c         write(*,'(a,8i8)')'ipoint=',ipoint,(indi(i),i=1,ndims)
         call routine(inc,ipoint,indi,ndims,iused,u)
         indi(1)=indi(1)+inc
         ipoint=ipoint+inc
         goto 101
      endif

 201  continue
c Reached the end.

      end
c******************************************************************
      subroutine rbroutine(inc,ipoint,indi,ndims,iused,u)
c Red/black routine, process squares alternately.
      integer ipoint,inc
      integer indi(ndims),iused(ndims)
      real u(*)

      parameter (mdims=10)
      integer ind1(mdims)
      data icount/0/
      data ind1/mdims*1000/

c iused is not used. To silence warnings spuriously use it:
      ic=iused(1)
c      write(*,'(a,11i8)')'ipoint,indi',ipoint,(indi(i),i=1,ndims)
c      write(*,*)'myipoint=',ipoint
c We build in correction of the increment here for red-black
c to change the parity. We have to remember the previous indis.
      ic=0
      do i=ndims,2,-1
         if(indi(i).gt.ind1(i))then
c we've carried
            ic=ic+(i-1)
c            write(*,*)'increment',ic
         endif
c Remember prior.
         ind1(i)=indi(i)
      enddo
      if(mod(ic,2).ne.0) then
         iaj=(1-2*mod(indi(1),2))
         indi(1)=indi(1)+iaj
         ipoint=ipoint+iaj
      endif

c Here is where we process the square chosen. Examples:
c sequence
c      u(ipoint+1)=mod(icount,9)+1
c height
      u(ipoint+1)=indi(3)+1
      icount=icount+1
      inc=2
      end
c************************************************************************
c Shell for skipping over the boundary of ndims-dimensional volume.
      subroutine bdyroutine(inc,ipoint,indi,ndims,iused,u)
      integer ipoint,inc
      integer indi(ndims),iused(ndims)
      real u(*)

      parameter (mdims=10)
c Structure vector needed for finding adjacent u values.
      integer iLs(mdims+1)
      common /iLscom/iLs

c Algorithm: if on a boundary face of dimension >1, steps of 1 (dim 1).
c Otherwise steps of iused(1)-1 or 1 on first or last (of dim 1).
      inc=1
      do n=ndims,2,-1
         if(indi(n).eq.0)then
c On boundary face 0 of dimension n>1. Break.
c This is where we put boundary setting for n>1
c Example:
            u(ipoint+1)=0.
c            if(n.eq.3)then
c Second derivative is zero:
c               u(ipoint+1)=2.*u(ipoint+1+iLs(n))-u(ipoint+1+2*iLs(n)
c First derivative is zero:
c               u(ipoint+1)=u(ipoint+1+iLs(n))
c            endif
            goto 101
         elseif(indi(n).eq.iused(n)-1)then
c On boundary face iused(n) of dimension n>1. Break.
c Example:
            u(ipoint+1)=0.
c 
c            if(n.eq.3) u(ipoint+1)=u(ipoint+1-iLs(n))
            goto 101
         endif
      enddo
c     We are not on any higher boundary.
c This is where the boundary setting is done for n=1
c Example:
      u(ipoint+1)=0.
      if(indi(n).eq.0)then
         inc=iused(1)-1
      elseif(indi(n).eq.iused(n)-1)then
         inc=1
      else
         write(*,*)'BDY Error. We should not be here',
     $        n,ipoint,indi
         stop
      endif
 101  continue
c      write(*,*)'indi,inc,iused,ipoint',indi,inc,iused,ipoint

      end
c************************************************************************
c***********************************************************************
      subroutine mditerarg(routine,ndims,ifull,iused,ipin,t,u,v,w)
c multidimensional iteration. For dimensions ndims, iterate over the
c array whose full dimensions are ifull(ndims) used iused(ndims).
c This version allows four (array) arguments: t,u,v,w to be passed.
c At each iteration, the routine (which should be declared external)
c
c         routine(inc,ipoint,indi,ndims,iused,t,u,v,w) is called,
c whose arguments are
c      integer ipin,inc
c      integer indi(ndims),iused(ndims)
c      real t(*),u(*),v(*),w(*)
c
c which accepts a pointer to the address in the array structure: ipoint
c referred to the full dimensions of u, so u(1+ipoint) is processed.  It
c returns the next increment in units of the lowest dimension: inc. The
c offsets for each dimension are passed in: indi(ndims) And can (but
c need not) be adjusted in routine. If they are, ipoint must also be
c adjusted equivalently in routine. Again: in routine, ipoint addresses
c u(ifull).  After the return of routine, ipoint and indi(1) are
c incremented by inc and wrapping is handled.
c
c Wrapping occurs at the
c _iused_ dimension length. In other words, if the incremented indi(id)
c passes the used length iused(id), then the next higher dimension is
c incremented by 1, and indi(id) decremented by iused(id). Thus inc is
c the increment in the lowest dimension index, referenced to the _used_
c dimensions, but ipoint is wrapped accounting for the knowledge of
c ifull so it points to the correct next position within the array. u is
c the data that is passed to the routine. On entry, ipin is the
c starting offset storage position from the start of u, but in this case
c (only), relative to the iused array structure. An error will occur if
c ipoint.gt.iused(1) on entry.  Iteration does not extend past the used
c bounds of u.  Therefore a shifted subarray should be treated by
c specifying ipoint=0 and giving the address of the subarray as u.
c
c Normal usage is to pass ipoint=0 to mditerate. Then the routine is called
c for the whole _used_ array at nodes spaced by inc in the 1st dimension.

      integer ndims
      integer ifull(ndims)
      integer iused(ndims)
      external routine
      real t(*),u(*),v(*),w(*)

      integer mdims
      parameter (mdims=10)
c Effective index in dimension, c-style (zero based)
      integer indi(mdims),indinp(mdims)
c Structure vector
      integer iLs(mdims+1)
      common /iLscom/iLs
      
c      write(*,*)'ndims',ndims,' ifull',ifull,' iused',iused
c Offset.
      call offsetexpand(ndims,ifull,ipin,indi)
      iLs(1)=1
      do n=1,ndims
         indinp(n)=indi(n)
         iLs(n+1)=iLs(n)*ifull(n)
         if(indinp(n)+iused(n).gt.ifull(n))then
            write(*,*) 'MDITERATE Error. Dimension',n,' ioff + iused',
     $           indinp(n),iused(n),' .gt. ifull',ifull(n)
         endif
      enddo
      ipoint=ipin

      inc=1
c      icount=0
      n=1
c Iteration over the multidimensional array
 101  continue
c      write(*,'(''('',i1,i4,'') '',$)')n,indi(n)
      if(indi(n)-indinp(n).gt.iused(n)-1)then
c     Overflow. Subtract off enough (inm) of next dimension
c     and move ipoint to appropriate position in full array.
         inm=0
 102     inm=inm+1
         ipoint=ipoint+iLs(n+1)-iused(n)*iLs(n)
         indi(n)=indi(n)-iused(n)
         if(indi(n)-indinp(n).gt.iused(n)-1)goto 102
c Increment the next level.
         n=n+1
         if(n.gt.ndims)goto 201
         indi(n)=indi(n)+inm
         goto 101
      elseif(n.gt.1)then
c We've carried and cleared an increment.
c Return stepwise to base level
         n=n-1
         goto 101
      else
c We're at the base level and have succeeded in incrementing.
c Do whatever we need to and increment indi(1) and ipoint
c         write(*,'(a,8i8)')'ipoint=',ipoint,(indi(i),i=1,ndims)
         call routine(inc,ipoint,indi,ndims,iused,t,u,v,w)
         indi(1)=indi(1)+inc
         ipoint=ipoint+inc
         goto 101
      endif

 201  continue
c Reached the end.

      end
c******************************************************************
      subroutine mditerset(u,ndims,ifull,iused,ipin,v)
c Iterate over the array u, setting its value equal to (scalar) v.
c Uses the mechanisms of mditerate. But inc is always 1.
      integer ndims
      integer ifull(ndims)
      integer iused(ndims)
      real u(*),v

      integer mdims
      parameter (mdims=10)
c Effective index in dimension, c-style (zero based)
      integer indi(mdims),indinp(mdims)
c Structure vector
      integer iLs(mdims+1)
      common /iLscom/iLs
      
c      write(*,*)'ndims',ndims,' ifull',ifull,' iused',iused
c Offset.
      call offsetexpand(ndims,ifull,ipin,indi)
      iLs(1)=1
      do n=1,ndims
         indinp(n)=indi(n)
         iLs(n+1)=iLs(n)*ifull(n)
         if(indinp(n)+iused(n).gt.ifull(n))then
            write(*,*) 'MDITERSET Error. Dimension',n,' ioff + iused',
     $           indinp(n),iused(n),' .gt. ifull',ifull(n)
         endif
      enddo
      ipoint=ipin

      inc=1
      n=1
c Iteration over the multidimensional array
 101  continue
c      write(*,'(''('',i1,i4,'') '',$)')n,indi(n)
      if(indi(n)-indinp(n).gt.iused(n)-1)then
c     Overflow. Subtract off enough (inm) of next dimension
c     and move ipoint to appropriate position in full array.
         inm=0
 102     inm=inm+1
         ipoint=ipoint+iLs(n+1)-iused(n)*iLs(n)
         indi(n)=indi(n)-iused(n)
         if(indi(n)-indinp(n).gt.iused(n)-1)goto 102
c Increment the next level.
         n=n+1
         if(n.gt.ndims)goto 201
         indi(n)=indi(n)+inm
         goto 101
      elseif(n.gt.1)then
c We've carried and cleared an increment.
c Return stepwise to base level
         n=n-1
         goto 101
      else
c We're at the base level and have succeeded in incrementing.
c Do whatever we need to and increment indi(1) and ipoint
         u(ipoint+1)=v
         indi(1)=indi(1)+inc
         ipoint=ipoint+inc
         goto 101
      endif

 201  continue
c Reached the end.

      end
c******************************************************************
      subroutine offsetexpand(ndims,ifull,index,indi)
c On entry index is the zero-based pointer to the position in the 
c array whose full dimensions are ifull(ndims). 
c On exit indi contains the corresponding (ndims)-dimensional offsets.
      integer ndims
      integer ifull(ndims),indi(ndims)
      integer index

      ind=index
      do i=1,ndims-1
         ind2=ind/ifull(i)
         indi(i)=ind-ind2*ifull(i)
         ind=ind2
      enddo
      indi(ndims)=ind
      if(ind.gt.ifull(3)) write(*,*)'indexexpand index too big',index
      end
c******************************************************************
      subroutine mditeradd(u,ndims,ifull,iused,ipin,v)
c Iterate over the array u, adding v(array) to it.
c Uses the mechanisms of mditerate. But inc is always 1.
c Normally, the starting pointer is ipin=0 for the full array.
      integer ndims
      integer ifull(ndims)
      integer iused(ndims)
      real u(*),v(*)

      integer mdims
      parameter (mdims=10)
c Effective index in dimension, c-style (zero based)
      integer indi(mdims),indinp(mdims)
c Structure vector
      integer iLs(mdims+1)
      common /iLscom/iLs
      
      call offsetexpand(ndims,ifull,ipin,indi)
      iLs(1)=1
      do n=1,ndims
         indinp(n)=indi(n)
         iLs(n+1)=iLs(n)*ifull(n)
         if(indinp(n)+iused(n).gt.ifull(n))then
            write(*,*) 'MDITERADD Error. Dimension',n,' ioff + iused',
     $           indinp(n),iused(n),' .gt. ifull',ifull(n)
         endif
      enddo
      ipoint=ipin

      inc=1
      n=1
c Iteration over the multidimensional array
 101  continue
c      write(*,'(''('',i1,i4,'') '',$)')n,indi(n)
      if(indi(n)-indinp(n).gt.iused(n)-1)then
c     Overflow. Subtract off enough (inm) of next dimension
c     and move ipoint to appropriate position in full array.
         inm=0
 102     inm=inm+1
         ipoint=ipoint+iLs(n+1)-iused(n)*iLs(n)
         indi(n)=indi(n)-iused(n)
         if(indi(n)-indinp(n).gt.iused(n)-1)goto 102
c Increment the next level.
         n=n+1
         if(n.gt.ndims)goto 201
         indi(n)=indi(n)+inm
         goto 101
      elseif(n.gt.1)then
c We've carried and cleared an increment.
c Return stepwise to base level
         n=n-1
         goto 101
      else
c We're at the base level and have succeeded in incrementing.
c Do whatever we need to and increment indi(1) and ipoint
         u(ipoint+1)=u(ipoint+1)+v(ipoint+1)
         indi(1)=indi(1)+inc
         ipoint=ipoint+inc
         goto 101
      endif

 201  continue
c Reached the end.
      end
c******************************************************************
      
