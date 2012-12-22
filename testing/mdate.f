c Obsolete routine replaced by mditerarg.
c***********************************************************************
      subroutine mditerate(routine,ndims,ifull,iused,ipin,u)
c multidimensional iteration. For dimensions ndims, iterate over the
c array whose full dimensions are ifull(ndims) used iused(ndims).
c At each iteration, 
c         routine(inc,ipoint,indi,ndims,iLs,iused,u) is called,
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
c      common /iLscom/iLs
      
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
         call routine(inc,ipoint,indi,ndims,iLs,iused,u)
         indi(1)=indi(1)+inc
         ipoint=ipoint+inc
         goto 101
      endif

 201  continue
c Reached the end.

      end
