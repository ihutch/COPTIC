c************************************************************************
c Shell for skipping over the boundary of ndims-dimensional volume.
      subroutine bdyroutine(inc,ipoint,indi,ndims,iused,u)
      integer ipoint,inc
      integer indi(ndims),iused(ndims)
      real u(*)

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
c***********************************************************************
      subroutine mditerarg(routine,ndims,ifull,iused,ipin,t,u,v,w)
c multidimensional iteration. For dimensions ndims, iterate over the
c array whose full dimensions are ifull(ndims) used iused(ndims).
c This version allows four (array) arguments: t,u,v,w to be passed.
c At each iteration, the routine (which should be declared external)
c
c         routine(inc,ipoint,indi,ndims,iLs,iused,t,u,v,w) is called,
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
      parameter (mdims=3)
c Effective index in dimension, c-style (zero based)
      integer indi(mdims),indinp(mdims)
c Structure vector
      integer iLs(mdims+1)
      
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
c        write(*,'(3i3,$)')(iused(ii),ii=1,3)
c        write(*,'(3i3,i7,i5)')(indi(i),i=1,ndims),ipoint,inc
         call routine(inc,ipoint,indi,ndims,iLs,iused,t,u,v,w)
         indi(1)=indi(1)+inc
         ipoint=ipoint+inc
         goto 101
      endif

 201  continue
c Reached the end.

      end
c*****************************************************************
c Use the generalized scalar mult call to set values
      subroutine mditerset(u,ndims,ifull,iused,ipin,v)
      call mditermults(u,ndims,ifull,iused,ipin,0.,v)
      end
c******************************************************************
c*****************************************************************
      integer function indexcontract(ndims,ifull,indi)
c On entry indi(ndims) contains the (ndims)-dimensional offsets
c to position in array whose full dimensions are ifull(ndims)
c Returns the corresponding linear offset index within that array.
      integer ndims
      integer ifull(ndims),indi(ndims)
      indexcontract=indi(ndims)
      do id=ndims-1,1,-1
         indexcontract=indi(id)+indexcontract*ifull(id)
      enddo
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
      if(ind.gt.ifull(ndims)) write(*,*)'offsetexpand index too big'
     $     ,index
      end
c******************************************************************
c New code for the rationalized mditerator.
c******************************************************************
      integer function mditerator(ndims,iview,indi,isw,iaux)
c Multidimensional iterator for ndims dimensions over an abstract
c view iview of an array. Zero-based offsets.
c On entry:
c    iview(3,ndims) contains the view of the array in the form
c       istart1,iend1,istride1;istart2,iend2,...,istridendims.
c    iaux(ndims) maybe contains extra information.
c    isw =0 iterate, else setup iview
c         1 set iview(istart,j)=iaux(j)
c         2 set iview(iend,j)=iaux(j)-1
c         3 set iview(istride,j)=iaux(j)
c         4 set iview(*,j) = 0,iaux(j)-1,1
c    indi(ndims) is the current offset-array (set to 0 if isw=4)
c On exit: indi is the new index-array. Optionally iview is [re]set.
c    Return value is 0 if iteration incomplete, else positive.
c Example of typical usage:
c      icomplete=mditerator(ndims,iview,indi,4,iused)
c 1    u(1+indexcontract(ndims,ifull,indi))=mod(1+indi(3),10)
c      if(mditerator(ndims,iview,indi,0,iused).eq.0)goto 1
c      
      integer ndims,iview(3,ndims),indi(ndims),isw,iaux(ndims)
      integer istart,iend,istride
      parameter (istart=1,iend=2,istride=3)

      mditerator=isw
      if(isw.eq.0)then
c-----------------------------------------------------
c Iterate
c Start at dimension 1.
         n=1
c Increment:
         indi(n)=indi(n)+iview(istride,n)
 101     continue
c      write(*,'(''('',i1,i4,'') '',$)')n,indi(n)
         if((indi(n)-iview(iend,n))*iview(istride,n).gt.0)then
c     Overflow. Subtract off enough (inm) of next dimension
            inm=0
 102        inm=inm+1
            indi(n)=indi(n)-(iview(iend,n)-iview(istart,n)+1)
c            write(*,*)inm,indi(n),iview(3,n),iview(2,n)
            if((indi(n)-iview(iend,n))*iview(istride,n).gt.0)goto 102
c Move to the next dimension:
            n=n+1
c Overflow of ndims dimension completes iteration:
            if(n.gt.ndims)then 
               mditerator=1
               return
            endif
c Increment this dimension with carry and restart:
            indi(n)=indi(n)+inm*iview(istride,n)
            goto 101
         elseif(n.gt.1)then
c We've carried and cleared an increment.
c Return stepwise to base level
            n=n-1
            goto 101
         endif
c Here when we have successfully incremented all levels.
c-------------------------------------------------------
      elseif(isw.eq.4)then
c Total setup
         do i=1,ndims
            iview(istart,i)=0
            iview(iend,i)=iaux(i)-1
            iview(istride,i)=1
            indi(i)=0
         enddo
      else
c Individual setup.
         do i=1,ndims
            iview(isw,i)=iaux(i)
            if(isw.eq.2)iview(isw,i)=iview(isw,i)-1
         enddo
      endif
      end
c*****************************************************************
      subroutine mditermults(u,ndims,ifull,iused,ipin,v,w)
c Iterate over the array u, multiplying it by v(scalar)
c and adding w(scalar): u = u*v+w.
c Reimplementation using mditerator.
c Normally, the starting pointer is ipin=0 for the full array.
      integer ndims,ifull(ndims),iused(ndims)
      real u(*),v,w
      integer mdims
      parameter (mdims=3)
c Effective index in dimension, c-style (zero based)
      integer indi(mdims),iview(3,mdims)
c This data statement serves to silence ftnchek. The first mditerator
c call actually initializes the iview and indi.
      data iview/mdims*0,mdims*0,mdims*0/
      data indi/mdims*0/
      
c Set to iterate the whole used array:
      icomplete=mditerator(ndims,iview,indi,4,iused)
c Reset starting indices from input pointer, in case non-zero.
      if(ipin.ne.0)call offsetexpand(ndims,ifull,ipin,indi)
 1      ii=1+indexcontract(ndims,ifull,indi)
        u(ii)=u(ii)*v+w
      if(mditerator(ndims,iview,indi,0,iused).eq.0)goto 1

      end
c*****************************************************************
      subroutine mditeradd(u,ndims,ifull,iused,ipin,v)
c Iterate over the array u, adding v(array) to it: u(i)=u(i)+v(i)
c Reimplementation using mditerator.
c Normally, the starting pointer is ipin=0 for the full array.
      integer ndims,ifull(ndims),iused(ndims)
      real u(*),v(*)
      integer mdims
      parameter (mdims=3)
c Effective index in dimension, c-style (zero based)
      integer indi(mdims),iview(3,mdims)
c This data statement serves to silence ftnchek. The first mditerator
c call actually initializes the iview and indi.
      data iview/mdims*0,mdims*0,mdims*0/
      data indi/mdims*0/
      
c Set to iterate the whole used array:
      icomplete=mditerator(ndims,iview,indi,4,iused)
c Reset starting indices from input pointer, in case non-zero.
      if(ipin.ne.0)call offsetexpand(ndims,ifull,ipin,indi)
 1      ii=1+indexcontract(ndims,ifull,indi)
        u(ii)=u(ii)+v(ii)
      if(mditerator(ndims,iview,indi,0,iused).eq.0)goto 1

      end
c*****************************************************************
      subroutine mditerave(u,ndims,ifull,iused,ipin,uave,istepave)
c Iterate over the array u, to produce uave in the form
c ((istepave-1)*uave + u)/istepave

c Normally, the starting pointer is ipin=0 for the full array.
      integer ndims,ifull(ndims),iused(ndims),istepave
      real u(*),uave(*)
      integer mdims
      parameter (mdims=3)
c Effective index in dimension, c-style (zero based)
      integer indi(mdims),iview(3,mdims)
c This data statement serves to silence ftnchek. The first mditerator
c call actually initializes the iview and indi.
      data iview/mdims*0,mdims*0,mdims*0/
      data indi/mdims*0/
      
      if(mdims.lt.ndims)then
         write(*,*)'Dimension error in mditerave',mdims,ndims
         stop
      endif
c Set to iterate the whole used array:
      icomplete=mditerator(ndims,iview,indi,4,iused)
c Reset starting indices from input pointer, in case non-zero.
      if(ipin.ne.0)call offsetexpand(ndims,ifull,ipin,indi)
 1      ii=1+indexcontract(ndims,ifull,indi)
        uave(ii)=(uave(ii)*(istepave-1.)+u(ii))/float(istepave)
      if(mditerator(ndims,iview,indi,0,iused).eq.0)goto 1

      end

c************************************************************************
c************************************************************************
      subroutine mditerarg1(routine,ndims,ifull,iused,ipin,t,u,v,w)
c [Reimplemented using mditerator. Runs about 50% slower for fast routines]
c multidimensional iteration. For dimensions ndims, iterate over the
c array whose full dimensions are ifull(ndims) used iused(ndims).
c This version allows four (array) arguments: t,u,v,w to be passed.
c At each iteration, the routine (which should be declared external)
c
c         routine(inc,ipoint,indi,ndims,iLs,iused,t,u,v,w) is called,
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
c Wrapping occurs at the _iused_ dimension length relative to the
c starting indices.  In other words, if the incremented
c indi(id)-indinp(id)+1 exceeds the used length iused(id), then the next
c higher dimension is incremented by 1, and indi(id) decremented by
c iused(id). Thus inc is the increment in the lowest dimension index,
c referenced to the _used_ dimensions, but ipoint is wrapped accounting
c for the knowledge of ifull so it points to the correct next position
c within the array. u is the data that is passed to the routine. On
c entry, ipin is the starting offset storage position from the start of
c u, but in this case (only), relative to the iused array structure. An
c error will occur if ipoint.gt.iused(1) on entry.  Iteration does not
c extend past the used bounds of u.  Therefore a shifted subarray should
c be treated by specifying ipoint=0 and giving the address of the
c subarray as u.
c
c Normal usage is to pass ipoint=0 to mditerate. Then the routine is called
c for the whole _used_ array at nodes spaced by inc in the 1st dimension.
      implicit none
      integer ndims,ipin
      integer ifull(ndims)
      integer iused(ndims)
      external routine
      real t(*),u(*),v(*),w(*)

      integer mdims,i
      parameter (mdims=3)
c Effective index in dimension, c-style (zero based)
      integer indi(mdims),iview(3,mdims)
      integer ipoint,inc,icomplete,indexcontract,mditerator
c This data statement serves to silence ftnchek. The first mditerator
c call actually initializes the iview and indi.
      data iview/mdims*0,mdims*0,mdims*0/
      data indi/mdims*0/

c Structure vector
      integer iLs(mdims+1)

      inc=1
c Set to iterate the whole used array:
      icomplete=mditerator(ndims,iview,indi,4,iused)
c Reset starting indices from input pointer, in case non-zero.
      call offsetexpand(ndims,ifull,ipin,indi)
      iLs(1)=1
      do i=1,ndims
c Structure vector must be set.
         iLs(i+1)=iLs(i)*ifull(i)
         if(indi(i)+iused(i).gt.ifull(i))then
            write(*,*) 'MDITERATE Error. Dimension',i,' indi + iused',
     $           indi(i),iused(i),' .gt. ifull',ifull(i)
            stop
         endif
c Iteration from starting index
         iview(1,i)=indi(i)
c With length equal to iused.
         iview(2,i)=iused(i)+indi(i)-1
      enddo
c      write(*,*)iview
c ipoint update cannot be put in mditerator, because it does not have
c the ifull information.
 1      ipoint=indexcontract(ndims,ifull,indi)
c        write(*,'(3i3,i7,i5)')(indi(i),i=1,ndims),ipoint,inc
        call routine(inc,ipoint,indi,ndims,iLs,iused,t,u,v,w)
c        write(*,*)indi,ipoint,inc,t(1+ipoint)
c Set the iterator view according to the increment of the first
c dimension. 
        iview(3,1)=inc
      if(mditerator(ndims,iview,indi,0,iused).eq.0)goto 1

      end
c************************************************************************
      subroutine mditerarg2(routine,ndims,ifull,iused,ipin,t,u,v,w)
c [Reimplemented using mditerator. With inline coding replacement]
c This was to explore further the timing questions but not completed.
c multidimensional iteration. For dimensions ndims, iterate over the
c array whose full dimensions are ifull(ndims) used iused(ndims).
c This version allows four (array) arguments: t,u,v,w to be passed.
c At each iteration, the routine (which should be declared external)
c
c         routine(inc,ipoint,indi,ndims,iLs,iused,t,u,v,w) is called,
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
c Wrapping occurs at the _iused_ dimension length relative to the
c starting indices.  In other words, if the incremented
c indi(id)-indinp(id)+1 exceeds the used length iused(id), then the next
c higher dimension is incremented by 1, and indi(id) decremented by
c iused(id). Thus inc is the increment in the lowest dimension index,
c referenced to the _used_ dimensions, but ipoint is wrapped accounting
c for the knowledge of ifull so it points to the correct next position
c within the array. u is the data that is passed to the routine. On
c entry, ipin is the starting offset storage position from the start of
c u, but in this case (only), relative to the iused array structure. An
c error will occur if ipoint.gt.iused(1) on entry.  Iteration does not
c extend past the used bounds of u.  Therefore a shifted subarray should
c be treated by specifying ipoint=0 and giving the address of the
c subarray as u.
c
c Normal usage is to pass ipoint=0 to mditerate. Then the routine is called
c for the whole _used_ array at nodes spaced by inc in the 1st dimension.
      implicit none
      integer ndims,ipin
      integer ifull(ndims)
      integer iused(ndims)
      external routine
      real t(*),u(*),v(*),w(*)

      integer mdims,i
      parameter (mdims=3)
c Effective index in dimension, c-style (zero based)
      integer indi(mdims),iview(3,mdims)
      integer ipoint,inc,icomplete,indexcontract,mditerator
      integer istart,iend,istride
      parameter (istart=1,iend=2,istride=3)
      integer n,inm
c This data statement serves to silence ftnchek. The first mditerator
c call actually initializes the iview and indi.
      data iview/mdims*0,mdims*0,mdims*0/
      data indi/mdims*0/

c Structure vector
      integer iLs(mdims+1)

      inc=1
c Set to iterate the whole used array:
      icomplete=mditerator(ndims,iview,indi,4,iused)
c Reset starting indices from input pointer, in case non-zero.
      call offsetexpand(ndims,ifull,ipin,indi)
      iLs(1)=1
      do i=1,ndims
c Structure vector must be set.
         iLs(i+1)=iLs(i)*ifull(i)
         if(indi(i)+iused(i).gt.ifull(i))then
            write(*,*) 'MDITERATE Error. Dimension',i,' indi + iused',
     $           indi(i),iused(i),' .gt. ifull',ifull(i)
            stop
         endif
c Iteration from starting index
         iview(1,i)=indi(i)
c With length equal to iused.
         iview(2,i)=iused(i)+indi(i)-1
      enddo
c      write(*,*)iview
 1      ipoint=indexcontract(ndims,ifull,indi)
c        write(*,'(3i3,i7,i5)')(indi(i),i=1,ndims),ipoint,inc
        call routine(inc,ipoint,indi,ndims,iLs,iused,t,u,v,w)
c        write(*,*)indi,ipoint,inc,t(1+ipoint)
c Set the iterator view according to the increment of the first
c dimension. Moving this outside the loop disallows changes by routine.
        iview(istride,1)=inc
        if(.false.)then
c Inline version of mditerator:
c-----------------------------------------------------
c Iterate
c Start at dimension 1.
         n=1
c Increment:
         indi(n)=indi(n)+iview(istride,n)
 101     continue
c      write(*,'(''('',i1,i4,'') '',$)')n,indi(n)
         if((indi(n)-iview(iend,n))*iview(istride,n).gt.0)then
c     Overflow. Subtract off enough (inm) of next dimension
            inm=0
 102        inm=inm+1
            indi(n)=indi(n)-(iview(iend,n)-iview(istart,n)+1)
c            write(*,*)inm,indi(n),iview(3,n),iview(2,n)
            if((indi(n)-iview(iend,n))*iview(istride,n).gt.0)goto 102
c Move to the next dimension:
            n=n+1
c Overflow of ndims dimension completes iteration:
            if(n.gt.ndims)then 
               goto 1
            endif
c Increment this dimension with carry and restart:
            indi(n)=indi(n)+inm*iview(istride,n)
            goto 101
         elseif(n.gt.1)then
c We've carried and cleared an increment.
c Return stepwise to base level
            n=n-1
            goto 101
         endif
c Here when we have successfully incremented all levels.
c-------------------------------------------------------
         else
      if(mditerator(ndims,iview,indi,0,iused).eq.0)goto 1
         endif

      end
      
c******************************************************************
c******************************************************************
c   Obsolete below here.
c******************************************************************
c******************************************************************

      subroutine mditermults0(u,ndims,ifull,iused,ipin,v,w)
c Iterate over the array u, multiplying it by v(scalar)
c and adding w(scalar): u = u*v+w.
c Uses the mechanisms of mditerate. But inc is always 1.
c Normally, the starting pointer is ipin=0 for the full array.
      integer ndims
      integer ifull(ndims)
      integer iused(ndims)
      real u(*),v,w

      integer mdims
      parameter (mdims=3)
c Effective index in dimension, c-style (zero based)
      integer indi(mdims),indinp(mdims)
c Structure vector
      integer iLs(mdims+1)
      
      call offsetexpand(ndims,ifull,ipin,indi)
      iLs(1)=1
      do n=1,ndims
         indinp(n)=indi(n)
         iLs(n+1)=iLs(n)*ifull(n)
         if(indinp(n)+iused(n).gt.ifull(n))then
            write(*,*) 'MDITERMULTS Error. Dimension',n,' ioff + iused',
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
         u(ipoint+1)=u(ipoint+1)*v+w
         indi(1)=indi(1)+inc
         ipoint=ipoint+inc
         goto 101
      endif

 201  continue
c Reached the end.
      end
c******************************************************************
      subroutine mditeradd0(u,ndims,ifull,iused,ipin,v)
c Iterate over the array u, adding v(array) to it.
c Uses the mechanisms of mditerate. But inc is always 1.
c Normally, the starting pointer is ipin=0 for the full array.
      integer ndims
      integer ifull(ndims)
      integer iused(ndims)
      real u(*),v(*)

      integer mdims
      parameter (mdims=3)
c Effective index in dimension, c-style (zero based)
      integer indi(mdims),indinp(mdims)
c Structure vector
      integer iLs(mdims+1)
    
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
