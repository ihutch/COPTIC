!************************************************************************
! Shell for skipping over the boundary of ndims-dimensional volume.
      subroutine bdyroutine(inc,ipoint,indi,ndims,iused,u)
      integer ipoint,inc
      integer indi(ndims),iused(ndims)
      real u(*)

! Algorithm: if on a boundary face of dimension >1, steps of 1 (dim 1).
! Otherwise steps of iused(1)-1 or 1 on first or last (of dim 1).
      inc=1
      do n=ndims,2,-1
         if(indi(n).eq.0)then
! On boundary face 0 of dimension n>1. Break.
! This is where we put boundary setting for n>1
! Example:
            u(ipoint+1)=0.
!            if(n.eq.3)then
! Second derivative is zero:
!               u(ipoint+1)=2.*u(ipoint+1+iLs(n))-u(ipoint+1+2*iLs(n)
! First derivative is zero:
!               u(ipoint+1)=u(ipoint+1+iLs(n))
!            endif
            goto 101
         elseif(indi(n).eq.iused(n)-1)then
! On boundary face iused(n) of dimension n>1. Break.
! Example:
            u(ipoint+1)=0.
! 
!            if(n.eq.3) u(ipoint+1)=u(ipoint+1-iLs(n))
            goto 101
         endif
      enddo
!     We are not on any higher boundary.
! This is where the boundary setting is done for n=1
! Example:
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
!      write(*,*)'indi,inc,iused,ipoint',indi,inc,iused,ipoint

      end
!***********************************************************************
      subroutine mditerarg(routine,ndims,ifull,iused,ipin,t,u,v,w,x)
! multidimensional iteration. For dimensions ndims, iterate over the
! array whose full dimensions are ifull(ndims) used iused(ndims).
! This version allows five (array) arguments: t,u,v,w,x to be passed.
! At each iteration, the routine (which should be declared external)
!
!         routine(inc,ipoint,indi,ndims,iLs,iused,t,u,v,w,x) is called,
! whose arguments are
!      integer ipin,inc
!      integer indi(ndims),iused(ndims)
!      real t(*),u(*),v(*),w(*),x(*)
!
! which accepts a pointer to the address in the array structure: ipoint
! referred to the full dimensions of u, so u(1+ipoint) is processed.  It
! returns the next increment in units of the lowest dimension: inc. The
! offsets for each dimension are passed in: indi(ndims) And can (but
! need not) be adjusted in routine. If they are, ipoint must also be
! adjusted equivalently in routine. Again: in routine, ipoint addresses
! u(ifull).  After the return of routine, ipoint and indi(1) are
! incremented by inc and wrapping is handled.
!
! Wrapping occurs at the
! _iused_ dimension length. In other words, if the incremented indi(id)
! passes the used length iused(id), then the next higher dimension is
! incremented by 1, and indi(id) decremented by iused(id). Thus inc is
! the increment in the lowest dimension index, referenced to the _used_
! dimensions, but ipoint is wrapped accounting for the knowledge of
! ifull so it points to the correct next position within the array. u is
! the data that is passed to the routine. On entry, ipin is the
! starting offset storage position from the start of u, but in this case
! (only), relative to the iused array structure. An error will occur if
! ipoint.gt.iused(1) on entry.  Iteration does not extend past the used
! bounds of u.  Therefore a shifted subarray should be treated by
! specifying ipoint=0 and giving the address of the subarray as u.
!
! Normal usage is to pass ipoint=0 to mditerate. Then the routine is called
! for the whole _used_ array at nodes spaced by inc in the 1st dimension.

      integer ndims
      integer ifull(ndims)
      integer iused(ndims)
      external routine
      real t(*),u(*),v(*),w(*)

      integer mdims
      parameter (mdims=3)
! Effective index in dimension, c-style (zero based)
      integer indi(mdims),indinp(mdims)
! Structure vector
      integer iLs(mdims+1)
      
!      write(*,*)'ndims',ndims,' ifull',ifull,' iused',iused
! Offset.
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
!      icount=0
      n=1
! Iteration over the multidimensional array
 101  continue
!      write(*,'(''('',i1,i4,'') '',$)')n,indi(n)
      if(indi(n)-indinp(n).gt.iused(n)-1)then
!     Overflow. Subtract off enough (inm) of next dimension
!     and move ipoint to appropriate position in full array.
         inm=0
 102     inm=inm+1
         ipoint=ipoint+iLs(n+1)-iused(n)*iLs(n)
         indi(n)=indi(n)-iused(n)
         if(indi(n)-indinp(n).gt.iused(n)-1)goto 102
! Increment the next level.
         n=n+1
         if(n.gt.ndims)goto 201
         indi(n)=indi(n)+inm
         goto 101
      elseif(n.gt.1)then
! We've carried and cleared an increment.
! Return stepwise to base level
         n=n-1
         goto 101
      else
! We're at the base level and have succeeded in incrementing.
! Do whatever we need to and increment indi(1) and ipoint
!        write(*,'(3i3,$)')(iused(ii),ii=1,3)
!        write(*,'(3i3,i7,i5)')(indi(i),i=1,ndims),ipoint,inc
         call routine(inc,ipoint,indi,ndims,iLs,iused,t,u,v,w,x)
         indi(1)=indi(1)+inc
         ipoint=ipoint+inc
         goto 101
      endif

 201  continue
! Reached the end.

      end
!*****************************************************************
      integer function indexcontract(ndims,ifull,indi)
! On entry indi(ndims) contains the (ndims)-dimensional offsets
! to position in array whose full dimensions are ifull(ndims)
! Returns the corresponding linear offset index within that array.
      integer ndims
      integer ifull(ndims),indi(ndims)
      indexcontract=indi(ndims)
      do id=ndims-1,1,-1
         indexcontract=indi(id)+indexcontract*ifull(id)
      enddo
      end
!******************************************************************
      subroutine offsetexpand(ndims,ifull,index,indi)
! On entry index is the zero-based pointer to the position in the 
! array whose full dimensions are ifull(ndims). 
! On exit indi contains the corresponding (ndims)-dimensional offsets.
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
     $     ,index,'  ndims=',ndims,' ifull=',ifull
      end
!******************************************************************
! New code for the rationalized mditerator.
!******************************************************************
      integer function mditerator(ndims,iview,indi,isw,iaux)
! Multidimensional iterator for ndims dimensions over an abstract
! view iview of an array. Zero-based offsets.
! On entry:
!    iview(3,ndims) contains the view of the array in the form
!       istart1,iend1,istride1;istart2,iend2,...,istridendims.
!    iaux(ndims) maybe contains extra information.
!    isw =0 iterate, else setup iview
!         1 set iview(istart,j)=iaux(j)
!         2 set iview(iend,j)=iaux(j)-1
!         3 set iview(istride,j)=iaux(j)
!         4 set iview(*,j) = 0,iaux(j)-1,1
!    indi(ndims) is the current offset-array (set to 0 if isw=4)
! On exit: indi is the new index-array. Optionally iview is [re]set.
!    Return value is 0 if iteration incomplete, else positive.
! Example of typical usage:
!      icomplete=mditerator(ndims,iview,indi,4,iused)
! 1    u(1+indexcontract(ndims,ifull,indi))=mod(1+indi(3),10)
!      if(mditerator(ndims,iview,indi,0,iused).eq.0)goto 1
!      
      integer ndims,iview(3,ndims),indi(ndims),isw,iaux(ndims)
      integer istart,iend,istride
      parameter (istart=1,iend=2,istride=3)

      mditerator=isw
      if(isw.eq.0)then
!-----------------------------------------------------
! Iterate
! Start at dimension 1.
         n=1
! Increment:
         indi(n)=indi(n)+iview(istride,n)
 101     continue
!      write(*,'(''('',i1,i4,'') '',$)')n,indi(n)
         if((indi(n)-iview(iend,n))*iview(istride,n).gt.0)then
!     Overflow. Subtract off enough (inm) of next dimension
            inm=0
 102        inm=inm+1
            indi(n)=indi(n)-(iview(iend,n)-iview(istart,n)+1)
!            write(*,*)inm,indi(n),iview(3,n),iview(2,n)
            if((indi(n)-iview(iend,n))*iview(istride,n).gt.0)goto 102
! Move to the next dimension:
            n=n+1
! Overflow of ndims dimension completes iteration:
            if(n.gt.ndims)then 
               mditerator=1
               return
            endif
! Increment this dimension with carry and restart:
            indi(n)=indi(n)+inm*iview(istride,n)
            goto 101
         elseif(n.gt.1)then
! We've carried and cleared an increment.
! Return stepwise to base level
            n=n-1
            goto 101
         endif
! Here when we have successfully incremented all levels.
!-------------------------------------------------------
      elseif(isw.eq.4)then
! Total setup
         do i=1,ndims
            iview(istart,i)=0
            iview(iend,i)=iaux(i)-1
            iview(istride,i)=1
            indi(i)=0
         enddo
      else
! Individual setup.
         do i=1,ndims
            iview(isw,i)=iaux(i)
            if(isw.eq.2)iview(isw,i)=iview(isw,i)-1
         enddo
      endif
      end
!*****************************************************************
      integer function mditer(ndims,ifull,iused,index)
! Uses the mditerator for full array without having to supply extra
! arrays. On first call iused is the used dimensions.  Works by
! returning iused as the current index until the array is completed,
! then it is reset. index is the compact index position into the
! array. It handles automatically the initialization and so does not
! have to be called differently for initialization.  
!
! Example of typical usage: 
!   icomplete=mditer(ndims,ifull,iused,index)
! 1  u(index)=mod(1+iused(3),10)
!   if(mditer(ndims,ifull,iused,index).eq.0)goto 1
!      
      implicit none
      integer ndims,index
      integer ifull(ndims),iused(ndims)
      integer isw,ndimsmax,i,nzeros
      parameter (ndimsmax=5,nzeros=3*ndimsmax)
      integer iview(3,ndimsmax),indexcontract,mditerator
! iview is explicitly zeroed to avoid uninitialized warnings.
      data isw/1/iview/nzeros*0/
      external indexcontract
      if(isw.ne.0)then
         if(.not.ndims.le.ndimsmax)Stop 'Too many dimensions in mditer'
         i=mditerator(ndims,iview,iused,4,iused)
! Now iview(2,*)=iused(*)-1
         mditer=1
         isw=0
      else
         mditer=mditerator(ndims,iview,iused,0,iused)
         if(mditer.ne.0)then
! Restore iused and return to initial state.
            do i=1,ndims
               iused(i)=iview(2,i)+1
            enddo
            isw=1
         endif
      endif
      index=indexcontract(ndims,ifull,iused)+1
      end

!*****************************************************************
      subroutine mditerset(u,ndims,ifull,iused,ipin,v)
! Iterate over the array u, setting u(i)=v scalar.
! Reimplementation using mditerator.
! Normally, the starting pointer is ipin=0 for the full array.
      integer ndims,ifull(ndims),iused(ndims)
      real u(*),v
      integer mdims
      parameter (mdims=3)
! Effective index in dimension, c-style (zero based)
      integer indi(mdims),iview(3,mdims)
! This data statement serves to silence ftnchek. The first mditerator
! call actually initializes the iview and indi.
      data iview/mdims*0,mdims*0,mdims*0/
      data indi/mdims*0/
      
! Set to iterate the whole used array:
      icomplete=mditerator(ndims,iview,indi,4,iused)
! Reset starting indices from input pointer, in case non-zero.
      if(ipin.ne.0)call offsetexpand(ndims,ifull,ipin,indi)
 1      ii=1+indexcontract(ndims,ifull,indi)
        u(ii)=v
      if(mditerator(ndims,iview,indi,0,iused).eq.0)goto 1

      end
!******************************************************************
!*****************************************************************
      subroutine mditermults(u,ndims,ifull,iused,ipin,v,w)
! Iterate over the array u, multiplying it by v(scalar)
! and adding w(scalar): u = u*v+w.
! Reimplementation using mditerator.
! Normally, the starting pointer is ipin=0 for the full array.
      integer ndims,ifull(ndims),iused(ndims)
      real u(*),v,w
      integer mdims
      parameter (mdims=3)
! Effective index in dimension, c-style (zero based)
      integer indi(mdims),iview(3,mdims)
! This data statement serves to silence ftnchek. The first mditerator
! call actually initializes the iview and indi.
      data iview/mdims*0,mdims*0,mdims*0/
      data indi/mdims*0/
      
! Set to iterate the whole used array:
      icomplete=mditerator(ndims,iview,indi,4,iused)
! Reset starting indices from input pointer, in case non-zero.
      if(ipin.ne.0)call offsetexpand(ndims,ifull,ipin,indi)
 1      ii=1+indexcontract(ndims,ifull,indi)
        u(ii)=u(ii)*v+w
      if(mditerator(ndims,iview,indi,0,iused).eq.0)goto 1

      end
!*****************************************************************
      subroutine mditeradd(u,ndims,ifull,iused,ipin,v)
! Iterate over the array u, adding v(array) to it: u(i)=u(i)+v(i)
! Reimplementation using mditerator.
! Normally, the starting pointer is ipin=0 for the full array.
      integer ndims,ifull(ndims),iused(ndims)
      real u(*),v(*)
      integer mdims
      parameter (mdims=3)
! Effective index in dimension, c-style (zero based)
      integer indi(mdims),iview(3,mdims)
! This data statement serves to silence ftnchek. The first mditerator
! call actually initializes the iview and indi.
      data iview/mdims*0,mdims*0,mdims*0/
      data indi/mdims*0/
      
! Set to iterate the whole used array:
      icomplete=mditerator(ndims,iview,indi,4,iused)
! Reset starting indices from input pointer, in case non-zero.
      if(ipin.ne.0)call offsetexpand(ndims,ifull,ipin,indi)
 1      ii=1+indexcontract(ndims,ifull,indi)
        u(ii)=u(ii)+v(ii)
      if(mditerator(ndims,iview,indi,0,iused).eq.0)goto 1

      end
!*****************************************************************
      subroutine mditerave(u,ndims,ifull,iused,ipin,uave,istepave)
! Iterate over the array u, to produce uave in the form
! ((istepave-1)*uave + u)/istepave

! Normally, the starting pointer is ipin=0 for the full array.
      integer ndims,ifull(ndims),iused(ndims),istepave
      real u(*),uave(*)
      integer mdims
      parameter (mdims=3)
! Effective index in dimension, c-style (zero based)
      integer indi(mdims),iview(3,mdims)
! This data statement serves to silence ftnchek. The first mditerator
! call actually initializes the iview and indi.
      data iview/mdims*0,mdims*0,mdims*0/
      data indi/mdims*0/
      
      if(mdims.lt.ndims)then
         write(*,*)'Dimension error in mditerave',mdims,ndims
         stop
      endif
! Set to iterate the whole used array:
      icomplete=mditerator(ndims,iview,indi,4,iused)
! Reset starting indices from input pointer, in case non-zero.
      if(ipin.ne.0)call offsetexpand(ndims,ifull,ipin,indi)
 1      ii=1+indexcontract(ndims,ifull,indi)
        uave(ii)=(uave(ii)*(istepave-1.)+u(ii))/float(istepave)
      if(mditerator(ndims,iview,indi,0,iused).eq.0)goto 1

      end

!************************************************************************
!************************************************************************
      subroutine mditerarg1(routine,ndims,ifull,iused,ipin,t,u,v,w)
! [Reimplemented using mditerator. Runs about 50% slower for fast routines]
! multidimensional iteration. For dimensions ndims, iterate over the
! array whose full dimensions are ifull(ndims) used iused(ndims).
! This version allows four (array) arguments: t,u,v,w to be passed.
! At each iteration, the routine (which should be declared external)
!
!         routine(inc,ipoint,indi,ndims,iLs,iused,t,u,v,w) is called,
! whose arguments are
!      integer ipin,inc
!      integer indi(ndims),iused(ndims)
!      real t(*),u(*),v(*),w(*)
!
! which accepts a pointer to the address in the array structure: ipoint
! referred to the full dimensions of u, so u(1+ipoint) is processed.  It
! returns the next increment in units of the lowest dimension: inc. The
! offsets for each dimension are passed in: indi(ndims) And can (but
! need not) be adjusted in routine. If they are, ipoint must also be
! adjusted equivalently in routine. Again: in routine, ipoint addresses
! u(ifull).  After the return of routine, ipoint and indi(1) are
! incremented by inc and wrapping is handled.
!
! Wrapping occurs at the _iused_ dimension length relative to the
! starting indices.  In other words, if the incremented
! indi(id)-indinp(id)+1 exceeds the used length iused(id), then the next
! higher dimension is incremented by 1, and indi(id) decremented by
! iused(id). Thus inc is the increment in the lowest dimension index,
! referenced to the _used_ dimensions, but ipoint is wrapped accounting
! for the knowledge of ifull so it points to the correct next position
! within the array. u is the data that is passed to the routine. On
! entry, ipin is the starting offset storage position from the start of
! u, but in this case (only), relative to the iused array structure. An
! error will occur if ipoint.gt.iused(1) on entry.  Iteration does not
! extend past the used bounds of u.  Therefore a shifted subarray should
! be treated by specifying ipoint=0 and giving the address of the
! subarray as u.
!
! Normal usage is to pass ipoint=0 to mditerate. Then the routine is called
! for the whole _used_ array at nodes spaced by inc in the 1st dimension.
      implicit none
      integer ndims,ipin
      integer ifull(ndims)
      integer iused(ndims)
      external routine
      real t(*),u(*),v(*),w(*)

      integer mdims,i
      parameter (mdims=3)
! Effective index in dimension, c-style (zero based)
      integer indi(mdims),iview(3,mdims)
      integer ipoint,inc,icomplete,indexcontract,mditerator
! This data statement serves to silence ftnchek. The first mditerator
! call actually initializes the iview and indi.
      data iview/mdims*0,mdims*0,mdims*0/
      data indi/mdims*0/

! Structure vector
      integer iLs(mdims+1)

      inc=1
! Set to iterate the whole used array:
      icomplete=mditerator(ndims,iview,indi,4,iused)
! Reset starting indices from input pointer, in case non-zero.
      call offsetexpand(ndims,ifull,ipin,indi)
      iLs(1)=1
      do i=1,ndims
! Structure vector must be set.
         iLs(i+1)=iLs(i)*ifull(i)
         if(indi(i)+iused(i).gt.ifull(i))then
            write(*,*) 'MDITERATE Error. Dimension',i,' indi + iused',
     $           indi(i),iused(i),' .gt. ifull',ifull(i)
            stop
         endif
! Iteration from starting index
         iview(1,i)=indi(i)
! With length equal to iused.
         iview(2,i)=iused(i)+indi(i)-1
      enddo
!      write(*,*)iview
! ipoint update cannot be put in mditerator, because it does not have
! the ifull information.
 1      ipoint=indexcontract(ndims,ifull,indi)
!        write(*,'(3i3,i7,i5)')(indi(i),i=1,ndims),ipoint,inc
        call routine(inc,ipoint,indi,ndims,iLs,iused,t,u,v,w)
!        write(*,*)indi,ipoint,inc,t(1+ipoint)
! Set the iterator view according to the increment of the first
! dimension. 
        iview(3,1)=inc
      if(mditerator(ndims,iview,indi,0,iused).eq.0)goto 1

      end
!******************************************************************
!******************************************************************
!   Obsolete below here.
!******************************************************************
!******************************************************************
