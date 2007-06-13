c***********************************************************************
      subroutine mditerate(ndims,ifull,iused,routine,u,ipoint)
c multidimensional iteration. For dimensions ndims, iterate over the
c array whose full dimensions are ifull(ndims) used iused(ndims).
c At each iteration, 
c         routine(inc,ipoint,indi,ndims,iused,u) is called
c which accepts a pointer to the address in the array structure: ipoint
c and returns the next increment in units of the lowest dimension: inc
c The index for each dimension is passed in: indi(ndims)
c And can be adjusted, but ipoint must also be adjusted equivalently.
c Increments wrap at the used dimension length. Pointer uses ifull.
c u is the data which is passed to the routine
c ipoint is the offset of the starting point from the start of u,
c referenced to u(iused()) as if it were the full array.
c Iteration does not extend past the used bounds of u.
c Therefore a shifted subarray should be treated by specifying ipoint=0
c and giving the address of the subarray as u. 
      integer ndims
      integer ifull(ndims)
      integer iused(ndims)
      external routine
      real u(*)

      integer mdims
      parameter (mdims=10)
c Effective index in dimension, c-style (zero based)
      integer indi(mdims)
c Structure vector
      integer iLs(mdims+1)
      common /iLscom/iLs
      
c      write(*,*)'ndims',ndims,' ifull',ifull,' iused',iused
      iLs(1)=1
      do n=1,ndims
         indi(n)=0
         iLs(n+1)=iLs(n)*ifull(n)
      enddo
     
c Offset. 
      indi(1)=ipoint

      inc=1
      icount=0
      n=1
c Iteration over the multidimensional array
 101  continue
c      write(*,'(''('',i1,i4,'') '',$)')n,indi(n)
      if(indi(n).gt.iused(n)-1)then
c     Overflow. Subtract off enough (inm) of next dimension
c     and move ipoint to appropriate position.
         inm=0
 102     inm=inm+1
         ipoint=ipoint+iLs(n+1)-iused(n)*iLs(n)
         indi(n)=indi(n)-iused(n)
         if(indi(n).gt.iused(n)-1)goto 102
c Increment the next level.
 103     n=n+1
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
