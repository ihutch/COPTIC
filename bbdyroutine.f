c**********************************************************************
      subroutine bdyset(ndims,ifull,iuds,cij,u,q)
c     Null version
      end
c**********************************************************************
      subroutine bdyset3(ndims,ifull,iuds,cij,u,q)
      integer ndims,ifull(ndims),iuds(ndims)
      real cij(*),u(*),q(*)
      external bdy3slope
c set the derivative to zero on boundaries 3.
      ipoint=0
      call mditerate(ndims,ifull,iuds,bdy3slope,u,ipoint)

c set the second derivative to zero on max j
c      do i=2,ni-1
c         u(i,nj)=relax*(2.*u(i,nj-1)-u(i,nj-2)) +(1.-relax)*u(i,nj)
c      enddo
      end
c**********************************************************************
c     L(u) + f(u) = q(x,y,...), 
c     where L is a second order elliptical differential operator 
c     represented by a difference stencil of specified coefficients,
c     f is some additional function, and q is the "charge density".
c f is exp(u) here for Boltzmann electrons and densities normalized
c to unity at infinity.
      real function faddu(u,fprime)
      real u,fprime
      fprime=exp(u)
      faddu=fprime
      end
c************************************************************************
      subroutine bdy3slope(inc,ipoint,indi,ndims,iused,u)
c Version of bdyroutine that sets derivative=0 on 3-boundary.
      integer ipoint,inc
      integer indi(ndims),iused(ndims)
      real u(*)

      parameter (mdims=10)
c Structure vector needed for finding adjacent u values.
c Can't be passed here because of mditerate argument conventions.
      integer iLs(mdims+1)
      common /iLscom/iLs

c Algorithm: if on a boundary face of dimension >1, steps of 1 (dim 1).
c Otherwise steps of iused(1)-1 or 1 on first or last (of dim 1).
      inc=1
      do n=ndims,2,-1
         if(indi(n).eq.0)then
c On boundary face 0 of dimension n>1. Break.
c This is where we put boundary setting for n>1
            u(ipoint+1)=0.
            if(n.eq.3)then
c Second derivative is zero:
c               u(ipoint+1)=2.*u(ipoint+1+iLs(n))-u(ipoint+1+2*iLs(n)
c First derivative is zero:
               u(ipoint+1)=u(ipoint+1+iLs(n))
            endif
            goto 101
         elseif(indi(n).eq.iused(n)-1)then
c On boundary face iused(n) of dimension n>1. Break.
            u(ipoint+1)=0.
            if(n.eq.3) u(ipoint+1)=u(ipoint+1-iLs(n))
            goto 101
         endif
      enddo
c     We are not on any higher boundary.
c This is where the boundary setting is done for n=1
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
