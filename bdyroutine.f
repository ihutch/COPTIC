c**********************************************************************
c This file gives examples of boundary setting for sormpi.
c The bdyset routine can be called anything, and name passed.
c Its arguments must of course be correct.
      subroutine bdyset(ndims,ifull,iuds,cij,u,q)
c      call bdysetnull
      call bdysetfree(ndims,ifull,iuds,cij,u,q)
      end
c**********************************************************************
      subroutine bdysetfree(ndims,ifull,iuds,cij,u,q)
      integer ndims,ifull(ndims),iuds(ndims)
      real cij(*),u(*),q(*)
c Specify external the boundary setting routine.
c If    Bit-0 of islp is not set, then use logarithmic derivative.
c else  use Mach slope condition (higher bits relevant).
      external bdyslopeDh,bdyslopescreen,bdymach
      include 'plascom.f'
      include 'meshcom.f'
      include 'slpcom.f'
      include 'myidcom.f'
      integer ifirst
      data ifirst/1/
      ipoint=0
c      write(*,*)'islp=',islp
      if(ibits(islp,0,1).eq.0)then
c Normal log phi-derivative=-1 :
c         slpD=-1.
c Adaptive boundary condition. Only approximate for non-spheres. 
c Direct logarithmic gradient setting.
         slpD=-(1.+rs*sqrt(1.+1./Ti)/debyelen)
c         write(*,*)'slpD=',slpD,Ti,debyelen
         call mditerate(bdyslopeDh,ndims,ifull,iuds,ipoint,u)
c Explicit screening uses buggy bdyslopescreen. Obsolete.
c         slpD=debyelen/sqrt(1.+1./Ti)
c         call mditerate(bdyslopescreen,ndims,ifull,iuds,ipoint,u)
      else
c      write(*,*)'islp=',islp,vd,debyelen
c Use Mach boundary condition on slope only.
c To make Face 3 phi=0. Set bit-3 of islp =8
c Mach boundary condition for drift vd. M-value in slpD.
c Drift angles larger than about 1.5 cause instabilities.
         if(slpD.gt.1.5)then
            slpD=min(vd,1.5)
            if(myid.eq.0)write(*,*)
     $           'Mach BC slpD too large. Reset to',slpD
         endif
         if(ibits(islp,13,1).eq.1)then
c Second derivative setting on trailing edge.
            dx=xn(ixnp(4))-xn(ixnp(4)-1)
            if(vd*debyelen.gt.0)then
               dxk2=(dx/(vd*debyelen))**2
               if(ifirst.eq.1 .and. myid.eq.0)then
                  write(*,*)'2nd deriv BC, dxk2=',dxk2
                  if(dxk2.gt.1)
     $                 write(*,*)'dxk2 too large',dxk2,' Limit to 1.'
c Turn off alternative                  islp=islp-8192
               endif
               ifirst=0
               if(dxk2.gt.1.)dxk2=1.
            else
               if(myid.eq.0)write(*,*)
     $              'Cannot use second derivative when vd=',vd,
     $              ' debylen=',debyelen,' Switching off.'
               islp=islp-8192
            endif
         endif
         call mditerate(bdymach,ndims,ifull,iuds,ipoint,u)
      endif
      end
c************************************************************************
      subroutine bdyslopeDh(inc,ipoint,indi,ndims,iused,u)
c Version of bdyroutine that sets logarithmic 'radial' gradient
c equal to D=slpD
c BC is du/dr=D u/r     in the form   (ub-u0)=  D*(ub+u0)*f/(1-f)
c where f = Sum_j[(xb_j+x0_j)dx_j]/(2*rm^2), dx=xb-x0
c Thus ub=u0(1-f-D.f)/(1-f+D.f) using radii from position (0,0,..)
      integer ipoint,inc
      integer indi(ndims),iused(ndims)
      real u(*)

c Structure vector needed for finding adjacent u values.
c Can't be passed here because of mditerate argument conventions.
      parameter (mdims=10)
      integer iLs(mdims+1)
      common /iLscom/iLs

      include 'meshcom.f'
      include 'slpcom.f'
      D=slpD
c Algorithm: take steps of 1 in all cases except when on a lower
c boundary face of dimension 1 (and not other faces).  There the step is
c iused(1)-1.
c------------------------------------------------------------------
c Between here and ^^^ is boundary setting. Adjust upper and lower.
      r2=0.
      fac=0.
      ipd=ipoint
      do n=1,ndims
c BC is du/dr=D u/r     in the form   (ub-u0)=  D*(ub+u0)*f/(1-f)
c where f = Sum_j[(xb_j+x0_j)dx_j]/(2*rm^2), dx=xb-x0
c Thus ub=u0(1-f-D.f)/(1-f+D.f)
c Here we are using radii from position (0,0,..)
         x=xn(ixnp(n)+1+indi(n))
         r2=r2+x*x
         if(indi(n).eq.0)then
c On lower boundary face
            dx=xn(ixnp(n)+1)-xn(ixnp(n)+2)
            ipd=ipd+iLs(n)
            fac=fac+x*dx
            if(n.eq.1)then
c The exception in step. Do not change!
               inc=iused(1)-1
            else
               inc=1
            endif
         elseif(indi(n).eq.iused(n)-1)then
c On upper boundary face
            dx=xn(ixnp(n)+1+indi(n))-xn(ixnp(n)+indi(n))
            ipd=ipd-iLs(n)
            fac=fac+x*dx
            inc=1
         endif
      enddo
      if(ipd.eq.ipoint)then
         write(*,*)'BDY function error; we should not be here'
         stop
      else
         fac=fac/(2.*r2)
         u(ipoint+1)=u(ipd+1) *(1.-fac+D*fac)/(1.-fac-D*fac)
      endif
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c      write(*,*)'indi,inc,iused,ipoint',indi,inc,iused,ipoint
      end
c************************************************************************
      subroutine bdymach(inc,ipoint,indi,ndims,iused,u)
c Version of bdyroutine that sets the BC on the x, y boundaries
c as being  du/dr + M du/dz =0. (r the cylindrical radius)
c The z-mach number, M, is the value of slpD in slpcom.
c Default z-boundary conditions: du/dz=0
c islp indicates other BC choices as follows:
c    Bit-3  (8): set BC at lower-z boundary u=0.
c    Bit-13 (8192): set BC at upper-z boundary u''= -k^2u, k=1/M.lambdaD
c dxk2 is in slpcom.
      integer ipoint,inc
      integer indi(ndims),iused(ndims)
      real u(*)
c Structure vector needed for finding adjacent u values.
c Can't be passed here because of mditerate argument conventions.
      parameter (mdims=10)
      integer iLs(mdims+1)
      common /iLscom/iLs
c Value of mach number passed in common.
      include 'meshcom.f'
      include 'slpcom.f'
      DM=slpD
c Algorithm: take steps of 1 in all cases except when on a lower
c boundary face of dimension 1 (and not other faces).  There the step is
c iused(1)-1.
c------------------------------------------------------------------
c Between here and ^^^ is boundary setting. Adjust upper and lower.
      r2=0.
      fac=0.
      ipd=ipoint
c The dimensions up to the last are treated the same.
      do n=1,ndims-1
c Here we are using radii from position (0,0,..)
         x=xn(ixnp(n)+1+indi(n))
         r2=r2+x*x
         if(indi(n).eq.0)then
c On lower boundary face
            dx=xn(ixnp(n)+1)-xn(ixnp(n)+2)
            ipd=ipd+iLs(n)
            fac=fac+x*dx
            if(n.eq.1)then
c The exception in step. Do not change!
               inc=iused(1)-1
            else
               inc=1
            endif
         elseif(indi(n).eq.iused(n)-1)then
c On upper boundary face
            dx=xn(ixnp(n)+1+indi(n))-xn(ixnp(n)+indi(n))
            ipd=ipd-iLs(n)
            fac=fac+x*dx
            inc=1
         endif
      enddo
c Finally we do the last dimension, z.
      n=ndims
      if(ipd.eq.ipoint)then
c We only are here if we are on z-faces, and not on x or y faces. 
c Fall back to du/dz=0.
         if(indi(ndims).eq.0)then
            u(ipoint+1)=u(ipoint+1+iLs(n))
         else
            if(ibits(islp,13,1).ne.0)then
               u(ipoint+1)=(2.-dxk2)*u(ipoint+1-iLs(n))
     $              -u(ipoint+1-2*iLs(n))
            else
               u(ipoint+1)=u(ipoint+1-iLs(n))
            endif
         endif
         inc=1
      else
c fac contains r.dr, r2 contains r^2
c u_p=u_d + (-M*(dx.x + dy.y)/r + dz).(du/dz)
         if(indi(n).eq.0)then
c On lower z-boundary face (as well as x or y)
            dz=xn(ixnp(n)+1)-xn(ixnp(n)+2)
            ipd=ipd+iLs(n)            
            inc=1
c z-derivative calculated from value 3
            dubydz=(u(ipd+1)-u(ipd+iLs(n)+1))
     $           /(xn(ixnp(n)+2)-xn(ixnp(n)+3))
         elseif(indi(n).eq.iused(n)-1)then
c On upper z-boundary face (as well as x or y)
            dz=xn(ixnp(n)+1+indi(n))-xn(ixnp(n)+indi(n))
            ipd=ipd-iLs(n)
            inc=1
c z-derivative calculated from value iused-3.
            dubydz=(u(ipd+1)-u(ipd-iLs(n)+1))
     $           /(xn(ixnp(n)+indi(n))-xn(ixnp(n)+indi(n)-1))
         else
c Not on z-face
            dz=0.
            dubydz=(u(ipd+1+iLs(n))-u(ipd+1-iLs(n)))
     $           /(xn(ixnp(n)+indi(n)+2)-xn(ixnp(n)+indi(n)))
         endif
c Set u-value by extrapolation using dubydz. 
c Since the radial(xy) derivative is equal to -M*du/dz,
c -M*(fac/r)*dubydz is the radial difference, and add z-difference.
         u(ipoint+1)=u(ipd+1)+
     $        (-DM*fac/sqrt(r2) + dz)*dubydz
      endif
c Special cases:
      if(ibits(islp,3,1).ne.0)then
c Bit 3(+1) (=8) set, put lower z-face to zero.
         if(indi(n).eq.0)u(ipoint+1)=0.
      endif
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c      write(*,*)'indi,inc,iused,ipoint',indi,inc,iused,ipoint
      end
c**********************************************************************
c The logic of the following might miss some corners.
c************************************************************************
      subroutine bdyslopescreen(inc,ipoint,indi,ndims,iused,u)
c Version of bdyroutine that sets logarithmic 'radial' gradient
c equal to that for a Debye screened potential: -(1+r/lambdas)
c slpD needs to be set to lamdas prior to call.
      integer ipoint,inc
      integer indi(ndims),iused(ndims)
      real u(*)

c Structure vector needed for finding adjacent u values.
c Can't be passed here because of mditerate argument conventions.
      parameter (mdims=10)
      integer iLs(mdims+1)
      common /iLscom/iLs

      real x(mdims)
      include 'meshcom.f'
      include 'slpcom.f'

      r2=r2indi(ndims,indi,x)
      r1=sqrt(r2)
      D=slpD

c Algorithm: take steps of 1 in all cases except
c when on a lower boundary face of dimension 1. 
c There the step is iused(1)-1.
      inc=1
      do n=ndims,1,-1
c------------------------------------------------------------------
c Between here and ^^^ is boundary setting. Adjust upper and lower.
         if(indi(n).eq.0)then
c The exception in step. Do not change!:
            if(n.eq.1)inc=iused(1)-1
c On lower boundary face
c du/dn= -(1.+r/D) r.n u/ r^2
            dx=xn(ixnp(n)+2)-xn(ixnp(n)+1)
            fac=(1.+0.5*dx/x(n))*r2*2./(-(1.+r1/D)*x(n)*dx)
c            write(*,*)n,dx,xn(ixnp(n)+1),fac,r2
c Centered BC:
c (u2-u1)/dx=(u2+u1)/2 *S*(x1+x2)/2 /rm^2  put fac=2*rm^2/S xm dx
c  u2(fac-1.)=u1(fac+1)
            u(ipoint+1)=u(ipoint+1+iLs(n))*(fac-1.)/(1.+fac)
            goto 101
         elseif(indi(n).eq.iused(n)-1)then
c On upper boundary face
c du/dn=-(1.+r/D) r.n u/ r^2
            dx=-(xn(ixnp(n)+1+indi(n))-xn(ixnp(n)+indi(n)))
            fac=(1.+0.5*dx/x(n))*r2*2./(-(1.+r1/D)*x(n)*dx)
            u(ipoint+1)=u(ipoint+1-iLs(n))*(fac-1.)/(1.+fac)
            goto 101
         endif
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      enddo
      write(*,*)'BDY Error. We should not be here',n,ipoint,indi
      stop
 101  continue
c      write(*,*)'indi,inc,iused,ipoint',indi,inc,iused,ipoint
      end
c**********************************************************************
c Return the radius squared and the vector position of point indexed as
c indi(ndims), relative to the center of the coordinate ranges,
c using information in meshcom.
      function r2indi(ndims,indi,x)
      integer indi(ndims)
      real x(ndims)
      include 'meshcom.f'
      if(ndims.ne.ndims_mesh)then
         write(*,*)'rindi dimension mismatch',ndims,ndims_mesh
         stop
      endif
      r2indi=0.
      do i=1,ndims
         x(i)=xn(ixnp(i)+indi(i)+1)-(xn(ixnp(i)+1)+xn(ixnp(i+1)))/2.
         r2indi=r2indi+x(i)**2
      enddo
      end
c********************************************************************
c********************************************************************
c Obsolete or Currently unused routines:
c**********************************************************************
      subroutine bdyset3sl(ndims,ifull,iuds,cij,u,q)
      integer ndims,ifull(ndims),iuds(ndims)
      real cij(*),u(*),q(*)
c Specify external the boundary setting routine.
      external bdy3slope 
c sets the derivative to zero on boundaries 3.
      ipoint=0
      call mditerate(bdy3slope,ndims,ifull,iuds,ipoint,u)
      end
c************************************************************************
      subroutine bdy3slope(inc,ipoint,indi,ndims,iused,u)
c Version of bdyroutine that sets derivative=0 on 3-boundary.
      integer ipoint,inc
      integer indi(ndims),iused(ndims)
      real u(*)

c Structure vector needed for finding adjacent u values.
c Can't be passed here because of mditerate argument conventions.
      parameter (mdims=10)
      integer iLs(mdims+1)
      common /iLscom/iLs

c Algorithm: take steps of 1 in all cases except
c when on a lower boundary face of dimension 1. 
c There the step is iused(1)-1.
      inc=1
      do n=ndims,1,-1
c------------------------------------------------------------------
c Between here and ^^^ is boundary setting. Adjust upper and lower.
         if(indi(n).eq.0)then
c The exception in step. Do not change!:
            if(n.eq.1)inc=iused(1)-1
c On lower boundary face
            u(ipoint+1)=0.
            if(n.eq.3)then
c First derivative is zero:
               u(ipoint+1)=u(ipoint+1+iLs(n))
            endif
c Second derivative is zero:
c               u(ipoint+1)=2.*u(ipoint+1+iLs(n))-u(ipoint+1+2*iLs(n))
            goto 101
         elseif(indi(n).eq.iused(n)-1)then
c On upper boundary face
            u(ipoint+1)=0.
            if(n.eq.3) u(ipoint+1)=u(ipoint+1-iLs(n))
            goto 101
         endif
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      enddo
      write(*,*)'BDY3s Error. We should not be here',n,ipoint,indi
      stop
 101  continue
c      write(*,*)'indi,inc,iused,ipoint',indi,inc,iused,ipoint
      end
c**********************************************************************
      subroutine bdyslope0(ndims,ifull,iuds,cij,u,q)
      integer ndims,ifull(ndims),iuds(ndims)
      real cij(*),u(*),q(*)
c Specify external the boundary setting routine.
      external bdyslopeDh,bdyslopescreen,bdymach
      include 'slpcom.f'
      ipoint=0
      islp=0
      slpD=0.
      call mditerate(bdyslopeDh,ndims,ifull,iuds,ipoint,u)
      end
c**********************************************************************
      subroutine bdysetnull()
c   (ndims,ifull,iuds,cij,u,q)
c     Null version
      end
c************************************************************************

c********************************************************************
