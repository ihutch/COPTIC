!**********************************************************************
! This file gives examples of boundary setting for sormpi.
! The bdyset routine can be called anything, and name passed.
! Its arguments must of course be correct.
      subroutine bdyset(mdims,ifull,iuds,cij,u,q)
      real u(*),cij(*),q(*)
      include 'ndimsdecl.f'
      include 'facebcom.f'
! Specify external the boundary setting routine.
      external bdyface 
      if(LF)then
! Rectangular face setting
!         write(*,*)'Calling bdyface setting',LF
         ipoint=0
         call mditerarg(bdyface,mdims,ifull,iuds,ipoint,u)
      else
         call bdysetfree(mdims,ifull,iuds,cij,u,q)
!      call bdysetnull
      endif
      end
!**********************************************************************
      subroutine bdysetfree(mdims,ifull,iuds,cij,u,q)
!      integer mdims,ifull(mdims),iuds(mdims)
      real u(*),cij(*),q(*)
! Specify external the boundary setting routine.
! If    Bit-0 of islp is not set, then use logarithmic derivative.
! else  use Mach slope condition (higher bits relevant).
      external bdyslopeDh,bdyslopescreen,bdymach
      include 'ndimsdecl.f'
      include 'plascom.f'
      include 'meshcom.f'
      include 'slpcom.f'
      include 'myidcom.f'
      integer ifirst
      data ifirst/1/
! Silence warnings about unused arguments
      use=cij(1)
      use=q(1)
      ipoint=0
!      write(*,*)'islp=',islp
      if(ibits(islp,0,1).eq.0)then
! Normal log phi-derivative=-1 :
!         slpD=-1.
! Adaptive boundary condition. Only approximate for non-spheres. 
! Direct logarithmic gradient setting.
         if(debyelen.eq.0)then
            slpD=-1.
         else
            slpD=-(1.+rs*sqrt(1.+1./Ti)/debyelen)
         endif
!         write(*,*)'slpD=',slpD,Ti,debyelen
         call mditerarg(bdyslopeDh,ndims,ifull,iuds,ipoint,u)
! Explicit screening uses buggy bdyslopescreen. Obsolete.
!         slpD=debyelen/sqrt(1.+1./Ti)
!         call mditerarg(bdyslopescreen,ndims,ifull,iuds,ipoint,u)
      else
!      write(*,*)'islp=',islp,vd,debyelen
! Use Mach boundary condition on slope only.
! To make Face 3 phi=0. Set bit-3 of islp =8
! Mach boundary condition for drift vd. M-value in slpD.
! Drift angles larger than about 1.5 cause instabilities.
         if(slpD.gt.1.5)then
            slpD=min(vd,1.5)
            if(myid.eq.0)write(*,*)
     $           'Mach BC slpD too large. Reset to',slpD
         endif
         if(ibits(islp,13,1).eq.1)then
! Second derivative setting on trailing edge.
            dx=xn(ixnp(4))-xn(ixnp(4)-1)
            if(vd*debyelen.gt.0)then
               dxk2=(dx/(vd*debyelen))**2
               if(ifirst.eq.1 .and. myid.eq.0)then
                  write(*,*)'2nd deriv BC, dxk2=',dxk2
                  if(dxk2.gt.1)
     $                 write(*,*)'dxk2 too large',dxk2,' Limit to 1.'
! Turn off alternative                  islp=islp-8192
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
         call mditerarg(bdymach,ndims,ifull,iuds,ipoint,u)
      endif
      end
!************************************************************************
      subroutine bdyslopeDh(inc,ipoint,indi,mdims,iLs,iused,u)
! Version of bdyroutine that sets logarithmic 'radial' gradient
! equal to D=slpD
! BC is du/dr=D u/r     in the form   (ub-u0)=  D*(ub+u0)*f/(1-f)
! where f = Sum_j[(xb_j+x0_j)dx_j]/(2*rm^2), dx=xb-x0
! Thus ub=u0(1-f-D.f)/(1-f+D.f) using radii from position (0,0,..)
      integer ipoint,inc
      integer indi(mdims),iused(mdims)
      real u(*)

! Structure vector needed for finding adjacent u values.
      integer iLs(mdims+1)
      include 'ndimsdecl.f'      
      include 'meshcom.f'
      include 'slpcom.f'
      D=slpD
!      write(*,*)'bdyslope',inc,ipoint,indi,iused
! Algorithm: take steps of 1 in all cases except when on a lower
! boundary face of dimension 1 (and not other faces).  There the step is
! iused(1)-1.
!------------------------------------------------------------------
! Between here and ^^^ is boundary setting. Adjust upper and lower.
      r2=0.
      fac=0.
      ipd=ipoint
      do n=1,ndims
! BC is du/dr=D u/r     in the form   (ub-u0)=  D*(ub+u0)*f/(1-f)
! where f = Sum_j[(xb_j+x0_j)dx_j]/(2*rm^2), dx=xb-x0
! Thus ub=u0(1-f-D.f)/(1-f+D.f)
! Here we are using radii from position (0,0,..)
         x=xn(ixnp(n)+1+indi(n))
         r2=r2+x*x
         if(indi(n).eq.0)then
! On lower boundary face
            dx=xn(ixnp(n)+1)-xn(ixnp(n)+2)
            ipd=ipd+iLs(n)
            fac=fac+x*dx
            if(n.eq.1)then
! The exception in step. Do not change!
               inc=iused(1)-1
            else
               inc=1
            endif
         elseif(indi(n).eq.iused(n)-1)then
! On upper boundary face
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
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!      write(*,*)'indi,inc,iused,ipoint',indi,inc,iused,ipoint
      end
!************************************************************************
      subroutine bdymach(inc,ipoint,indi,mdims,iLs,iused,u)
! Version of bdyroutine that sets the BC on the x, y boundaries
! as being  du/dr + M du/dz =0. (r the cylindrical radius)
! The z-mach number, M, is the value of slpD in slpcom.
! Default z-boundary conditions: du/dz=0
! islp indicates other BC choices as follows:
!    Bit-3  (8): set BC at lower-z boundary u=0.
!    Bit-13 (8192): set BC at upper-z boundary u''= -k^2u, k=1/M.lambdaD
! dxk2 is in slpcom.
      integer ipoint,inc
      integer indi(mdims),iused(mdims)
      real u(*)
! Structure vector needed for finding adjacent u values.
      integer iLs(mdims+1)
! Value of mach number passed in common.
      include 'ndimsdecl.f'
      include 'meshcom.f'
      include 'slpcom.f'
      DM=slpD
! Algorithm: take steps of 1 in all cases except when on a lower
! boundary face of dimension 1 (and not other faces).  There the step is
! iused(1)-1.
!------------------------------------------------------------------
! Between here and ^^^ is boundary setting. Adjust upper and lower.
      r2=0.
      fac=0.
      ipd=ipoint
! The dimensions up to the last are treated the same.
      do n=1,ndims-1
! Here we are using radii from position (0,0,..)
         x=xn(ixnp(n)+1+indi(n))
         r2=r2+x*x
         if(indi(n).eq.0)then
! On lower boundary face
            dx=xn(ixnp(n)+1)-xn(ixnp(n)+2)
            ipd=ipd+iLs(n)
            fac=fac+x*dx
            if(n.eq.1)then
! The exception in step. Do not change!
               inc=iused(1)-1
            else
               inc=1
            endif
         elseif(indi(n).eq.iused(n)-1)then
! On upper boundary face
            dx=xn(ixnp(n)+1+indi(n))-xn(ixnp(n)+indi(n))
            ipd=ipd-iLs(n)
            fac=fac+x*dx
            inc=1
         endif
      enddo
! Finally we do the last dimension, z.
      n=ndims
      if(ipd.eq.ipoint)then
! We only are here if we are on z-faces, and not on x or y faces. 
! Fall back to du/dz=0.
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
! fac contains r.dr, r2 contains r^2
! u_p=u_d + (-M*(dx.x + dy.y)/r + dz).(du/dz)
         if(indi(n).eq.0)then
! On lower z-boundary face (as well as x or y)
            dz=xn(ixnp(n)+1)-xn(ixnp(n)+2)
            ipd=ipd+iLs(n)            
            inc=1
! z-derivative calculated from value 3
            dubydz=(u(ipd+1)-u(ipd+iLs(n)+1))
     $           /(xn(ixnp(n)+2)-xn(ixnp(n)+3))
         elseif(indi(n).eq.iused(n)-1)then
! On upper z-boundary face (as well as x or y)
            dz=xn(ixnp(n)+1+indi(n))-xn(ixnp(n)+indi(n))
            ipd=ipd-iLs(n)
            inc=1
! z-derivative calculated from value iused-3.
            dubydz=(u(ipd+1)-u(ipd-iLs(n)+1))
     $           /(xn(ixnp(n)+indi(n))-xn(ixnp(n)+indi(n)-1))
         else
! Not on z-face
            dz=0.
            dubydz=(u(ipd+1+iLs(n))-u(ipd+1-iLs(n)))
     $           /(xn(ixnp(n)+indi(n)+2)-xn(ixnp(n)+indi(n)))
         endif
! Set u-value by extrapolation using dubydz. 
! Since the radial(xy) derivative is equal to -M*du/dz,
! -M*(fac/r)*dubydz is the radial difference, and add z-difference.
         u(ipoint+1)=u(ipd+1)+
     $        (-DM*fac/sqrt(r2) + dz)*dubydz
      endif
! Special cases:
      if(ibits(islp,3,1).ne.0)then
! Bit 3(+1) (=8) set, put lower z-face to zero.
         if(indi(n).eq.0)u(ipoint+1)=0.
      endif
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!      write(*,*)'indi,inc,iused,ipoint',indi,inc,iused,ipoint
      end
!**********************************************************************
! The logic of the following might miss some corners.
!************************************************************************
      subroutine bdyslopescreen(inc,ipoint,indi,mdims,iLs,iused,u)
! Version of bdyroutine that sets logarithmic 'radial' gradient
! equal to that for a Debye screened potential: -(1+r/lambdas)
! slpD needs to be set to lamdas prior to call.
      integer ipoint,inc
      integer indi(mdims),iused(mdims)
      real u(*)

! Structure vector needed for finding adjacent u values.
      integer iLs(mdims+1)
      include 'ndimsdecl.f'
      include 'meshcom.f'
      include 'slpcom.f'
      real x(ndims)

      r2=r2indi(ndims,indi,x)
      r1=sqrt(r2)
      D=slpD

! Algorithm: take steps of 1 in all cases except
! when on a lower boundary face of dimension 1. 
! There the step is iused(1)-1.
      inc=1
      do n=ndims,1,-1
!------------------------------------------------------------------
! Between here and ^^^ is boundary setting. Adjust upper and lower.
         if(indi(n).eq.0)then
! The exception in step. Do not change!:
            if(n.eq.1)inc=iused(1)-1
! On lower boundary face
! du/dn= -(1.+r/D) r.n u/ r^2
            dx=xn(ixnp(n)+2)-xn(ixnp(n)+1)
            fac=(1.+0.5*dx/x(n))*r2*2./(-(1.+r1/D)*x(n)*dx)
!            write(*,*)n,dx,xn(ixnp(n)+1),fac,r2
! Centered BC:
! (u2-u1)/dx=(u2+u1)/2 *S*(x1+x2)/2 /rm^2  put fac=2*rm^2/S xm dx
!  u2(fac-1.)=u1(fac+1)
            u(ipoint+1)=u(ipoint+1+iLs(n))*(fac-1.)/(1.+fac)
            goto 101
         elseif(indi(n).eq.iused(n)-1)then
! On upper boundary face
! du/dn=-(1.+r/D) r.n u/ r^2
            dx=-(xn(ixnp(n)+1+indi(n))-xn(ixnp(n)+indi(n)))
            fac=(1.+0.5*dx/x(n))*r2*2./(-(1.+r1/D)*x(n)*dx)
            u(ipoint+1)=u(ipoint+1-iLs(n))*(fac-1.)/(1.+fac)
            goto 101
         endif
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      enddo
      write(*,*)'BDY Error. We should not be here',n,ipoint,indi
      stop
 101  continue
!      write(*,*)'indi,inc,iused,ipoint',indi,inc,iused,ipoint
      end
!**********************************************************************
! Return the radius squared and the vector position of point indexed as
! indi(ndims), relative to the center of the coordinate ranges,
! using information in meshcom.
      function r2indi(mdims,indi,x)
      integer indi(mdims)
      real x(mdims)
      include 'ndimsdecl.f'
      include 'meshcom.f'
      if(ndims.ne.ndims)then
         write(*,*)'rindi dimension mismatch',ndims,ndims
         stop
      endif
      r2indi=0.
      do i=1,ndims
         x(i)=xn(ixnp(i)+indi(i)+1)-(xn(ixnp(i)+1)+xn(ixnp(i+1)))/2.
         r2indi=r2indi+x(i)**2
      enddo
      end
!********************************************************************
      subroutine bdyface(inc,ipoint,indi,mdims,iLs,iused,u)
! Boundary routine for general setting in terms of face.
! Coefficients (2ndims) in facebcom. AF phi + BF dphi/dn + CF =0
! CF may vary linearly with position with coefs CxyzF(3,6)
! If B=0, A is assumed =1. Otherwise precalculated coeffs are used
! AmBF= (A/2-B/dn), ApBF= (A/2+B/dn) where dn is the outward mesh step. 
      integer ipoint,inc
      integer indi(mdims),iused(mdims)
      real u(*)
      include 'ndimsdecl.f'
      include 'meshcom.f'
      include 'facebcom.f'
! Structure vector needed for finding adjacent u values.
      integer iLs(ndims+1)

      integer iupper,idn,id
! Algorithm: take steps of 1 in all cases except
! when on a lower boundary face of dimension 1. 
! There the step is iused(1)-1.
      inc=1
      do n=ndims,1,-1
!------------------------------------------------------------------
! Between here and ^^^ is boundary setting. Adjust upper and lower.
         if(indi(n).eq.0)then
! The exception in step. Do not change!:
            if(n.eq.1)inc=iused(1)-1
! On lower boundary face 1,2,3
            iupper=-1
            idn=n
         elseif(indi(n).eq.iused(n)-1)then
! On upper boundary face 4,5,6
            iupper=1
            idn=n+ndims
         else
! Not on face in this dimension.
            goto 102
         endif
         if(LPF(n))then
! Here we set the boundary because otherwise the final array won't have
! it set. But we ought not to during iterations because it's unnecessary.
! During iterations we use the new bdyshare routine.
            u(ipoint+1)=u(ipoint+1-iupper*(iused(n)-2)*iLs(n))
         else
! Only if we are not on a periodic face:
            if(LCF(idn))then
! Variable C. Calculate:
               C=C0F(idn)
               do id=1,ndims
                  C=C+CxyzF(id,idn)*xn(ixnp(id)+indi(id)+1)
               enddo
            else
! Simple short cut.
               C=C0F(idn)
            endif
            if(BF(idn).eq.0)then
! Assume A=1.
               u(ipoint+1)=-C
            else
               u(ipoint+1)=-(C+AmBF(idn)*u(ipoint+1-iupper*iLs(n)))
     $              /ApBF(idn)
            endif
         endif
         goto 101
 102     continue
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      enddo
      write(*,*)'BDYface Error. We should not be here',n,ipoint,indi
      stop
 101  continue
!      write(*,*)'indi,inc,iused,ipoint',indi,inc,iused,ipoint
      end
!********************************************************************
      subroutine bdyfaceinit(idn,CFpin)
! Initialize (if idn<1). Turn off this type of BC if idn=0.
! Or update the face boundary conditions if 0<idn<=ndims.
! For index idn = 1,2,3 (lower) 4,5,6 (upper) face.
! You can't literally equivalence them but 
!      real Ain,Bin,C0in,Cxyzin(ndims) == CFpin(6)
      integer idn
      include 'ndimsdecl.f'
      include 'meshcom.f'
      real CFpin(3+ndims)
      include 'facebcom.f'
      integer id

      if(idn.le.-2)then
! Print out settings:
         write(*,*)' AF      BF       C0F     CxyzF'
     $        ,'                  AmBF     ApBF  LF LCF LPF'
         do id=1,2*ndims
            write(*,'(8f8.4,7L3)')AF(id),BF(id),C0F(id) ,(CxyzF(ii,id)
     $           ,ii=1,3),AmBF(id),ApBf(id),LF,LCF(id)
     $           ,LPF(mod(id-1,3)+1)
         enddo
      elseif(idn.le.0)then
! Initialize A=-1, B=0, C=0
         do id=1,2*ndims
            AF(id)=1.
            BF(id)=0.
            C0F(id)=0.
! Uniform C:
            LCF(id)=.false.
! And switch off this type of BC:
            if(idn.eq.0)LF=.false.
         enddo
      elseif(idn.le.2*ndims)then
         LF=.true.
         AF(idn)=CFpin(1)
         BF(idn)=CFpin(2)
         if(BF(idn).eq.0.and.AF(idn).eq.0)AF(idn)=1.
         C0F(idn)=CFpin(3)
         do id=1,ndims
            CxyzF(id,idn)=CFpin(3+id)
            if(CxyzF(id,idn).ne.0)LCF(idn)=.true.
         enddo
         if(BF(idn).ne.0.)then
! Calculate the extra coefficients.
            if(idn.gt.ndims)then
               id=idn-ndims
               dn=xn(ixnp(id+1))-xn(ixnp(id+1)-1)
!               write(*,*)'FACEPOS Upper',dn,xn(ixnp(id+1)),xn(ixnp(id+1)
!     $              -1)
            else
               id=idn
               dn=xn(ixnp(id)+1)-xn(ixnp(id)+2)
!               write(*,*)'FACEPOS Lower',dn,xn(ixnp(id)+1),xn(ixnp(id)
!     $              +2)
            endif
! Outward facing normal differential: dn is obtained by taking abs()
            dn=abs(dn)
            if(dn.eq.0)stop 'bdyface init ERROR. dn=0'
! phi_b = -[C+(A/2-B/dn)phi_i]/[A/2+B/dn]
            AmBF(idn)=AF(idn)/2.-BF(idn)/dn
            ApBF(idn)=AF(idn)/2.+BF(idn)/dn
            if(ApBF(idn).eq.0.)then
               write(*,*)'bdyface coefficient singularity.'
     $              ,idn,dn
               write(*,*)'There''s a problem with those values.'
     $              ,' They must be changed.'
               stop
            endif
         endif
      else
         write(*,*)'bdyfaceinit ERROR. Incorrect index',idn
      endif
!      write(*,*)'BDYFACEINIT',idn,CFpin
      end

!********************************************************************
!********************************************************************
! Obsolete or Currently unused routines:
!**********************************************************************
      subroutine bdyset3sl(ndims,ifull,iuds,cij,u,q)
      integer ndims,ifull(ndims),iuds(ndims)
      real cij(*),u(*),q(*)
! Specify external the boundary setting routine.
      external bdy3slope 
! sets the derivative to zero on boundaries 3.
! Silence warnings about unused arguments
      use=cij(1)
      use=q(1)
      ipoint=0
      call mditerarg(bdy3slope,ndims,ifull,iuds,ipoint,u)
      end
!************************************************************************
      subroutine bdy3slope(inc,ipoint,indi,ndims,iLs,iused,u)
! Version of bdyroutine that sets derivative=0 on 3-boundary.
      integer ipoint,inc
      integer indi(ndims),iused(ndims)
      real u(*)

! Structure vector needed for finding adjacent u values.
      integer iLs(ndims+1)

! Algorithm: take steps of 1 in all cases except
! when on a lower boundary face of dimension 1. 
! There the step is iused(1)-1.
      inc=1
      do n=ndims,1,-1
!------------------------------------------------------------------
! Between here and ^^^ is boundary setting. Adjust upper and lower.
         if(indi(n).eq.0)then
! The exception in step. Do not change!:
            if(n.eq.1)inc=iused(1)-1
! On lower boundary face
            u(ipoint+1)=0.
            if(n.eq.3)then
! First derivative is zero:
               u(ipoint+1)=u(ipoint+1+iLs(n))
            endif
! Second derivative is zero:
!               u(ipoint+1)=2.*u(ipoint+1+iLs(n))-u(ipoint+1+2*iLs(n))
            goto 101
         elseif(indi(n).eq.iused(n)-1)then
! On upper boundary face
            u(ipoint+1)=0.
            if(n.eq.3) u(ipoint+1)=u(ipoint+1-iLs(n))
            goto 101
         endif
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      enddo
      write(*,*)'BDY3s Error. We should not be here',n,ipoint,indi
      stop
 101  continue
!      write(*,*)'indi,inc,iused,ipoint',indi,inc,iused,ipoint
      end
!**********************************************************************
      subroutine bdyslope0(ndims,ifull,iuds,cij,u,q)
! Set the boundary normal slope to zero.
      integer ndims,ifull(ndims),iuds(ndims)
      real cij(*),u(*),q(*)
! Specify external the boundary setting routine.
      external bdyslopeDh,bdyslopescreen,bdymach
      include 'slpcom.f'
! Silence warnings about unused arguments
      use=cij(1)
      use=q(1)
      ipoint=0
      islp=0
      slpD=0.
      call mditerarg(bdyslopeDh,ndims,ifull,iuds,ipoint,u)
      end
!**********************************************************************
      subroutine bdysetnull()
!   (ndims,ifull,iuds,cij,u,q)
!     Null version
      end
!************************************************************************

!********************************************************************
