!**********************************************************************
      subroutine diaginject(x)
! Accumulate diagnostics of reinjections, etc.
! Position and velocity of reinjection.
      real x(*)
! Potential of reinjection
!      real phi
! Number of launches needed to get this reinjection
!      integer ilaunch
! Diagnostic storage etc. Here we are assuming a 3-D problem which is
! required by this reinjection scheme. 2-D position on reinjection
! surface.      
      parameter (pi=3.1415927)
      parameter (eps=.0001)
      include 'reincom.f'
      logical lfirst
      data lfirst/.true./
      save
!----------------------------
! Initializations:
      if(lfirst)then
! Set up cell arrays (for plotting) and zero counts.
         do i=1,ndth
! cos(theta) center cells (use nint to choose) (-1 - +1)
            reincth(i)=-1.+((i-1)*2.+1.)/ndth
            cthtot(i)=0.
         enddo
         do i=1,ndpsi
! psi angle center cells (-pi - +pi)
            reinpsi(i)=pi*(-1.+((i-1)*2.+1.)/ndpsi)
            do j=1,ndth
               reinpos(j,i)=0.
            enddo
            psitot(i)=0.
         enddo
! fix any rounding errors at top and bottom of range:
         dcth=2.*(1.+eps)/ndth
         dpsi=2.*(1.+eps)*pi/ndpsi
         lfirst=.false.
! Velocity
         vrange=6.
         dvf=2.*vrange/ndth
         do i=1,ndth
            vfv(i)=(-vrange+((i-1)+.5)*dvf)
            do j=1,3
               fv(i,j)=0.
            enddo
            sv(i)=0.
            vs(i)=vrange*(i-.5)/ndth
         enddo
      endif
!------------------------------
! Particle position
      r=sqrt(x(1)**2+x(2)**2+x(3)**2)
      cth=x(3)/r
      psi=atan2(x(2),x(1))
! Find array position.
      icth=nint((cth+(1.+eps))/dcth+.5)
      ipsi=nint((psi+(1.+eps)*pi)/dpsi+.5)

      if(icth.lt.1 .or. icth.gt.ndth
     $     .or. ipsi.lt.1 .or. ipsi.gt.ndpsi)then
         write(*,*)'diaginj error',icth,ipsi,cth,psi,r,(x(j),j=1,3)
         stop
      endif
! Increment counts.
      reinpos(icth,ipsi)=reinpos(icth,ipsi)+1
      cthtot(icth)=cthtot(icth)+1
      psitot(ipsi)=psitot(ipsi)+1
!      write(*,*)icth,ipsi,cthtot(icth),psitot(ipsi)

! Particle velocity
      sp=0.
      do j=1,3
         iv=nint((x(3+j)+vrange)/dvf+.5)
         if(iv.lt.1)iv=1
         if(iv.gt.ndth)iv=ndth
         fv(iv,j)=fv(iv,j)+1
         sp=sp+x(3+j)**2
      enddo
      sp=sqrt(sp)
      iv=nint(ndth*sp/vrange+0.5000001)
      if(iv.gt.ndth)then 
!         write(*,*)'Reindiag Velocity range exceeded.',sp,iv
         iv=ndth
      elseif(iv.lt.0)then
!         write(*,*)'Reindiag Velocity range exceeded.',sp,iv
         iv=0
      endif
      sv(iv)=sv(iv)+1
      end
!**********************************************************************
      subroutine plotinject(Ti)
! Diagnostic storage etc. Here we are assuming a 3-D problem which is
! required by this reinjection scheme. 2-D position on reinjection
! surface.      

      parameter (pi=3.1415927)
      include 'reincom.f'

      real gaussian(ndth)
! Plot some diagnostic information from reinjection.

!      write(*,*)cthtot
!      write(*,*)psitot
      call automark(reincth,cthtot,ndth,1)
      call axlabels('cos(!Aq!@)','count')
      call pltend()
      call automark(reinpsi,psitot,ndth,1)
      call axlabels('angle psi','count')
      call pltend()

      counts=0.
      do i=1,ndth
         counts=counts+cthtot(i)
      enddo

! Actually the flux is not a gaussian. So this isn't yet useful
      do i=1,ndth
         gaussian(i)=exp(-vfv(i)**2/(2*Ti))*counts*
     $        (vfv(2)-vfv(1))/sqrt(2.*pi*Ti)
      enddo
      write(*,*)Ti,gaussian(1)

      call multiframe(3,1,0)
      call automark(vfv,fv(1,1),ndth,1)
      call boxtitle('!A1!@ x, !A2!@ y, !A3!@ z')
      call axlabels(' ','x count')
      call automark(vfv,fv(1,2),ndth,2)
      call axlabels(' ','y count')
      call automark(vfv,fv(1,3),ndth,3)
!      call polymark(vfv,fv(1,2),ndth,2)
!      call polymark(vfv,fv(1,3),ndth,3)
!      call polyline(vfv,gaussian,ndth)
      call axlabels('velocity distributions','z count')
      call pltend()
      call multiframe(0,0,0)

      call automark(vs,sv,ndth,1)
      call axlabels('speed distribution','count')
      call pltend()

      end
