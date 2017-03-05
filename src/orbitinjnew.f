! This version of injection is not ready for multiple species.
! It also does not do fractional injection numbers.

      subroutine reinject(xr,ilaunch)
      real xr(*)
      integer ilaunch
      include 'reincom.f'
      logical lfirst
      data lfirst/.true./
! This call passes address of the particle slot to reiject.      
      dt=1.
      ibcr=0
      icolntype=0
      if(lfirst)then
!         write(*,*)'Initializing reinjection'
         call injinit(icolntype,ibcr)
         lfirst=.false.
      endif
      call thisreinject(1,dt,icolntype,ibcr,xr,ilaunch)
!      write(*,*)'Returning from reinject',(xr(j),j=1,6)
      end
!***********************************************************************
! General version allows choice of reinjection scheme.
!***********************************************************************
      subroutine thisreinject(i,dt,icolntype,bcr,xp,ilaunch)
      implicit none
      integer i,icolntype,ilaunch
      real dt
      integer bcr 
      real xp(6,*)
      if(bcr.ne.0) then
!         call maxreinject(i,dt,bcr)
      elseif(icolntype.eq.1.or.icolntype.eq.5) then
! Injection from fv distribution at the boundary.
!         call fvreinject(i,dt,icolntype)
      elseif(icolntype.eq.2.or.icolntype.eq.6)then
! Injection from a general gyrotropic distribution at infinity
!         call ogenreinject(i,dt)
      else
! Injection from a shifted maxwellian at infinity
         call oreinject(i,dt,xp,ilaunch)
      endif
      end
!***********************************************************************
      subroutine injinit(icolntype,bcr)
      implicit none
      integer icolntype
      integer bcr

      if(bcr.ne.0) then
! Injection from a maxwellian at boundary?
!         call maxinjinit(bcr)
      elseif(icolntype.eq.1.or.icolntype.eq.5) then
! Injection from fv distribution at the boundary.
!         call fvinjinit(icolntype)
      elseif(icolntype.eq.2.or.icolntype.eq.6)then
! Injection from a general gyrotropic distribution at infinity
!         call ogeninjinit(icolntype)
      else
! Injection from a shifted maxwellian at infinity
         call oinjinit()
      endif
      end
!***********************************************************************
!***********************************************************************
! Other versions are in other source files.
      subroutine oreinject(i,dt,xp,ilaunch)
      implicit none
      integer i
      real dt
      real xp(6,*)
      integer ilaunch
      include 'ndimsdecl.f'
      include 'plascom.f'
      include 'reincom.f'
!      real averein,debyelen,pi,Ti,vd
!      real Vcom
!      integer nr
!      integer nrfull,nrused
!      integer nthsize
!      parameter (nthsize=2)
!      integer nrsize
!      parameter (nrsize=2)
!      integer nvel,ntrapre
!      parameter (nvel=2)
!      real r(nrsize)
!      real th(nthsize),phi(1,nthsize)
!      real getphihere
! End of problematic variables (fix).

! Common data:
      real eup
      parameter (eup=1.e-10)
      external pu
!      logical istrapped
! Testing
      real vdist(nvel)
      real tdist(nthsize)
      real crdist(nthsize),cidist(nthsize)
      common/rtest/crdist,cidist,tdist,vdist

! Making all variables declared:
      real a1,brc,brcsq,ceta,chium2,cosal,sinal,crt,czt,dv,eta
      real expuu2,p2,rcyl,rp,seta,srt,szt,u,ua,ua1,ua3
      real Uc,uu2,vscale,y,y1,zt,alcos,alsin
      integer icr,idum,ierr,iv
      real rri

! In this routine we work in velocity units relative to ion thermal till end.
      vscale=sqrt(2.*Ti)
      idum=1
      ilaunch=0
! Silence warnings.
      brcsq=0.
      ierr=0
      rp=0.
 1    continue
      ilaunch=ilaunch+1
      if(ilaunch.gt.100)then
         write(*,*)'ilaunch excessive. averein=',averein,' brcsq=',
     $        brcsq,' ierr=',ierr,' rp=',rp
         stop
      endif
! Pick normal velocity from cumulative Pu
      call ranlux(y1,1)
      call finvtfunc(pu,nvel,y1,u)
      iv=int(u)
      dv=u-iv
      u=dv*Vcom(iv+1)+(1.-dv)*Vcom(iv)
      if(.not.dv.le.1)write(*,*)'Error in u calculation',y1,u,iv,dv
      vdist(iv)=vdist(iv)+1.
! Pick angle from cumulative Pc.
      call ranlux(y,1)
! Here the drift velocity is scaled to the ion temperature.
      Uc=vd/vscale
      uu2=2.*Uc*u
      if(uu2.gt.50.) then
         crt=1.+alog(y)/uu2
      elseif(uu2.lt.-50.) then
         crt=-1.+alog(1-y)/uu2
      elseif(abs(uu2).lt.1.e-5)then
         crt=2.*y -1.
      else
         expuu2=exp(uu2)
! This expression is evaluated very inaccurately if expuu2 is nearly 1.
! That is why such a large limit on abs(uu2) is adopted.
         crt=alog(y*expuu2 + (1-y)/expuu2)/uu2
! The following do not do any better at solving this problem.
!         crt=alog( (y*expuu2 + (1-y)/expuu2)**(1./uu2))
!         crt=-1.+alog(1.+(expuu2**2-1.)*y)/uu2
      endif
      if(.not. abs(crt).le.1)then
!         write(*,*)'Strange crt:',crt,y,expuu2,uu2
! It seems impossible to avoid occasional strange results when uu2 is small.
         crt=2.*y-1.
      endif
! Testing angular distribution.
      icr=int((1.+crt)*0.5*((nthsize-1)-1) + 1.5)
      crdist(icr)=crdist(icr)+1.
! End of testing distribution monitor.
      srt=sqrt(1.- crt**2)
! Pick angle zt of poloidal impact and angle eta of impact parameter
      call ranlux(zt,1)
      zt=zt*2.*pi
      czt=cos(zt)
      szt=sin(zt)
      call ranlux(eta,1)
      eta=eta*2.*pi
      ceta=cos(eta)
      seta=sin(eta)
! Choose impact parameter, preventing overflow.
      chium2=-averein/Ti/(u+eup)**2
      if(chium2.le.-1.) then
!         write(*,*)'Impossible chium2=',chium2,' averein=', averein,
!     $        ' u=',u,' iv=',iv
         goto 1
!         stop
      endif
!      if(.not.lfixedn)chium2=0.
      call ranlux(brcsq,1)
      brcsq=brcsq*(1.+chium2)
! Reject a particle that will not reach boundary.
      if(brcsq.lt.0.) then
         goto 1
      endif
      brc=sqrt(brcsq)
! Get cosine and sine of impact angle relative to distant position.
! Based on integration.
      p2=brcsq*2.*Ti*u**2
      ierr=0
      if(debyelen.gt..001)then
! Orbit integration angle calculation.
! There is an overflow with this at zero debyelen. Ought to be properly fixed.
         call alphaint(p2,brcsq,cosal,sinal,ierr)
         if(ierr.ne.0)goto 1
!      write(*,'(4f9.4)')cosal-alcos(brc,chium2),sinal-alsin(brc,chium2)
!     Now ilaunch is the number of launches at infinity it took to get
!     one that reached the boundary.
      else
! Alternative based on analytic orbit calculation.
! Used for low debyelen, but really assumes negligible boundary potential.
!         call alcossin(brc,chium2,cosal,sinal)
         cosal=alcos(brc,chium2)
         sinal=alsin(brc,chium2)
      endif
! Install reinjection position
      a1=crt*ceta*sinal+srt*cosal
! Ensure that we are just inside the rs sphere
      rri=0.99999*rs
      xp(1,i)=rri*(czt*a1+szt*seta*sinal)
      xp(2,i)=rri*(-szt*a1+czt*seta*sinal)
      xp(3,i)=rri*(-srt*ceta*sinal + crt*cosal)
! Injection velocity components normalized in the rotated frame:
      ua1=-brc*cosal -sqrt(1.+chium2-brcsq)*sinal
      ua3=brc*sinal - sqrt(1.+chium2-brcsq)*cosal
      ua=crt*ceta*ua1+srt*ua3
! Install reinjection velocity in Te-scaled units
      u=u*vscale
      xp(4,i)=u*(czt*ua+szt*seta*ua1)
      xp(5,i)=u*(-szt*ua+czt*seta*ua1)
      xp(6,i)=u*(-srt*ceta*ua1 + crt*ua3)
! If we don't recalculate rp, then we don't trap NANs in the random choices.
      rcyl=xp(1,i)**2+xp(2,i)**2
      rp=rcyl+xp(3,i)**2
!      rp=rs
!      write(*,*)'oreinject',rp
! Reject particles that are already outside. But this is more to do
! with detecting NaNs.
      if(.not.rp.lt.rs*rs)then
!      if(.not.rp.le.r(nr)*r(nr))then
         write(*,'(a,4f8.4)')'Relaunch rp,xp(1,i),xp(2,i),xp(3,i)',
     $        rp,xp(1,i),xp(2,i),xp(3,i)
         write(*,'(a,7f8.4)')'vscale,u,brc,chium2,brcsq,sinal,cosal',
     $        vscale,u,brc,chium2,brcsq,sinal,cosal
         goto 1
      endif
      end
!********************************************************************
! Initialize the distributions describing reinjected particles
      subroutine oinjinit()
      implicit none
! Common data:
      include 'ndimsdecl.f'
      include 'plascom.f'
      include 'reincom.f'

! To fix exclusion we need
!      real Vcom(1),pu1(1),pu2(1)
!      real Ti,vd,pi
!      integer myid,nvel


! Making variables declared
      integer i
      real expu0,u0,Uc,uminus,uplus,vspread,erfcc

! Here the drift velocity is scaled to the ion temperature.
! And U's are in units of sqrt(2T/m), unlike vd.
      Uc=abs(vd)/sqrt(2.*Ti)
! Range of velocities (times (Ti/m_i)^(1/2)) permitted for injection.
      vspread=5.+abs(Uc)

! Can't use these formulas for Uc exactly equal to zero.
      if(abs(Uc).lt.1.e-4)then
         if(Uc.lt.1.e-20) Uc=1.e-20
         do i=1,nvel
            u0= vspread*(i-1.)/(nvel-1.)
            Vcom(i)=u0
            expu0=exp(-u0**2)
            pu2(i)=2.*Uc*expu0
            pu1(i)=0.5*(4.*u0**2*Uc + 2.*Uc)*expu0
     $           +(Uc**2 +0.5)*pu2(i)
         enddo
      else
         do i=1,nvel
            u0= vspread*(i-1.)/(nvel-1.)
            Vcom(i)=u0
            uplus=u0+Uc
            uminus=u0-Uc
            pu2(i)=0.5*sqrt(pi)*(erfcc(uminus)-erfcc(uplus))
            pu1(i)=0.5*(-uminus*exp(-uplus**2)+uplus*exp(-uminus**2))
     $           +(Uc**2 +0.5)*pu2(i)
         enddo
      endif
      end
!***********************************************************************
!***********************************************************************
! Calculate the cumulative probability for velocity index iu such that
!         u= vspread*(iu-1.)/(nvel-1.)   as per injinit
      real function pu(iu)
      implicit none
      real pudenom
      integer iu
!     averein is the average potential of reinjected particles, which is
!     used as an estimate of the potential at the reinjection boundary.
!     It is expressed in units of Te so needs to be scaled to Ti.
      include 'ndimsdecl.f'
      include 'plascom.f'
      include 'reincom.f'

! needed from common:
!      real averein,Ti,pu1(1),pu2(1)

      pudenom=pu1(1)-pu2(1)*averein/Ti
      pu=1.- (pu1(iu)-pu2(iu)*averein/Ti)/pudenom
      end
!********************************************************************
! Given a monotonic (increasing?) 
! function Q(x) on a 1-D grid x=1..nq, solve Q(x)=y for x.
! That is, invert Q to give x=Q^-1(y).
      subroutine finvtfunc(Q,nq,y,x)
! Somehow this breaks the passing of a function reference.
      implicit none
!      real external Q
      real Q
      external Q
      integer nq
      real y,x
!
      integer iqr,iql,iqx
      real Qx,Qr,Ql,Qd
      Ql=Q(1)
      Qr=Q(nq)
      iql=1
      iqr=nq
      if((y-Ql)*(y-Qr).gt.0.) then
! Value is outside the range.
         x=0
         return
      endif
 200  if(iqr-iql.eq.1)goto 210
      iqx=(iqr+iql)/2
      Qx=Q(iqx)
!      write(*,*)y,Ql,Qx,Qr,iql,iqr
! Formerly .lt. which is an error.
      if((Qx-y)*(Qr-y).le.0.) then
         Ql=Qx
         iql=iqx
      else
         Qr=Qx
         iqr=iqx
      endif
      goto 200
 210  continue
! Now iql and iqr, Ql and Qr bracket Q
!      x=(y-Ql)/(Qr-Ql)+iql
! Trap errors caused by flat sections.
      Qd=Qr-Ql
      if(Qd.eq.0.)then
         x=(iql+iqr)/2.
      else
         x=(y-Ql)/Qd+iql
      endif
      end
!**********************************************************************
! Inverse square law (phi\propto 1/r) injection functions:
!**********************************************************************
      subroutine alcossin(s,c,cosal,sinal)
      real s,c,cosal,sinal,alcos,alsin
      cosal=alcos(s,c)
      sinal=alsin(s,c)
      end
!**********************************************************************
      real function alcos(s,c)
      real s,c,r
      if(abs(s).le.1.e-12*abs(c))then
         alcos=-1.
         return
      else
         r=c/(2.*s)
         alcos=-(sqrt(1.+c-s**2)-(s-r)*r)/(1+r**2)
      endif
      end
!**********************************************************************
      real function alsin(s,c)
      real s,c,r
      if(s.le.1.e-12*c)then
         alsin=0.
         return
      else
         r=c/(2.*s)
         alsin=(sqrt(1.+c-s**2)*r+(s-r))/(1+r**2)
      endif
      end
!**********************************************************************
!**********************************************************************
! Return the angle cosine and sine for impact of a particle calculated
! by integrating the angle formula for a general central force.
! We integrate in the variable q=1/r from 0 to 1, provided that we don't
! encounter a barrier. If we do encounter one, we return an error.
!
! This version designed to work with uneven spacing if
! necessary, and trapezoidal integration. And uses init2ext
      subroutine alphaint(p2,b2,cosal,sinal,ierr)
      implicit none
!     The angular momentum and impact parameter squared
!     in units such that injection radius is 1.
      real p2, b2
!     The angle values returned.
      real cosal, sinal
!     The error return signal if non-zero.
      integer ierr
! Common data:
      include 'ndimsdecl.f'
      include 'plascom.f'
! Not needed:
      include 'reincom.f'

!     Need real adeficit,averein
!     Need debyelen


      integer iqsteps
! This choice ensures iqsteps is large enough to accommodate a profile
! read in for processing with orbitint, but might not be the best choice for
! regular use.
!      parameter (iqsteps=nrsize+1)
      parameter (iqsteps=100+1)
      real phibye(iqsteps),phiei(iqsteps)
      real qp(iqsteps)
!      real pp(iqsteps)
      real alpha,b2i,d1,d2,p2i2,sa,xlambda
      integer i,iqs
      logical uninitialized
      data uninitialized/.true./
      save
!      real extpot
! Statement function
! Inverse square law potential.
!      extpot(q)=averein*q
! Not currently in use.

! First time through, do the initialization.
      if(uninitialized)then
         iqs=iqsteps
         uninitialized=.false.
! Screening length accounting for both ions and electrons.
         xlambda=debyelen/sqrt(1.+1./Ti)
         if(xlambda.lt.1.e-6)xlambda=1.e-6

!$$$         if(diagrho(1).eq.999.)then
!$$$c Special case to use the read in data.
!$$$c Some common data used in abnormal ways.
!$$$            qp(1)=0.
!$$$            phibye(1)=0.
!$$$            if(ninner.gt.iqsteps)stop 'Alphaint Too many r-cells read.'
!$$$            do i=1,ninner
!$$$               qp(i+1)=1./rcc(ninner+1-i)
!$$$               phibye(i+1)=diagphi(ninner+1-i)
!$$$            enddo
!$$$            iqs=ninner+1
!$$$            if(qp(iqs).ne.1.)then
!$$$               iqs=iqs+1
!$$$               qp(iqs)=1.
!$$$c Extrapolate linearly in q.
!$$$               phibye(iqs)=phibye(iqs-1)+
!$$$     $              (phibye(iqs-1)-phibye(iqs-2))*
!$$$     $              (qp(iqs)-qp(iqs-1))/(qp(iqs-1)-qp(iqs-2))
!$$$            endif
!$$$c The read-in case just uses phibye not phiei. And scales phibye(iqs) to 1
!$$$c putting the absolute edge value into averein.
!$$$            averein=phibye(iqs)
!$$$            adeficit=0.
!$$$            do i=1,iqs
!$$$               phibye(i)=phibye(i)/phibye(iqs)
!$$$            enddo
!$$$c End of read-in potential special case.
!$$$         else
! Specify the q-array. It goes from 0 to 1.
            do i=1,iqs
               qp(i)=((i-1.)/(iqs-1.))
            enddo
            call initext(iqs,qp,phibye,phiei,rs,xlambda)
!$$$         endif
! Diagnostics
!         write(*,'(3f8.4)')(qp(j),phibye(j),phiei(j),j=1,iqs)
!         write(*,*)'averein,adeficit',averein,adeficit

!         if(averein.ne.0)then
!     When used by SCEPTIC averein is zero the first time, so this will
!     not be called.
! No longer true when restarting code is used.
!            do i=1,iqs
!               pp(i)=averein*phibye(i)-adeficit*phiei(i)
!            enddo
!            write(*,*)'adeficit=',adeficit
!            call autoplot(qp,pp,iqs)
!            call axlabels('q','potential')
!            call pltend()
!         endif
      endif
! End of initialization section.
! Do the integration for the orbit.
      ierr=0
      b2i=1./b2
      p2i2=2./p2
      sa=b2i - qp(1)**2 - p2i2*(averein*phibye(1)-adeficit*phiei(1))
      d1=(1./sqrt(sa))
!         write(*,*)d1,d2,b2i,p2i2,qp(1),phibye(1),phiei(1)
! Inverse square case.
!      d1=1./sqrt(b2i)
      alpha=0.
! Trapezoidal rule.
      do i=2,iqs
         d2=d1
         sa=b2i - qp(i)**2 - p2i2*(averein*phibye(i)-adeficit*phiei(i))
! Inverse square case.
!         sa=b2i - qp(i)**2 - p2i2*extpot(qp(i))
         if(sa .le. 0.) goto 2
         d1=(1./sqrt(sa))
         alpha=alpha+(qp(i)-qp(i-1))*(d1+d2)*.5
      enddo
!      write(*,*)'alpha=',alpha
! Negative sign for definition of alpha relative to the forward direction.
      cosal=-cos(alpha)
      sinal=sin(alpha)
      if(.not.sinal.le.1.)then
         write(*,*)'averein,adeficit,alpha',averein,adeficit,alpha
      endif
      return
 2    ierr=i
!      write(*,101)i,b2,p2,averein,adeficit
! 101  format('Barrier: i=',i3,' b2=',f8.1,
!     $     ' p2=',f8.4,' averein=',f8.4,' adeficit=',f8.4)
      end
!********************************************************************
      subroutine avereinset(aver)
! Set the averein value
      include 'reincom.f'
      averein=aver
      end
!********************************************************************
      subroutine adeficitset(adefic)
      include 'reincom.f'
      adeficit=adefic
      end
!*******************************************************************
      subroutine plascomset(d,T,v,r,p)
      include 'ndimsdecl.f'
      include 'plascom.f'
      debyelen=d
      Ti=T
      vd=v
      rs=r
      phip=p
      end
!*********************************************************************
      subroutine rhoinfcalc(dtin)
! Obtain the rhoinf to be used in calculating the electron shielding,
! based upon the number and average potential of the reinjections.
      include 'ndimsdecl.f'
      include 'plascom.f'
! No time-averaging for now.
! Use particle information for initializing.
      include 'partcom.f'

! Moved nrein and phirein reductions to psumreduce.
      if(nrein.ne.0)then
! Calculate rhoinf from nrein if there are enough.
         chi=max(-phirein/Ti,-0.5)
         riest=(nrein/dtin) /
     $        (sqrt(Ti)*
     $        smaxflux(vd/sqrt(2.*Ti),chi)
     $        *rs**2 )
         rhoinf=riest
!         write(*,*)'nrein,dtin,Ti,vd,phirein,chi,rs,numprocs=',
!     $        nrein,dtin,Ti,vd,phirein,chi,rs,numprocs
      else
         if(rhoinf.lt.1.e-4)then
! Approximate initialization
            rhoinf=numprocs*n_part/(4.*3.1415926*rs**3/3.)
            write(*,*)'Rhoinf in rhoinfcalc approximated as',rhoinf
     $           ,numprocs,n_part,rs
         endif
! Else just leave it alone.
      endif
!      write(*,*)'Ending rhoinfcalc',rhoinf,nrein,n_part
      end
!*********************************************************************
      subroutine ninjcalc(dtin)
! Given ripernode, decide the number of reinjections per step ninjcomp
! for average edge potential.
      include 'ndimsdecl.f'
      include 'plascom.f'
! No time-averaging for now.
! Particle information
      include 'partcom.f'

      pinjcompa(1)=0.
      if(n_part.ne.0)return
! Calculate ninjcomp from ripernode
      chi=min(-phirein/Ti,0.5)
      ninjcomp=nint(ripernode*dtin* 
     $        (sqrt(Ti)*
     $        smaxflux(vd/sqrt(2.*Ti),chi)
     $        *rs**2 ))
      nrein=ninjcomp*numprocs
      if(ninjcomp.le.0)ninjcomp=1
      n_part=int(ripernode*4.*3.1415926*rs**3/3.)
      rhoinf=ripernode*numprocs
      if(n_part.gt.n_partmax)then
         write(*,*)'ERROR. Too many particles required.'
         write(*,101)rhoinf,n_part,n_partmax
 101     format('rhoinf=',f8.2,'  needs n_part=',i8
     $        ,'  which exceeds n_partmax=',i8)
         stop
      endif
      end

!********************************************************************
      real function smaxflux(uc,chi)
!     Return the total flux to a unit radius sphere from a unit density
!     maxwellian distribution shifted by velocity
      real uc
!     normalized to sqrt(2T/m), in a spherically symmetric potential
!     having a value on the sphere normalized to Ti of minus
      real chi
      real eps,pi
      data eps/1.e-3/pi/3.1415927/
      erf=1.-erfcc(uc)
      sqpi=sqrt(pi)
      if(abs(uc).lt.eps) then
         erfbyu=(2./sqpi)*(1.-uc**2 /3.)
      else
         erfbyu=erf/uc
      endif
      smaxflux=pi*sqrt(2.)*(uc*erf +(0.5+chi)*erfbyu + exp(-uc**2)/sqpi)
      end
!**********************************************************************
! Set Boundaries that need treatment specific to problem.
      subroutine geominit(myid)
      include 'ndimsdecl.f'
      include '3dcom.f'
      include 'plascom.f'

! Inner boundary setting ----------------- Specific to this problem. 
! First object is sphere of radius rc and potential phi.
      rc=obj_geom(oradius,1)
      phip=-obj_geom(oabc+2,1)/obj_geom(oabc,1)
! Outer boundary setting -----------------
! Don't do this if there's no second object because things then break.
! Second object is bounding sphere of radius rs.
      rs=obj_geom(oradius,2)
! Override the boundary condition of object 2 with an OML condition.
      xlambda=debyelen/sqrt(1.+1./Ti)
      a=1./xlambda+1./rs
      b=1.
      x=rs/xlambda
      adeficit=((1.-2.*phip/Ti)*rc**2/(4.*debyelen**2))
! IHH approximation to exp(x)E1(x) valid to 0.5% for positive x.
      c= (adeficit/rs)*(alog(1.+1./x)-.56/(1.+4.1*x+0.9*x*x))
      if(myid.eq.0)write(*,*)'Outer boundary a,b,c'
     $     ,a,b,c
      call objsetabc(2,a,b,c)
      call adeficitset(adeficit)
      end
!********************************************************************
      subroutine cavereinset(phi)
! Null
      include 'ndimsdecl.f'
      include 'reincom.f'
      include 'partcom.f'
       caverein=crelax*phi+(1.-crelax)*caverein
!      write(*,*)'CAvereinset',caverein
      end
!********************************************************************
!********************************************************************
      subroutine colreinit(myid,ispecies)
! Initialize and normalize the cdistflux factors from the
! Collisional distribution data, presumed already calculated by pinit.
! Based upon the ipartperiod settings.
      include 'ndimsdecl.f'
      include 'cdistcom.f'
      include 'partcom.f'
      if(ncdists(ispecies).eq.0)return
      ctot=0.
      do i=1,ndims
!         if(myid.eq.0)write(*,*)ipartperiod(i),cdistflux(i)
         if(ipartperiod(i).ge.3)cdistfluxs(i,ispecies)=0.
! Normalize the cdistflux and cdistcum to face area.
         cdistfluxs(i,ispecies)=cdistfluxs(i,ispecies)*fcarea(i)
         ctot=ctot+cdistfluxs(i,ispecies)
      enddo
      if(ctot.ne.0.)then
         cdistcums(1,ispecies)=0.
         do i=1,ndims
            cdistfluxs(i,ispecies)=cdistfluxs(i,ispecies)/ctot
            cdistcums(i+1,ispecies)=cdistcums(i,ispecies)
     $           +cdistfluxs(i,ispecies)
         enddo
! Avoid rounding problems.
         cdistcums(ndims+1,ispecies)=1.
      else
         if(myid.eq.0)write(*,*
     $        )'PROBLEM. colreinject: No reinjection faces'
      endif

! Now evaluate the cumulative distribution in 3 normal-directions.
      do id=1,ndims
         fxvcols(1,id,ispecies)=0.
         do i=1,ncdists(ispecies)
            fxvcols(i+1,id,ispecies)=fxvcols(i,id,ispecies)
     $           +abs(vcols(id,i,ispecies))
         enddo
      enddo
! There's a resolution issue in that a million steps can hardly
! be resolved by single precision. However, it is unlikely that
! substantial statistical distortion will occur. If it did we could
! use double precision.
!      if(myid.eq.0)write(*,'(a,7f8.4)')' colreinit completed',cdistcum
!     $     ,cdistflux
!      if(myid.eq.0)write(*,'(a,3f16.4)')' fxvcol:',(fxvcol(ncdist+1,j),j
!     $     =1,3)
      end
!********************************************************************
      real function ranlenposition(id)
! Return a random fractional position in the coordinate direction id,
! accounting for the density scale length if present
      include 'ndimsdecl.f'
      include 'meshcom.f'
      include 'creincom.f'
      include 'plascom.f'
      real expsa(ndims),expsi(ndims)
      logical lfirst
      data lfirst/.true./
      save lfirst,expsa,expsi

      if(lfirst)then
         do i=1,ndims
            g=gn(i)
            s0=gp0(i)
            si=min(xmeshstart(i),xmeshend(i))-s0
            sa=max(xmeshstart(i),xmeshend(i))-s0
            expsa(i)=exp(g*sa)
            expsi(i)=exp(g*si)
!            write(*,*)i,expsa(i),expsi(i),g
         enddo
         lfirst=.false.
      endif
      g=gn(id)
      call ranlux(P,1)
      if(abs(g).ne.0)then
         sp=gp0(id)+alog(P*expsa(id)+(1.-P)*expsi(id))/g
      else
         sp=(1.-P)*xmeshstart(id)+P*xmeshend(id)
      endif
      ranlenposition=sp*0.999999+.0000005
      if(.false.)then
         write(*,*)' ranlenpos',id,P,sp,g,expsi,expsa
         write(*,*)'Ranlenpos error',id,P,sp,g,ranlenposition
         write(*,*)expsi,expsa
      endif
      end
!********************************************************************
