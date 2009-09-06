      subroutine reinject(xr,ilaunch)
      real xr(*)
      integer ilaunch
      include 'reincom.f'
      logical lfirst
      data lfirst/.true./
c This call passes address of the particle slot to reiject.      
      dt=1.
      ibcr=0
      icolntype=0
      if(lfirst)then
         write(*,*)'Initializing reinjection'
         call injinit(icolntype,ibcr)
         lfirst=.false.
      endif
      call thisreinject(1,dt,icolntype,ibcr,xr,ilaunch)
c      write(*,*)'Returning from reinject',(xr(j),j=1,6)
      end
c***********************************************************************
c General version allows choice of reinjection scheme.
c***********************************************************************
      subroutine thisreinject(i,dt,icolntype,bcr,xp,ilaunch)
      implicit none
      integer i,icolntype,ilaunch
      real dt
      integer bcr 
      real xp(6,*)
      if(bcr.ne.0) then
c         call maxreinject(i,dt,bcr)
      elseif(icolntype.eq.1.or.icolntype.eq.5) then
c Injection from fv distribution at the boundary.
c         call fvreinject(i,dt,icolntype)
      elseif(icolntype.eq.2.or.icolntype.eq.6)then
c Injection from a general gyrotropic distribution at infinity
c         call ogenreinject(i,dt)
      else
c Injection from a shifted maxwellian at infinity
         call oreinject(i,dt,xp,ilaunch)
      endif
      end
c***********************************************************************
      subroutine injinit(icolntype,bcr)
      implicit none
      integer icolntype
      integer bcr

      if(bcr.ne.0) then
c Injection from a maxwellian at boundary?
c         call maxinjinit(bcr)
      elseif(icolntype.eq.1.or.icolntype.eq.5) then
c Injection from fv distribution at the boundary.
c         call fvinjinit(icolntype)
      elseif(icolntype.eq.2.or.icolntype.eq.6)then
c Injection from a general gyrotropic distribution at infinity
c         call ogeninjinit(icolntype)
      else
c Injection from a shifted maxwellian at infinity
         call oinjinit()
      endif
      end
c***********************************************************************
c***********************************************************************
c Other versions are in other source files.
      subroutine oreinject(i,dt,xp,ilaunch)
      implicit none
      integer i
      real dt
      real xp(6,*)
      integer ilaunch
      include 'plascom.f'
      include 'rancom.f'
      include 'reincom.f'
c      real averein,debyelen,pi,Ti,vd
c      real Vcom
c      integer nr
c      integer nrfull,nrused
c      integer nthsize
c      parameter (nthsize=2)
c      integer nrsize
c      parameter (nrsize=2)
c      integer nvel,ntrapre
c      parameter (nvel=2)
c      real r(nrsize)
c      real th(nthsize),phi(1,nthsize)
c      real getphihere
c End of problematic variables (fix).

c Common data:
      real eup
      parameter (eup=1.e-10)
      external pu
c      logical istrapped
c Testing
      real vdist(nvel)
      real tdist(nthsize)
      real crdist(nthsize),cidist(nthsize)
      common/rtest/crdist,cidist,tdist,vdist

c Making all variables declared:
      real a1,brc,brcsq,ceta,chium2,cosal,sinal,crt,czt,dv,eta
      real expuu2,p2,rcyl,rp,seta,srt,szt,u,ua,ua1,ua3
      real Uc,uu2,vscale,y,y1,zt,ran1,alcos,alsin
      integer icr,idum,ierr,iv
      real rri

c In this routine we work in velocity units relative to ion thermal till end.
      vscale=sqrt(2.*Ti)
      idum=1
      ilaunch=0
c Silence warnings.
      brcsq=0.
      ierr=0
      rp=0.
 1    continue
      ilaunch=ilaunch+1
      if(ilaunch.gt.1000)then
         write(*,*)'ilaunch excessive. averein=',averein,' brcsq=',
     $        brcsq,' ierr=',ierr,' rp=',rp
         stop
      endif
c Pick normal velocity from cumulative Pu
      y1=ran1(idum)
      call finvtfunc(pu,nvel,y1,u)
      iv=int(u)
      dv=u-iv
      u=dv*Vcom(iv+1)+(1.-dv)*Vcom(iv)
      if(.not.dv.le.1)write(*,*)'Error in u calculation',y1,u,iv,dv
      vdist(iv)=vdist(iv)+1.
c Pick angle from cumulative Pc.
      y=ran1(idum)
c Here the drift velocity is scaled to the ion temperature.
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
c This expression is evaluated very inaccurately if expuu2 is nearly 1.
c That is why such a large limit on abs(uu2) is adopted.
         crt=alog(y*expuu2 + (1-y)/expuu2)/uu2
c The following do not do any better at solving this problem.
c         crt=alog( (y*expuu2 + (1-y)/expuu2)**(1./uu2))
c         crt=-1.+alog(1.+(expuu2**2-1.)*y)/uu2
      endif
      if(.not. abs(crt).le.1)then
c         write(*,*)'Strange crt:',crt,y,expuu2,uu2
c It seems impossible to avoid occasional strange results when uu2 is small.
         crt=2.*y-1.
      endif
c Testing angular distribution.
      icr=int((1.+crt)*0.5*((nthsize-1)-1) + 1.5)
      crdist(icr)=crdist(icr)+1.
c End of testing distribution monitor.
      srt=sqrt(1.- crt**2)
c Pick angle zt of poloidal impact and angle eta of impact parameter
      zt=ran1(idum)*2.*pi
      czt=cos(zt)
      szt=sin(zt)
      eta=ran1(idum)*2.*pi
      ceta=cos(eta)
      seta=sin(eta)
c Choose impact parameter, preventing overflow.
      chium2=-averein/Ti/(u+eup)**2
      if(chium2.le.-1.) then
c         write(*,*)'Impossible chium2=',chium2,' averein=', averein,
c     $        ' u=',u,' iv=',iv
         goto 1
c         stop
      endif
c      if(.not.lfixedn)chium2=0.
      brcsq=ran1(idum)*(1.+chium2)
c Reject a particle that will not reach boundary.
      if(brcsq.lt.0.) then
         goto 1
      endif
      brc=sqrt(brcsq)
c Get cosine and sine of impact angle relative to distant position.
c Based on integration.
      p2=brcsq*2.*Ti*u**2
      ierr=0
      if(debyelen.gt..001)then
c Orbit integration angle calculation.
c There is an overflow with this at zero debyelen. Ought to be properly fixed.
         call alphaint(p2,brcsq,cosal,sinal,ierr)
         if(ierr.ne.0)goto 1
c      write(*,'(4f9.4)')cosal-alcos(brc,chium2),sinal-alsin(brc,chium2)
c     Now ilaunch is the number of launches at infinity it took to get
c     one that reached the boundary.
      else
c Alternative based on analytic orbit calculation.
c Used for low debyelen, but really assumes negligible boundary potential.
c         call alcossin(brc,chium2,cosal,sinal)
         cosal=alcos(brc,chium2)
         sinal=alsin(brc,chium2)
      endif
c Install reinjection position
      a1=crt*ceta*sinal+srt*cosal
c Ensure that we are just inside the rs sphere
      rri=0.99999*rs
      xp(1,i)=rri*(czt*a1+szt*seta*sinal)
      xp(2,i)=rri*(-szt*a1+czt*seta*sinal)
      xp(3,i)=rri*(-srt*ceta*sinal + crt*cosal)
c Injection velocity components normalized in the rotated frame:
      ua1=-brc*cosal -sqrt(1.+chium2-brcsq)*sinal
      ua3=brc*sinal - sqrt(1.+chium2-brcsq)*cosal
      ua=crt*ceta*ua1+srt*ua3
c Install reinjection velocity in Te-scaled units
      u=u*vscale
      xp(4,i)=u*(czt*ua+szt*seta*ua1)
      xp(5,i)=u*(-szt*ua+czt*seta*ua1)
      xp(6,i)=u*(-srt*ceta*ua1 + crt*ua3)
c If we don't recalculate rp, then we don't trap NANs in the random choices.
      rcyl=xp(1,i)**2+xp(2,i)**2
      rp=rcyl+xp(3,i)**2
c      rp=rs
c      write(*,*)'oreinject',rp
c Reject particles that are already outside. But this is more to do
c with detecting NANS.
      if(.not.rp.lt.rs*rs)then
c      if(.not.rp.le.r(nr)*r(nr))then
         write(*,*)'Relaunch',rp,xp(1,i),xp(2,i),xp(3,i)
         write(*,*)vscale,u,brc,chium2,brcsq,sinal,cosal
         goto 1
      endif
      end
c********************************************************************
c Initialize the distributions describing reinjected particles
      subroutine oinjinit()
      implicit none
c Common data:
      include 'plascom.f'
      include 'rancom.f'
      include 'reincom.f'

c To fix exclusion we need
c      real Vcom(1),pu1(1),pu2(1)
c      real Ti,vd,pi
c      integer myid,nvel


c Making variables declared
      integer i
      real expu0,u0,Uc,uminus,uplus,vspread,erfcc

c Here the drift velocity is scaled to the ion temperature.
c And U's are in units of sqrt(2T/m), unlike vd.
      Uc=abs(vd)/sqrt(2.*Ti)
c Range of velocities (times (Ti/m_i)^(1/2)) permitted for injection.
      vspread=5.+abs(Uc)

c Can't use these formulas for Uc exactly equal to zero.
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
      call srand(myid+1)
      end
c***********************************************************************
c***********************************************************************
c Calculate the cumulative probability for velocity index iu such that
c         u= vspread*(iu-1.)/(nvel-1.)   as per injinit
      real function pu(iu)
      implicit none
      real pudenom
      integer iu
c     averein is the average potential of reinjected particles, which is
c     used as an estimate of the potential at the reinjection boundary.
c     It is expressed in units of Te so needs to be scaled to Ti.
      include 'plascom.f'
      include 'rancom.f'
      include 'reincom.f'

c needed from common:
c      real averein,Ti,pu1(1),pu2(1)

      pudenom=pu1(1)-pu2(1)*averein/Ti
      pu=1.- (pu1(iu)-pu2(iu)*averein/Ti)/pudenom
      end
c********************************************************************
c Given a monotonic (increasing?) 
c function Q(x) on a 1-D grid x=1..nq, solve Q(x)=y for x.
c That is, invert Q to give x=Q^-1(y).
      subroutine finvtfunc(Q,nq,y,x)
c Somehow this breaks the passing of a function reference.
      implicit none
c      real external Q
      real Q
      external Q
      integer nq
      real y,x
c
      integer iqr,iql,iqx
      real Qx,Qr,Ql,Qd
      Ql=Q(1)
      Qr=Q(nq)
      iql=1
      iqr=nq
      if((y-Ql)*(y-Qr).gt.0.) then
c Value is outside the range.
         x=0
         return
      endif
 200  if(iqr-iql.eq.1)goto 210
      iqx=(iqr+iql)/2
      Qx=Q(iqx)
c      write(*,*)y,Ql,Qx,Qr,iql,iqr
c Formerly .lt. which is an error.
      if((Qx-y)*(Qr-y).le.0.) then
         Ql=Qx
         iql=iqx
      else
         Qr=Qx
         iqr=iqx
      endif
      goto 200
 210  continue
c Now iql and iqr, Ql and Qr bracket Q
c      x=(y-Ql)/(Qr-Ql)+iql
c Trap errors caused by flat sections.
      Qd=Qr-Ql
      if(Qd.eq.0.)then
         x=(iql+iqr)/2.
      else
         x=(y-Ql)/Qd+iql
      endif
      end
c**********************************************************************
C Inverse square law (phi\propto 1/r) injection functions:
c**********************************************************************
      subroutine alcossin(s,c,cosal,sinal)
      real s,c,cosal,sinal,alcos,alsin
      cosal=alcos(s,c)
      sinal=alsin(s,c)
      end
c**********************************************************************
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
c**********************************************************************
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
c**********************************************************************
c**********************************************************************
c Return the angle cosine and sine for impact of a particle calculated
c by integrating the angle formula for a general central force.
c We integrate in the variable q=1/r from 0 to 1, provided that we don't
c encounter a barrier. If we do encounter one, we return an error.
c
c This version designed to work with uneven spacing if
c necessary, and trapezoidal integration. And uses init2ext
      subroutine alphaint(p2,b2,cosal,sinal,ierr)
      implicit none
c     The angular momentum and impact parameter squared
c     in units such that injection radius is 1.
      real p2, b2
c     The angle values returned.
      real cosal, sinal
c     The error return signal if non-zero.
      integer ierr
c Common data:
      include 'plascom.f'
c Not needed:
c      include 'rancom.f'
      include 'reincom.f'

c     Need real adeficit,averein
c     Need debyelen


      integer iqsteps
c This choice ensures iqsteps is large enough to accommodate a profile
c read in for processing with orbitint, but might not be the best choice for
c regular use.
c      parameter (iqsteps=nrsize+1)
      parameter (iqsteps=100+1)
      real phibye(iqsteps),phiei(iqsteps)
      real qp(iqsteps),pp(iqsteps)
      real alpha,b2i,d1,d2,p2i2,sa,xlambda
      integer i,iqs
      logical uninitialized
      data uninitialized/.true./
      save
c      real extpot
c Statement function
c Inverse square law potential.
c      extpot(q)=averein*q
c Not currently in use.

c First time through, do the initialization.
      if(uninitialized)then
         iqs=iqsteps
         uninitialized=.false.
c Screening length accounting for both ions and electrons.
         xlambda=debyelen/sqrt(1.+1./Ti)
         if(xlambda.lt.1.e-6)xlambda=1.e-6

C$$$         if(diagrho(1).eq.999.)then
C$$$c Special case to use the read in data.
C$$$c Some common data used in abnormal ways.
C$$$            qp(1)=0.
C$$$            phibye(1)=0.
C$$$            if(ninner.gt.iqsteps)stop 'Alphaint Too many r-cells read.'
C$$$            do i=1,ninner
C$$$               qp(i+1)=1./rcc(ninner+1-i)
C$$$               phibye(i+1)=diagphi(ninner+1-i)
C$$$            enddo
C$$$            iqs=ninner+1
C$$$            if(qp(iqs).ne.1.)then
C$$$               iqs=iqs+1
C$$$               qp(iqs)=1.
C$$$c Extrapolate linearly in q.
C$$$               phibye(iqs)=phibye(iqs-1)+
C$$$     $              (phibye(iqs-1)-phibye(iqs-2))*
C$$$     $              (qp(iqs)-qp(iqs-1))/(qp(iqs-1)-qp(iqs-2))
C$$$            endif
C$$$c The read-in case just uses phibye not phiei. And scales phibye(iqs) to 1
C$$$c putting the absolute edge value into averein.
C$$$            averein=phibye(iqs)
C$$$            adeficit=0.
C$$$            do i=1,iqs
C$$$               phibye(i)=phibye(i)/phibye(iqs)
C$$$            enddo
C$$$c End of read-in potential special case.
C$$$         else
c Specify the q-array. It goes from 0 to 1.
            do i=1,iqs
               qp(i)=((i-1.)/(iqs-1.))
            enddo
c            write(*,*)iqs,r(nr),xlambda,averein,adeficit
c     Since phi is specified in real space units, we tell the initialization
c     function what the rmax really is, and it does the transformation.
c            call initext(iqs,qp,phibye,phiei,r(nr),xlambda)
            call initext(iqs,qp,phibye,phiei,rs,xlambda)
C$$$         endif
c Diagnostics
c         write(*,'(3f8.4)')(qp(j),phibye(j),phiei(j),j=1,iqs)
c         write(*,*)'averein,adeficit',averein,adeficit
         if(averein.ne.0)then
c     When used by SCEPTIC averein is zero the first time, so this will
c     not be called.
            do i=1,iqs
               pp(i)=averein*phibye(i)-adeficit*phiei(i)
            enddo
c            write(*,*)'adeficit=',adeficit
            call autoplot(qp,pp,iqs)
            call axlabels('q','potential')
            call pltend()
         endif
      endif
c End of initialization section.
c Do the integration for the orbit.
      ierr=0
      b2i=1./b2
      p2i2=2./p2
      sa=b2i - qp(1)**2 - p2i2*(averein*phibye(1)-adeficit*phiei(1))
      d1=(1./sqrt(sa))
c Inverse square case.
c      d1=1./sqrt(b2i)
      alpha=0.
c Trapezoidal rule.
      do i=2,iqs
         d2=d1
         sa=b2i - qp(i)**2 - p2i2*(averein*phibye(i)-adeficit*phiei(i))
c Inverse square case.
c         sa=b2i - qp(i)**2 - p2i2*extpot(qp(i))
         if(sa .le. 0.) goto 2
         d1=(1./sqrt(sa))
         alpha=alpha+(qp(i)-qp(i-1))*(d1+d2)*.5
      enddo
c      write(*,*)'alpha=',alpha
c Negative sign for definition of alpha relative to the forward direction.
      cosal=-cos(alpha)
      sinal=sin(alpha)
c      write(*,*)'averein,adeficit',averein,adeficit,alpha
      return
 2    ierr=i
c      write(*,101)i,b2,p2,averein,adeficit
c 101  format('Barrier: i=',i3,' b2=',f8.1,
c     $     ' p2=',f8.4,' averein=',f8.4,' adeficit=',f8.4)
      end
c********************************************************************
      subroutine avereinset(aver)
c Set the averein value
      include 'reincom.f'
      averein=aver
      end
c********************************************************************
      subroutine adeficitset(adefic)
      include 'reincom.f'
      adeficit=adefic
      end
c*******************************************************************
      subroutine plascomset(d,T,v,r,p)
      include 'plascom.f'
      debyelen=d
      Ti=T
      vd=v
      rs=r
      phip=p
      end
c*********************************************************************
      subroutine rhoinfcalc(dtin)
c Obtain the rhoinf to be used in calculating the electron shielding,
c based upon the number and average potential of the reinjections.
      include 'plascom.f'
c No time-averaging for now.
c Use particle information for initializing.
      include 'partcom.f'

c Moved nrein and phirein reductions to psumreduce.
      if(nrein.ne.0)then
c Calculate rhoinf from nrein if there are enough.
         chi=max(-phirein/Ti,-0.5)
         riest=(nrein/dtin) /
     $        (sqrt(Ti)*
     $        smaxflux(vd/sqrt(2.*Ti),chi)
     $        *rs**2 )
         rhoinf=riest
c         write(*,*)'nrein,dtin,Ti,vd,phirein,chi,rs,numprocs=',
c     $        nrein,dtin,Ti,vd,phirein,chi,rs,numprocs
      else
         if(rhoinf.lt.1.e-4)then
c Approximate initialization
            rhoinf=numprocs*n_part/(4.*3.1415926*rs**3/3.)
            write(*,*)'Rhoinf in rhoinfcalc approximated as',rhoinf
     $           ,numprocs,n_part,rs
         endif
c Else just leave it alone.
      endif
      end
c*********************************************************************
      subroutine nreincalc(dtin)
c Given rhoinf, decide the number of reinjections per step ninjcomp
c for average edge potential.
      include 'plascom.f'
c No time-averaging for now.
c Particle information
      include 'partcom.f'

c Calculate ninjcomp from rhoinf
      chi=min(-phirein/Ti,0.5)
      ninjcomp=nint(rhoinf*dtin* 
     $        (sqrt(Ti)*
     $        smaxflux(vd/sqrt(2.*Ti),chi)
     $        *rs**2 ))
      nrein=ninjcomp
      if(ninjcomp.le.0)ninjcomp=1
      n_part=int(rhoinf*4.*3.1415926*rs**3/3./numprocs)
      if(n_part.gt.n_partmax)then
         write(*,*)'ERROR. Too many particles required.'
         write(*,101)rhoinf,n_part,n_partmax
 101     format('rhoinf=',f8.2,'  needs n_part=',i8
     $        ,'  which exceeds n_partmax=',i8)
         stop
      endif
      end

c********************************************************************
      real function smaxflux(uc,chi)
c     Return the total flux to a unit radius sphere from a unit density
c     maxwellian distribution shifted by velocity
      real uc
c     normalized to sqrt(2T/m), in a spherically symmetric potential
c     having a value on the sphere normalized to Ti of minus
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
c**********************************************************************
c Set Boundaries that need treatment specific to problem.
      subroutine geominit(myid)
c Inner boundary setting ----------------- Specific to this problem. 
c First object is sphere of radius rc and potential phi.
      include '3dcom.f'
      include 'plascom.f'

      rc=obj_geom(5,1)
      phip=-obj_geom(10,1)/obj_geom(8,1)
c Outer boundary setting -----------------
c Second object is bounding sphere of radius rs.
      rs=obj_geom(5,2)
c But use a tad more for the mesh size
      rsmesh=obj_geom(5,2)*1.00001
c Override the boundary condition of object 2 with an OML condition.
      xlambda=debyelen/sqrt(1.+1./Ti)
      a=1./xlambda+1./rs
      b=1.
      x=rs/xlambda
      adeficit=((1.-2.*phip/Ti)*rc**2/(4.*debyelen**2))
c IHH approximation to exp(x)E1(x) valid to 0.5% for positive x.
      c= (adeficit/rs)*(alog(1.+1./x)-.56/(1.+4.1*x+0.9*x*x))
      if(myid.eq.0)write(*,*)'Outer boundary a,b,c'
     $     ,a,b,c
      call objsetabc(2,a,b,c)
      call adeficitset(adeficit)
      end
