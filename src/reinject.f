! This file contains the reinjection code which is specific to a
! particular (type of) boundary. 
! It must provide the routines:
!    reinject(xr,nrein,caverein) 
!       which reinijects a particle in position/velocity slot xr(2*ndims)
!       and returns the number of launches. (Pass xp(1,i)).
!       This routine should initialize itself appropriately.
!    rhoinfcalc(dt)     
!        which should return the value of rhoinfinity through common.
!    ninjcalc(dt)
!        which decides how many reinjections per step ninjcomp.
! Any other code needed for these routines ought to be in here too.
! The majority of the needed data is passed in partcom.f
! 
! This case is a spherical boundary.
! It is not ready for multiple species or fractional injections.
!***********************************************************************
      subroutine reinject(xr,ilaunch,caverein)
      parameter (mdims=3)
      real xr(3*mdims)
!      integer i
      include 'ndimsdecl.f'
! Common data for these routines:
      include 'reincom.f'
! Plasma common data
      include 'plascom.f'
      logical lfirst
      save lfirst
      data lfirst/.true./

      if(lfirst) then
         call injinit()
         lfirst=.false.
      endif

      vscale=sqrt(Ti)
      vdi=vd/vscale
      cerr=0.
      idum=1
 1    continue
      call ranlux(y,1)
! Pick angle from cumulative Q.
      call invtfunc(Qcom,nQth,y,x)
      ic1=x
      if(x.lt.1. .or. x.ge.float(nQth))then
         write(*,*)  'REINJECT Q-Error'
         write(*,*)'y,x,nQth=',y,x,nQth
         write(*,*)'Qcom=',Qcom
         goto 1
      endif
!      if(ic1.ge.nQth)ic1=ic1-1
      ic2=ic1+1
      dc=x-ic1
      call ranlux(y,1)
! Pick normal velocity from cumulative G.
      call invtfunc(Gcom(1,ic1),nvel,y,v1)
      call invtfunc(Gcom(1,ic2),nvel,y,v2)
      vr=dc*v2+(1.-dc)*v1
      if(vr.lt.1. .or. vr.ge.nvel) then
         write(*,*) 'REINJECT V-Error'
         write(*,*) y,v1,v2,ic1,ic2,nvel
         goto 1
      endif
      iv=vr
      dv=vr-iv
      vr=dv*Vcom(iv+1)+(1.-dv)*Vcom(iv)
! New angle interpolation.
      ct=1.-2.*(x-1.)/(nQth-1)
      ic1=x
      ic2=ic1+1
      dc=x-ic1
! ct is cosine of the angle of the velocity -- opposite to the radius.      
      st=sqrt(1.- ct**2)
! Now we have cosine theta=c and normal velocity normalized to v_ti.
! Theta and phi velocities are (shifted) Maxwellians but we are working
! in units of vti.
      vt=gasdev(0)- st*vdi
      vp=gasdev(0)
! All velocities now.
      call ranlux(p,1)
      p=2.*pi*p
      cp=cos(p)
      sp=sin(p)
!      write(*,*)ct,st,cp,sp
! If velocity is normalized to sqrt(Te/mi), and Ti is Ti/Te really,
! then a distribution with standard deviation sqrt(Ti/mi) is obtained
! from the unit variance random distribution multiplied by sqrt(Ti)=vscale.
      xr(6)=(vr*ct - vt*st)*vscale
      xr(5)=((vr*st+ vt*ct)*sp + vp*cp)*vscale
      xr(4)=((vr*st+ vt*ct)*cp - vp*sp)*vscale

      rfrac=0.9999
      xr(3)=-rfrac*rs*ct
      xr(2)=-rfrac*rs*st*sp
      xr(1)=-rfrac*rs*st*cp

!      write(*,'(''Reinject vr,vt,vp='',3f8.3)') vr,vt,vp
!      write(*,'(''Reinject i,xr=''i6,6f8.3)')i,(xr(kk),kk=1,6)

!      vv2=vt**2 + vr**2 + vp**2
! Only one launch attempt here.
      ilaunch=1
! Direct ic1 usage
      end
!********************************************************************
! Initialize the distributions describing reinjected particles
      subroutine injinit()
! Common data:
      include 'ndimsdecl.f'
      include 'reincom.f'
      include 'plascom.f'
      real gam(nQth)

! Range of velocities (times (Ti/m_i)^(1/2)) permitted for injection.
      vspread=5.+abs(vd)/sqrt(Ti)
! Random interpolates
      sq2pi=1./sqrt(2.*pi)
      sq2=1./sqrt(2.)
      Qcom(1)=0.
      dqp=0.
      do i=1,nQth
! Qth is the cosine angle of the ith angle interpolation position. 
! We used to used th(i). This is the equivalent definition.
         Qth=1.-2.*(i-1.)/(nQth-1.)
! Here the drift velocity is scaled to the ion temperature.
         vdr=vd*Qth/sqrt(Ti)
         dqn= sq2pi*exp(-0.5*vdr**2) + .5*vdr*erfcc(-sq2*vdr)
         if(i.gt.1) Qcom(i)=Qcom(i-1) +dqp +dqn
! Gamma is the total flux over all velocities at this angle. 
! But it is never used 
         gam(i)=dqn
         dqp=dqn
! At this angle,
         do j=1,nvel
! Contruct the Cumulative distribution of radial velocity
! On the mesh Vcom
            Vcom(j)=vspread*(j-1.)/(nvel-1.)
! The cumulative distribution is Gcom, 
            Gcom(j,i)=(dqn - sq2pi*exp(-0.5*(Vcom(j)-vdr)**2)
     $           - .5*vdr*erfcc(sq2*(Vcom(j)-vdr)) )/dqn
         enddo
         Gcom(1,i)=0.
         Gcom(nvel,i)=1.
      enddo
      do i=2,nQth
         Qcom(i)=Qcom(i)/Qcom(nQth)
      enddo
! Now Gcom(*,i) is the cumulative distribution of radial velocity at cos(Qth)
! normalized to the ion thermal velocity, not sqrt(T_e/m_i).
! And Qcom() is the cumulative distribution in cosine angles Qth
     
! 501  format(a,11f8.4)
      end
!*********************************************************************
      subroutine rhoinfcalc(dtin)
! Return the rhoinf to be used in calculating the electron shielding,
! based upon the number and average potential of the reinjections.
      include 'ndimsdecl.f'
      include 'plascom.f'
! No time-averaging for now.
! Use particle information for initializing.
      include 'partcom.f'

! Moved nrein and phirein reductions to psumreduce.

      if(nrein.ne.0)then
! Calculate rhoinf from nrein if there are enough.
         chi=min(-phirein/Ti,0.5)
         riest=(nrein/dtin) /
     $        (sqrt(Ti)*
     $        smaxflux(vd/sqrt(2.*Ti),chi)
     $        *rs**2 )
         rhoinf=riest
      else
         if(rhoinf.lt.1.e-4)then
! Approximate initialization
            rhoinf=numprocs*n_part/(4.*3.1415926*rs**3/3.)
            write(*,*)'Rhoinf in rhoinfcalc approximated as',rhoinf
     $           ,numprocs,n_part,rs
         endif
! Else just leave it alone.
      endif
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
         write(*,101)ripernode,n_part,n_partmax
 101     format('ripernode=',f8.2,'  needs n_part=',i8
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
!********************************************************************
      subroutine avereinset()
      return
      end
!********************************************************************
      subroutine geominit()
      return
      end
!********************************************************************
      subroutine cavereinset(phi)
! Null
      include 'reincom.f'
      include 'ndimsdecl.f'
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
