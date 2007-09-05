c***********************************************************************
      subroutine reinject(i,xr)
      parameter (mdims=3)
      real xr(3*mdims)
      integer i
c Common data:
      include 'rancom.f'
      common /myidcom/myid,nprocs
c Hack of required information probably needs to be passed or common.
c      real Ti,vd,pi,rs
c      parameter (Ti=1.,vd=0.,pi=3.1415927,rs=.45)
      include 'plascom.f'
      logical lfirst
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
      y=ran0(idum)
c Pick angle from cumulative Q.
      call invtfunc(Qcom,nQth,y,x)
      ic1=x
      if(x.lt.1. .or. x.ge.float(nQth))then
         write(*,*)  'REINJECT Q-Error'
         write(*,*)'y,x,nQth=',y,x,nQth
         write(*,*)'Qcom=',Qcom
         goto 1
      endif
c      if(ic1.ge.nQth)ic1=ic1-1
      ic2=ic1+1
      dc=x-ic1
      y=ran0(idum)
c Pick normal velocity from cumulative G.
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
c New angle interpolation.
      ct=1.-2.*(x-1.)/(nQth-1)
      ic1=x
      ic2=ic1+1
      dc=x-ic1
c Old version used th() directly.
c      ct=th(ic1)*(1.-dc)+th(ic2)*dc
c      write(*,*)'ic1,ic2,dc,ct',ic1,ic2,dc,ct
c ct is cosine of the angle of the velocity -- opposite to the radius.      
      st=sqrt(1.- ct**2)
c Now we have cosine theta=c and normal velocity normalized to v_ti.
c Theta and phi velocities are (shifted) Maxwellians but we are working
c in units of vti.
      vt=gasdev(idum)- st*vdi
      vp=gasdev(idum)
c All velocities now.
      p=2.*pi*ran0(idum)
      cp=cos(p)
      sp=sin(p)
c      write(*,*)ct,st,cp,sp
c If velocity is normalized to sqrt(Te/mi), and Ti is Ti/Te really,
c then a distribution with standard deviation sqrt(Ti/mi) is obtained
c from the unit variance random distribution multiplied by sqrt(Ti)=vscale.
      xr(6)=(vr*ct - vt*st)*vscale
      xr(5)=((vr*st+ vt*ct)*sp + vp*cp)*vscale
      xr(4)=((vr*st+ vt*ct)*cp - vp*sp)*vscale

      rfrac=0.9999
      xr(3)=-rfrac*rs*ct
      xr(2)=-rfrac*rs*st*sp
      xr(1)=-rfrac*rs*st*cp

c      write(*,'(''Reinject vr,vt,vp='',3f8.3)') vr,vt,vp
c      write(*,'(''Reinject i,xr=''i6,6f8.3)')i,(xr(kk),kk=1,6)

c      vv2=vt**2 + vr**2 + vp**2
c Do the outer flux accumulation.
      nrein=nrein+1
c Direct ic1 usage
      end
c********************************************************************
c Initialize the distributions describing reinjected particles
      subroutine injinit()
c Common data:
      include 'rancom.f'
      common /myidcom/myid,nprocs
      real gam(nQth)
c      character*1 work(nvel,nth)
c Hack of required information probably needs to be passed or common.
c      real Ti,vd,pi
c      parameter (Ti=1.,vd=0.,pi=3.1415927)
      include 'plascom.f'

c Range of velocities (times (Ti/m_i)^(1/2)) permitted for injection.
      vspread=5.+abs(vd)/sqrt(Ti)
c Random interpolates
      sq2pi=1./sqrt(2.*pi)
      sq2=1./sqrt(2.)
      Qcom(1)=0.
      dqp=0.
      do i=1,nQth
c Qth is the cosine angle of the ith angle interpolation position. 
c We used to used th(i). This is the equivalent definition.
         Qth=1.-2.*(i-1.)/(nQth-1.)
c Here the drift velocity is scaled to the ion temperature.
         vdr=vd*Qth/sqrt(Ti)
         dqn= sq2pi*exp(-0.5*vdr**2) + .5*vdr*erfcc(-sq2*vdr)
         if(i.gt.1) Qcom(i)=Qcom(i-1) +dqp +dqn
c Gamma is the total flux over all velocities at this angle. 
c But it is never used 
         gam(i)=dqn
         dqp=dqn
c At this angle,
         do j=1,nvel
c Contruct the Cumulative distribution of radial velocity
c On the mesh Vcom
            Vcom(j)=vspread*(j-1.)/(nvel-1.)
c The cumulative distribution is Gcom, 
            Gcom(j,i)=(dqn - sq2pi*exp(-0.5*(Vcom(j)-vdr)**2)
     $           - .5*vdr*erfcc(sq2*(Vcom(j)-vdr)) )/dqn
         enddo
         Gcom(1,i)=0.
         Gcom(nvel,i)=1.
      enddo
      do i=2,nQth
         Qcom(i)=Qcom(i)/Qcom(nQth)
      enddo
c Now Gcom(*,i) is the cumulative distribution of radial velocity at cos(Qth)
c normalized to the ion thermal velocity, not sqrt(T_e/m_i).
c And Qcom() is the cumulative distribution in cosine angles Qth
c      work(1,1)=' '
     
 501  format(a,11f8.4)
c      write(*,*)'Initializing Random Numbers:', myid
      call srand(myid)
      end

