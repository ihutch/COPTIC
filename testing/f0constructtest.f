      integer nphi,nu
      parameter (nphi=100,nu=100)
      real psi,um,xmax
      real phiarray(0:NPHI),us(0:NPHI),xofphi(0:NPHI)
      real den(0:NPHI),denuntrap(0:NPHI),dentrap(0:NPHI)
      real tilden(0:NPHI)
      real f0(-2*nphi:2*nphi),u0(-2*nphi:2*nphi)
      real u(-nu:nu),f(-nu:nu),cumf(-nu:nu)
      integer nptsmax
      parameter (nptsmax=500)
      real xp(nptsmax),yp(nptsmax)

c    1: plot distribution tests. 2: run 1M calls for timing.
      idebug=1

      psi=1.
      um=0.0
c Hole (decay) length (from runBGKint)
      coshlen=4.+psi/2
c The flattop length toplen. Negligible for large negative values      
      tl=-1./psi
      xmax=1.3*findxofphi(psi/(NPHI),psi,coshlen,tl,0.,50.,7)
      call f0Construct(nphi,psi,um,xmax,coshlen,tl,
     $     phiarray,us,xofphi,den,denuntrap,dentrap,tilden,f0,u0)

c Initialize u-range
      umax=3.
      du=umax/nu
      do i=-nu,nu
         u(i)=i*du
      enddo


      if(idebug.eq.1)then
         call autoplot(xofphi,phiarray,nphi+1)
         call axlabels('x','phi')
         call pltend

      call autoplot(u0(0),f0(0),nphi+1)
      call axlabels('u0','f0 trapped')
      call pltend

      call autoplot(u0(-2*nphi),f0(-2*nphi),4*nphi+1)
      call axlabels('u0','f0')
      call pltend

      nline=10
c Plot distributions.
      call pltinit(-umax,umax,0.,.8)
      call axis
      call axis2
      call axlabels('u','f(u)')
      do i=0,nline
         phi=i*psi/nline
         call GetDistribAtPhi(psi,um,nphi,f0,u0,phi,nu,u,f,cumf)
         call color(i+1)
         call polyline(u,f,2*nu+1)
      enddo
      call pltend

c Plot cumulative distribs
      call pltinit(-umax,umax,0.,1.2)
      call axis
      call axis2
      call axlabels('u','!AJ!@!uu!u f(u)du normalized')
!      nline=5
      do i=0,nline
         phi=i*psi/nline
         call GetDistribAtPhi(psi,um,nphi,f0,u0,phi,nu,u,f,cumf)
         call color(i+1)
         call polyline(u,cumf/cumf(nu),2*nu+1)
      enddo
      call pltend
      endif

      if(idebug.eq.2)then
! This test takes 2.7 seconds for nphi=50, nu=50, 5.4s for nu=100.
      do i=0,1000000
         call ranlux(r,1)
         phi=r*psi
         call GetDistribAtPhi(psi,um,nphi,f0,u0,phi,nu,u,f,cumf)
! This adds negligible time:
         call ranlux(r,1)
         cf=r*cumf(nu)
         iupos=interp(cumf,2*nphi+1,cf,upos) 
      enddo
      endif

      if(idebug.eq.3)then
! Test findxofran
         npts=30
         nbi=15
         xmax=10.
         xmin=-xmax
         do i=1,npts
            xp(i)=-xmax+2.*xmax*float(i)/npts
            yp(i)=(xp(i)-xmin
     $           +derivphiofx(xp(i),psi,coshlen,t1))/(xmax-xmin)
         enddo
         call autoplot(xp,yp,npts)
         call pltend
         write(*,*)'Testing findxofran for xmax=',xmax
         do i=1,npts
            ranf=float(i)/npts-1.e-6
            x=findxofran(ranf,psi,coshlen,tl,-xmax,xmax,nbi)
            write(*,*)i,ranf,x
            yp(i)=yp(i)
            xp(i)=xp(i)
         enddo
         call autoplot(xp,yp,npts)
         call axlabels('findxofran','ranf')
         call pltend
      endif

      end
