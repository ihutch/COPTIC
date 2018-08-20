! Quiet Particle initialization.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This version is intended for situations with no objects and no collisions.
! The uniformity at low-k is improved by depositing uniformly by qblocks,
! with random placement within the qblock. When there are fewer particles
! remaining than the total number of qblocks, the qblock size is increased
! reducing the total number of qblocks and the process iterated till all
! particles are placed.

! The number of qblocks per dimension is kept less than nqblksmax,
! but greater than zero. And for single-cell sides it is unity.
! It is initialized as nqbset but constrained by those limits.

      subroutine qinit(wavespec)
      implicit none
! Common data:
      include 'ndimsdecl.f'
      include 'partcom.f'
      include 'plascom.f'
      include 'myidcom.f'
      include 'meshcom.f'
      real wavespec(2*ndims+1)
      external ranlenposition,gasdev
! Factor by which we leave additional space for increased number of 
! particles if not fixed:
      real slotsurplus
      real thetamax
      parameter (thetamax=1.)
! Local variables to satisfy implicit none
      integer i,i1,islotmax,ispecies,ko,ntries,islot
!      integer k
      real theta,ranlenposition,gasdev
      integer nqbset(ndims),nqblks(ndims),mlen
      
      holetoplen=holepow
! nqbset parameters adjusted only from 1 parameter: nqblkmax [-pi...]
      do i=1,ndims
         mlen=ixnp(i+1)-ixnp(i)
         if(mlen-2.eq.1)then   ! Single cell sides need no quieting.
            nqbset(i)=1
         else
            nqbset(i)=max(1,nqblkmax)           ! Quiet short sides.
!         nqbset(i)=max(1,min(nqblkmax,mlen-2)) ! Don't quiet short.
         endif
         nqblks(i)=nqbset(i)
      enddo

! No special orbits.
      norbits=0
      i1=1

! Point to the bottom of the particle stack for start of species 1.
      iicparta(1)=1
      do ispecies=1,nspecies
! Conveniently here initialize distribution numbers.
         if(ninjcompa(ispecies).gt.0)then
            slotsurplus=1.1  ! This is safer at 1.3 than 1.1
         else
            slotsurplus=1.
         endif
         ntries=0
! Set tperp to a tiny number if species is infinitely magnetized.
         theta=Bt*eoverms(ispecies)*dt
         if(abs(theta).gt.thetamax)Tperps(ispecies)=1.e-24
! This version requires separable velocity distributions.
! Perhaps some safety checks for this here.
! Scale the number of particles for higher slots to the ions(?)
         nparta(ispecies)=nparta(1)/numratioa(ispecies)
         islotmax=nparta(ispecies)+iicparta(ispecies)-1
         if((islotmax+1)*slotsurplus.gt.n_partmax)then
            write(*,*)'Too few particle slots',n_partmax
     $           ,' max, to accommodate',islotmax
            write(*,*)'multiplied by slotsurplus factor'
     $           ,slotsurplus,int(islotmax*slotsurplus)
            stop
         endif
! ----------------------------- Actual Particle Setting ----------
            islot=iicparta(ispecies)+i1-1
            call fillqblocksize(nqblks,islot,islotmax,ispecies)
!            write(*,'(i8,7f8.4)')(i,(x_part(k,i),k=1,6),x_part(iflag,i)
!     $           ,i=1,100)
!------------------------------- End of Actual Particle Setting --
! Clean up.
!! The maximum used slot for this species
         iocparta(ispecies)=islotmax
! Start of next slot-set may give a gap for overflow.
         iicparta(ispecies+1)=int((islotmax+1)*slotsurplus)
! Zero the overflow slots' flag
         do i=iocparta(ispecies)+1,iicparta(ispecies+1)-1
            x_part(iflag,i)=0
         enddo
         if(myid.eq.0)then
            write(*,101)ispecies,nprocs,
     $           iocparta(ispecies)-iicparta(ispecies)+1,ntries
 101        format(' Initialized species',i2,i4,'x',i7
     $           ,' ntries=',i7,$)
            if(nspecies.gt.0)write(*,'(a,3i8)')' Slots'
     $           ,iicparta(ispecies),iocparta(ispecies)
     $           ,ninjcompa(ispecies)
         endif
! Initialize orbit tracking
         do ko=1,norbits
            iorbitlen(ko)=0
         enddo         
! Don't shift for special particle the subsequent sections.
         i1=1
      enddo             ! End of species iteration

! Set flag of unused slots to 0
      do i=iicparta(ispecies),n_partmax
         x_part(iflag,i)=0
      enddo
! Allow the last species to fill the array by setting ghost init slot:
      iicparta(ispecies)=n_partmax+1
 ! Apply initial wave displacement
      if(.not.lnotallp)call wavedisplace(wavespec)
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine fillqblocksize(nqblks,islot,islotmax,ispecies)
! Fill the qblocks of current size with as many complete sets of particles 
! as can fill no more than islotmax, leaving a remainder. 
! Then reduce the number of blocks and iterate till the remainder is zero.
! Increment islot for each particle added. ntries tallying not needed.
! If ldouble, then inject 2 particles with opposite velocities but same
! position to null out the initial current density.
      include 'ndimsdecl.f'      
      include 'meshcom.f'
      include 'myidcom.f'
      integer nqblks(ndims),islot,islotmax,ispecies

      integer iview(3,ndims),indi(ndims)     !Iterator
      integer nqbt
      logical ldouble
      data ldouble/.true./ ! Whether to place 2 counter-velocity particles.

      nremain=islotmax-islot+1
      nper=1           ! Number injected per placeqbk. 1 or 2.
      if(ldouble)nper=2
      nqbt=nper
      do i=1,ndims
         nqbt=nqbt*nqblks(i)
      enddo
      if(myid.eq.0)then
         write(*,*)'qinit:ispecies  nfills     islot       nqblks ...'
     $        ,'             distributing'
      endif
 2    nfills=int(nremain/nqbt)      ! assuming 1 particle per fill per qblock.
      if(myid.eq.0.and.nremain/float(islotmax).gt.0.01)
     $     write(*,'(8i10)')ispecies,nfills,islot,nqblks,nfills*nqbt
      nremain=nremain-nfills*nqbt
      do i=1,nfills
! Fill using iterator
         icomplete=mditerator(ndims,iview,indi,4,nqblks)
 1       continue
! Place one particle in block indi.
           call placeqblk(nqblks,indi,islot,ispecies,ldouble)
           islot=islot+1
         if(mditerator(ndims,iview,indi,0,nqblks).eq.0)goto 1
      enddo
      if(nremain.le.nper-1)return
! Else reduce nqblks by 2, and iterate
      nqbt=nper
      do i=1,ndims
         nqblks(i)=max(1,nqblks(i)/2)
         nqbt=nqbt*nqblks(i)
      enddo
      goto 2

      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine placeqblk(nqblks,indi,islot,ispecies,ldouble)
! Place a particle randomly in the block addressed by indi
!
! For initializing particles with a trapped hole region
! The 1-d parallel direction is that of B if it is coordinate-aligned,
! or else dimension 1 (x).
!
! Hole parameters [switch-mnemonic,defaults] are 
! holepsi, the peak hole potential [psi,0] 
! holeum, the drift speed of the hole relative to the distribution [v,0]
! holelen, the parallel length of the hole [l,4*Debyelen]
! holetoplen, the stretch parameter for the hole top [tl,-1/psi]
! holerad, the transverse hole radius [r,0]
! holespeed, the rest-frame hole speed [derived: =holeum+vds]

      include 'ndimsdecl.f'      
      include 'meshcom.f'
      include 'partcom.f'
      include 'plascom.f'
      include 'myidcom.f'
      integer nqblks(ndims),indi(ndims),islot,ispecies
      logical linmesh,ldouble
      integer nbi
      parameter (nbi=16)
! Hole-related parameters:
      real phimin
      parameter (phimin=.01)
      integer nphi,nu,nup
      parameter (nphi=50,nu=100,nup=2*nu+1)
      real psi,um,xmax
      real phiarray(0:NPHI),us(0:NPHI),xofphi(0:NPHI)
      real den(0:NPHI),denuntrap(0:NPHI),dentrap(0:NPHI)
      real tilden(0:NPHI)
      real f0(-2*nphi:2*nphi),u0(-2*nphi:2*nphi),cumf(-nu:nu)
!      real f1(-2*nphi:2*nphi),u1(-2*nphi:2*nphi),cumf1(-nu:nu)
      real u(-nu:nu),f(-nu:nu),du
      real tisq,tisq2,tisqperp,umax,coshlen
      integer lastspecies,id
      data lastspecies/0/
      save

!----------------------------------------
      if(ispecies.ne.lastspecies)then    !Initialization
! These values must be set even for zero hole depth.
         idebug=0
         do i=1,ndims
            if(Bfield(i).ne.0.)goto 2
         enddo
         i=1
 2       id=i
         tisq=sqrt(Ts(ispecies)*abs(eoverms(ispecies)))
         tisq2=sqrt(2.)*tisq
         tisqperp=sqrt(Tperps(ispecies)*abs(eoverms(ispecies)))

         if(holepsi.ne.0.)then ! Hole Initialization
! Calculate holespeed using vds component in projection dimension.
! Holeum is minus f drift speed relative to hole so holeum=-vds+holespeed:
            holespeed=holeum+vds(ispecies)*vdrift(id)
            if(myid.eq.0.and.ispecies.eq.hspecies)
     $           write(*,'(2a,i2,a,f7.3,a,f7.3)')' Hole particle'
     $           ,' initialization. Trapping id',id,'. Speed',holespeed
     $           ,' psi',holepsi
c Initialize u-range
            umax=4.
            du=umax/nu
            do i=-nu,nu
               u(i)=i*du
            enddo
            psi=holepsi
    ! Here holeum is in standard velocity units, but um is in tisq2 units.
            um=holeum/tisq2
c Hole (decay) length
            coshlen=holelen     ! Now set in cmdline +psi/2
c The flattop length holetoplen. Negligible for large negative values.     
            xmax=1.3*findxofphi(psi/(NPHI),psi,coshlen,holetoplen,0.,50.
     $           ,7)
            call f0Construct(nphi,psi,um,xmax,coshlen,holetoplen,
     $           phiarray,us,xofphi,den,denuntrap,dentrap,tilden,f0,u0)
            if(idebug.eq.1)then
               call autoplot(u0(-2*nphi),f0(-2*nphi),4*nphi+1)
               call axlabels('u0','f0')
               call pltend
            endif
         endif
      endif
      lastspecies=ispecies
!----------------------------------------
 1    r2=0.
      iend=ndims
      if(holepsi.ne.0)iend=ndims-1
      do ii=1,iend            ! Transverse or Total Position set.
         i=mod(id+ii-1,ndims)+1
         call ranlux(ran,1)
         fp=(indi(i)+ran)/float(nqblks(i))
         fp=max(.000001,min(.999999,fp))
         x_part(i,islot)=(1.-fp)*xmeshstart(i)+fp*xmeshend(i)
         if(holerad.ne.0)r2=r2+x_part(i,islot)**2
      enddo
      psiradfac=1.
      if(holepsi.ne.0.)then
         call ranlux(ran,1)
         fp=(indi(id)+ran)/float(nqblks(id))
         fp=max(.000001,min(.999999,fp))
         if(holerad.ne.0)psiradfac=exp(-r2/holerad**2)
! Hole density non-uniformity: transverse local value of peak potential.
         x_part(id,islot)=findxofran(fp,psiradfac*psi,coshlen,holetoplen
     $        ,xmeshstart(id),xmeshend(id),nbi)
      endif
      phi=psiradfac*phiofx(x_part(id,islot),psi,coshlen,holetoplen)

!                            ! Velocities
      do i=1,ndims
         if(i.eq.id)then
! Hole normal direction.
            if(abs(psi).lt.phimin.or.fp.lt.phimin.or.1-fp.lt.phimin)then
! In non-hole region. Shortcut to external distrib.
               x_part(ndims+i,islot)=tisq*gasdev(myid)
     $              +vds(ispecies)*vdrift(i)
            else
! Use cumf interpolation. We use the central psi value, but local phi.
               call GetDistribAtPhi(psi,um,nphi,f0,u0,phi,nu,u,f,cumf)
               call ranlux(fp,1)
               fp=fp*cumf(nu)
               ixp=interp(cumf,nup,fp,p)
               if(ixp.gt.0.and.ixp.lt.nup)then
                  x_part(ndims+i,islot)=
     $                 tisq2*((1-p+ixp)*u(-nu-1+ixp)+(p-ixp)*u(-nu+ixp))
     $                 +holespeed
               else
                  write(*,*)'placeqblk interpolation error',ixp,p
                  stop
               endif
            endif
         else
            x_part(ndims+i,islot)=tisqperp*gasdev(myid)
     $           +vds(ispecies)*vdrift(i)
         endif
      enddo
! The previous timestep length.
      x_part(idtp,islot)=0.
! Initialize the mesh fraction data in x_part.
      call partlocate(x_part(1,islot),ixp,xfrac,iregion,linmesh)
! This test rejects particles exactly on mesh boundary:
      if(.not.linmesh)then
         write(*,*)'Particle',islot,' outside mesh',x_part(1,islot)
         goto 1
      endif
      x_part(iflag,islot)=1

      if(ldouble)then ! Double particle initialization
         islot=islot+1
         do i=1,ndims
            x_part(i,islot)=x_part(i,islot-1)
            x_part(i+ndims,islot)=-x_part(i+ndims,islot-1)
            x_part(i+2*ndims,islot)=x_part(i+2*ndims,islot-1)
            x_part(iflag,islot)=1                        
         enddo
         x_part(idtp,islot)=0.         
      endif

      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Given the distribution function in the hole frame, which moves 
! relative to the external Maxwellian at speed um, at phi=psi: 
! f0 at an array over velocities
! u0, obtain the distribution function and its integral giving cumulative
! probability distribution at a lower potential phi over given u-array.
! The +- limits of u**2 must be > psi. (Still in hole frame)
      subroutine GetDistribAtPhi(psi,um,nphi,f0,u0,phi,nu,u,f,cumf)
      integer nphi,nu
      real f0(-2*nphi:2*nphi),u0(-2*nphi:2*nphi),psi
      real u(-nu:nu),f(-nu:nu),cumf(-nu:nu)
      real um
      
      if(phi.gt.psi.or.psi.ge.min(u(-nu)**2,u(nu)**2))then
         write(*,*)'GetDistriAtPhi error. psi   phi  u(nu)  u0max'
         write(*,*)psi,phi,u(nu),u0(2*nphi)
         stop
      endif
      sqpi=sqrt(3.1415926)
! Fill the f and cumf arrays.
      do i=-nu,nu
         if(u(i)**2-phi.gt.0)then  
! Passing particles must be done directly because their slope singularity
! at the separatrix undermines accuracy of reverse interpolation.
            uinf=sign(sqrt(u(i)**2-phi),u(i))
            f(i)=exp(-(uinf+um)**2)/sqpi
         else
! Trapped: Interpolate f0,u0 to find f0i=f0(u0i)
            u0i=sign(sqrt(u(i)**2+psi-phi),u(i))
!               write(*,*)i,ipos,pos,u0i,u(i)**2,psi,phi
            ipos=interp(u0,4*nphi+1,u0i,pos)
            if(ipos.le.0)then
               write(*,*)i,ipos,pos,u0i,u(i)**2,psi,phi
               stop 'GetDistribAtPhi interpolation error'
            endif
            ip=ipos-2*nphi-1    ! refer to -2*nphi:2*nphi
            f(i)=(1.-pos+ipos)*f0(ip) + (pos-ipos)*f0(ip+1)
         endif
         if(i.le.-nu)then
            cumf(-nu)=0.
         else
            cumf(i)=cumf(i-1)+0.5*(f(i-1)+f(i))*(u(i)-u(i-1))
         endif
      enddo

      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Construct the distribution f0 at u0, using BGKint.
      subroutine f0Construct(nphi,psi,um,xmax,coshlen,tl,
     $     phiarray,us,xofphi,den,denuntrap,dentrap,tilden,f0,u0)

      integer nphi
      real psi,um,xmax
      real phiarray(0:NPHI),us(0:NPHI),xofphi(0:NPHI)
      real den(0:NPHI),denuntrap(0:NPHI),dentrap(0:NPHI)
      real tilden(0:NPHI)
      real f0(-2*nphi:2*nphi),u0(-2*nphi:2*nphi)

      du=(4.-sqrt(psi))/nphi
      sqpi=sqrt(3.1415926)
c Call in BGKint such a way as to put into the negative trapped section
c because it returns f0 u0 in reverse order.
      call BGKint(nphi,psi,-um,xmax,coshlen,tl,
     $     phiarray,us,xofphi,den,denuntrap,dentrap,tilden
     $     ,f0(-nphi),u0(-nphi))
c Copy across to the positive trapped section.
      do i=1,nphi
c         write(*,*)i,u0(-i),f0(-i)
         u0(i)=u0(-i)
         f0(i)=f0(-i)
         u0(-i)=-u0(-i)
      enddo
c Fill in the untrapped distribution (maybe not needed?)
      do i=1,nphi
         ui=sqrt(psi)+i*du
         uinf=sqrt(abs(ui**2-psi))   ! prevent rounding to negative.
         u0(nphi+i) = ui
         f0(nphi+i) =exp(-(uinf+um)**2)/sqpi
         u0(-nphi-i)=-ui
         f0(-nphi-i)=exp(-(uinf-um)**2)/sqpi
      enddo
c Now f0(-2*nphi:2*nphi) is the entire distribution at phi=psi.

      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c**********************************************************************
c BGKint solves the integral equation to find the trapped
c distribution function for an electron hole. 
c The potential shape is specified via a function phiofx(x) which can
c be made whatever one wishes.
c The untrapped density is given by function untrappeddensimple(phi,um)
c Units of x are debyelengths, of potential Te/e, of time omega_p^{-1}
c But u is v/sqrt(2), i.e. normalized to sqrt(2Te/me).

      subroutine BGKint(nphi,psi,um,xmax,coshlen,tl,
     $     phi,us,xofphi,den,denuntrap,dentrap,tilden,f,u0)
c In:
c nphi is the number of phi (i.e. u^2) positions.
c psi the maximum, zero the minimum. um the maxwellian shift,
c coshlen and t1 are parameters passed to the routine phiofx 
c for calculating: phi, the potential (grid). us is the sqrt of phi. 
c xofphi is the corresponding (positive) postion. 
c Out:
c den is total density, denuntrap is the untrapped electron density,
c dentrap the trapped density, tilden the difference between trapped 
c density and a flat-top. 
c f(u0) is the distribution function, 
c u0=sqrt(psi-i*phistep)  =speed at x=0 when total energy is -phi. 
c But u0 runs from u0(0)=u_s to u0(nphi)=0. (i.e. reverse order)
      integer nphi
      real psi,um,xmax
      real phi(0:NPHI),us(0:NPHI),xofphi(0:NPHI)
      real den(0:NPHI),denuntrap(0:NPHI),dentrap(0:NPHI)
      real tilden(0:NPHI)
      real f(0:nphi),u0(0:nphi)

      real pi
      parameter (pi=3.1415926)
      parameter (nbi=20)
      real delx

      phistep=psi/(NPHI)
      delx=4.*xmax/nphi
      sphistep=sqrt(phistep)
      flatf=exp(-um**2)/sqrt(pi)
c 
      do i=0,NPHI
         phi(i)=i*phistep
      enddo
      do i=0,NPHI
c Find the xofphi by bisection.
         xofphi(i)=findxofphi(phi(i),psi,coshlen,tl,0.,xmax,nbi)
c Calculate the total density -d^2\phi/dx^2 as a function of potential,
c at the nodes.
         xc=xofphi(i)
         den(i)=1.+(phiofx(xc+delx,psi,coshlen,tl)
     $        +phiofx(xc-delx,psi,coshlen,tl)
     $        -2.*phiofx(xc,psi,coshlen,tl))/delx**2
         us(i)=sqrt(phi(i))
c Get the untrapped electron density at this potential and drift.
         denuntrap(i)=untrappeddensimple(phi(i),um)
         dentrap(i)=den(i)-denuntrap(i)
c Density difference c.f. flat:
         tilden(i)=dentrap(i)-2*us(i)*flatf
      enddo

c f(u) = (1/pi) \int_0^{psi-u^2} dn/d\phi d\phi/sqrt(\psi-u^2-phi).
c u^2=psi-i*phistep, phi=j*phistep, so sqrt -> (i-j)*psistep.
      u0(0)=sqrt(psi)
      f(0)=flatf
      do i=1,NPHI
         u0(i)=sqrt(psi-i*phistep)
         fi=0.
         do j=1,i
c We calculate the dndphi based upon \tilde f's density rather than on the
c total density, because this avoids big errors near the separatrix.
            dndphi=(tilden(j)-tilden(j-1))/phistep
            fi=fi+dndphi*2.*sphistep*(sqrt(i-j+1.)-sqrt(float(i-j)))
         enddo
         f(i)=fi/pi+flatf
      enddo

      end
c*********************************************************************
c Flattened sech^4 potential function.
      real function phiofx(x,psi,coshlen,toplen)
      real x,psi,coshlen,toplen
      et=exp(-toplen)
      xo=x/coshlen
      if(.not.abs(xo).le.10.)then
         phiofx=0.
      else
         phiofx=psi*(1.+et)
     $        /(1.+et*cosh(xo)**4)
      endif
      end
c*********************************************************************
c Derivative of phiofx
      real function derivphiofx(x,psi,coshlen,toplen)
      real x,psi,coshlen,toplen
      et=exp(-toplen)
      xo=x/coshlen
      if(.not.abs(xo).le.10.)then
         derivphiofx=0.
      else
         derivphiofx=-4.*psi/coshlen *et*(1.+et)*sinh(xo)*cosh(xo)**3
     $        /(1+et*cosh(xo)**4)**2
      endif
      end
c********************************************************************
      real function findxofphi(phiv,psi,coshlen,tl,xmin,xmax,nbi)
c Solve by bisection phiv=phiofx(xc,psi,coshlen,tl) and return xc.
c The intial x-range is [0,xmax]. Up to nbi bisections are allowed.
c This version uses a function, with three extra parameters
c psi and coshlen and tl. 
      real phiv,psi,coshlen,tl,xmin,xmax
      integer nbi
      real phiofx
      external phiofx
      xa=xmin
      xb=xmax
      xc=0
c This extra tiny value prevents rounding errors from causing outside
c range potentials.
      fva=phiofx(xa,psi,coshlen,tl)-phiv+1.e-7
c         fvb=phiofx(xb,psi,coshlen,tl)-phiv
c Allow a value beyond the xmax range so long as phiv is non-negative.
      fvb=-phiv
      if(fva*fvb.gt.0.)then
         write(*,*)xa,fva,xb,fvb,phiv,phiofx(xa,psi,coshlen,tl)
     $        ,phiofx(xb,psi,coshlen,tl)
         stop 'Potential outside range'
      endif
      do j=1,nbi
         xc=0.5*(xa+xb)
         fvc=phiofx(xc,psi,coshlen,tl)-phiv
         if(fvc.eq.0.)then
            goto 1
         elseif(sign(1.,fvc).eq.sign(1.,fva))then
            xa=xc
            fva=fvc
         else
            xb=xc
            fvb=fvc
         endif
      enddo
 1    continue
      findxofphi=xc
      end
c********************************************************************
      real function findxofran(ranf,psi,coshlen,tl,xmin,xmax,nbi)
c Find by bisection the x-position given the range fraction ranf s.t.
c   ranf=(xc-xmin)/(xmax-xmin)+derivphiofx(xc...) 
c and return xc. This expression is the cumulative probability.
c Up to nbi bisections are allowed.
c Shortcuts are taken if potential is less than a range of relevance.
      real ranf,psi,coshlen,tl,xmin,xmax
      integer nbi
      real phimin
      parameter (phimin=.01)
      real phiofx
      external phiofx
      if(abs(psi).lt.phimin.or.ranf.lt.phimin.or.1-ranf.lt.phimin)then
 ! Effectively linear shortcut.
         findxofran=xmin+ranf*(xmax-xmin)
         return
      endif
      xa=xmin
      xb=xmax
      xc=0.
      fva=(xa-xmin+derivphiofx(xa,psi,coshlen,tl))/(xmax-xmin)-ranf
      fvb=(xb-xmin+derivphiofx(xb,psi,coshlen,tl))/(xmax-xmin)-ranf
c Allow a value beyond the xmax range so long as ranf is non-negative.
      if(fva*fvb.gt.0.)then
         write(*,*)xa,fva,xb,fvb,phiv,phiofx(xa,psi,coshlen,tl)
     $        ,phiofx(xb,psi,coshlen,tl)
         stop 'Potential outside range'
      endif
      do j=1,nbi
         xc=0.5*(xa+xb)
         fvc=(xc-xmin+derivphiofx(xc,psi,coshlen,tl))/(xmax-xmin)-ranf
         if(fvc.eq.0.)then
            goto 1
         elseif(sign(1.,fvc).eq.sign(1.,fva))then
            xa=xc
            fva=fvc
         else
            xb=xc
            fvb=fvc
         endif
      enddo
 1    continue
      findxofran=xc
      end
c*********************************************************************
c Untrapped density function for a shifted Maxwellian
c  untrappedden(\phi,u_m) = {2\over \sqrt{\pi}}\int_{\sqrt{\phi}}^\infty
c     \exp(-u^2+\phi-u_m^2)\cosh(2u_m\sqrt{u^2-\phi}) du.
c   = {1\over\sqrt{\pi}} \int_{\sqrt{\phi}}^\infty \sum_{\pm}
c             \exp(-[\pm\sqrt{u^2-\phi}-u_m]^2)du
c written as f = \int g du.

c This is the density relative to the background density of untrapped
c particles having a Maxwellian background shifted by a speed um,
c measured in units of sqrt(2T/m), at a potential energy -phi
c measured in units of T/e.

c The integration is carried out on uniform grid, initially with 
c velocity spacing dui. The spacing is halved until relative difference
c from the prior integral is less than df.
c [Some sort of Gaussian integration would probably be faster but tricky
c because of the infinite range.]
c umswitch is currently unused.
      real function untrappeddensimple(phi,um)
      real phi,um

      parameter (dui=.04,df=1.e-5,np=10000,umswitch=4.)

c silence warnings not really needed.
      f=0.
      fp=0.
      fm=0.
      v=0.
      gm1=0.
      f1=0
      g1=1.
      g1p=1.
      g1m=1.
      gp=0.
      gm=0.

      um2=um**2
      du=dui
      nstepit=7

      if(phi.lt.0)stop 'Negative phi in untrappedden call not allowed'
c Iterate over steps sizes
      do k=1,nstepit
c Integrate
         do i=0,np
            u=sqrt(phi)+i*du
            u2=u**2
            if(i.eq.0)then
               g=exp(-um2)
               gp=g
               gm=g
               g1=g
               if(g.eq.0)stop 'untrappedden error um2 overflow'
               f=0.
               fp=0.
               fm=0.
            else
c This is simply more reliable:
               gp=exp(-(sqrt(u2-phi)+abs(um))**2)
               gm=exp(-(sqrt(u2-phi)-abs(um))**2)
               g=0.5*(gp+gm)
c This implicitly multiplies by 2:               
               f=f+(u-v)*(g+gm1)
               fp=fp+(u-v)*(gp+g1p)*0.5
               fm=fm+(u-v)*(gm+g1m)*0.5
            endif
            if(.not.g.ge.0 .or. .not.g.lt.1.e30)then
c Error trap
               write(*,*)'u,um,phi,g error',u,um,phi,g
c               stop
               g=0.
            endif
            f=fp+fm
            gm1=g
            g1p=gp
            g1m=gm
            v2=u2
            v=u
c If new contributions are negligible, break
            if(g.lt.g1*df)goto 1
            if(gp.lt.g1*df.and.gm.lt.g1*df)
     $           write(*,*)'g,gp,gm,g1',g,gp,gm,g1
         enddo
         write(*,*)'untrappedden exhausted number of steps np',i
 1       continue
c If converged, break
         if(abs(f1-f).lt.df)goto 2
         f1=f
         du=du/2.
      enddo
      write(*,*)'untrappedden exhausted step iterations',k,f1,f,du
      write(*,*)'Steps, du, u, f, um',i,du,u,f,um
      write(*,*)u,g,f
 2    continue

      untrappeddensimple=f/sqrt(3.1415926)

      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Introduce a plane-wave displacement into the initialized particles
! Only for use with a uniform domain in which particle displacement
! can be consistently treated as periodic: xi= dispvec exp(i n.x/L)
! API: Is by the input array wavespec(nwspec) of which the elements
!      1    Whether anything is to be done (zero means no).
!      2:1+ndims          Mode numbers in 3 coordinate directions.
!      2+ndims:1+2*ndims  Magnitude and direction of dispvec
! Call only after qinit to enforce meaningful logic.
      subroutine wavedisplace(wavespec)
      implicit none
      include 'ndimsdecl.f'
      include 'meshcom.f'
      include 'partcom.f'
      integer nwspec
      parameter (nwspec=2*ndims+1)
      real wavespec(nwspec)
      real kw(ndims)
      logical linmesh
      real cosphase,phase,xfrac(ndims)
      integer i,id,iregion,ixp(ndims),js

!      write(*,*)'Applying initial wave',wavespec
      if(wavespec(1).eq.0)return
      do id=1,ndims    ! Ensure precise periodicity
         kw(id)=nint(wavespec(1+id))*2.*3.1415926
     $        /(xmeshend(id)-xmeshstart(id))
         write(*,*)wavespec(1+id),xmeshend(id),xmeshstart(id),kw(id)
      enddo
      do js=1,nspecies
         do i=iicparta(js),iocparta(js)
            if(x_part(iflag,i).ne.0)then
               phase=0.                ! calculate phase
               do id=1,ndims
                  phase=phase+kw(id)*x_part(id,i)
               enddo
               cosphase=cos(phase)
               if(.not.abs(cosphase).le.1)write(*,*)'cosphase',cosphase
               do id=1,ndims           ! Displace particles
                  x_part(id,i)=x_part(id,i)
     $                 +wavespec(1+ndims+id)*cosphase
               enddo
! Initialize the mesh fraction data in x_part. This also checks periodicity.
               call partlocate(x_part(1,i),ixp,xfrac,iregion,linmesh)
               if(.not.linmesh)then
                  write(*,*)'wavedisplace linmesh ERROR',ixp,xfrac
     $                 ,iregion,linmesh
               endif
            endif
         enddo
      enddo
      

      end
