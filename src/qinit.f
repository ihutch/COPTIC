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

! No special orbits.
      norbits=0
      i1=1

! Point to the bottom of the particle stack for start of species 1.
      iicparta(1)=1
      do ispecies=1,nspecies
! nqbset parameters adjusted only from 1 parameter: nqblkmax [-pi...]
         do i=1,ndims
            mlen=ixnp(i+1)-ixnp(i)
            if(mlen-2.eq.1)then ! Single cell sides need no quieting.
               nqbset(i)=1
            else
               nqbset(i)=max(1,nqblkmax) ! Quiet short sides.
!         nqbset(i)=max(1,min(nqblkmax,mlen-2)) ! Don't quiet short.
            endif
            nqblks(i)=nqbset(i)
         enddo
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
! Scale the number of particles for higher slots to first(? NOT)
         nparta(ispecies)=nparta(ispecies)/numratioa(ispecies)
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
! Apply initial wave displacement, now regardless of actual periodicity
!      if(.not.lnotallp)call wavedisplace(wavespec)
      call wavedisplace(wavespec)
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
      include 'plascom.f'
      include 'partcom.f'
      include 'fvcom.f'
      integer nqblks(ndims),islot,islotmax,ispecies

      integer iview(3,ndims),indi(ndims)     !Iterator
      integer nqbt,iworkspecies
      logical ldouble ! Whether to place 2 counter-velocity particles.
      data ldouble/.true./iworkspecies/0/

      if(iworkspecies.ne.ispecies)then 
         ldouble=.true. ! Try to do double setting for each species.
         iworkspecies=ispecies
! Double particle initialization is possible only if the distribution
! has reflectional symmetry, which requires vdrift and vhole to be equal.
! At the moment we disallow double if any drift is non-zero.
! Double helps by forcing the electron drift to exactly zero.
! So it might help to generalize this to allow ldouble with vds!=0.
         if(vds(ispecies).ne.0.or.holeum.ne.0..or.holerad.ne.0.
     $        .or.nc(ispecies).ne.1)then
            ldouble=.false.
            if(myid.eq.0)then
               write(*,*)'Finite drift. No double particle placement.'
            endif
         endif
      endif
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
! Fill using iterator
      icomplete=mditerator(ndims,iview,indi,4,nqblks)
 1     continue
       do i=1,nfills
! Place one particle in block indi.
         call placeqblk(nqblks,indi,islot,ispecies,ldouble)
         islot=islot+1
       enddo
       if(mditerator(ndims,iview,indi,0,nqblks).eq.0)goto 1
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
      include 'fvcom.f'
      integer nqblks(ndims),indi(ndims),islot,ispecies
      logical linmesh,ldouble
      integer nbi,ixp(ndims)
      parameter (nbi=16)
! Hole-related parameters:
      real phimin,phiprev
      parameter (phimin=.01)
      integer nphi,nu,nup
      parameter (nphi=50,nu=100,nup=2*nu+1)
      real psi,um,xmax
      real phiarray(0:NPHI),us(0:NPHI),xofphi(0:NPHI)
      real den(0:NPHI),denuntrap(0:NPHI),dentrap(0:NPHI)
      real tilden(0:NPHI),denion(0:nphi)
      real pden(-nphi-1:nphi+1),cump(-nphi-1:nphi+1)
      real xpden(-nphi-1:nphi+1)
      real f0(-2*nphi:2*nphi),u0(-2*nphi:2*nphi),cumf(-nu:nu)
!      real f1(-2*nphi:2*nphi),u1(-2*nphi:2*nphi),cumf1(-nu:nu)
      real u(-nu:nu),f(-nu:nu),du
      real tisq,tisq2,tisqperp,umax,coshlen,kp2,kp2find
      real rhov(ndims),xfrac(ndims)
      character*30 string
      integer lastspecies,id
      data lastspecies/0/
      save

!----------------------------------------
      if(ispecies.ne.lastspecies)then 
! First call Initialization of this species.
! These values must be set even for zero hole depth.
         phiprev=-1.   ! Never skip the first fvhill call.
         idebug=0
!         idebug=2      ! Plot distribution setup
!         idebug=1      ! Write out setup arrays.
         do i=1,ndims
            if(Bfield(i).ne.0.)goto 2
         enddo
         i=1
 2       id=i
! The tisq numbers are really thermal velocities accounting for mass.
         tisq=sqrt(Ts(ispecies)*abs(eoverms(ispecies)))
         tisq2=sqrt(2.)*tisq
         tisqperp=sqrt(Tperps(ispecies)*abs(eoverms(ispecies)))

         if(holepsi.ne.0..and.ispecies.eq.hspecies)then ! Hole Initialization
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
c Transverse k^2 
            kp2=0.
            if(holerad.ne.0)kp2=4./holerad**2
            kp2find=kp2  ! Central value. Maybe to be changed later.
c The flattop length holetoplen. Negligible for large negative values.     
            xmax=12.
            call f0Construct(nphi,psi,um,xmax,coshlen,holetoplen,
     $           phiarray,us,xofphi,den,denuntrap,dentrap,tilden,f0,u0
     $           ,kp2,denion)
            if(idebug.eq.1.and.myid.eq.0)then
               write(*,*)'xofphi  phiarray  den  denion  denuntrap',
     $              ' dentrap  tilden'
               write(*,'(7f8.4)')(xofphi(ik),phiarray(ik),den(ik)
     $              ,denion(ik),denuntrap(ik),dentrap(ik),tilden(ik)
     $              ,ik=0 ,nphi)
            endif
            if(idebug.eq.2.and.myid.eq.0)then
               call multiframe(2,1,3)
               call autoplot(sqrt(2.)*u0,f0/sqrt(2.),4
     $              *nphi+1)
               call axis2
               call axlabels('v!ds!d/(T!ds!d/m!ds!d)!u1/2!u',
     $              'f(v)!AX!@(T!ds!d/m!ds!d)!u1/2!u')
               call legendline(.7,.9,0,' electrons(1)')
               if(nc(2).ne.0)call legendline(.05,.9,258,' ions, s=')
               call winset(.true.)
               do is=2,3
                  if(nc(is).ne.0)then
                     call dashset(is)
                     call vofvinit(is)
                     call fvhill(nc(is),dcc(1,is),vsc(1,is),vtc(1,is),0.
     $                    ,delphi,0.,1,.true.,nofv,vofv,fofv,Pfofv)
                     call polyline(vofv,fofv,nofv)
                     call iwrite(is,iwidth,string)
                     call legendline(.05,1.-.1*is,0,' '
     $                    //string(1:iwidth))
                  endif
               enddo
               call dashset(0)
               call autoplot(xofphi,den,nphi+1)
               call legendline(.7,.8,0,' n!de!d (adj)')
               call axis2
               call axlabels('x','n')
               call dashset(1)
               if(nc(2).ne.0)call legendline(.7,.9,0,' n!di!d')
               call polyline(xofphi,denion,nphi+1)
               call dashset(2)
               call polyline(xofphi,denuntrap+dentrap,nphi+1)
               call legendline(.7,.7,0,' n!det!d+n!dep!d')
               call dashset(0)
               call polyline([0.,xofphi(0)],[1.,1.],2)
               call pltend
               call multiframe(0,0,0)
            endif
         endif
      endif
      lastspecies=ispecies
!----------------------------------------
 1    r2=0.
! Gyro velocity setting and gyroradius (drift added later)
      rhov(id)=0.
      i1=mod(id,ndims)+1   ! Transverse directions.
      i2=mod(id+1,ndims)+1
      do i=i1,i2
            vp=tisqperp*gasdev(myid) 
            x_part(ndims+i,islot)=vp
      enddo
! v=x^qB/m => v^qB/m=-x(qB/m)^2; so gyro-vector rhov=-v^Bt/(eoverm*Bt^2)
      rhov(i1)=-x_part(ndims+i2,islot)/(eoverms(ispecies)*Bt) 
      rhov(i2)= x_part(ndims+i1,islot)/(eoverms(ispecies)*Bt)
!----------------------------------------
! Position setting, all dims (reset x_id later if necessary)
      do i=1,ndims
         call ranlux(ran,1)
         fp=(indi(i)+ran)/float(nqblks(i)) ! The fraction of total mesh
         fp=max(.000001,min(.999999,fp))
         x_part(i,islot)=(1.-fp)*xmeshstart(i)+fp*xmeshend(i)
         if(i.ne.id)then
            r2=r2+x_part(i,islot)**2  ! Perpendicular particle radius^2
! Alternative gyrocenter radius^2: makes peak lower.
!            r2=r2+(x_part(i,islot)-rhov(i))**2
         else
            fpid=fp   ! Save the random position for parallel
         endif
      enddo
!----------------------------------------
      psiradfac=1.
      if(holepsi.ne.0)then     ! Maybe Reset normal position
         if(ispecies.eq.hspecies)then
            kp2find=0.
            if(holerad.ne.0)then ! Account for transverse variation
               psiradfac=exp(-r2/holerad**2)
               kp2find=(4./holerad**2)*(1-r2/holerad**2)
            endif
! Hole parallel nonuniformity at transverse local value of peak potential.
            if(.not.lfv)then  ! Set electrons for Ion distrib uniform
               x_part(id,islot)=findxofran(fpid,psiradfac*psi,coshlen
     $              ,holetoplen,xmeshstart(id),xmeshend(id),nbi,kp2find)
            else        ! Set electrons accounting for ions non-uniform
               isw=0
               foundx=findxofcum(fpid,nphi,den,pden,cump,xpden,xofphi
     $              ,isw,im,xf,xmeshstart(id),xmeshend(id))
               x_part(id,islot)=foundx
            endif
         elseif(lfv)then   ! Ion distribution nonuniform. (MultiGaussian).
            isw=0
            foundx=findxofcum(fpid,nphi,denion,pden,cump,xpden,xofphi
     $           ,isw,im,xf,xmeshstart(id),xmeshend(id))
            x_part(id,islot)=foundx
         endif
      endif
!----------------------------------------
! Hole normal (i.e. parallel) velocity setting
      phi=psiradfac*phiofx(x_part(id,islot),psi,coshlen,holetoplen)
      if(ispecies.ne.hspecies)then  ! Ions
         if(nc(ispecies).ne.0)then  ! MultiGausian
            if(phiprev.eq.-1.and.myid.eq.0)then
               write(*,'(a,i2,a,i2)')'MultiGaussian Species:',ispecies
     $              ,' components:',nc(ispecies)
            endif
            if(abs(phi-phiprev).gt.phimin)then ! Construct Pofv
               isigma=int(sign(1.,x_part(id,islot))) ! Strange result.
               isigma=1
               delphi=0.
               call vofvinit(ispecies)
               call fvhill(nc(ispecies),dcc(1,ispecies),vsc(1,ispecies)
     $              ,vtc(1,ispecies),holepsi,delphi,phi,isigma,.true.
     $              ,nofv,vofv,fofv,Pfofv)
               isw=0     ! Need to recalculate cumulative distrib.
               phiprev=phi
            endif
            call ranlux(fp,1)
            fp=fp*Pfofv(nofv)
            fps=fp
            ixpx=interp(Pfofv,nofv,fp,p)
            if(ixpx.gt.0.and.ixpx.lt.nofv)then
! Multigauss velocity is in units of sqrt(Tr/ms) convert to sqrt(Tr/m1)
!               x_part(ndims+id,islot)=tisq*((1-p+ixpx)*vofv(-1+ixpx)+(p
               x_part(ndims+id,islot)=sqrt(abs(eoverms(ispecies)))*((1-p
     $             +ixpx)*vofv(-1+ixpx)+(p-ixpx) *vofv(ixpx)) +holespeed
            else
               write(*,'(a,i5,7f8.4)')'placeqblk MultiGauss intp error'
     $              ,ixpx,fps,p,fp,phiprev
!               write(*,'(10f8.4)')Pofv
               stop
            endif
            
         else                       ! Shifted Single Gaussian
               x_part(ndims+id,islot)=tisq*(gasdev(myid)
     $              +vds(ispecies)*vdrifts(id,ispecies))
         endif
      elseif(abs(phi).lt.phimin.or.fp.lt.phimin.or.1-fp.lt.phimin)then
! Hole species but non-hole region. Shortcut to external distrib.
            x_part(ndims+id,islot)=tisq*(gasdev(myid)
     $           +vds(ispecies)*vdrifts(id,ispecies))
      else
! Use cumf interpolation. We use the central psi value, but local phi.
         call GetDistribAtPhi(psi,um,nphi,f0,u0,phi,nu,u,f,cumf)
         call ranlux(fp,1)
         fp=fp*cumf(nu)
         ixpx=interp(cumf,nup,fp,p)
         if(ixpx.gt.0.and.ixpx.lt.nup)then
            x_part(ndims+id,islot)=
     $           tisq2*((1-p+ixpx)*u(-nu-1+ixpx)+(p-ixpx)*u(-nu+ixpx))
     $           +holespeed
         else
            write(*,*)'placeqblk interpolation error',ixpx,p
            stop
         endif

      endif
! Add transverse ExB velocity  Was only on holespecies. Now on all.
      if(holerad.ne.0)then
         Etot=(2./holerad**2)*sqrt(r2)*phi ! Gradient of Gaussian
         E1=Etot*x_part(i1,islot)/sqrt(r2) ! Project transverse
         E2=Etot*x_part(i2,islot)/sqrt(r2)
         x_part(ndims+i1,islot)=x_part(ndims+i1,islot)+E2/Bt            
         x_part(ndims+i2,islot)=x_part(ndims+i2,islot)-E1/Bt            
      endif
! Add background transverse drift.
      do i=mod(id,ndims)+1,mod(id+1,ndims)+1
            x_part(ndims+i,islot)=x_part(ndims+i,islot)            
     $           +tisq*vds(ispecies)*vdrifts(i,ispecies)
      enddo
! End of parallel velocity initialization.
!----------------------------------------
! The previous timestep length.
      x_part(idtp,islot)=0.
! Initialize the mesh fraction data in x_part.
      call partlocate(x_part(1,islot),ixp,xfrac,iregion,linmesh)
! This test rejects particles exactly on mesh boundary:
! It is broken when -fbounds-check is used!!!!!
      if(.not.linmesh)then 
         write(*,*)'qinit Particle',islot,' outside mesh',x_part(1
     $        ,islot),linmesh
!         stop
         goto 1
      endif
      x_part(iflag,islot)=1

      if(ldouble)then ! Double particle initialization
         islot=islot+1
         do i=1,ndims
            x_part(i,islot)=x_part(i,islot-1)
! Reverse velocity accounting for drift possible only when there
! is no drift of the hole relative to the background.
!            x_part(i+ndims,islot)=
!     $           2.*vds(ispecies)*vdrift(i)-x_part(i+ndims,islot-1)
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
      subroutine f0Construct(nphi,psi,um,xmax,coshlen,tl, phiarray,us
     $     ,xofphi,den,denuntrap,dentrap,tilden,f0,u0,kp2,denion)
      include 'ndimsdecl.f'
      include 'fvcom.f'   
      include 'partcom.f' ! Needed for nspecies.

      integer nphi
      real psi,um,xmax,kp2
      real phiarray(0:NPHI),us(0:NPHI),xofphi(0:NPHI),denion(0:nphi)
      real den(0:NPHI),denuntrap(0:NPHI),dentrap(0:NPHI)
      real tilden(0:NPHI)
      real f0(-2*nphi:2*nphi),u0(-2*nphi:2*nphi)
      external denionfun

      du=(4.-sqrt(psi))/nphi
      sqpi=sqrt(3.1415926)
c Call BGKint in such a way as to put into the negative trapped section
c because it returns f0 u0 in reverse order.
      if(lfv)then
!         write(*,*)'               f0construct non-Maxwellian version'
         call BGKintnew(nphi,psi,-um,xmax,coshlen,tl,
     $        phiarray,us,xofphi,den,denuntrap,dentrap,tilden
     $        ,f0(-nphi),u0(-nphi),kp2,denionfun,denion)
      else  ! Zero ion response version.
!         write(*,*)'               f0construct ions uniform'
         call BGKint(nphi,psi,-um,xmax,coshlen,tl,
     $        phiarray,us,xofphi,den,denuntrap,dentrap,tilden
     $        ,f0(-nphi),u0(-nphi),kp2)
      endif
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
c**********************************************************************
c BGKintnew solves the integral equation to find the trapped and hence
c entire distribution function for an electron hole.
c The potential shape is specified via a function phiofx(x) which can
c be made whatever one wishes. Arrays are constructed from nphi,psi.
c Ion density is specified by external function denionfun(phi)
c Units of x are debyelengths, of potential Te/e, of time omega_p^{-1}
c But u is v/sqrt(2), i.e. normalized to sqrt(2Te/me).

      subroutine BGKintnew(nphi,psi,um,xmax,coshlen,tl,phi,us,xofphi
     $     ,den,denuntrap,dentrap,tilden,f,u0,kp2,denionfun,denion)
c In:
c  (0:nphi) is the array length of phi (i.e. u^2) positions.
c  psi the maximum, zero the minimum potential.
c  um the maxwellian shift in units of sqrt(2Te/me).

c  coshlen and t1 are parameters passed to the routine phiofx 
c  for calculating: phi, the potential (grid). 
c  denionfun: the function of (phi,isigma) that returns ion density
c Out:
c  phi is the equally spaced potential array, running from 0 to psi.
c  us is the sqrt of phi, the separatrix velocity.
c  xofphi is the corresponding (positive) position.
c  den is total electron density, denuntrap is the untrapped e-density,
c  dentrap the trapped e-density, tilden the difference between trapped 
c  density and a flat-top. 
c  denion is the ion density
c  f(u) is the electron distribution function, as a function of u at x=0:
c  u0=sqrt(psi-i*phistep) 
c Calls:
c  function denionfun(phi,isigma)
c  function phiofx(x,psi,...)
c  function findxofphi(phi,psi,...)
c  function untrapden(phi,um)
      integer NPHI
      real psi,um,xmax
      real phi(0:NPHI),us(0:NPHI),xofphi(0:NPHI),denion(0:nphi)
      real den(0:NPHI),denuntrap(0:NPHI),dentrap(0:NPHI),tilden(0:NPHI)
      real f(0:NPHI),u0(0:NPHI),kp2
      real pi
      parameter (pi=3.1415926)
      parameter (nbi=20)
      real delx
      external denionfun

      phistep=psi/(NPHI)
!      delx=4.*xmax/nphi
      delx=4.*xmax/nphi
      sphistep=sqrt(phistep)
      flatf=exp(-um**2)/sqrt(pi)
      isigma=1 ! Ignoring asymmetries for now.
c 
      denionint=0.
      denint=0.
      do i=0,NPHI
         phi(i)=i*phistep
      enddo
      do i=0,NPHI
c Find the xofphi by bisection.
         xofphi(i)=findxofphi(phi(i),psi,coshlen,tl,0.,xmax,nbi)
c Calculate the total density -d^2\phi/dx^2 as a function of potential,
c at the nodes.
         xc=xofphi(i)
c ne=ni+d^2\phi/dx^2
         denion(i)=denionfun(phi(i),isigma)
         den(i)=denion(i)+(phiofx(xc+delx,psi,coshlen,tl)
     $        +phiofx(xc-delx,psi,coshlen,tl)
     $        -2.*phiofx(xc,psi,coshlen,tl))/delx**2
     $        -2.*kp2*phi(i)  ! Transverse k term set from holerad.
         us(i)=sqrt(phi(i))
c Get the untrapped electron density at this potential and drift.
         denuntrap(i)=untrappeddensimple(phi(i),um)
         dentrap(i)=den(i)-denuntrap(i)
c Density difference c.f. flat:
         tilden(i)=dentrap(i)-2*us(i)*flatf
!         if(i.eq.nphi)xofphi(i)=0.
         if(i.gt.0)then
            denint=denint-0.5*(den(i)+den(i-1))*(xofphi(i)-xofphi(i-1))
            denionint=denionint-0.5*(denion(i)+denion(i-1))*(xofphi(i)
     $           -xofphi(i-1))
         endif
      enddo
! This section ensures ni=ne in the rest of the plasma.
!      write(*,*)'denint,denionint',denint,denionint
      endcont=abs(xofphi(0)-xofphi(1))*0.5
      adjfac=(denionint-endcont)/(denint-endcont*den(0))
      do i=1,nphi  ! Rescale den to denion
         den(i)=den(i)*adjfac
      enddo
      den(0)=1.
c f(u) = (1/pi) \int_0^{psi-u^2} dn/d\phi d\phi/sqrt(\psi-u^2-phi).
c u^2=psi-i*phistep, phi=j*phistep, so sqrt -> (i-j)*psistep.
      u0(0)=sqrt(psi)
      f(0)=flatf
      do i=1,NPHI
         u0(i)=sqrt((nphi-i)*phistep)
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
      real function denionfun(phi,isigma)
      include 'ndimsdecl.f'
      include 'fvcom.f'
      include 'partcom.f'
      logical linit,ldummy
      data linit/.false./

      ldummy=.true.
      do i=1,nspecies ! Might not be necessary if not called.
         if(i.ne.hspecies.and.nc(i).ne.0)ldummy=.false.
      enddo
      if(ldummy)then
c Dummy function that just returns 1.
         denionfun=1.
         return
      else
! Not now needed. vofvinit is called in coptic at start.
!         if(.not.linit)then
!            call vofvinit(2)
!         endif
         denionfun=0.
         delphi=0
!         isigma=1               ! For now ignore asymmetry.
         do i=1,nspecies
            if(i.ne.hspecies.and.nc(i).ne.0)then
               call vofvinit(i)
               call fvhill(nc(i),dcc(1,i),vsc(1,i),vtc(1,i),holepsi
     $              ,delphi,phi,isigma,.true.,nofv,vofv,fofv,Pfofv)
               denionfun=denionfun+Pfofv(nofv)
!               write(*,*)'Species',i,' phi',phi,' denionfun',denionfun
            endif
         enddo
      endif
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine vofvinit(isp)
! Initialize velocity range supposing hole is stationary. (vh=0.)
      include 'ndimsdecl.f'
      include 'fvcom.f'
      include 'partcom.f'
      integer isp
      vh=0.
      vmax=4.4*maxval(vtc(:,isp))
     $   +max(abs(maxval(vsc(:,isp))-vh),abs(minval(vsc(:,isp))-vh))
!      write(*,*)'isp=',isp,' vmax=',vmax
      do i=1,nofv
         vofv(i)=-vmax+2.*vmax*(i-1.)/(nofv-1.)
      enddo
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate velocity distributions and densities for shifted Maxwellian
! distant distributions on a single-humped repelling potential hill of
! maximum potential phimax, at potential phi, on side isigma of the hill.
! Also return \int fdv in Pfofv. 
! On entry:
! nc is the number of Maxwellian components 
! dc(nc), vs(nc), vt(nc) is each component's density, vshift, sqrt(T/M)
! phimax is hill's peak potential relative to mean distant potential
! delphi is potential difference phi(+inf)-phi(-inf) across the hill. 
! phi is the potential at which to find the distribution.
! isigma is the side of the hill on which this potential lies.
! vofv(nofv) is the monotonic array of velocities to use.
! All velocities are relative to the hill's rest frame.
! On exit:
! fofv(nofv) is the total velocity distribution at vofv
! Pfofv(nofv) is the integral of fofv dv. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine fvhill(nc,dc,vs,vt,phimax,delphi,phi,isigma,local,nofv
     $     ,vofv,fofv,Pfofv)
      real :: dc(nc),vs(nc),vt(nc),phimax,phi,delphi
      integer :: isigma, nofv
      real, dimension(nofv) :: vofv, fofv, Pfofv
      logical :: local
      
      if(isigma.ne.1.and.isigma.ne.-1)stop 'isigma must be +-1'
!     Decide the meaning of f(v) when delphi !=0.  
      if(local)then             ! f(v) is relative to different phi_infty
         phix2=2.*(phi    -isigma*delphi/2.) ! 2*(phi   -phiinf)
         phimx2=2.*(phimax-isigma*delphi/2.) ! 2*(phimax-phiinf)
         delphix2=2.*isigma*delphi  
         if(phix2.lt.0)then
            write(*,*)'Impossible phi (<phiinf)',phi,isigma,delphi
            stop
         endif
      else                      ! f(v)= f( sigma*sqrt(v^2+2*phi) )
         phix2=2.*phi
         phimx2=2.*phimax
         delphix2=0.
      endif
      if(phimx2.lt.phix2)then   ! Move the threshold far away.
         vthresh=40.*(vofv(nofv)-vofv(1))
      else
         vthresh=isigma*sqrt(phimx2-phix2)
      endif
      Pfofv=0.
      fofv=0.
      vdiffm=vofv(1)-vthresh
      fim=0.
      pam=0.
      pa=0.
      sqtwopi=sqrt(2.*3.1415926)
      do i=1,nofv
         vsign=sign(1.,vofv(i))
         vdiff=vofv(i)-vthresh
         enx2=phix2+vofv(i)**2
         if(int(vsign).ne.isigma)then ! Moving inward
            vinf=vsign*sqrt(max(0.,enx2))
         elseif((phix2+vofv(i)**2).gt.phimx2)then ! Passing
            vinf=vsign*sqrt(max(0.,delphix2+enx2))
         else                   ! Reflected
            vinf=-vsign*sqrt(max(0.,enx2))
         endif
         fi=0.
         do j=1,nc              ! f(v)=finf(vinf)
            fi=fi+dc(j)*exp(-0.5*(vinf-vs(j))**2/vt(j)**2)
     $           /(vt(j)*sqtwopi)
         enddo
         fofv(i)=fi
         if(i.gt.1)pa=pam+ (abs(vdiffm)*fim+abs(vdiff)*fi)
     $        /(abs(vdiffm)+abs(vdiff))*(vofv(i)-vofv(i-1))
!        Actually using this interpolation always smooths ftrapped glitches.
!        But I am not 100% certain it is always justified.
         Pfofv(i)=Pfofv(i)+pa
         if(.not.Pfofv(i).lt.1.e30)then
            write(*,*)'fvhill Bad Pfofv',i,Pfofv(i),pam,pa,vdiffm,vdiff
     $           ,vofv(i),vofv(i-1)
            stop
         endif
         pam=pa
         vdiffm=vdiff
         fim=fi
      enddo
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c**********************************************************************
c BGKint solves the integral equation to find the trapped
c distribution function for an electron hole. 
c The potential shape is specified via a function phiofx(x) which can
c be made whatever one wishes.
c The untrapped density is given by function untrappeddensimple(phi,um)
c Units of x are debyelengths, of potential Te/e, of time omega_p^{-1}
c But u is v/sqrt(2), i.e. normalized to sqrt(2Te/me).

      subroutine BGKint(nphi,psi,um,xmax,coshlen,tl,phi,us,xofphi,den
     $     ,denuntrap,dentrap,tilden,f,u0,kp2)
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
      real psi,um,xmax,kp2
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
c Calculate the total density d^2\phi/dx^2 as a function of potential,
c at the nodes.
         xc=xofphi(i)
         den(i)=1.+(phiofx(xc+delx,psi,coshlen,tl)
     $        +phiofx(xc-delx,psi,coshlen,tl)
     $        -2.*phiofx(xc,psi,coshlen,tl))/delx**2
     $        -2.*kp2*phi(i)  ! Transverse k term set from holerad.
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
      include 'ndimsdecl.f'
      include 'partcom.f'
c Factor by which to flatten f at separatrix; gfac=0.6 is flat.
      gfac=holegfac
      et=exp(min(-toplen,10.))
      xo=x/coshlen
      if(.not.abs(xo).le.10.)then
         phiofx=0.
      else
         phishape=(1.+et)/(1.+et*cosh(xo)**4)
         phiofx=psi*phishape*(1.+(phishape-1.)*gfac)
      endif
      end
c*********************************************************************
c Derivative of phiofx
      real function derivphiofx(x,psi,coshlen,toplen)
      real x,psi,coshlen,toplen
      include 'ndimsdecl.f'
      include 'partcom.f'
      gfac=holegfac
      et=exp(min(-toplen,10.))
      xo=x/coshlen
      if(.not.abs(xo).le.10.)then
         derivphiofx=0.
      else
         coshxo=cosh(xo)
         phishape=(1.+et)/(1.+et*coshxo**4)
         derivphiofshape=-4./coshlen *et*(1.+et)*sinh(xo)*coshxo**3
     $        /(1+et*coshxo**4)**2
         derivphiofx=psi*derivphiofshape*(1.+(2.*phishape-1.)*gfac)
      endif
      end
c********************************************************************
c Integral of phiofx, assumed symmetric, using trapezoidal interpolation.
c It uses external phiofx function
      real function intphiofx(x,psi,coshlen,toplen)
      real x,psi,coshlen,toplen
      integer nx
      parameter (nx=200)
      real phiint(nx),cl,tl,xmax
      external phiofx
      data cl/0/tl/0/xmax/0/
      save phiint,cl,tl,xmax
      if(coshlen.ne.cl.or.toplen.ne.tl)then ! (Re)Initialize
         cl=coshlen
         tl=toplen
         xmax=tl+3.5*cl  ! Choose xmax range.
         phiint(1)=0.
         phiprev=1.
         phiint(1)=0.
         dx=xmax/(nx-1.)
         do i=2,nx      ! Precalculate \int_0^|x| phi normalized (shape).
            xi=(i-1.)*dx
            phishape=phiofx(xi,1.,coshlen,toplen)
            phiint(i)=phiint(i-1)+(phishape+phiprev)*0.5*dx
            phiprev=phishape
!            write(*,*)i,phiint(i),phishape
         enddo
       endif
      fx=(nx-1.)*abs(x)/xmax
      ix=min(int(fx)+1,nx-1)
      fx=fx-ix+1
      intphiofx=sign(((1-fx)*phiint(ix)+fx*phiint(ix+1))*psi,x)
!      write(*,*)ix,x,fx,intphiofx,phiint(ix),phiint(ix+1)
! Integral is from 0 to x: and antisymmetric.
      end
c********************************************************************
      real function findxofphi(phiv,psi,coshlen,tl,xmin,xmax,nbi)
c Solve by bisection phiv=phiofx(xc,psi,coshlen,tl) and return xc.
c The intial x-range is [xmin,xmax]. Up to nbi bisections are allowed.
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
      real function findxofcum(ranf,nphi,den,pden,cump,xpden,xofphi,isw
     $     ,im,xf,xmin,xmax)
c Replacement for findxofran that uses a specified den(phi).
c and NOT the potential profile implied by psi, coshlen, tl.
      real ranf,xf,xmin,xmax
      integer nphi,isw,im
      real den(0:nphi),xofphi(0:nphi)
      real pden(-nphi-1:nphi+1),cump(-nphi-1:nphi+1),xpden(-nphi-1:nphi
     $     +1)

c Local variables
c      integer np
c      parameter (np=100)
c      real p(np),x(np),cum(np),r
      
      if(isw.eq.0)then   ! Set up arrays.
         do i=-nphi,nphi
            pden(i)=den(nphi-abs(i))
            xpden(i)=sign(xofphi(nphi-abs(i)),float(i))
!            write(*,*)'i,pden(i),xpden(i)',i,pden(i),xpden(i)
         enddo
         pden(-nphi-1)=1.
         xpden(-nphi-1)=xmin
         pden(nphi+1)=1.
         xpden(nphi+1)=xmax
      endif
      findxofcum=xofcum(2*nphi+3,pden,xpden,cump,ranf,isw,im,xf)
      end
c********************************************************************
      real function findxofran(ranf,psi,coshlen,tl,xmin,xmax,nbi,kp2)
c Find by bisection the x-position given the range fraction ranf s.t.
c   ranf=(xc-xmin+derivphiofx(xc...))/(xmax-xmin) 
c and return xc. This expression is the cumulative probability because
c P\propto \int n_e dx = \int d/dx(\phi')+1 dx=\phi'+x+ const.
c The "block" boundaries are non-uninformly spaced when a hole is present
c This routine knows nothing about blocks per se, but scales density via
c non-uniform mapping of ranf to x-position.
c Up to nbi bisections are allowed.
c Shortcuts are taken if potential is less than a range of relevance.
c Version 2 with additional parameter kp2 for 3-D holes.
      real ranf,psi,coshlen,tl,xmin,xmax,kp2
      integer nbi
      real phimin
      parameter (phimin=.01)
      real derivphiofx,intphiofx ! These two functions specify phiofx
      external derivphiofx,intphiofx
      if(abs(psi).lt.phimin.or.ranf.lt.phimin.or.1-ranf.lt.phimin)then
 ! Effectively linear shortcut.
         findxofran=xmin+ranf*(xmax-xmin)
         return
      endif
      xa=xmin
      xb=xmax
      xc=0.
      if(kp2.eq.0)then
         fva=(xa-xmin+derivphiofx(xa,psi,coshlen,tl))/(xmax-xmin)-ranf
         fvb=(xb-xmin+derivphiofx(xb,psi,coshlen,tl))/(xmax-xmin)-ranf
      else
         fva=(xa-xmin+derivphiofx(xa,psi,coshlen,tl)
     $        -kp2*(intphiofx(xa,psi,coshlen,tl)
     $           -intphiofx(xmin,psi,coshlen,tl)))
     $        /(xmax-xmin-kp2*(intphiofx(xmax,psi,coshlen,tl)
     $           -intphiofx(xmin,psi,coshlen,tl)))
     $        -ranf
         fvb=(xb-xmin+derivphiofx(xb,psi,coshlen,tl)
     $        -kp2*(intphiofx(xb,psi,coshlen,tl)
     $           -intphiofx(xmin,psi,coshlen,tl)))
     $        /(xmax-xmin-kp2*(intphiofx(xmax,psi,coshlen,tl)
     $           -intphiofx(xmin,psi,coshlen,tl)))
     $        -ranf
      endif
c Allow a value beyond the xmax range so long as ranf is non-negative.
      if(fva*fvb.gt.0.)then
         write(*,*)xa,fva,xb,fvb,phiv,phiofx(xa,psi,coshlen,tl)
     $        ,phiofx(xb,psi,coshlen,tl)
         stop 'Potential outside range'
      endif
      do j=1,nbi
         xc=0.5*(xa+xb)
         if(kp2.eq.0)then
           fvc=(xc-xmin+derivphiofx(xc,psi,coshlen,tl))/(xmax-xmin)-ranf
         else
            fvc=(xc-xmin+derivphiofx(xc,psi,coshlen,tl)
     $           -kp2*(intphiofx(xc,psi,coshlen,tl)
     $           -intphiofx(xmin,psi,coshlen,tl)))
     $           /(xmax-xmin-kp2*(intphiofx(xmax,psi,coshlen,tl)
     $           -intphiofx(xmin,psi,coshlen,tl)))
     $           -ranf
         endif
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
      integer i,id,iregion,js
      integer ixp(ndims)
      integer savedipartperiod

!      write(*,*)'Applying initial wave',wavespec
      if(wavespec(1).eq.0)return
      do id=1,ndims    ! Ensure precise periodicity
         kw(id)=nint(wavespec(1+id))*2.*3.1415926
     $        /(xmeshend(id)-xmeshstart(id))
!         write(*,*)wavespec(1+id),xmeshend(id),xmeshstart(id),kw(id)
      enddo
! Cope with vreset particles.
      savedipartperiod=ipartperiod(1)
      ipartperiod(1)=4
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
! Initialize the mesh fraction data in x_part.
! This also relocates periodically into the mesh because ipartperiod=4.
               call partlocate(x_part(1,i),ixp,xfrac,iregion,linmesh)
               if(.not.linmesh)then
                  write(*,*)'wavedisplace linmesh ERROR',ixp,xfrac
     $                 ,iregion,linmesh
               endif
            endif
         enddo
      enddo
      ipartperiod(1)=savedipartperiod
      end
c*********************************************************************
      subroutine testintphiofx
      integer nx
      parameter (nx=100)
      real phiint(nx),xp(nx),intphiofx,kp2,dint(nx)
      external intphiofx
      coshlen=4.
      toplen=1.
      kp2=0.
      psi=1.
      xmin=-5.*coshlen
      xmax=5.*coshlen
      dx=(xmax-xmin)/(nx-1.)
      dint(1)=0.
      prev=0.
      do i=1,nx
         xp(i)=xmin+(i-1)*dx
         phiint(i)=intphiofx(xp(i),psi,coshlen,toplen)
! Direct integration for verification:
         if(i.gt.1)then
            phi=phiofx(xp(i),psi,coshlen,toplen)
            dint(i)=dint(i-1)+(prev+phi)*0.5*dx
            prev=phi
         endif
      enddo
      offset=0.5*(dint(1)+dint(nx))
      dint=dint - offset
      call autoplot(xp,phiint,nx)
      call axlabels('x','!AJf!@dx')
      call color(4)
      call dashset(1)
      call polyline(xp,dint,nx)
      call dashset(0)
      call color(15)
      call pltend
      end
c*********************************************************************
c Optional main test.
c      call testintphiofx
c      end
c*********************************************************************
