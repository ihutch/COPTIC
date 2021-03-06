! pinit is the only routine in this file that is called from other
! files.  Hence it can safely be replaced wholesale by another file.
! However, colreinit in cartreinject.f presumes collisional distribution
! data has already been calculated by routines here. 
!***********************************************************************
! Initializing particles.
      subroutine pinit(sprior)
      implicit none
      real sprior
! Common data:
      include 'ndimsdecl.f'
      include 'partcom.f'
      include 'plascom.f'
      include 'myidcom.f'
      include 'meshcom.f'
      include 'colncom.f'
      include 'cdistcom.f'
      external linregion,ranlenposition,gasdev
      logical linregion
! Local dummy variables for partlocate.
      real xfrac(ndims)
      integer ixp(ndims)
      logical linmesh
! Factor by which we leave additional space for increased number of 
! particles if not fixed:
      real slotsurplus
      real thetamax
      parameter (thetamax=1.)
! Local variables to satisfy implicit none
      integer i,i1,islotmax,ispecies,j,k,k1,k2,ko,nclim,nonzero,ntries
      integer iregion
      real theta,tisq,vzave,ranlenposition,gasdev
!-----------------------------------------------------------------
      npassthrough=0
! A special orbit. Preserved only for ispecies=1.
! Pinit resets x_part. So set it for the special first particle.
      x_part(1,1)=2.
      x_part(2,1)=0.
      x_part(3,1)=0.
! Prior half-step radial velocity
      x_part(4,1)=0.5*sprior*(abs(phip)/x_part(1,1)**2)
!      x_part(4,1)=0.5*dt*(abs(phip)/x_part(1,1)**2)
!      x_part(4,1)=0.25*dt*(abs(phip)/x_part(1,1)**2)
! Tangential velocity of circular orbit at r=?.
      x_part(5,1)=sqrt(abs(phip/x_part(1,1))-x_part(4,1)**2)
      x_part(6,1)=0.
      x_part(iflag,1)=1
      i1=2
      call partlocate(x_part(1,1),ixp,xfrac,iregion,linmesh)
      if(.not.linmesh.or..not.linregion(ibool_part,ndims,x_part(1
     $     ,1)))then
         if(myid.eq.0.)then
            write(*,*)'WARNING Special particle-1 outside region.'
     $           ,' Resetting from coordinates'
            write(*,'(9f7.3)')(x_part(i,1),i=1,9)
         endif
         i1=1
      endif
!-----------------------------------------------------------------
! Point to the bottom of the particle stack for start of species 1.
      iicparta(1)=1
! Decide whether the external distribution is separable.
      do ispecies=1,nspecies
! Conveniently here initialize distribution numbers.
         ncdists(ispecies)=0
         if(ninjcompa(ispecies).gt.0)then
! This is safer at 1.3 than 1.1
            slotsurplus=1.1
         else
            slotsurplus=1.
         endif
         ntries=0
         Eneutral=0.
! Set tperp to a tiny number if species is infinitely magnetized.
         theta=Bt*eoverms(ispecies)*dt
         if(abs(theta).gt.thetamax)Tperps(ispecies)=1.e-24
         notseparable(ispecies)=0
         if(colntime.ne.0..and.ispecies.eq.1)then
! At this point vperp refers to the perp part of flow, set by cmdline.
! Initialize the reinjection particles and discover the full Eneutral.
            call colninit(0,myid)
            call colreinit(myid,ispecies)
! Finalize the Eneutral fraction and related settings.
            Eneutral=Enfrac*Eneutral
            if(Eneutral.ne.0.and.colpow.ne.0)notseparable(ispecies)=1
! Add on the orthogonal EnxB drift, so as to make vperp the velocity of
! the frame of reference in which the background E-field is truly zero:
            do k=1,ndims
               k1=mod(k,ndims)+1
               k2=mod(k+1,ndims)+1
               vperp(k)=vperp(k)+(Eneutral/Bt)
     $              *(vdrift(k1)*Bfield(k2)-vdrift(k2)*Bfield(k1))
            enddo
         elseif(Tperps(ispecies).ne.Ts(ispecies)
     $           )then
! Collisionless anisotropic distribution is notseparable.
! Infinite B case picked up a bug but it is fixed so it can use this too.
!     $        .and.Tperps(ispecies).gt.1.e-24)then
            notseparable(ispecies)=2
            if(myid.eq.0)write(*,'(a,2f8.4,a,i3)')
     $           ' Non-separable anisotropic T.'
     $           ,Ts(ispecies),Tperps(ispecies),' species',ispecies
         else
! Count the number of non-zero vdrift components. If it is more than one
! then vdrift is not along a coordinate axis => nonseparable.
            nonzero=0
            do k=1,ndims
               if(vdrifts(k,ispecies).ne.0)nonzero=nonzero+1
            enddo
            if(nonzero.gt.1)then
               notseparable(ispecies)=3
               if(myid.eq.0)write(*,*)'Non-separable oblique vdrift'
!     $           ,(vdrifts(k,ispecies),k=1,ndims),' species',ispecies
            endif
! Test
!               notseparable(ispecies)=3
         endif

         if(notseparable(ispecies).ge.2)then
!----------------------------------------
! This is where the collisionless distribution is really set.
! More complicated distributions should redo this setting differently.
! How many prior particles to save:
!            ncdists(ispecies)=0
! How many extra to add:
            nclim=ncdistmax
!            write(*,*)'Calling maxwellstats',nclim,ispecies
            call maxwellstats(nclim,vzave,ispecies)
            call colreinit(myid,ispecies)
         endif
         if(myid.eq.0.and.ldistshow.and.notseparable(ispecies).ne.0)
     $        call colndistshow(ispecies)
!      write(*,*)i1,' iicparta(1)',iicparta(1),iicparta(ispecies)
! Scale the number of particles for higher slots to the ions(?)
         nparta(ispecies)=nparta(1)/numratioa(ispecies)
!     $        *sqrt(abs(eoverms(1)/eoverms(ispecies)))
         islotmax=nparta(ispecies)+iicparta(ispecies)-1
         if((islotmax+1)*slotsurplus.gt.n_partmax)then
            write(*,*)'Too few particle slots',n_partmax
     $           ,' max, to accommodate',islotmax
            write(*,*)'multiplied by slotsurplus factor'
     $           ,slotsurplus,int(islotmax*slotsurplus)
            stop
         endif
         tisq=sqrt(Ts(ispecies)*abs(eoverms(ispecies)))
! ----------------------------- Actual Particle Setting ----------
         do i=iicparta(ispecies)+i1-1,islotmax
            x_part(iflag,i)=1
 1          continue
            ntries=ntries+1
! New position choice including density gradients.
            do j=1,ndims
               x_part(j,i)=ranlenposition(j)
            enddo
            if(.not.linregion(ibool_part,ndims,x_part(1,i)))then
!            write(*,*)'Initialization of',i,' wrong region',inewregion
!     $           ,(x_part(kk,i),kk=1,3)
!     If we are not in the plasma region, try again if we are doing fixed
!     particle number else just set this slot empty.
               if(ninjcompa(ispecies).ne.0)then
                  x_part(iflag,i)=0
               else
                  goto 1
               endif
            endif
            if(notseparable(ispecies).eq.0. .and. Eneutral.eq.0)then
! Shifted Gaussians.
               x_part(4,i)=tisq*gasdev(myid) + vds(ispecies)*vdrift(1)
               x_part(5,i)=tisq*gasdev(myid) + vds(ispecies)*vdrift(2)
               x_part(6,i)=tisq*gasdev(myid) + vds(ispecies)*vdrift(3)
            else
               call colvget(x_part(4,i),ispecies)
            endif

! One idea that might improve the initialization would be to add the 
! (electric) potential energy here, by scaling up the total velocity,
! but keeping the direction the same. 

! However, to make this accurate, one would have to have solved the
! initial potential better than the simple Laplace solve currently
! used. A solution of the electron-shielded case might be nearer the
! mark, but it would not be correct when T_i<T_e. It is perhaps
! advantageous because it constitutes a minimal estimate of the
! shielding, which is probably what one wants.

! One advantage of using the Laplace or electron-shielded case is that
! it enhances the kinetic energy. Thus if the potential rises towards
! zero during initial evolution, one has indeed avoided spuriously
! trapped ions. 

! Asymmetries are of course not accounted for. The particles in the
! close wake are in fact going to have mostly negative velocities. They
! are therefore completely unmodelled in this initialization, which
! fills in positive drifting ions even where they could not be present.

! The previous timestep length.
            x_part(idtp,i)=0.
! Initialize the mesh fraction data in x_part.
            call partlocate(x_part(1,i),ixp,xfrac,iregion,linmesh)
! This test rejects particles exactly on mesh boundary:
            if(.not.linmesh)goto 1
         enddo
!------------------------------- End of Actual Particle Setting --
! The maximum used slot for this species
         iocparta(ispecies)=i-1
! Start of next slot-set may give a gap for overflow.
         iicparta(ispecies+1)=int(i*slotsurplus)
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
      enddo
! Set flag of unused slots to 0
      do i=iicparta(ispecies),n_partmax
         x_part(iflag,i)=0
      enddo
! Allow the last species to fill the array:
      iicparta(ispecies)=n_partmax+1
      end
!***********************************************************************
! Return a random velocity from the (precalculated) distribution heap.
      subroutine colvget(v,ispecies)
      implicit none
      include 'ndimsdecl.f'
      include 'cdistcom.f'
      real v(ndims)
      integer ispecies
      real ra
      integer nc,i

      call ranlux(ra,1)
      nc=int(ncdists(ispecies)*ra)
      do i=1,ndims
         v(i)=vcols(i,nc,ispecies)
      enddo
      end
!*********************************************************************
!*********************************************************************
! Wrapper of Version to do fullscale neutral statistics.
      subroutine colninit(ncin,myid)
      implicit none
      integer ncin,myid
      real relax,vzave
      include 'ndimsdecl.f'
      include 'plascom.f'
      include 'colncom.f'
      include 'cdistcom.f'
      integer nclim,k

      if(Enfrac.eq.0.)then
! Initialize just from the neutral distribution.
         Eneutral=0.
         if(myid.eq.0)write(*,*)'Enfrac=0. Not initializing cdist.'
      else
! The number of samples.
         if(ncin.eq.0)then
            nclim=ncdistmax
         else
            nclim=ncin
         endif
!         write(*,*)'nclim=',nclim
! Initial guess at what whole Eneutral really needs to be to give the
! drift requested. Also set relaxation rate for iterations.
! Eneutral is based only on species-1 drifts.
         Eneutral=(vds(1)-vneutral)/colntime*(1-0.5*colpow)
         if(myid.eq.0)write(*,'(a,$)')' Colninit velocity: '
         relax=1.
         if(colpow.gt.0)relax=0.5
         do k=1,3
! Iterate over guesses at Eneutral. 
            call colnstats(nclim,vzave) 
            vzave=vzave/nclim
            if(myid.eq.0)write(*,'(f6.4,'', '',$)')vzave
            if(colpow.eq.0.)goto 4
            Eneutral=Eneutral*(relax*((vds(1)-vneutral)/vzave-1.)+1.)
!            write(*,*)Eneutral
         enddo
! End of Eneutral iteration.
 4       continue
         if(myid.eq.0)write(*,'(i7,a,f8.4)')nclim,' Colninit Eneutral='
     $        ,Eneutral
      endif
      end
!*********************************************************************
      subroutine colnstats(nclim,vzave)
! Routine for initializing particles when collisions and drifts
! driven by Eneutral are present. An array of velocities that 
! represents the background distribution of ions is calculated by
! advancing from random birth velocity representative of the neutral
! distribution. The collision frequency is proportional to velocity
! to the power colpow. We use the "Null Collision Method" of Skullerud
! to weight the distribution of the neutrals when the collision frequency
! is velocity dependent. 

! The velocities are deposited in cdistcom variable vcol(i,j) for j=1,nclim.
! In addition, the sum of all the absolute velocity cpts in cdistflux(i).
! It decides weight of 3 coordinate directions for reinjection.
! vzave returns the sum of the vdrift-components of the velocities.

      implicit none
      integer nclim
      real vzave
      include 'ndimsdecl.f'
      include 'cdistcom.f'
      include 'plascom.f'
      include 'colncom.f'
      real ctprime
      real gasdev
      external gasdev

! Local variables
      integer i
      real tisq,dt0,torb,ttic,v2
      real v(ndims),vn(ndims)
      real theta,stheta,ctheta,dvpar,delta
      real colnu,rd
      parameter (delta=0.2)

      tisq=sqrt(Ti)

! Decide the duration of the time tic dt0, and the duration of the
! possible-collision time ctprime.  Take one or other to be a smallish
! fraction of the collision time if the collision time is variable, so
! that sufficient resolution is present to do a reasonable integral of a
! variable collision frequency.  Taking dt0 small relative to colntime
! makes the perpendicular velocity and orbit length representation
! statistically worse. But having substantially it bigger than ctprime
! increases the computational effort of the initialization.  The actual
! collision time is colntime*(v^2/Ti)^(colpow/2.)  So if colpow=0. the
! collision time is colntime independent of v.  
      if(colpow.eq.0.)then
!         if(.false.)then
         dt0=1.*colntime
         ctprime=.99*colntime
      else
         dt0=.1*colntime
         ctprime=.08*colntime
      endif
!         write(*,*)'colntime,dt0,ctprime',colntime,dt0,ctprime

! One might imagine a need to ensure dt does not happen to coincide with
! a low order rational multiple of the cyclotron period. But there is no 
! sign in output of this need. Probably the random collision length does
! enough to randomize even a bad coincidence.
! Iterate over adjustments to Eneutral.
      do i=1,ndims
         cdistflux(i)=0.
      enddo
      ncdist=0
      vzave=0.
      ttic=dt0
!--------------------------------------------
! Start of a new orbit
! 2    continue 
! Inject from neutral distribution
      v2=delta*Ti
      do i=1,ndims
         v(i)=tisq*gasdev(0)
         v(i)=v(i)+vneutral*vdrift(i)
         v2=v2+v(i)**2
      enddo
! torb is the orbit time in velocity-scaled units. ttic in unscaled.
      call ranlux(torb,1)
      torb=-alog(torb+1.e-15)*ctprime

!      write(*,*)'ctprime=',ctprime
!      write(*,*)ttic,torb,Eneutral,v,Ti,Bt
 3    if(ttic.le.torb)then
! If ttic<torb, advance to tic.
! Subtract the timestep from time to the orbit end in scaled units
         torb=torb-ttic
         v2=delta*Ti
         ncdist=ncdist+1
         if(Bt.eq.0.)then
! B-field free case
            do i=1,ndims
!               if(i.eq.ndims)v(i)=v(i)+Eneutral*ttic
               v(i)=v(i)+Eneutral*vdrift(i)*ttic
               v2=v2+v(i)**2
               vcol(i,ncdist)=v(i)
               cdistflux(i)=cdistflux(i)+abs(v(i))
            enddo
!            write(*,*)'ncdist,v(3)',ncdist,v(3),ttic,torb
         else
! Finite B-field case
            theta=Bt*ttic
! Rotation is counterclockwise for ions.
            stheta=sin(-theta)
            ctheta=cos(theta)
! Subtract off the perpendicular velocity of frame in which E is zero.
            do i=1,ndims
               v(i)=v(i)-vperp(i)
            enddo
! Rotate the velocity
            call rotate3(v,stheta,ctheta,Bfield)
! Acceleration along the B-direction:
            dvpar=0.
            do i=1,ndims
               dvpar=dvpar+Eneutral*vdrift(i)*Bfield(i)*ttic
            enddo
! And add back the drift velocity, and the parallel acceleration.
            do i=1,ndims
               v(i)=v(i)+vperp(i)+dvpar*Bfield(i)
               v2=v2+v(i)**2
               vcol(i,ncdist)=v(i)
               cdistflux(i)=cdistflux(i)+abs(v(i))
            enddo
!               write(*,*)'ncdist=',ncdist
!               write(*,*)'theta,v',theta,v
         endif
!         vzave=vzave+vcol(ndims,ncdist)
         do i=1,ndims
            vzave=vzave+vcol(i,ncdist)*vdrift(i)
         enddo
! Now we are at the end of a tic. Set the next tic length
         ttic=dt0
! If we haven't exhausted the sample number, start next tic
         if(ncdist.lt.nclim)goto 3
      else
! If ttic>torb, advance to next orbit without storing a sample
! Here, rather than just restarting the orbit with the neutral
! distribution we need to do the null collision calculation and goto 3
! not 2. 
! Advance to a possible collision.
         if(Bt.eq.0.)then
!            v(ndims)=v(ndims)+Eneutral*torb 
            do i=1,ndims
               v(i)=v(i)+Eneutral*vdrift(i)*torb
            enddo
         else
! Finite B-field case
            theta=Bt*torb
! Rotation is counterclockwise for ions.
            stheta=sin(-theta)
            ctheta=cos(theta)
! Subtract off the perpendicular velocity of frame in which E is zero.
            do i=1,ndims
               v(i)=v(i)-vperp(i)
            enddo
! Rotate the velocity
            call rotate3(v,stheta,ctheta,Bfield)
! Acceleration along the B-direction:
            dvpar=0.
            do i=1,ndims
               dvpar=dvpar+Eneutral*vdrift(i)*Bfield(i)*torb
            enddo
! And add back the drift velocity, and the parallel acceleration.
            do i=1,ndims
               v(i)=v(i)+vperp(i)+dvpar*Bfield(i)
            enddo
         endif
! Possible collision: Choose a neutral velocity and calculate (v-vn)^2
         v2=0.
         do i=1,ndims
            vn(i)=tisq*gasdev(0)
            vn(i)=vn(i)+vneutral*vdrift(i)
            v2=v2+(v(i)-vn(i))**2
         enddo
! Calculate the ratio of the actual collision time at this relative
! velocity to the possible-collision time:
         colnu=((v2+delta*Ti)/(Ti+delta*v2))**(colpow/2.)
     $        *colntime/ctprime
! Select a fraction of these collisions weighted by this ratio.
         call ranlux(rd,1)
         rd=rd*colnu
         if(rd.lt.1.)then
! This collision occurred. 
            do i=1,ndims
               v(i)=vn(i)
            enddo
         else
! No collision. Just carry on.
         endif
         ttic=ttic-torb
! Choose next orbit time.
         call ranlux(torb,1)
         torb=-alog(torb+1.e-15)*ctprime
         if(.not. torb.lt.100.*ctprime)then
            write(*,*)'torb problem',torb,ctprime
         endif
         goto 3
      endif
!      write(*,*)'v(i)',v
      end
!****************************************************************
      subroutine colndistshow(ispecies)
      include 'ndimsdecl.f'
      include 'cdistcom.f'
      integer nxbin,nybin
      parameter (nxbin=50,nybin=50,idim1=1,idim2=3)
      real fvxvy(nxbin,nybin),vxa(nxbin,nybin),vya(nxbin,nybin)
      real work(0:nxbin+1,0:nybin+1)
      real vlimit(2,ndims)

! Defaults
      do id=1,ndims
         vlimit(1,id)=1.
         vlimit(2,id)=-1.
      enddo
! Adjust vlimits
!      do j=1,ncdist-1 I did not understand the -1.
      do j=1,ncdist
         vx=vcols(idim1,j,ispecies)
         if(vx.gt.vlimit(1,1))vlimit(1,1)=vx
         if(vx.lt.vlimit(2,1))vlimit(2,1)=vx
         vy=vcols(idim2,j,ispecies)
         if(vy.gt.vlimit(1,2))vlimit(1,2)=vy
         if(vy.lt.vlimit(2,2))vlimit(2,2)=vy
      enddo
! Velocity arrays
      do i=1,nxbin
      do j=1,nybin
         vxa(i,j)=vlimit(1,1)+(vlimit(2,1)-vlimit(1,1))*(i-1.)/(nxbin
     $        -1.)
         vya(i,j)=vlimit(1,2)+(vlimit(2,2)-vlimit(1,2))*(j-1.)/(nybin
     $        -1.)
      enddo
      enddo

! Accumulate the fvxvy distribution
         do k=1,nybin
            do j=1,nxbin
               fvxvy(j,k)=0.
            enddo
         enddo
!         do j=1,ncdist-1  I don't understand the -1
         do j=1,ncdist
            vx=vcols(idim1,j,ispecies)
            vy=vcols(idim2,j,ispecies)
            nx=nint(nxbin*(vx-vlimit(1,1))/(vlimit(2,1)-vlimit(1,1)))
            ny=nint(nybin*(vy-vlimit(1,2))/(vlimit(2,2)-vlimit(1,2)))
            if(nx.ge.1.and.nx.le.nxbin.and.ny.ge.1.and.ny.le.nybin)
     $           then
               fvxvy(nx,ny)=fvxvy(nx,ny)+1.
            endif
         enddo
!      write(*,*)vxa
98     call pltinit(0.,1.,0.,1.)
!       Plot the surface. With axes (2-D). Web color 10, axis color 7.
        j=2 + 256*10 + 256*256*7
        call surf3d(vxa,vya,fvxvy,nxbin,nxbin,nybin,j,work)
        call color(15)
        call ax3labels('x','z','f ')
        if(ieye3d().ne.0) goto 98

        open(2,file='fvxvydist.dat')
        write(2,*)(fvxvy(nxbin/2,j),j=1,nybin)
        close(2)
      end

!****************************************************************
      subroutine maxwellstats(nclim,vzave,ispecies)
! Routine for initializing particles with a general anisotropic and
! shifted maxwellian distribution. An array of velocities that
! represents the background distribution of ions is calculated The
! velocities are deposited in cdistcom variable vcols(i,j,ispecies) for
! j=1,nclim.  In addition, the sum of all the absolute velocity cpts in
! cdistflux(i).  It decides weight of 3 coordinate directions for
! reinjection.  vzave returns the sum of the vdrift-components of the
! velocities.

! If Bt is non-zero then probably the distribution is gyrotropic, and
! thus the direction of any temperature asymmetry should be
! perpendicular or parallel to B. B is the preferred axis.  If Bt is
! zero, then only the velocity direction cosines vdrift determines a
! preferred axis. 

! Ts will be taken as representing the temperature(s) in the direction of 
! the preferred axis. Tperps are then the temperatures in the perpendicular
! direction (isotropic for now).

! The triplet of unit vectors is 1-parallel (to B or v); 2-perpendicular to
! 1 and to x, unless 1 happens to be in the x-direction, in which case 2
! is perpendicular to 1 and y; 3-perpendicular to 1 and 2.

      implicit none
      integer nclim,ispecies
      real vzave

      include 'ndimsdecl.f'
      include 'cdistcom.f'
      include 'plascom.f'
      real vdirs(ndims,ndims+1,nspeciesmax),rmag
      real vns(nspeciesmax),vnp(nspeciesmax)
      integer i,j
      real gasdev
      external gasdev
      save vns,vnp,vdirs

!      write(*,*)'Bt,Bfield',Bt,Bfield

      if(Bt.eq.0.)vpars(ispecies)=vds(ispecies)
! Define the ndims directions 1 parallel and (ndims-1) perpendicular.
      do i=1,ndims
         cdistfluxs(i,ispecies)=0.
         if(Bt.ne.0)then
! Parallel is B-direction 
            vdirs(i,1,ispecies)=Bfield(i)
         else
! vdrifts defines parallel direction
            vdirs(i,1,ispecies)=vdrifts(i,ispecies)
         endif
         vdirs(i,ndims,ispecies)=0.
         vdirs(i,ndims+1,ispecies)=0.
      enddo
! Alternative vectors to generate perpendicular.
      vdirs(1,ndims,ispecies)=1.
      vdirs(2,ndims+1,ispecies)=1.
! Construct perpendicular vectors:
      call crossprod(vdirs(1,1,ispecies),vdirs(1,ndims
     $     ,ispecies),vdirs(1,2,ispecies),rmag)
! If we were unlucky and chose a direction that coincides, redo.
      if(rmag.le.0)
     $     call crossprod(vdirs(1,1,ispecies),vdirs(1,ndims+1
     $     ,ispecies),vdirs(1,2,ispecies),rmag)
      do i=1,ndims
         vdirs(i,2,ispecies)=vdirs(i,2,ispecies)/rmag
      enddo
      call crossprod(vdirs(1,1,ispecies),vdirs(1,2
     $     ,ispecies),vdirs(1,3,ispecies),rmag)
! Should already be normalized because 1st two orthonormal.
! Now we have 3 orthogonal unit vectors in vdirs(1:3,1:3,ispecies)
! Save the thermal velocity
      vns(ispecies)=sqrt(Ts(ispecies)*abs(eoverms(ispecies)))
      vnp(ispecies)=sqrt(Tperps(ispecies)*abs(eoverms(ispecies)))
!      write(*,'(3f10.4)'),((vdirs(i,j,ispecies),i=1,3),j=1,3)
!      write(*,*)'ispecies,vns,vnp=',ispecies,vns,vnp
      if(ncdists(ispecies)+nclim.gt.ncdistmax)then
         write(*,*)'Trying to initialize too many random particles',
     $        ncdists(ispecies),'+',nclim,' in pinit.'
      else
         do j=ncdists(ispecies)+1,ncdists(ispecies)+nclim
! Draw random parallel and perpendicular velocities, and add on the
! drift velocity. 
            do i=1,ndims
               vcols(i,j,ispecies)
     $              =vns(ispecies)*gasdev(0)*vdirs(i,1,ispecies)
     $              +vnp(ispecies)*gasdev(0)*vdirs(i,2,ispecies)
     $              +vnp(ispecies)*gasdev(0)*vdirs(i,3,ispecies)
     $                     +vdrifts(i,ispecies)*vds(ispecies)
               vzave=vzave+vcols(i,j,ispecies)*vdrift(i)
! Evaluate the total absolute fluxes in the three coordinate directions.
! This required cdistfluxs to be real*8 to avoid rounding errors:
               cdistfluxs(i,ispecies)=cdistfluxs(i,ispecies)
     $              +abs(vcols(i,j,ispecies))
            enddo
!            if(ispecies.eq.2.and.abs(vcols(2,j,ispecies)).gt.0.001
!     $           .and.j.le.10)then
!               write(*,*)'maxwell vynonzero',j,
!     $              (vcols(i,j,ispecies),i=1,3)
!            endif
         enddo
      endif
      ncdists(ispecies)=ncdists(ispecies)+nclim

!      write(*,*)'Completed Maxwellstats',nclim,ncdists(ispecies),vzave
      end


