c***********************************************************************
c Initializing particles.
      subroutine pinit(sprior)
c Common data:
      include 'ndimsdecl.f'
      include 'partcom.f'
      include 'plascom.f'
      include 'myidcom.f'
      include '3dcom.f'
      include 'meshcom.f'
      include 'colncom.f'
      external linregion,ranlenposition
      logical linregion
c Local dummy variables for partlocate.
      real xfrac(ndims)
      integer ixp(ndims)
      logical linmesh
c Factor by which we leave additional space for increased number of 
c particles if not fixed:
      real slotsurplus
c-----------------------------------------------------------------
c A special orbit. Preserved only for ispecies=1.
c Pinit resets x_part. So set it for the special first particle.
      x_part(1,1)=2.
      x_part(2,1)=0.
      x_part(3,1)=0.
c Prior half-step radial velocity
      x_part(4,1)=0.5*sprior*(abs(phip)/x_part(1,1)**2)
c      x_part(4,1)=0.5*dt*(abs(phip)/x_part(1,1)**2)
c      x_part(4,1)=0.25*dt*(abs(phip)/x_part(1,1)**2)
c Tangential velocity of circular orbit at r=?.
      x_part(5,1)=sqrt(abs(phip/x_part(1,1))-x_part(4,1)**2)
      x_part(6,1)=0.
      x_part(iflag,1)=1
      i1=2
      call partlocate(x_part(1,1),ixp,xfrac,iregion,linmesh)
      if(.not.linmesh.or..not.linregion(ibool_part,ndims,x_part(1
     $     ,1)))then
         if(myid.eq.0.)then
            write(*,*)'WARNING Special particle-1 outside region.'
            write(*,*)'COORDINATES',(x_part(i,1),i=1,9),' Resetting.'
         endif
         i1=1
      endif
c-----------------------------------------------------------------
c Point to the bottom of the particle stack for start of species 1.
      iicparta(1)=1
      do ispecies=1,nspecies
         if(ninjcompa(ispecies).gt.1)then
            slotsurplus=1.3
         else
            slotsurplus=1.
         endif
         ntries=0
         Eneutral=0.
         if(colntime.ne.0..and.ispecies.eq.1)then
c At this point vperp refers to the perp part of vd, set by cmdline.
c Initialize the reinjection particles and discover the full Eneutral.
            call colninit(0,myid)
            if(myid.eq.0.and.ldistshow)call colndistshow()
            call colreinit(myid)
c Finalize the Eneutral fraction and related settings.
            Eneutral=Enfrac*Eneutral
c Add on the orthogonal EnxB drift, so as to make vperp the velocity of
c the frame of reference in which the background E-field is truly zero:
            do k=1,ndims
               k1=mod(k,ndims)+1
               k2=mod(k+1,ndims)+1
               vperp(k)=vperp(k)+(Eneutral/Bt)
     $              *(vdrift(k1)*Bfield(k2)-vdrift(k2)*Bfield(k1))
            enddo
         endif
         tisq=sqrt(Ts(ispecies)*abs(eoverms(ispecies)))
c      write(*,*)i1,' iicparta(1)',iicparta(1),iicparta(ispecies)
c Scale the number of particles for higher slots to the ions(?)
         nparta(ispecies)=nparta(1)
     $        *sqrt(abs(eoverms(1)/eoverms(ispecies)))
         do i=iicparta(ispecies)+i1-1,nparta(ispecies)
     $        +iicparta(ispecies)-1
            x_part(iflag,i)=1
 1          continue
            ntries=ntries+1
c New position choice including density gradients.
            do j=1,ndims
               x_part(j,i)=ranlenposition(j)
            enddo
c     If we are not in the plasma region, try again.
            if(.not.linregion(ibool_part,ndims,x_part(1,i)))then
c            write(*,*)'Initialization of',i,' wrong region',inewregion
c     $           ,(x_part(kk,i),kk=1,3)
               goto 1
            endif
            if(Eneutral.eq.0.)then
               x_part(4,i)=tisq*gasdev(myid) + vd*vdrift(1)
               x_part(5,i)=tisq*gasdev(myid) + vd*vdrift(2)
               x_part(6,i)=tisq*gasdev(myid) + vd*vdrift(3)
            else
               call colvget(x_part(4,i))
            endif

c One idea that might improve the initialization would be to add the 
c (electric) potential energy here, by scaling up the total velocity,
c but keeping the direction the same. 

c However, to make this accurate, one would have to have solved the
c initial potential better than the simple Laplace solve currently
c used. A solution of the electron-shielded case might be nearer the
c mark, but it would not be correct when T_i<T_e. It is perhaps
c advantageous because it constitutes a minimal estimate of the
c shielding, which is probably what one wants.

c One advantage of using the Laplace or electron-shielded case is that
c it enhances the kinetic energy. Thus if the potential rises towards
c zero during initial evolution, one has indeed avoided spuriously
c trapped ions. 

c Asymmetries are of course not accounted for. The particles in the
c close wake are in fact going to have mostly negative velocities. They
c are therefore completely unmodelled in this initialization, which
c fills in positive drifting ions even where they could not be present.

c The previous timestep length.
            x_part(idtp,i)=0.
c Initialize the mesh fraction data in x_part.
            call partlocate(x_part(1,i),ixp,xfrac,iregion,linmesh)
c This test rejects particles exactly on mesh boundary:
            if(.not.linmesh)goto 1
         enddo
c----------------------------------------
c The maximum used slot for this species
         iocparta(ispecies)=i-1
c Start of next slot-set may give a gap for overflow.
         iicparta(ispecies+1)=int(i*slotsurplus)
c Zero the overflow slots' flag
         do i=iocparta(ispecies)+1,iicparta(ispecies+1)-1
            x_part(iflag,i)=0
         enddo
         if(myid.eq.0)then
            if(nspecies.gt.0)write(*,*)'Slots',ispecies
     $           ,iicparta(ispecies),iocparta(ispecies)
     $           ,ninjcompa(ispecies)

            write(*,101)ispecies,nprocs,
     $           iocparta(ispecies)-iicparta(ispecies)+1,ntries
 101        format(' Initialized species',i2,i4,'x',i7
     $           ,' ntries=',i7,$)
         endif
c Initialize rhoinf:
         if(rhoinf.eq.0.)rhoinf=numprocs*n_part/(4.*pi*rs**3/3.)
c Initialize orbit tracking
         do ko=1,norbits
            iorbitlen(ko)=0
         enddo
         
c Don't shift for special particle the subsequent sections.
         i1=1
      enddo
c Set flag of unused slots to 0
      do i=iicparta(ispecies),n_partmax
         x_part(iflag,i)=0
      enddo
c Allow the last species to fill the array:
      iicparta(ispecies)=n_partmax+1
      end
c***********************************************************************
      subroutine locateinit()
c Initialize particles loaded from restart files so that they have the
c correct mesh-fraction coordinate in their upper values.
c If they are outside the mesh or particle region, drop them.
c Common data:
      include 'ndimsdecl.f'
      include 'partcom.f'
c      include 'myidcom.f'
      include '3dcom.f'
      include 'meshcom.f'
c Local dummy variables for partlocate.
      real xfrac(ndims)
      integer ixp(ndims)
      logical linmesh,linregion
      external linregion

      do i=1,ioc_part
         call partlocate(x_part(1,i),ixp,xfrac,iregion,linmesh)
c         write(*,*)linmesh,linregion(ibool_part,ndims,x_part(1,i))
         if(.not.linmesh .or. .not.
     $        linregion(ibool_part,ndims,x_part(1,i)))then
            x_part(iflag,i)=0
         endif
      enddo
      nrein=0
c      write(*,*)'Redetermined particle mesh locations (locateinit)'
      end
c***********************************************************************
c Return a random velocity from the (precalculated) distribution heap.
      subroutine colvget(v)
      implicit none
      include 'ndimsdecl.f'
      include 'cdistcom.f'
      real v(ndims)
      integer nc,i
      real ran1
      external ran1

      nc=ncdist*ran1(1)
      do i=1,ndims
         v(i)=v_col(i,nc)
      enddo
      end
c*********************************************************************
c*********************************************************************
c Wrapper of Version to do fullscale neutral statistics.
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
c Initialize just from the neutral distribution.
         Eneutral=0.
         if(myid.eq.0)write(*,*)'Enfrac=0. Not initializing cdist.'
      else
c The number of samples.
         if(ncin.eq.0)then
            nclim=ncdistmax
         else
            nclim=ncin
         endif
c         write(*,*)'nclim=',nclim
c Initial guess at what whole Eneutral really needs to be to give the
c drift requested. Also set relaxation rate for iterations.
         Eneutral=(vd-vneutral)/colntime*(1-0.5*colpow)
         if(myid.eq.0)write(*,'(a,$)')' Colninit velocity: '
         relax=1.
         if(colpow.gt.0)relax=0.5
         do k=1,3
c Iterate over guesses at Eneutral. 
            call colnstats(nclim,vzave) 
            vzave=vzave/nclim
            if(myid.eq.0)write(*,'(f6.4,'', '',$)')vzave
            if(colpow.eq.0.)goto 4
            Eneutral=Eneutral*(relax*((vd-vneutral)/vzave-1.)+1.)
c            write(*,*)Eneutral
         enddo
c End of Eneutral iteration.
 4       continue
         if(myid.eq.0)write(*,'(i7,a,f8.4)')nclim,' Colninit Eneutral='
     $        ,Eneutral
      endif
      end
c*********************************************************************
      subroutine colnstats(nclim,vzave)
c Routine for initializing particles when collisions and drifts
c driven by Eneutral are present. An array of velocities that 
c represents the background distribution of ions is calculated by
c advancing from random birth velocity representative of the neutral
c distribution. The collision frequency is proportional to velocity
c to the power colpow. We use the "Null Collision Method" of Skullerud
c to weight the distribution of the neutrals when the collision frequency
c is velocity dependent. 

c The velocities are deposited in cdistcom variable v_col(i,j) for j=1,nclim.
c In addition, the sum of all the velocities in cdistflux(i).
c vzave returns the sum of the vdrift-components of the velocities.

      implicit none
      integer nclim
      real vzave
      include 'ndimsdecl.f'
      include 'cdistcom.f'
      include 'plascom.f'
      include 'colncom.f'
      real gasdev,ran1,ctprime
      external gasdev, ran1

c Local variables
      integer i
      real tisq,dt0,torb,ttic,v2
      real v(ndims),vn(ndims)
      real theta,stheta,ctheta,dvpar,delta
      real colnu,rd
      parameter (delta=0.2)

      tisq=sqrt(Ti)

c Decide the duration of the time tic dt0, and the duration of the
c possible-collision time ctprime.  Take one or other to be a smallish
c fraction of the collision time if the collision time is variable, so
c that sufficient resolution is present to do a reasonable integral of a
c variable collision frequency.  Taking dt0 small relative to colntime
c makes the perpendicular velocity and orbit length representation
c statistically worse. But having substantially it bigger than ctprime
c increases the computational effort of the initialization.  The actual
c collision time is colntime*(v^2/Ti)^(colpow/2.)  So if colpow=0. the
c collision time is colntime independent of v.  
      if(colpow.eq.0.)then
c         if(.false.)then
         dt0=1.*colntime
         ctprime=.99*colntime
      else
         dt0=.1*colntime
         ctprime=.08*colntime
      endif
c         write(*,*)'colntime,dt0,ctprime',colntime,dt0,ctprime

c One might imagine a need to ensure dt does not happen to coincide with
c a low order rational multiple of the cyclotron period. But there is no 
c sign in output of this need. Probably the random collision length does
c enough to randomize even a bad coincidence.
c Iterate over adjustments to Eneutral.
      do i=1,ndims
         cdistflux(i)=0.
      enddo
      ncdist=0
      vzave=0.
      ttic=dt0
c--------------------------------------------
c Start of a new orbit
 2    continue 
c Inject from neutral distribution
      v2=delta*Ti
      do i=1,ndims
         v(i)=tisq*gasdev(i)
c         if(i.eq.ndims)v(i)=v(i)+vneutral
         v(i)=v(i)+vneutral*vdrift(i)
         v2=v2+v(i)**2
      enddo
c torb is the orbit time in velocity-scaled units. ttic in unscaled.
      torb=-alog(ran1(1)+1.e-15)*ctprime

c      write(*,*)'ctprime=',ctprime
c      write(*,*)ttic,torb,Eneutral,v,Ti,Bt
 3    if(ttic.le.torb)then
c If ttic<torb, advance to tic.
c Subtract the timestep from time to the orbit end in scaled units
         torb=torb-ttic
         v2=delta*Ti
         ncdist=ncdist+1
         if(Bt.eq.0.)then
c B-field free case
            do i=1,ndims
c               if(i.eq.ndims)v(i)=v(i)+Eneutral*ttic
               v(i)=v(i)+Eneutral*vdrift(i)*ttic
               v2=v2+v(i)**2
               v_col(i,ncdist)=v(i)
               cdistflux(i)=cdistflux(i)+abs(v(i))
            enddo
c            write(*,*)'ncdist,v(3)',ncdist,v(3),ttic,torb
         else
c Finite B-field case
            theta=Bt*ttic
c Rotation is counterclockwise for ions.
            stheta=sin(-theta)
            ctheta=cos(theta)
c Subtract off the perpendicular velocity of frame in which E is zero.
            do i=1,ndims
               v(i)=v(i)-vperp(i)
            enddo
c Rotate the velocity
            call rotate3(v,stheta,ctheta,Bfield)
c Acceleration along the B-direction:
            dvpar=0.
            do i=1,ndims
               dvpar=dvpar+Eneutral*vdrift(i)*Bfield(i)*ttic
            enddo
c And add back the drift velocity, and the parallel acceleration.
            do i=1,ndims
               v(i)=v(i)+vperp(i)+dvpar*Bfield(i)
               v2=v2+v(i)**2
               v_col(i,ncdist)=v(i)
               cdistflux(i)=cdistflux(i)+abs(v(i))
            enddo
c               write(*,*)'ncdist=',ncdist
c               write(*,*)'theta,v',theta,v
         endif
c         vzave=vzave+v_col(ndims,ncdist)
         do i=1,ndims
            vzave=vzave+v_col(i,ncdist)*vdrift(i)
         enddo
c Now we are at the end of a tic. Set the next tic length
         ttic=dt0
c If we haven't exhausted the sample number, start next tic
         if(ncdist.lt.nclim)goto 3
      else
c If ttic>torb, advance to next orbit without storing a sample
c Here, rather than just restarting the orbit with the neutral
c distribution we need to do the null collision calculation and goto 3
c not 2. 
c Advance to a possible collision.
         if(Bt.eq.0.)then
c            v(ndims)=v(ndims)+Eneutral*torb 
            do i=1,ndims
               v(i)=v(i)+Eneutral*vdrift(i)*torb
            enddo
         else
c Finite B-field case
            theta=Bt*torb
c Rotation is counterclockwise for ions.
            stheta=sin(-theta)
            ctheta=cos(theta)
c Subtract off the perpendicular velocity of frame in which E is zero.
            do i=1,ndims
               v(i)=v(i)-vperp(i)
            enddo
c Rotate the velocity
            call rotate3(v,stheta,ctheta,Bfield)
c Acceleration along the B-direction:
            dvpar=0.
            do i=1,ndims
               dvpar=dvpar+Eneutral*vdrift(i)*Bfield(i)*torb
            enddo
c And add back the drift velocity, and the parallel acceleration.
            do i=1,ndims
               v(i)=v(i)+vperp(i)+dvpar*Bfield(i)
            enddo
         endif
c Possible collision: Choose a neutral velocity and calculate (v-vd)^2
         v2=0.
         do i=1,ndims
            vn(i)=tisq*gasdev(i)
c            if(i.eq.ndims)vn(i)=vn(i)+vneutral
            vn(i)=vn(I)+vneutral*vdrift(i)
            v2=v2+(v(i)-vn(i))**2
         enddo
c Calculate the ratio of the actual collision time at this relative
c velocity to the possible-collision time:
         colnu=((v2+delta*Ti)/(Ti+delta*v2))**(colpow/2.)
     $        *colntime/ctprime
c Select a fraction of these collisions weighted by this ratio.
         rd=ran1(1)*colnu
         if(rd.lt.1.)then
c This collision occurred. 
            do i=1,ndims
               v(i)=vn(i)
            enddo
         else
c No collision. Just carry on.
         endif
         ttic=ttic-torb
c Choose next orbit time.
         torb=-alog(ran1(1)+1.e-15)*ctprime
         if(.not. torb.lt.100.*ctprime)then
            write(*,*)'torb problem',torb,ctprime
         endif
         goto 3
      endif
c      write(*,*)'v(i)',v
      end
c****************************************************************
      subroutine colndistshow()
      include 'ndimsdecl.f'
      include 'cdistcom.f'
      integer nxbin,nybin
      parameter (nxbin=50,nybin=50,idim1=1,idim2=3)
      real fvxvy(nxbin,nybin),vxa(nxbin,nybin),vya(nxbin,nybin)
      real work(0:nxbin+1,0:nybin+1)
      real vlimit(2,ndims)

c Defaults
      do id=1,ndims
         vlimit(1,id)=3.
         vlimit(2,id)=-2.
      enddo
c Adjust vlimits
      do j=1,ncdist-1
         vx=v_col(idim1,j)
         if(vx.gt.vlimit(1,1))vlimit(1,1)=vx
         if(vx.lt.vlimit(2,1))vlimit(2,1)=vx
         vy=v_col(idim2,j)
         if(vy.gt.vlimit(1,2))vlimit(1,2)=vy
         if(vy.lt.vlimit(2,2))vlimit(2,2)=vy
      enddo
c Velocity arrays
      do i=1,nxbin
      do j=1,nybin
         vxa(i,j)=vlimit(1,1)+(vlimit(2,1)-vlimit(1,1))*(i-1.)/(nxbin
     $        -1.)
         vya(i,j)=vlimit(1,2)+(vlimit(2,2)-vlimit(1,2))*(j-1.)/(nybin
     $        -1.)
      enddo
      enddo

c Accumulate the fvxvy distribution
         do k=1,nybin
            do j=1,nxbin
               fvxvy(j,k)=0.
            enddo
         enddo
         do j=1,ncdist-1
            vx=v_col(idim1,j)
            vy=v_col(idim2,j)
            nx=nint(nxbin*(vx-vlimit(1,1))/(vlimit(2,1)-vlimit(1,1)))
            ny=nint(nybin*(vy-vlimit(1,2))/(vlimit(2,2)-vlimit(1,2)))
            if(nx.ge.1.and.nx.le.nxbin.and.ny.ge.1.and.ny.le.nybin)
     $           then
               fvxvy(nx,ny)=fvxvy(nx,ny)+1.
            endif
         enddo
c      write(*,*)vxa
98     call pltinit(0.,1.,0.,1.)
c       Plot the surface. With axes (2-D). Web color 10, axis color 7.
        j=2 + 256*10 + 256*256*7
        call surf3d(vxa,vya,fvxvy,nxbin,nxbin,nybin,j,work)
        call color(15)
        if(ieye3d().ne.0) goto 98

        open(2,file='fvxvydist.dat')
        write(2,*)(fvxvy(nxbin/2,j),j=1,nybin)
        close(2)
      end

c****************************************************************
