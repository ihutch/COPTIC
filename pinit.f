c***********************************************************************
c Initializing particles.
      subroutine pinit(sprior)
c Common data:
      include 'partcom.f'
      include 'plascom.f'
      include 'myidcom.f'
      include '3dcom.f'
      include 'meshcom.f'
      include 'colncom.f'
      external linregion
      logical linregion
c Local dummy variables for partlocate.
      real xfrac(ndims_mesh)
      integer ixp(ndims_mesh)
      logical linmesh

c-----------------------------------------------------------------
c A special orbit.
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
      if_part(1)=1
      i1=1
      call partlocate(1,ixp,xfrac,iregion,linmesh)
      if(.not.linmesh.and.myid.eq.0)then
         write(*,*)'WARNING Special particle-1 outside region.'
         write(*,*)'Coordinates',(x_part(i,1),i=1,9),' Resetting.'
      endif
      if(linmesh)i1=2
c-----------------------------------------------------------------
      ntries=0
c      ntrapped=0
c The maximum used slot is the same as the number of particles initially
      ioc_part=n_part
      call colninit(0)
      call colreinit()
c     We initialize the 'true' particles'
      tisq=sqrt(Ti)
      do i=i1,n_part
         if_part(i)=1
 1       continue
c         write(*,'(i8,$)')i
         ntries=ntries+1
         do j=1,ndims_mesh
            x_part(j,i)=xmeshstart(j)+
     $        ran1(myid)*(xmeshend(j)-xmeshstart(j))
         enddo
c     If we are not in the plasma region, try again.
         if(.not.linregion(ibool_part,npdim,x_part(1,i)))then
c            call partlocate(1,ixp,xfrac,inewregion,linmesh)
c            write(*,*)'Initialization of',i,' wrong region',inewregion
c     $           ,(x_part(kk,i),kk=1,3)
            goto 1
         endif
c Here is the section for which we want to fashion an alternative:
         if(Eneutral.eq.0.)then
c         if(.true.)then
            x_part(4,i)=tisq*gasdev(myid)
            x_part(5,i)=tisq*gasdev(myid)
            x_part(6,i)=tisq*gasdev(myid) + vd
         else
            call colvget(x_part(4,i))
         endif

c One idea that might improve the initialization would be to add the 
c (electric) potential energy here, by scaling up the total velocity,
c but keeping the direction the same. 

c However, to make this accurate, one would have to have solved the
c initial potential better than the simple Laplace solve currently
c used. A solution of the electron-shielded case might be nearer the
c mark, but it would not be correct when Ti<Te. It is perhaps
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

         dtprec(i)=0.
c Initialize the mesh fraction data in x_part.
         call partlocate(i,ixp,xfrac,iregion,linmesh)
c This test rejects particles exactly on mesh boundary:
         if(.not.linmesh)goto 1
      enddo
c Set flag of unused slots to 0
      do i=n_part+1,n_partmax
         if_part(i)=0
      enddo
      if(myid.eq.0)then
       write(*,'('' Initialized '',i3,''x n='',i7,'' ntries='',i7,$)')
     $     nprocs,n_part,ntries
      endif
c Initialize rhoinf:
      if(rhoinf.eq.0.)rhoinf=numprocs*n_part/(4.*pi*rs**3/3.)
c Initialize orbit tracking
      do ko=1,norbits
         iorbitlen(ko)=0
      enddo
      end
c***********************************************************************
      subroutine locateinit()
c Initialize particles loaded from restart files so that they have the
c correct mesh-fraction coordinate in their upper values.
c If they are outside the mesh or particle region, drop them.
c Common data:
      include 'partcom.f'
c      include 'myidcom.f'
      include '3dcom.f'
      include 'meshcom.f'
c Local dummy variables for partlocate.
      real xfrac(ndims_mesh)
      integer ixp(ndims_mesh)
      logical linmesh,linregion
      external linregion

      do i=1,ioc_part
         call partlocate(i,ixp,xfrac,iregion,linmesh)
c         write(*,*)linmesh,linregion(ibool_part,ndims_mesh,x_part(1,i))
         if(.not.linmesh .or. .not.
     $        linregion(ibool_part,ndims_mesh,x_part(1,i)))if_part(i)=0
      enddo
      nrein=0
c      write(*,*)'Redetermined particle mesh locations (locateinit)'
      end
c***********************************************************************
c Return a random velocity from the (precalculated) distribution heap.
      subroutine colvget(v)
      implicit none
      include 'cdistcom.f'
      real v(nc_ndims)
      integer nc,i
      real ran1
      external ran1

      nc=ncdist*ran1(1)
      do i=1,nc_ndims
         v(i)=v_col(i,nc)
      enddo
      end
c*********************************************************************
      subroutine colninit(ncin)
c Routine for initializing particles when collisions and drifts
c driven by Eneutral are present. An array of velocities that 
c represents the background distribution of ions is calculated by
c advancing from random birth velocity representative of the neutral
c distribution. The collision frequency is proportional to velocity
c to the power colpow.

      implicit none
      integer ncin
      include 'cdistcom.f'
      include 'plascom.f'
      include 'colncom.f'
      real gasdev,ran1
      external gasdev, ran1

c Local variables
      integer nclim,i,k
      real tisq,dt0,torb,ttic,v2,vzave,orbtic
      real v(nc_ndims)
c The number of samples.
      if(ncin.eq.0)then
         nclim=ncdistmax
      else
         nclim=ncin
      endif
      tisq=sqrt(Ti)
      if(Eneutral.eq.0.)then
c No Eneutral. Initialize just from the neutral distribution.
         dt0=0.
c For now, don't initialize
         return
      else
c Decide the duration of the time tic. Take it to be a smallish fraction
c of the collision time so that sufficient resolution is present to do a
c reasonable integral of a variable collision frequency, but not too
c small otherwise the perpendicular velocity and orbit length
c representation becomes statistically worse.
c The actual collision time is coltime*(v^2/Ti)^(colpow/2.)
c So if colpow=0. the collision time is coltime independent of v.
         if(colpow.eq.0.)then
            dt0=.5*colntime
         else
            dt0=0.03*colntime
         endif
      endif
      do i=1,nc_ndims
         cdistflux(i)=0.
      enddo

c Initial guess at what Eneutral really needs to be to give the
c drift requested. Based on Ti=0.1.
      Eneutral=Eneutral*(1.+(2.5*colpow+2.)*colpow)
      write(*,'(a,$)')' Colinit velocity: '

c Iterate over adjustments to Eneutral.
      do k=1,4
      ncdist=0
      vzave=0.
      ttic=dt0
c--------------------------------------------
c Start of a new orbit
 2    continue 
c torb is the orbit time in velocity-scaled units. ttic in unscaled.
      torb=-alog(ran1(1))*colntime
c Inject from neutral distribution
      v2=0.
      do i=1,nc_ndims
         v(i)=tisq*gasdev(i)
         if(i.eq.nc_ndims)v(i)=v(i)+vneutral
         v2=v2+v(i)**2
      enddo
c The torb consumption rate depends on velocity.
c For example, if cross-section is proportional to velocity to the 
c power -p, then dynamic colnfreq \propto v^{1.-p}. =v^2(1.-p)/2.
c Constant cross-section pow=colpow*0.5=0.5; constant nu colpow=0.
c The ratio of orbit consumption rate to tic consumption rate is:
      orbtic=(v2/Ti)**(colpow/2.)

 3    if(ttic*orbtic.le.torb)then
c If ttic*orbtic<torb, advance to tic.
c Subtract the timestep from time to the orbit end in scaled units
         torb=torb-ttic*orbtic
         v2=0.
         if(Bt.eq.0.)then
c B-field free case
            ncdist=ncdist+1
            do i=1,nc_ndims
               if(i.eq.nc_ndims)v(i)=v(i)+Eneutral*ttic
               v2=v2+v(i)**2
               v_col(i,ncdist)=v(i)
               cdistflux(i)=cdistflux(i)+abs(v(i))
            enddo
         endif
         vzave=vzave+v_col(nc_ndims,ncdist)
         orbtic=(v2/Ti)**(colpow/2.)
c Now we are at the end of a tic. Set the next tic length
         ttic=dt0
c If we haven't exhausted the sample number, start next tic
         if(ncdist.lt.nclim)goto 3
      else
c If ttic*orbtic>torb, advance to next orbit without storing a sample
         ttic=ttic-torb/orbtic
         goto 2
      endif
c--------------------------------------------
c      write(*,*)'End of colninit',dt0,colntime,colpow
c     $     ,vzave/nclim,Eneutral

      write(*,'(f6.4,'', '',$)')vzave/nclim
      if(colpow.eq.0.)goto 4
      Eneutral=Eneutral*(1.5*((vd-vneutral)*nclim/vzave-1.)+1.)
      enddo
 4    continue
      write(*,'(a,f8.4)')' Eneutral=',Eneutral

      end

c The Eneutral=(vd-vn)/colntime. This defines what colntime must be
c for a variable collision frequency. Suppose nu=A v^p, then the mean
c velocity is given by \int 
