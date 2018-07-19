! Quiet Particle initialization.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This version is intended for situations with no objects and no collisions.
! The uniformity at low-k is improved by depositing uniformly by qblocks,
! with random placement within the qblock. When there are fewer particles
! remaining than the total number of qblocks, the qblock size is increased
! reducing the total number of qblocks and the process iterated till all
! particles are placed.

! The number of qblocks per dimension is kept less than the number of
! cells-2 and 30, but greater than zero. It is initialized as nqbset but
! constrained by those limits.

      subroutine qinit()
      implicit none
! Common data:
      include 'ndimsdecl.f'
      include 'partcom.f'
      include 'plascom.f'
      include 'myidcom.f'
      include 'meshcom.f'
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

! Till nqbset parameter made adjustable, set it algorithmically.
! Later we'll remove the nqbset setting.
      do i=1,ndims
         mlen=ixnp(i+1)-ixnp(i)
         nqbset(i)=max(1,min(30,mlen-2))
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
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine fillqblocksize(nqblks,islot,islotmax,ispecies)
! Fill the qblocks of current size with as many complete sets of particles 
! as possible up to nps, leaving a remainder. Increment islot for each
! particle added.
      include 'ndimsdecl.f'      
      include 'meshcom.f'
      include 'myidcom.f'
      integer nqblks(ndims),islot,islotmax,ispecies

      integer iview(3,ndims),indi(ndims)     !Iterator
      integer nqbt

      if(myid.eq.0)then
!        write(*,*)'fillqblocksize',nqblks,islot,islotmax
         write(*,*)'        i    nfills     islot       nqblks'
      endif
      nremain=islotmax-islot+1
      nqbt=1
      do i=1,ndims
         nqbt=nqbt*nqblks(i)
      enddo
 2    nfills=int(nremain/nqbt)      ! assuming 1 particle per fill per qblock.
      nremain=nremain-nfills*nqbt
      do i=1,nfills
! Fill using iterator
         icomplete=mditerator(ndims,iview,indi,4,nqblks)
 1       continue
! Place one particle in block indi.
         call placeqblk(nqblks,indi,islot,ispecies)         
         islot=islot+1
         if(mditerator(ndims,iview,indi,0,nqblks).eq.0)goto 1
         if(mod(i,max(1,nfills/30)).eq.0.and.myid.eq.0)
     $        write(*,'(8i10)')i,nfills,islot,nqblks
      enddo

      if(nremain.le.0)return
! Else reduce nqblks by 2, and redo
      nqbt=1
      do i=1,ndims
         nqblks(i)=max(1,nqblks(i)/2)
         nqbt=nqbt*nqblks(i)
      enddo
      goto 2

      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine placeqblk(nqblks,indi,islot,ispecies)
! Place a particle randomly in the block addressed by indi
! This version just uniform distribution.
      include 'ndimsdecl.f'      
      include 'meshcom.f'
      include 'partcom.f'
      include 'plascom.f'
      include 'myidcom.f'
      integer nqblks(ndims),indi(ndims),islot,ispecies
      logical linmesh

 1    do i=1,ndims              ! Position
         call ranlux(ran,1)
         fp=(indi(i)+ran)/float(nqblks(i))
         fp=max(.000001,min(.999999,fp))
! Here's where to adjust spatial distribution.
         x_part(i,islot)=(1.-fp)*xmeshstart(i)+fp*xmeshend(i)
      enddo

      tisq=sqrt(Ts(ispecies)*abs(eoverms(ispecies)))
! Shifted Gaussians.                     ! Velocities
      x_part(4,islot)=tisq*gasdev(myid) + vds(ispecies)*vdrift(1)
      x_part(5,islot)=tisq*gasdev(myid) + vds(ispecies)*vdrift(2)
      x_part(6,islot)=tisq*gasdev(myid) + vds(ispecies)*vdrift(3)
!      write(*,*)islot,(x_part(i+3,islot),i=1,3)
!      write(*,*)ispecies,nspecies,vds(ispecies),vdrift(3)
! The previous timestep length.
      x_part(idtp,islot)=0.
! Initialize the mesh fraction data in x_part.
      call partlocate(x_part(1,islot),ixp,xfrac,iregion,linmesh)
! This test rejects particles exactly on mesh boundary:
      if(.not.linmesh)goto 1
      x_part(iflag,islot)=1

      end
