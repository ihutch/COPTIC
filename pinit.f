c***********************************************************************
c Initializing particles.
      subroutine pinit()
c Common data:
      include 'partcom.f'
      include 'plascom.f'
      include 'myidcom.f'
      include '3dcom.f'
      external linregion
      logical linregion

      ntries=0
c      ntrapped=0
      rmax=0.99999*rs
c      rmax2=rmax*rmax
c The maximum used slot is the same as the number of particles initially
      ioc_part=n_part
c     We initialize the 'true' particles'
      do i=1,n_part
         if_part(i)=1
 1       continue
         ntries=ntries+1
         x_part(1,i)=rmax*(2.*ran1(myid)-1.)
         x_part(2,i)=rmax*(2.*ran1(myid)-1.)
         x_part(3,i)=rmax*(2.*ran1(myid)-1.)
c     If we are not in the plasma region, try again.
         if(.not.linregion(ibool_part,npdim,x_part(1,i)))then
c            write(*,*)'Injection of',i,' wrong region',inewregion
c     $           ,(x_part(kk,i),kk=1,3)
            goto 1
         endif
         Ti0=Ti
         tisq=sqrt(Ti0)
         x_part(4,i)=tisq*gasdev(myid)
         x_part(5,i)=tisq*gasdev(myid)
         x_part(6,i)=tisq*gasdev(myid) + vd
      enddo
c Set flag of unused slots to 0
      do i=n_part+1,n_partmax
         if_part(i)=0
      enddo
c      write(*,*)'Initialized ','id=',myid,
c     $     '  n=',n_part,'  ntries=',ntries
c Initialize rhoinf:
      if(rhoinf.eq.0.)rhoinf=numprocs*n_part/(4.*pi*rmax**3/3.)
c Initialize orbit tracking
      do ko=1,norbits
         iorbitlen(ko)=0
      enddo
      end
c***********************************************************************
