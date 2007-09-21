c***********************************************************************
c Initializing particles.
      subroutine pinit()
c Common data:
      include 'partcom.f'
      include 'plascom.f'

      ntries=0
c      ntrapped=0
      rmax=0.99999*rs
      rmax2=rmax*rmax
c     We initialize the 'true' particles'
      do i=1,n_part
         if_part(i)=1
 1       continue
         ntries=ntries+1
         x_part(1,i)=rmax*(2.*ran0(idum)-1.)
         x_part(2,i)=rmax*(2.*ran0(idum)-1.)
         x_part(3,i)=rmax*(2.*ran0(idum)-1.)
         inewregion=insideall(npdim,x_part(1,i))
c     If we are not in the plasma region, try again.
         if(inewregion.ne.iregion_part)then
c            write(*,*)'Injection of',i,' wrong region',inewregion,
c     $           ' instead of',iregion_part
c     $           ,(x_part(kk,i),kk=1,3)
            goto 1
         endif
         Ti0=Ti
         tisq=sqrt(Ti0)
         x_part(4,i)=tisq*gasdev(idum)
         x_part(5,i)=tisq*gasdev(idum)
         x_part(6,i)=tisq*gasdev(idum) + vd
      enddo
c Set flag of unused slots to 0
      do i=n_part+1,n_partmax
         if_part(i)=0
      enddo
      write(*,*)'Initialized ','id=',myid,
     $     '  n=',n_part,'  ntries=',ntries
c     $     ,'  ntrapped=',ntrapped
c Initialize rhoinf:
      if(rhoinf.eq.0.)rhoinf=numprocs*n_part/(4.*pi*rmax**3/3.)
c Initialize orbit tracking
      do ko=1,norbits
         iorbitlen(ko)=0
      enddo
      end
c***********************************************************************
