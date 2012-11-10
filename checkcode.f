c***********************************************************************
c This contains code that checks on the writing and reading back of
c code state. It is for debugging restarts.
c***********************************************************************
      subroutine checkuqcij(ifull,u,q,psum,volumes,cij)
      implicit none
      integer ifull(3)
      real u(ifull(1),ifull(2),ifull(3)),q(ifull(1),ifull(2),ifull(3))
      real psum(ifull(1),ifull(2),ifull(3)),volumes(ifull(1),ifull(2)
     $     ,ifull(3))
      real cij(7,ifull(1),ifull(2),ifull(3))
c Local duplicates.
      integer Li2
      parameter (Li2=32)
      real psum2(Li2,Li2,Li2),volumes2(Li2,Li2,Li2)
      real u2(Li2,Li2,Li2),q2(Li2,Li2,Li2)
      real cij2(7,Li2,Li2,Li2)

      integer i,j,k,l,ic
      logical linit,lend
      data linit/.false./lend/.false./
      save

      if(.not.linit)then
         linit=.true.
         open(41,file='uqcijnew',status='new',
     $        form='unformatted',err=110)
         open(40,file='uqcij',status='old',
     $        form='unformatted',err=101)
         goto 105
 101     write(*,*)
         write(*,*)'=====Could not open uqcij file. Just writing.'
         close(40)
         lend=.true.
 105     continue
      endif

      write(41,err=103)(((u(i,j,k),q(i,j,k),(cij(l,i,j,k),l=1,7),psum(i
     $     ,j,k),volumes(i,j,k),i=1,Li2),j=1,Li2),k=1,Li2)

      if(lend) return

      read(40,err=102,end=102)(((u2(i,j,k),q2(i,j,k),(cij2(l,i,j,k),l=1
     $     ,7),psum2(i,j,k),volumes2(i,j,k),i=1,Li2),j=1,Li2),k=1,Li2)

      write(*,*)
      write(*,*)'Successfully read:u2,q2,cij2,psum2,volumes2'

      ic=0
      do k=1,min(Li2,ifull(3))
         do j=1,min(Li2,ifull(2))
            do i=1,min(Li2,ifull(1))

      if(u(i,j,k).ne.u2(i,j,k))then
         write(*,*)'***** u difference',i,j,k,u(i,j,k),u2(i,j,k)
         ic=ic+1
      endif
      if(q(i,j,k).ne.q2(i,j,k))then
         write(*,*)'***** q difference',i,j,k,q(i,j,k),q2(i,j,k)
         ic=ic+1
      endif
      if(psum(i,j,k).ne.psum2(i,j,k))then
        write(*,*)'***** psum difference',i,j,k,psum(i,j,k),psum2(i,j,k)
         ic=ic+1
      endif
      if(volumes(i,j,k).ne.volumes2(i,j,k))then
         write(*,*)'***** volumes difference',i,j,k,
     $        volumes(i,j,k),volumes2(i,j,k)
         ic=ic+1
      endif

      do l=1,7
         if(cij(l,i,j,k).ne.cij2(l,i,j,k))then
            write(*,'(''***** cij difference'',4i3,2e14.6)')
     $           l,i,j,k,cij(l,i,j,k),cij2(l,i,j,k)
            ic=ic+1
         endif
      enddo
      if(ic.gt.12) goto 100
          enddo
        enddo
      enddo

 100  continue
      write(*,*)'====== Finished uqcijckeck'

      return

 102  write(*,*)'*****Unexpected End of uqcij file'
      lend=.true.
      return

 103  write(*,*)'======Could not write to uqcijnew'
      lend=.true.
      return

 110  write(*,*)
      write(*,*)'*****Could not open uqcijnew'
      return

      end
c***********************************************************************
      subroutine checkx()
      implicit none
c The storage used here defines how far our checking goes for bigger
c actual cases. Don't want to use too much, or too little.
      integer icp
      parameter (icp=100000)
      real x_part2(9,icp)
      integer n_part2,if_part2(icp),iregion_part2,ioc_part2
      integer nrein2,numprocs2,ninjcomp2
      real dt2,rhoinf2,phirein2
      logical ldiags2

      include 'partcom.f'
      integer i,j,ic,ii,jj
      logical linit,lend
      data linit/.false./lend/.false./
      save

      if(.not.linit)then
         linit=.true.
         open(51,file='xpartnew',status='new',
     $        form='unformatted',err=110)
         open(50,file='xpart',status='old',
     $        form='unformatted',err=101)
         goto 105
 101     write(*,*)'=====Could not open xpart file. Just writing.'
         close(50)
         lend=.true.
 105     continue
      endif

      write(51,err=103)ioc_part
      write(51,err=103)n_part
     $     ,((x_part(ii,jj),ii=1,9),if_part(jj),jj=1,min(ioc_part,icp))
     $     ,iregion_part,
     $     dt,ldiags,rhoinf,nrein,phirein,numprocs,ninjcomp
      if(lend) return

      read(50,err=102)ioc_part2
      read(50,err=102)n_part2
     $     ,((x_part2(ii,jj),ii=1,9),if_part2(jj),
     $     jj=1,min(ioc_part2,icp))
     $     ,iregion_part2,
     $     dt2,ldiags2,rhoinf2,nrein2,phirein2,numprocs2,ninjcomp2

      if(ioc_part.ne.ioc_part2.or.n_part.ne.n_part2.or.
     $     iregion_part.ne.iregion_part2) then
         write(*,*)'***** particle count difference',
     $        ioc_part,n_part,iregion_part,
     $        ioc_part2,n_part2,iregion_part2
      endif
         
      if(dt.ne.dt2.or.ldiags.neqv.ldiags2.or.rhoinf.ne.rhoinf2.or.
     $     nrein.ne.nrein2.or.phirein.ne.phirein2.or.
     $     numprocs.ne.numprocs2.or.ninjcomp.ne.ninjcomp2)then
         write(*,*)'***** dt etc difference',
     $     dt,ldiags,rhoinf,nrein,phirein,numprocs,ninjcomp,
     $     dt2,ldiags2,rhoinf2,nrein2,phirein2,numprocs2,ninjcomp2
      endif

      ic=0
      do j=1,min(ioc_part,icp)
         if(if_part(j).ne.if_part2(j))then
            write(*,*)'***** if_part difference',
     $           j,if_part(j),if_part2(j)
            ic=ic+1
         endif
         do i=1,9
            if(x_part(i,j).ne.x_part2(i,j))then
               write(*,*)'***** x_part difference',
     $              i,j,x_part(i,j),x_part2(i,j)
               ic=ic+1
            endif
            if(ic.gt.50) goto 100
         enddo
      enddo

 100  continue
      write(*,*)'====== Finished xpartcheck'

      return

 102  write(*,*)'*****End of xpart file'
      lend=.true.
      return

 103  write(*,*)'======Could not write to xpartnew'
      lend=.true.
      return

 110  write(*,*)'*****Could not open xpartnew'
      return

      end
c****************************************************************
      subroutine checkdelta(delta,deltaold)

      logical linit,lend
      data linit/.false./lend/.false./
      save

      if(.not.linit)then
         linit=.true.
         open(31,file='checknew',status='new',
     $        form='unformatted',err=110)
         open(30,file='checkdelta',status='old',
     $        form='unformatted',err=101)
      endif

      write(31,err=103)delta

      if(lend) return

      read(30,err=102,end=102)deltaold

      if(delta.eq.deltaold)then
         return
      else
         write(*,*)'***** Delta difference',deltaold,delta
      endif

      return

 101  write(*,*)'*****Could not open checkdelta file.'
      close(30)
      return

 102  write(*,*)'*****End of checkdelta file'
      lend=.true.
      return

 103  write(*,*)'Could not write to checknew'
      lend=.true.
      return

 110  write(*,*)'*****Could not open checknew'
      return

      end
