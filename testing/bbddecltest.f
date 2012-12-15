      integer ndimsdecl
      parameter (ndimsdecl=3)
      include 'bbdydecl.f'

      do i=1,ndimsdecl
         iLs(i)=1
         lperiod(i)=.false.
         ifull(i)=64
         iuds(i)=32
         icoords(i)=1
         myside(i)=24
         idims(i)=0
      enddo

      mycartid=9999
      myid=8
      icommcart=333
      myorig=1
      write(*,*)iLs
      write(*,*)ifull
      write(*,*)iuds
      write(*,*)idone
      write(*,*)ndimsdecl
      write(*,*)idims
      write(*,*)lperiod
      write(*,*)icoords
      write(*,*)iLcoords
      write(*,*)myside
      write(*,*)myorig
      write(*,*)icommcart
      write(*,*)mycartid
      write(*,*)myid

      call passed(iLs,ifull,iuds,idone,ndimsdecl,idims,lperiod
     $     ,icoords,iLcoords,myside,myorig,icommcart,mycartid,myid)

      end
c**************************************
      subroutine passed(iLs,ifull,iuds,idone,ndimsdecl,idims,lperiod
     $     ,icoords,iLcoords,myside,myorig,icommcart,mycartid,myid)
c      parameter (ndimsdecl=3)
      include 'bbdydecl.f'


      write(*,*)'passed::::::::::::::::::::'

      write(*,*)iLs
      write(*,*)ifull
      write(*,*)iuds
      write(*,*)idone
      write(*,*)ndimsdecl
      write(*,*)idims
      write(*,*)lperiod
      write(*,*)icoords
      write(*,*)iLcoords
      write(*,*)myside
      write(*,*)myorig
      write(*,*)icommcart
      write(*,*)mycartid
      write(*,*)myid


      end
