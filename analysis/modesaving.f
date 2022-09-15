c********************************************************************
      subroutine savemodes(nff,nuf,nfx,nux,nmodes,phimodes,time,xn)
! write to file the phimodes we have found
      integer nff,nuf,nfx,nux,nmodes
      complex phimodes(nff,nfx,nmodes)
      real time(nff),xn(nfx)
      character*30 string
      string='1 saved complex phimodes'
      open(8,file='savedmodes.dat',status='new',form='unformatted')
      write(8)string
      write(8)nuf,nux,nmodes
!      write(*,*)nuf,nux,nmodes,nff,nfx
      write(8)(((phimodes(i,j,m),i=1,nuf),j=1,nux),m=1,nmodes)
      write(8)(time(i),i=1,nuf)
      write(8)(xn(j),j=1,nux)
      close(8)
!      write(*,*)'Wrote phimode(10,5,2)',phimodes(10,5,2)
      end
c********************************************************************
      subroutine unsavemodes(nff,nuf,nfx,nux,nmodes,nmr,phimodes,time,xn   &
     &     ,ierr)
! read from file the phimodes saved.
! on exit ierr =0 if success, =1 if saved data unavailable.
      integer nff,nuf,nfx,nux,nmodes,ierr
      complex phimodes(nff,nfx,nmodes)
      real time(nff),xn(nfx)
      character*30 string
      ierr=0
      open(9,file='savedmodes.dat',status='old',form='unformatted'         &
     &     ,err=101)
      read(9)string
      write(*,*)'Reading back version ',string
      read(9)nuf,nux,nmr
      if(nuf.gt.nff.or.nux.gt.nfx.or.nmr.gt.nmodes)then
         write(*,*)'File used dimensions', nuf,nux,nmr
         write(*,*)'incompatible with',nff,nfx,nmodes,' declared.'
         stop
      endif
      write(*,*)'       nfiles        nx        nmodes   ',
     $     '   nfmax       nxmax'
      write(*,*)nuf,nux,nmr,nff,nfx
      read(9)(((phimodes(i,j,m),i=1,nuf),j=1,nux),m=1,nmr)
      read(9)(time(i),i=1,nuf)
      read(9)(xn(j),j=1,nux)
      close(9)
!      write(*,*)'Read phimode(10,5,2)',phimodes(10,5,2)
      return
 101  continue
      write(*,*)'No file savedmodes.dat. Creating modes'
      ierr=1
      end
