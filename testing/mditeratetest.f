c Testing:
c******************************************************************
      program mditeratetest
      integer ndims
      parameter (ndims=3)
      integer ifull(ndims)
      integer iused(ndims),iu2(ndims)
      external rbroutine,bdyroutine

      parameter (ifull1=100,ifull2=100,ifull3=100)
      parameter (if3=ifull1*ifull2*ifull3)
      real u(if3)


      ifull(1)=ifull1
      ifull(2)=ifull2
      ifull(3)=ifull3
      iused(1)=20
      iused(2)=10
      iused(3)=4
      iu2(1)=iused(1)-4
      iu2(2)=iused(2)-4
      iu2(3)=iused(3)

      do j=1,if3
            u(j)=0.
      enddo


      write(*,*)'**************** first call ************'
      mdims=3
      do i=1,1
         ipoint=mod(i+1,2)
         call mditerate(rbroutine,mdims,ifull,iused,ipoint,u)
      enddo
      iform=1
      uscale=1
      call udisplay(mdims,u,ifull,iused,iform,uscale)

      write(*,*)'******** second call **********'
      ipoint=0
      call mditerate(bdyroutine,mdims,ifull,iused,ipoint,u)
      ipoint=0
      uscale=1
      call udisplay(mdims,u,ifull,iused,iform,uscale)

      do j=1,if3
            u(j)=0
      enddo

      write(*,*)'******** Third Call Restricted **********'
      
      ipoint=0
      call mditerate(rbroutine,mdims,ifull,iu2,
     $     ipoint,u(3+2*ifull(1)))
      ipoint=0
      uscale=1
      call udisplay(mdims,u,ifull,iused,iform,uscale)

      end
c************************************************************************
c
c Test of 10x 100x100x100 takes 0.933s
c 100 takes 11.4s. Essentially identical for 10000x 100x10x10
c Basic speed for this large mesh is 10M mesh points per second.
c On laptop.
