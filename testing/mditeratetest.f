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
      iused(1)=99
      iused(2)=100
      iused(3)=100
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
         call mditerarg(rbroutine,mdims,ifull,iused,ipoint,u)
      enddo
c      stop

      iform=1
      uscale=1
      call udisplay(mdims,u,ifull,iused,iform,uscale)

      write(*,*)'******** second call **********'
      ipoint=0
      call mditerarg(bdyroutine,mdims,ifull,iused,ipoint,u)
      ipoint=0
      uscale=1
      call udisplay(mdims,u,ifull,iused,iform,uscale)

      do j=1,if3
            u(j)=0
      enddo

      write(*,*)'******** Third Call Restricted **********'
      
      ipoint=0
      call mditerarg(rbroutine,mdims,ifull,iu2,
     $     ipoint,u(3+2*ifull(1)))
      ipoint=0
      uscale=1
      call udisplay(mdims,u,ifull,iused,iform,uscale)

      end
c************************************************************************
c New mditerarg is slower, 8+s vs 5+s optimized. Profiling:
c 1000 x 99x100^2
c time   seconds   seconds    calls   s/call   s/call  name    
c 31.94      4.45     4.45 495000000     0.00     0.00  rbroutine_
c 25.25      7.96     3.51 495001000     0.00     0.00  mditerator_
c 21.08     10.89     2.93 495000000     0.00     0.00  indexcontract_
c 15.54     13.06     2.16     1000     0.00     0.01  mditerarg_
c Old mditerarg:
c time   seconds   seconds    calls   s/call   s/call  name    
c 50.78      3.57     3.57 495000000     0.00     0.00  rbroutine_
c 46.23      6.83     3.25     1000     0.00     0.01  mditerarg_
c Looks like having the extra calls is the main cost.
c Probably it is a poor idea actually to replace mditerarg.
