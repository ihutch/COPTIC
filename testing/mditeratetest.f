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
      iused(1)=10
      iused(2)=10
      iused(3)=9
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
c      do i=1,1000
      do i=1,1
         ipoint=mod(i+1,2)
c or mditerarg1 is the reimplemented version.
         call mditerarg(rbroutine,mdims,ifull,iused,ipoint,u,t,v,w)
      enddo
c      stop

      iform=1
      uscale=1
      call udisplay(mdims,u,ifull,iused,iform,uscale)

c      write(*,*)'******** second call **********'
c      ipoint=0
c      call mditerarg(bdyroutine,mdims,ifull,iused,ipoint,u,t,v,w)
c      ipoint=0
c      uscale=1
c      call udisplay(mdims,u,ifull,iused,iform,uscale)

      do j=1,if3
            u(j)=0
      enddo
      write(*,*)'******** Third Call Restricted **********'
      
      ipoint=0
      call mditerarg(rbroutine,mdims,ifull,iu2,
     $     ipoint,u(3+2*ifull(1)),t,v,w)
      ipoint=0
      uscale=1
      call udisplay(mdims,u,ifull,iused,iform,uscale)

      end
c************************************************************************
c New mditerarg is slower, 8+s vs 5+s optimized. Profiling:
c 1000 x 99x100^2
c  %   cumulative   self              self     total 
c time   seconds   seconds    calls   s/call   s/call  name    
c 31.94      4.45     4.45 495000000     0.00     0.00  rbroutine_
c 25.25      7.96     3.51 495001000     0.00     0.00  mditerator_
c 21.08     10.89     2.93 495000000     0.00     0.00  indexcontract_
c 15.54     13.06     2.16     1000     0.00     0.01  mditerarg_
c
c Old mditerarg:
c time   seconds   seconds    calls   s/call   s/call  name    
c 50.78      3.57     3.57 495000000     0.00     0.00  rbroutine_
c 46.23      6.83     3.25     1000     0.00     0.01  mditerarg_
c Looks like having the extra calls is the main cost.
c Probably it is a poor idea actually to replace mditerarg.

c******************************************************************
      subroutine rbroutine(inc,ipoint,indi,ndims,iLs,iused,u)
c Red/black routine, process squares alternately.
      integer ipoint,inc
      integer indi(ndims),iused(ndims),iLs(ndims+1)
      real u(*)

      parameter (mdims=10)
      integer ind1(mdims)
      data icount/0/
      data ind1/mdims*1000/

c iused is not used. To silence warnings spuriously use it:
      ic=iused(1)
c      write(*,'(a,11i8)')'ipoint,indi',ipoint,(indi(i),i=1,ndims)
c      write(*,*)'myipoint=',ipoint,'  iused',iused
c We build in correction of the increment here for red-black
c to change the parity. We have to remember the previous indis.
      ic=0
      do i=ndims,2,-1
         if(indi(i).gt.ind1(i))then
c we've carried
            ic=ic+(i-1)
c            write(*,*)'increment',ic
         endif
c Remember prior.
         ind1(i)=indi(i)
      enddo
      if(mod(ic,2).ne.0) then
         iaj=(1-2*mod(indi(1),2))
         indi(1)=indi(1)+iaj
         ipoint=ipoint+iaj
      endif

c Here is where we process the square chosen. Examples:
c sequence
c      u(ipoint+1)=mod(icount,9)+1
c height
      u(ipoint+1)=indi(3)+1
      icount=icount+1
      inc=2
      end
