      program pointtest
c Test the accis point plotting ability
      parameter (npts=5000,ny=20,v=.03)
      real x(npts),y(npts)

      call pfset(3)
      do kk=1,10
      call pltinit(0.,1.,-1.,1.)
      call axis()
      call winsetmargin(.true.,0.002)
      do j=1,ny
         call color(mod(j,16))
         do i=1,npts
            x(i)=(i-1.)/(npts-1.)+v*kk
            y(i)=sin(10.*(x(i)-v*kk))*exp(-x(i))+.2*j/ny
            call vecw(x(i),y(i),-1)
         enddo
      enddo
      call accisflush()
      call usleep(1000)
      enddo

!      call pltend()

      end

! Tests show glx driver takes 1.0 sec for 100,000 points. 9.6s for 1M.
!            vecx driver 0.6s 100k, 5.5s 1M.  Both for 10 plots.
! Thus we are drawing about 2M points per second.
! That is so fast as not to be worth playing around in X to fix. 
