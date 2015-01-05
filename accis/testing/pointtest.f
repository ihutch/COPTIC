      program pointtest
c Test the accis point plotting ability
      parameter (npts=300,ny=20)
      real x(npts),y(npts)

      call pfset(3)
      call pltinit(0.,1.,-1.,1.)
      call axis()

      do j=1,ny
         call color(mod(j,16))
         do i=1,npts
            x(i)=(i-1.)/(npts-1.)
            y(i)=sin(10.*x(i))*exp(-x(i))+.2*j/ny
            call vecw(x(i),y(i),-1)
         enddo
      enddo

      call pltend()

      end
