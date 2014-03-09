C Test of driver routines, flushing, and animation.

      integer npts
      parameter (npts=100)
      real x(100),y(100)

      do i=1,npts
         x(i)=i
         y(i)=sin(x(i)/10.)
      enddo

      call autoinit(x,y,npts)

      call axis()

      if(.false.)then
         call polyline(x,y,npts)
      else
         call vecw(x(1),y(1),0)
         do i=2,npts
            call vecw(x(i),y(i),1)
            call vecw(x(i),y(i),0)
            call accisflush()
            call usleep(5000)
         enddo
      endif
      

      call pltend()

      end
