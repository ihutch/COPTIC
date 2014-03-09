      program box
      integer length,long
      parameter (length=10)
      real x(length),y(length)
      data x/-1.,1.,1.,-1.,-1.,-.5,0.,.5,0.,-.5/
      data y/-1.,-1.,1.,1.,-1.,0.,.5,0.,-.5,0./
      call pfset(2)
      call pltinit(-1.,1.,-1.,1.)
      call color(5)
      call polyline(x,y,length)
      call pathfill()
c      do i=1,4
c         x(i)=0.5*x(i)
c         y(i)=0.5*y(i)
c      enddo
c      call color(2)
c      call polyline(x,y,length)
c      call pathfill()
      call pltend()
      end
