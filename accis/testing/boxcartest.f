      program boxcartest
      integer npts
      parameter (npts=200)
      real x(npts),y(npts),yave(npts)

      do i=1,npts
         x(i)=i
         y(i)=sin(x(i)/5.)
      enddo

      nb=9
      call boxcarave(npts,nb,y,yave)
c      call triangave(npts,nb,y,yave)

      call autoplot(x,y,npts)
      call color(2)
      call polyline(x,yave,npts)
      call pltend()

      call autoplot(x,y,npts)
      call color(3)
      call smoothline(x,y,npts,nb)
      call pltend()

      call autoplot(x,y,npts)
      call color(3)
      call smoothline(x,y,npts,-nb)
      call pltend()


      end

c      include 'boxcarave.f'
