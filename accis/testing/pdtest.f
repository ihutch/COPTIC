c     Test polydraw
      real x(10),y(10)
      external accircle,acx,acplus,acast,acgen,actrid

c Inverted triangle
      x(1)=0.
      y(1)=-.33*sqrt(3.)
      x(2)=0.5
      y(2)=.17*sqrt(3.)
      x(3)=-x(2)
      y(3)=y(2)
      x(4)=x(1)
      y(4)=y(1)
      call acgset(x,y,4,1)

      call pfset(2)
      do i=1,10
         x(i)=i/10.
         y(i)=x(i)**2
      enddo
c      call AXREGION(.31,.61,.1,.5) 
      call pltinit(0.,1.,0.,1.)
      call axis()
c      call polydraw(x,y,10,acx)
      call polydraw(x,y,10,acast)
      call polydraw(x,y,10,actrid)
      call polydraw(y,x,10,acgen)
      call polydraw(y,x,10,acplus)
      call charsize(.06,.06)
      call polydraw(x,y,10,accircle)
      call charsize(0.,0.)
      call polydraw(x,x,10,acx)
      call polymark(y,y,10,128+ichar('*'))
      call pltend()

      end
