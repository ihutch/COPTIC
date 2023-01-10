c Illustration/Test of the use of manautoinit for mixed initialization
c of plots between autoscaling and specified limits.
      integer length
      parameter (length=100)
      real x(length),y(length),err(length),ym(length)
      integer i,irows,icolumns,itype

c Make test arrays.
      do 2 i=1,length
         x(i)=float(i)*1.
         y(i)=1.*(1.4+sin(0.2*x(i)))
         err(i)=0.4*sin(x(i)*sin(x(i)))
         ym(i)=y(i)-0.5*err(i)
    2 continue

      isw=0
      call manautoinit(x,y,length,isw,sxmin,sxmax,symin,symax)
      call axis()
      call polyline(x,y,length)
      call pltend()

      isw=2
      sxmax=50.
      call manautoinit(x,y,length,isw,sxmin,sxmax,symin,symax)
      call axis()
      call polyline(x,y,length)
      call pltend()

      isw=5
      sxmin=-10.
      symin=-1.
      call manautoinit(x,y,length,isw,sxmin,sxmax,symin,symax)
      call axis()
      call polyline(x,y,length)
      call pltend()


      end
