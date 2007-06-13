c Plot convergence
      integer ndata
      parameter (ndata=5)
      integer ni(ndata)
      real ermax(ndata),ersd(ndata),dx(ndata)
      real xl(2),yl(2)
      
      data ni/10,20,40,80,100/
      data ermax/0.01724,0.004208,0.001210,0.0003077,0.0002338/
      data ersd/0.005111,0.0013607,0.0003119,9.302E-05,7.775E-05/

      do i=1,ndata
         dx(i)=1./(ni(i)-1.)
      enddo

      xl(1)=.01
      xl(2)=.1
      do i=1,2
         yl(i)=xl(i)**2
      enddo

      call pltinit(0.,1.,0.,1.)
      call scalewn(.01,.1,.0001,.01,.true.,.true.)
      call axis()
      call polyline(dx,ermax,ndata)
      call polymark(dx,ermax,ndata,1)
c      call lautoplot(dx,ermax,ndata,.true.,.true.)
      call axlabels('dx','Error')
      call color(3)
      call polyline(dx,ersd,ndata)
      call polymark(dx,ersd,ndata,2)
      call color(4)
      call polyline(xl,yl,2)
      call pltend()

      end
