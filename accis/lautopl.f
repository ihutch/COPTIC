c************************************************************************
c	Automatic setup of arrays with logarithmic scales.
      subroutine lautoinit(x,y,n,lx,ly)
      real x(1),y(1)
      integer n
      logical lx,ly
      include 'plotcom.h'
      real xmin,xmax,xc,ymin,ymax,yc,xfac,xdelta
      integer i,nxfac
      xmax=x(1)
      xmin=xmax
      ymax=y(1)
      ymin=ymax
      do 1 i=2,n
	 xc=x(i)
	 yc=y(i)
	 if(xc.lt.xmin)xmin=xc
	 if(xc.gt.xmax)xmax=xc
	 if(yc.lt.ymin)ymin=yc
	 if(yc.gt.ymax)ymax=yc
    1 continue
      if(.not.lx)then
	 call fitrange(xmin,xmax,ticnum,nxfac,xfac,xdelta,xmin,xmax)
      else
	 xmin=10.**(nint(log10(xmin)-0.49999))
	 xmax=10.**(nint(log10(xmax)+0.49999))
      endif
      if(.not.ly)then
	 call fitrange(ymin,ymax,ticnum,nxfac,xfac,xdelta,ymin,ymax)
      else
	 ymin=10.**(nint(log10(ymin)-0.49999))
	 ymax=10.**(nint(log10(ymax)+0.49999))
      endif
      call pltinit(xmin,xmax,ymin,ymax)
      call scalewn(xmin,xmax,ymin,ymax,lx,ly)
      return
      end
c********************************************************************
c         Automatic plotting of Arrays*/
      subroutine lautoplot(x,y,n,lx,ly)
      real x(1),y(1)
      integer n
      logical lx,ly
      call lautoinit(x,y,n,lx,ly)
      call axis()
      call polyline(x,y,n)
      return
      end
c********************************************************************
c         Automatic symbol plotting of Arrays*/
      subroutine lautomark(x,y,n,lx,ly,isym)
      real x(1),y(1)
      integer n,isym
      logical lx,ly
      call lautoinit(x,y,n,lx,ly)
      call axis()
      if(isym.gt.0)call polymark(x,y,n,isym)
      return
      end
