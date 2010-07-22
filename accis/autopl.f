c********************************************************************
c Automatic plotting of an array versus its index
      subroutine yautoplot(y,n)
      include 'plotcom.h'
      integer n
      real y(n)
      real xmin,xmax,ymin,ymax
      real xfac,xdelta
      integer nxfac
      call minmax(y,n,ymin,ymax)
c      call fitrange(1.,float(n),ticnum,nxfac,xfac,xdelta,xmin,xmax)
      call fitrange(ymin,ymax,ticnum,nxfac,xfac,xdelta,ymin,ymax)
c      call pltinit(xmin,xmax,ymin,ymax)
      call pltinit(1.,float(n),ymin,ymax)
      call axis()
      call ypolyline(y,n)
      end
c********************************************************************
c crude plotting of y versus its index. No dashed line capability.
      subroutine ypolyline(y,n)
      integer n,i
      real y(n)
      call vecw(1.,y(1),0)
      do 1 i=2,n
	 call vecw(float(i),y(i),1)
    1 continue
      end
c********************************************************************
c         Automatic plotting of Arrays*/
      subroutine autoplot(x, y, n)
      real x(n),y(n)
      integer n
      call autoinit(x,y,n)
      call axis()
      call polyline(x,y,n)
      return
      end
c******************************************************************
      subroutine autoinit(x,y,n)
      real x(n),y(n)
      integer n
      include 'plotcom.h'
      real xmin,xmax,ymin,ymax
      real xfac,xdelta
      integer nxfac
      call minmax(x,n,xmin,xmax)
      call minmax(y,n,ymin,ymax)
      call fitrange(xmin,xmax,ticnum,nxfac,xfac,xdelta,xmin,xmax)
      call fitrange(ymin,ymax,ticnum,nxfac,xfac,xdelta,ymin,ymax)
      call pltinit(xmin,xmax,ymin,ymax)
      end
c********************************************************************
c    Automatic symbol plotting of Arrays*/
      subroutine automark(x, y, n, isym)
      real x(n),y(n)
      integer n,isym
      call autoinit(x,y,n)
      call axis()
      if(isym.gt.0)call polymark(x,y,n,isym)
      return
      end
c******************************************************************
      subroutine auto3init(x,y,z,n)
      real x(n),y(n),z(n)
      integer n
      save
      include 'plotcom.h'
      real xmin(3),xmax(3),dx(3)
      real xfac,xdelta
      integer nxfac,i
      real x2,y2,z2
      call minmax(x,n,xmin(1),xmax(1))
      call minmax(y,n,xmin(2),xmax(2))
      call minmax(z,n,xmin(3),xmax(3))
      do i=1,3
         call fitrange(xmin(i),xmax(i),ticnum,nxfac,
     $        xfac,xdelta,xmin(i),xmax(i))
      enddo
C We assume we want to preserve metric of the geometry, so:
      imax=1
      do i=1,3
         dx(i)=xmax(i)-xmin(i)
         if(dx(i).ge.dx(imax)) imax=i
      enddo
      do i=1,3
         if(i.ne.imax) then
            xmid=(xmin(i)+xmax(i))*0.5
            xmin(i)=xmid+dx(imax)*(xmin(i)-xmid)/dx(i)
            xmax(i)=xmid+dx(imax)*(xmax(i)-xmid)/dx(i)
         endif
      enddo
c Need a way to determine eye position. Use defaults.
      call geteye(x2,y2,z2)
      call pltinit(0.,1.,0.,1.)
      call SCALE3(xmin(1),xmax(1),xmin(2),xmax(2),xmin(3),xmax(3))
      call trn32(0.,0.,0.,x2,y2,z2,1)
c determine icorner:
      xb=0
      yb=0
      zb=0
      if(x2.ge.0)xb=1
      if(y2.ge.0)yb=1
      if(z2.ge.0)zb=1
      icorner= (2*zb-1)*( (1 +3*yb) + (1 - 2*yb)*xb )
c Then decide on optimal axis labels. Not done now.
c      write(*,*)'icorner',icorner
      icorner=icorner
c      write(*,*) 'icorner=',icorner
c      call cubed(-mod(icorner+2,4))
c or draw full cube      
      call cubed(0)
c I don't know why the calling of cubed changes the behaviour of axproj.
      call axproj(icorner)
c      call axis()
      end
