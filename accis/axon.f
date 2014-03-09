      subroutine axon(x,y,z,iLx,nx,ny)
c Plot axonometrically z(x,y) dim (iLx\nx,ny).
      integer nx,ny,ud,iLx
      real z(iLx,ny),x(nx),y(ny),zmin,zmax,xt,yt,zt
      real wx2nx,wy2ny,xnreg,ynreg,ydelta,y1st,ylast,ypscl,xrm
      real xfrac,yfrac,xni,xna,yni,yna
      integer i,j,nyfac
      include 'plotcom.h'
c Read eye.
      call geteye(xfrac,yfrac,zmin)
      call hidinit(0.,1.)
      if(y(ny)-y(1).eq.0.) stop 'AXON: needs y(ny).ne.y(1).'
c Set axon scaling
      call trn32(xt,yt,zt,xfrac/(y(ny)-y(1)),yt,yfrac/(y(ny)-y(1)),2)
      ynreg=.7-.1
      xrm=.1-min(0.,xfrac)
      xnreg=.9-max(0.,xfrac)-xrm
      call axregion(xrm,xrm+xnreg,.1,.1+ynreg)
      call minmax2(z,iLx,nx,ny,zmin,zmax)
      if(abs(zmin/(zmax-zmin)).lt.0.1)zmin=0.
      zmax=zmin+(1.+yfrac)*(zmax-zmin)
      call pltinit(x(1),x(nx),zmin,zmax)
      call xaxis(0.,0.)
      if(xfrac.lt.0) then
	 naxpt=naxmax
	 call ticrev
      endif
      call yaxis(0.,0.)
      if(xfrac.lt.0)then
	 naxpt=naxmin
	 call ticrev
      endif
      call fitrange(y(1),y(ny),7,nyfac,ypscl,ydelta,y1st,ylast)
      if(nyfac.le.2.and.nyfac.ge.0)then
	 nyfac=0.
	 ypscl=1.
      endif
      if(xfrac.gt.0)then
	 i=nx
	 xni= wx2nx(x(i))
         xna= wx2nx(x(i)+xfrac*(x(nx)-x(1))/xnreg)
	 yni= wy2ny(zmin)
         yna= wy2ny(zmin+yfrac*(zmax-zmin)/ynreg)
	 zmin=y(1)/ypscl
	 zmax=y(ny)/ypscl
      else
	 i=1
	 xna= wx2nx(x(i))
         xni= wx2nx(x(i)+xfrac*(x(nx)-x(1))/xnreg)
	 yna= wy2ny(zmin)
         yni= wy2ny(zmin+yfrac*(zmax-zmin)/ynreg)
	 zmin=y(ny)/ypscl
	 zmax=y(1)/ypscl
      endif
      call gaxis(zmin,zmax,nyfac,y1st/ypscl,ydelta/ypscl,
     $	   xni,xna,yni,yna,.true.,.false.)
      do 3 j=1,ny
	 do 4 i=1,nx
	    ud=i-1
	    call trn32(wx2nx(x(i)),y(j)-y(1),wy2ny(z(i,j))
     $		 ,xt,yt,zt,0)
	    call hidvecn(xt,yt,ud)
    4	 continue
    3 continue
      call pltend
      end


