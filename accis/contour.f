
c_______________________________________________________________________

      subroutine CONTOUR(z,x,y,nx,ny,cl,nl)

c  Contour a function        z(1:nx,1:ny)
c  On a mesh of dimensions   nx times ny
c  with coordinates          x(nx,ny),y(nx,ny)
c  at levels                 cl(nl)
c  whose number is           nl

c  External routines called 
c  Polyline       to draw line segment(s).
c  pltinit
c  axis

	integer nx,ny
      real xp(5),yp(5),z(nx,ny),x(nx,ny),y(nx,ny),cl(nl)
      real zc(5),zd1,zd2
      integer i,j,i1,i2,j1,k,kl,kk,kb,kp,k0
      integer ic(5),jc(5)

c Indices for the unit cell.
      data ic/0,0,1,1,0/
      data jc/0,1,1,0,0/

c Initialize the plot.
      call pltinit(x(1,1),x(nx,1),y(1,1),y(1,ny))
      call axis()
      zd1=nint((x(nx,1)-x(1,1))/10.)
      zd2=nint((y(1,ny)-y(1,1))/10.)

c For a cell, find if contour goes through its boundary.
c If so, interpolate for the intersection points and draw.
      do 1000 i=1,nx-1
	 do 2000 j=1,ny-1
	    zc(1)=z(i,j)
	    zc(2)=z(i,j+1)
	    zc(3)=z(i+1,j+1)
	    zc(4)=z(i+1,j)
	    zc(5)=zc(1)
	    do 3000 n=1,nl
	       zl=cl(n)
	       kp=0
	       k0=0
	       do 2100 k=1,4
		  zd1=zc(k)-zl
		  zd2=zl-zc(k+1)
c		  if(zd1*zd2)2100,9,10
		  if(zd1*zd2.lt.0.)then
                     goto 2100
                  elseif(zd1*zd2.gt.0.)then
                     goto 10
                  endif
c    9		  if(zd2)2100,11,2100
    9		  if(zd2.ne.0)then 
                     goto 2100
                  else
                     goto 11
                   endif
   11		  if(zd1.eq.0) then
		  k0=k0+1
		  goto 2100
		  endif
   10		  kp=kp+1
		  i1=i+ic(k)
		  i2=i+ic(k+1)
		  j1=j+jc(k)
		  j2=j+jc(k+1)
		  xp(kp)=(x(i1,j1)*zd2+x(i2,j2)*zd1)/(zd1+zd2)
		  yp(kp)=(y(i1,j1)*zd2+y(i2,j2)*zd1)/(zd1+zd2)
 2100	       continue
	       if(kp.gt.1)then
		  if(kp.eq.2) then
c  If this case is standard, less than two double zero segments.
		     if(k0.lt.2) then
			call polyline(xp,yp,kp)

C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c Sections to deal with special cases.

		     elseif (k0.eq.2) then
c  Special case of three nodes having z=zl; two double zero segments.
c  Join the adjacent nodes.
			kl=0
			do 2300 kk=1,5
			   i1=i+ic(kk)
			   j1=j+jc(kk)
c If this is a zero
			   if((z(i1,j1)-zl).eq.0) then
c then we have found the next zero
			      kl=kl+1
			      xp(mod(kl+1,4)+1)=x(i1,j1)
			      yp(mod(kl+1,4)+1)=y(i1,j1)
c else this is a nonzero, and if we have already found a zero
			   elseif(kl.gt.0)then
c then skip and determine the bottom of the polyline.
			      kl=kl+1
			      kb=mod(kl+2,4)+1
			   endif
 2300			continue
			call polyline(xp(kb),yp(kb),3)
			   endif

		     elseif(kp.eq.4) then
c  Case of two line-segments in a cell.
c  Join the pairs of points that are closest.
			if(((xp(1)-xp(2))**2+(yp(1)-yp(2))**2
     &          +(xp(3)-xp(4))**2+(yp(3)-yp(4))**2).le.
     &		((xp(1)-xp(4))**2+(yp(1)-yp(4))**2
     &          +(xp(2)-xp(3))**2+(yp(2)-yp(3))**2)) then
			   call polyline(xp,yp,2)
			   call polyline(xp(3),yp(3),2)
			else
			   xp(5)=xp(1)
			   yp(5)=yp(1)
			   call polyline(xp(4),yp(4),2)
			   call polyline(xp(2),yp(2),2)
			endif

		     else
c  Here if kp=3
c  Special case of two vertices of the cell having z=zl.
			do 2200 k=1,4
			   zd1=zc(k)-zl
			   zd2=zl-zc(k+1)
			   if(zd1*zd2.gt.0) then
c		Project across the cell from the side intersection.
			    i1=i+ic(k)
			    i2=i+ic(k+1)
			    j1=j+jc(k)
			    j2=j+jc(k+1)
			    xp(1)=(x(i1,j1)*zd2+x(i2,j2)*zd1)/(zd1+zd2)
			    yp(1)=(y(i1,j1)*zd2+y(i2,j2)*zd1)/(zd1+zd2)
			    i1=i+ic(mod(k+1,4)+1)
			    i2=i+ic(mod(k+2,4)+1)
			    j1=j+jc(mod(k+1,4)+1)
			    j2=j+jc(mod(k+2,4)+1)
			    xp(2)=(x(i1,j1)*zd1+x(i2,j2)*zd2)/(zd1+zd2)
			    yp(2)=(y(i1,j1)*zd1+y(i2,j2)*zd2)/(zd1+zd2)
			    call polyline(xp,yp,2)
c		Join the opposing points.
			    xp(3)=x(i1,j1)
			    yp(3)=y(i1,j1)
			    xp(4)=x(i2,j2)
			    yp(4)=y(i2,j2)
			    call polyline(xp(3),yp(3),2)
			    goto 3100
			   endif
 2200			continue
 3100			continue
                     endif
c  End of special sections.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		     endif
 3000	    continue
 2000	 continue
 1000 continue
c      call pltend()
      return
      end

