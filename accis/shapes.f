! Return a rectangle with rounded corners using nround-point arcs
! The lower left and upper right vertices are in x1,y1; x2,y2
! The radii of the corners are rx, ry.
! Examples: 0,0,1,1,.2,.2,20,xp,yp          rounded square
!         : -1,-1,1,1,1,1,20                circle
!         : -2,-1,2,1,2,1,20                ellipse
!         : -2,-1,2,1,1,0.5,20              rectangle rounded by aspect
!         : -2,-1,2,1,1,1,20                circle ends with x-insert
!         : -2,-1,2,1,0,0,1                 plain rectangle
! Dimensions of output arrays xp,yp must be at least 4*nround+1
      subroutine rectround(x1,y1,x2,y2,rx,ry,nround,xp,yp)
      real x1,y1,x2,y2,rx,ry
      integer nround
      real xp(4*nround+1),yp(4*nround+1)
! The drawing path starts at the top of the left corner arc and proceeds
! anticlockwise.
      real xstep(4),ystep(4)
      data xstep/0.,1.,1.,0./ystep/0.,0.,1.,1./
      j=0
      if(1.999*rx.gt.abs(x2-x1).or.1.999*ry.gt.abs(y2-y1))then
         write(*,*)'Incompatible radii',rx,ry
      endif
      theta=-3.1415926
      do i=1,4
         xc=x1+(x2-x1)*xstep(i)+ rx*(1.-2*xstep(i))
         yc=y1+(y2-y1)*ystep(i)+ ry*(1.-2*ystep(i))
         do nr=1,nround
            j=j+1
            xp(j)=xc+rx*cos(theta)
            yp(j)=yc+ry*sin(theta)
            theta=theta+3.1415926/(2.*nround)
         enddo
      enddo
      j=j+1
      xp(j)=xp(1)
      yp(j)=yp(1)
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
c      integer nr,np
c      parameter  (nr=20,np=4*nr+1)
c      real xp(np),yp(np)
c      
c      x1=-2
c      y1=-1
c      x2=2.
c      y2=1.
c      rx=1.
c      ry=1.
c      call rectround(x1,y1,x2,y2,rx,ry,nr,xp,yp)
c      call pltinit(x1,x2,y1,x2)
c      call polyline(xp,yp,np)
c      call pltend
c
c      end
