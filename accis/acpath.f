c**************************************************************************
c Dummy acpathdoc must be replaced by explicitly linked version if
c path documentation is actually desired.
      subroutine acpathdoc(imax,xc,yc,cv)
      real xc(imax),yc(imax),cv
      integer imax
c API: after the construction of the whole contour path
c     acpathdoc is called with arguments
c     cv the contour value 
c     imax (in) the length of the xc, yc arrays if positive
c          or if negative an indicator
c     xc, yc, the position of points in path arrays.
c Typical usage might be to declare the new acpathdoc to have a common
c block into which xc,yc are entered to provide the contour as a polyline
c do i=1,imax; xpath(i)=xc(i); ypath(i)=yc(i);  enddo 
      end
