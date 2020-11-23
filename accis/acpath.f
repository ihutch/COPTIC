c**************************************************************************
c Dummy acpathdoc must be replaced by explicitly linked version if
c path documentation is actually desired. This version turns it off.
      subroutine acpathdoc(imax,xc,yc,cv)
      real xc(imax),yc(imax),cv
      integer imax
c API: after the construction of the whole contour path
c     acpathdoc is called with arguments
c     cv the contour value 
c     imax (in) the length of the xc, yc arrays if positive
c          or if zero an indicator of prior problems
c     xc, yc, the position of points in path arrays, world coords.
c Typical usage might be to declare the new acpathdoc to have a common
c block into which xc,yc are entered to provide the contour as a polyline
c This is a convenient common block to use in alternate acpathdocs:
      integer nacpmax,nacp
      parameter (nacpmax=1000)
      real xacpath(nacpmax),yacpath(nacpmax),accv
      common /acpathcom/nacp,xacpath,yacpath,accv
c Change the following sign to positive to activate setting.
c By default imax is always positive so -nacpmax turns setting off. 
c Adjust nacpmax to what you want then use the common block elsewhere.
      if(imax.le.-nacpmax)then
         do i=1,imax
            xacpath(i)=xc(i)
            yacpath(i)=yc(i)
         enddo
         accv=cv
         nacp=imax
      endif
      end
