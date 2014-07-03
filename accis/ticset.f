      subroutine ticset(xlen,ylen,xoff,yoff,nxw,nyw,nxp,nyp)
c Set various attribute of axis tics and labels.
      real xlen,ylen,xoff,yoff
      integer nxw,nxp,nyw,nyp
      include 'plotcom.h'
      if(xlen.eq.0..and.ylen.eq.0..and.xoff.eq.0..and.yoff.eq.0.)then
c Default tics
         xticlen=0.015
         yticlen=0.015
         xticoff=-.03
         yticoff=-.02
      else
c Set tics
         xticlen=xlen
         yticlen=ylen
         xticoff=xoff
         yticoff=yoff
      endif
      if((nxw.eq.0).and.(nxp.eq.0).and.(nyw.eq.0).and.(nyp.eq.0))then
c Default labels
         nxlabw=4
         nxlabp=1
         nylabw=4
         nylabp=1
      else
c Set labels
         nxlabw=nxw
         nxlabp=nxp
         nylabw=nyw
         nylabp=nyp
      endif
      end
c*****************************************************************/
      subroutine ticrev()
c Reverse the side of the axes that tics are drawn.
      include 'plotcom.h'
      xticoff=-xticoff
      yticoff=-yticoff
      xticlen=-xticlen
      yticlen=-yticlen
      end
c***********************************************************************
      subroutine axptset(px,py)
c Set the axis point in fractions of the axis region.
c Default is 0.,0.
      real px,py
      include 'plotcom.h'
      naxpt=px*(naxmax-naxmin) + naxmin
      naypt=py*(naymax-naymin) + naymin
      end
c***********************************************************************
      subroutine axis2()
      include 'plotcom.h'
      call axptset(1.,1.)
      call ticlabtog()
      call ticrev()
      call axis()
      call ticrev()
      call ticlabtog()
      call axptset(0.,0.)
      end
c***********************************************************************
      subroutine ticlabtog()
c Toggle tic labels on and off.
      include 'plotcom.h'
      integer ticwidth 
      save ticwidth
      if(nxlabw.le.0)then
         if(ticwidth.gt.0)then
            nxlabw=ticwidth
         else
            nxlabw=4
         endif
         ticwidth=0
      else
         ticwidth=nxlabw
         nxlabw=0
      endif
      end
c***********************************************************************
      subroutine ticnumset(numtics)
c Toggle tic labels on and off.
      include 'plotcom.h'
      integer numtics
      if(numtics .lt. 0) then 
         ticnum=6
      else
         ticnum=numtics
      endif
      end
