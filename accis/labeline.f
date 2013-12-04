C********************************************************************
      subroutine labeline(x,y,npts,label,nlabel)
c Draw a polyline with imbedded label of nlabel characters.
c 9 Aug 92.
c If nlabel eq -99, then set the interval between labels to the value of 
c x in normalized units, and the ipen (up down) to the value of y. 
c If x=0. set it to default, 0.3. If ipen.ne.0 then draw line.
      integer npts,nlabel
      real x(npts),y(npts)
      character*(*) label
      character llstr*10
      integer i,nlab1
      include 'plotcom.h'
c
      real nx,ny,cc,cs,theta,nxs,nys
      real vlen,dx,dy,cx,cy,plen,flen,dlen
      real wx2nx, wy2ny
      integer ncmax
      parameter (ncmax=40)
      real clen(ncmax)
      real linlen,curdist
c clen is the arc length in normalized units of the the ith line
c segment. curdist is the starting fractional part of the ith arc.
      real wstr
      integer j,lstr,iline
      data curdist/0./
      data linlen/0.3/iline/1/
      save linlen,iline

c      write(*,*)npts,label,nlabel
c Speed version with no label.
      if(nlabel.eq.0)then
         call polyline(x,y,npts)
	 return
      elseif(nlabel.eq.-99)then
         if(x(1).eq.0.) then
            linlen=0.3
            iline=1
         else
            linlen=x(1)
            if(y(1).eq.0)then
               iline=0
            else
               iline=1
            endif
            return
         endif
      endif
c Standard version. NLabel non-zero.
c store current character direction.
      cc=chrscos
      cs=chrssin
      lstr=1
      nlab1=min(nlabel,ncmax-1)+1
      call vecn(wx2nx(x(1)),wy2ny(y(1)),0)
      do 4 i=1,nlab1-1
c Fix July 99 for G77. Has to be terminated else wstr can't work.
	 llstr=label(i:i)//char(0)
	 clen(i)=wstr(llstr)
c	 write(*,*)llstr(1:1),' ',clen(i)
    4 continue
      clen(nlab1)=linlen
c Start with line.
      j=nlab1
c Draw line
      do 3 i=2,npts
c We shall bypass vecw and go straight to normal.
	 nx=wx2nx(x(i))
	 ny=wy2ny(y(i))
	 nxs=nx
	 nys=ny
c Lengths of total vector:
	 cx= crsrx
	 cy= crsry
	 dx=nx-cx
	 dy=ny-cy
	 vlen=sqrt(dx*dx+dy*dy)
c Partial length remaining:
	 plen=vlen
c	    if(vlen.eq.0)return
c Distance to end of segment(or character)
    1	 dlen=(clen(j)-curdist)
	 if(plen.gt.dlen)then
c Vector longer than this segment. Write segment and iterate.
	    curdist=0.
	    plen=plen-dlen
	    theta=atan2(dy,dx)
	    chrscos=cos(theta)
	    chrssin=sin(theta)
	    if(j.ne.nlab1)then
	       if(clen(j).gt.0)then
		  llstr(lstr:lstr)=label(j:j)
		  llstr(lstr+1:lstr+1)=char(0)
		  call drcstr(llstr)
		  lstr=1
	       else
c Store zero length characters. May be controls.
		  llstr(lstr:lstr)=label(j:j)
		  lstr=lstr+1
	       endif
	    else
	       flen=dlen/vlen
	       nx= crsrx+dx*flen
	       ny= crsry+dy*flen
               call vecn(nx,ny,iline)
	    endif
	    j=j+1
	    if(j.gt.nlab1)j=1
	    goto 1
	 elseif(j.eq.nlab1)then
c Vector ends before segment. If this is the line segment draw to end
c of vector and quit. Else do nothing.
	    curdist=plen+curdist
	    nx=nxs
	    ny=nys
            call vecn(nx,ny,iline)
	 endif
    3 continue
      chrscos=cc
      chrssin=cs
      end

c*********************************************************************
c      program tlabline
c      integer i,imax,nlab
c      parameter (imax=50)
c      real x(imax),y(imax),ymin,ymax
c      character label*(30)
c      external wx2nx,wy2ny
c
c      ymin=0
c      ymax=0
c      do 1 i=1,imax
c	 x(i)=i
c	 y(i)=sin(3.*i/float(imax))
c	 if(y(i).gt.ymax)ymax=y(i)
c	 if(y(i).lt.ymin)ymin=y(i)
c    1 continue
c    2 write(*,*)'Enter Label (CR) length'
c      read(*,'(a)')label
c      read(*,*)nlab
c      call pltinit(x(1),x(imax),ymin,ymax)
c      call axis
c      call labeline(x,y,imax,label,nlab)
c      call pltend()
c      goto 2
c      end



