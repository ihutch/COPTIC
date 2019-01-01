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
c dashlen is the arc length in normalized units of the the ith line
c segment. dashdist is the starting fractional part of the ith arc.
c Segments alternate pen down, pen up.
      logical ldash
      real dashlen,dashdist
      integer MASKNO,dashmask,jmask
      parameter (MASKNO=4)
      dimension dashmask(MASKNO),dashlen(MASKNO)
      common/dashline/ldash,dashlen,dashdist,dashmask,jmask
      integer cud

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
         llstr=label(i:i)//char(0)
         clen(i)=wstr(llstr)
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
c Partial length remaining to end of vector:
         plen=vlen
c           if(vlen.eq.0)return
c Distance to end of segment(or character)
    1    dlen=(clen(j)-curdist)
         if(plen.gt.dlen)then
c Vector longer than this label segment. Write segment and iterate.
            curdist=0.
            plen=plen-dlen
            theta=atan2(dy,dx)
            chrscos=cos(theta)
            chrssin=sin(theta)
            if(j.ne.nlab1)then
c Character writing:
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
               if(ldash)then
c Dashline distance to end of segment
c                  write(*,*)jmask,dashlen(jmask),dashdist,plen
 2                dlen=(dashlen(jmask)-dashdist)
                  if(plen.gt.dlen)then
c Vector longer than this dash segment. Draw segment and iterate.
                     plen=plen-dlen
                     flen=dlen/vlen
                     nx= crsrx+dx*flen
                     ny= crsry+dy*flen
                     cud=dashmask(jmask)
                     call vecn(nx,ny,cud)
                     dashdist=0
                     jmask=mod(jmask,MASKNO)+1
                     goto 2
                  else
c Vector ends before dash segment. Draw to end of vector and quit.
                     dashdist=plen+dashdist
c                     write(*,*)'dashdist',dashdist,plen,dlen
                     nx=crsrx+dx*plen/vlen
                     ny=crsry+dy*plen/vlen
                     cud=dashmask(jmask)
                     call vecn(nx,ny,cud)
                  endif
               else
c Replaces this undashed line:
                  flen=dlen/vlen
                  nx= crsrx+dx*flen
                  ny= crsry+dy*flen
                  call vecn(nx,ny,iline)
               endif
            endif
            j=j+1
            if(j.gt.nlab1)j=1
            goto 1

         elseif(j.eq.nlab1)then
c Vector ends before segment. If this is the line segment draw to end
c of vector and quit. Else do nothing.
            if(ldash)then
               curdist=plen+curdist
c Dashline distance to end of segment
 5             dlen=(dashlen(jmask)-dashdist)
               if(plen.gt.dlen)then
c Vector longer than this segment. Draw segment and iterate.
                  plen=plen-dlen
                  flen=dlen/vlen
                  nx= crsrx+dx*flen
                  ny= crsry+dy*flen
                  cud=dashmask(jmask)
                  call vecn(nx,ny,cud)
                  dashdist=0
                  jmask=mod(jmask,MASKNO)+1
                  goto 5
               else
c Vector ends before segment. Draw to end of vector and quit.
                  dashdist=plen+dashdist
                  nx=nxs
                  ny=nys
                  cud=dashmask(jmask)
                  call vecn(nx,ny,cud)
               endif
            else
c Replaces this:
               curdist=plen+curdist
               nx=nxs
               ny=nys
               call vecn(nx,ny,iline)
            endif
         endif
    3 continue
      chrscos=cc
      chrssin=cs
      end

c*********************************************************************
c$$$      program tlabline
c$$$      integer i,imax,nlab
c$$$      parameter (imax=20)
c$$$      real x(imax),y(imax),ymin,ymax
c$$$      character label*(30)
c$$$      external wx2nx,wy2ny
c$$$
c$$$      ymin=0
c$$$      ymax=0
c$$$      do 1 i=1,imax
c$$$        x(i)=i
c$$$        y(i)=sin(3.*i/float(imax))
c$$$        if(y(i).gt.ymax)ymax=y(i)
c$$$        if(y(i).lt.ymin)ymin=y(i)
c$$$    1 continue
c$$$    2 write(*,*)'Enter Label (CR)'
c$$$      read(*,'(a)')label
c$$$      nlab=lentrim(label)
c$$$      call pltinit(x(1),x(imax),ymin,ymax)
c$$$      call dashset(5)
c$$$      call axis
c$$$      call labeline(x,y,imax,label,nlab)
c$$$      call pltend()
c$$$      end



