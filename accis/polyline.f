C********************************************************************
      subroutine polyline(x,y,npts)
c Dashed line version.
      integer npts
      real x(npts),y(npts)
      integer i
      include 'plotcom.h'
c Dashed line code
      real nx,ny
      real vlen,dx,dy,cx,cy,plen,flen,dlen
      real wx2nx, wy2ny
      integer cud

c dashlen is the arc length in normalized units of the the ith line
c segment. dashdist is the starting fractional part of the ith arc.
c Segments alternate pen down, pen up.
      logical ldash
      real dashlen,dashdist
      integer MASKNO,dashmask,jmask
      parameter (MASKNO=4)
      dimension dashmask(MASKNO),dashlen(MASKNO)
      common/dashline/ldash,dashlen,dashdist,dashmask,jmask

      if(npts.le.0) return
      call vecw(x(1),y(1),0)
      do 3 i=2,npts
         if(.not.ldash) then
            call vecw(x(i),y(i),1)
         else
c We shall bypass vecw and go straight to normal.
            nx=wx2nx(x(i))
            ny=wy2ny(y(i))
c Lengths of total vector:
            cx=nx
            cy=ny
            dx=nx- crsrx
            dy=ny- crsry
            vlen=sqrt(DX*DX+DY*DY)
c Partial length remaining:
            plen=vlen
            if(vlen.eq.0)return
c Distance to end of segment
    1       dlen=(dashlen(jmask)-dashdist)
            if(plen.gt.dlen)then
c Vector longer than this segment. Draw segment and iterate.
               dashdist=0
               plen=plen-dlen
               flen=dlen/vlen
               nx= crsrx+dx*flen
               ny= crsry+dy*flen
               cud=dashmask(jmask)
c              call optvecn(nx,ny,cud)
               call vecn(nx,ny,cud)
               jmask=mod(jmask,MASKNO)+1
               goto 1
            else
c Vector ends before segment. Draw to end of vector and quit.
               dashdist=plen+dashdist
               nx=cx
               ny=cy
               cud=dashmask(jmask)
c              call optvecn(nx,ny,cud)
               call vecn(nx,ny,cud)
            endif
         endif
    3 continue
      end
C********************************************************************
      subroutine dashset(i)
c Set some line styles. 0 solid (dashing off).
      integer i
      external dashdata

      logical ldash
      real dashlen,dashdist
      integer MASKNO,dashmask,jmask
      parameter (MASKNO=4)
      dimension dashmask(MASKNO),dashlen(MASKNO)
      common/dashline/ldash,dashlen,dashdist,dashmask,jmask
      save
      if(i.eq.0)then
         ldash=.false.
         return
      elseif(i.eq.1)then
c Long dashes.
         dashlen(1)=.03
         dashlen(2)=.03
         dashlen(3)=.03
         dashlen(4)=.03
      elseif(i.eq.2)then
c Short dashes.
         dashlen(1)=.01
         dashlen(2)=.01
         dashlen(3)=.01
         dashlen(4)=.01
      elseif(i.eq.3)then
c Long/Short.
         dashlen(1)=.03
         dashlen(2)=.01
         dashlen(3)=.01
         dashlen(4)=.01
      elseif(i.eq.4)then
c 'Dots'.
         dashlen(1)=.002
         dashlen(2)=.01
         dashlen(3)=.002
         dashlen(4)=.01
      elseif(i.eq.5)then
c 'Medium/dot'.
         dashlen(1)=.02
         dashlen(2)=.01
         dashlen(3)=.002
         dashlen(4)=.01
      elseif(i.eq.6)then
c 'Long Dashes short breaks'.
         dashlen(1)=.03
         dashlen(2)=.01
         dashlen(3)=.03
         dashlen(4)=.01
      elseif(i.eq.7)then
c 'Medium Dashes short breaks'.
         dashlen(1)=.02
         dashlen(2)=.005
         dashlen(3)=.02
         dashlen(4)=.005
      elseif(i.eq.8)then
c 'Short Dashes shorter breaks'.
         dashlen(1)=.01
         dashlen(2)=.005
         dashlen(3)=.01
         dashlen(4)=.005
      elseif(i.eq.9)then
c 'Dot short'.
         dashlen(1)=.002
         dashlen(2)=.005
         dashlen(3)=.01
         dashlen(4)=.01
      elseif(i.eq.10)then
c 'Long/Short breaks'.
         dashlen(1)=.02
         dashlen(2)=.005
         dashlen(3)=.02
         dashlen(4)=.02
      elseif(i.eq.11)then
c 'Long/Medium long breaks'.
         dashlen(1)=.03
         dashlen(2)=.02
         dashlen(3)=.015
         dashlen(4)=.02
      elseif(i.eq.12)then
c 'Long with tiny breaks'.
         dashlen(1)=.03
         dashlen(2)=.005
         dashlen(3)=.03
         dashlen(4)=.005
      elseif(i.eq.13)then
c 'Dots with long/short breaks'.
         dashlen(1)=.003
         dashlen(2)=.005
         dashlen(3)=.003
         dashlen(4)=.01
      elseif(i.eq.14)then
c 'Short, dot '.
         dashlen(1)=.01
         dashlen(2)=.01
         dashlen(3)=.003
         dashlen(4)=.01
      elseif(i.eq.15)then
c 'Short short break'.
         dashlen(1)=.01
         dashlen(2)=.006
         dashlen(3)=.01
         dashlen(4)=.02
      endif
      dashdist=1.e-6
      jmask=1
      ldash=.true.
      return
      end
C********************************************************************
      subroutine polymark(x,y,nx,nmark)
      integer nx,nmark,i
      character*4 mark
      real x(*),y(*),xp,yp,wx2nx,wy2ny
      include 'plotcom.h'
      ipf=pfPS
c      if(nmark.eq.3)then
c         mark='!A3'//char(0)
c  That does not quite align well. Better not to mix thinkgs up.
      if(nmark.lt.10)then
         pfPS=0
         mark=char(nmark+176)//char(0)
      elseif(nmark.eq.10)then
         mark='+'//char(0)
      elseif(nmark.eq.11)then
         mark='!AX'//char(0)
      elseif(nmark.eq.12)then
         mark='!A*'//char(0)
      elseif(nmark.eq.13)then
         mark='!A-'//char(0)
      elseif(nmark.eq.15)then
         mark='!A'//char(48)//char(0)
      else
         mark=char(nmark)//char(0)
      endif
      do 1 i=1,nx
         xp=wx2nx(x(i))
         yp=wy2ny(y(i))
      if(nmark.eq.14)then
         call actrid(xp,yp)
      else
         call jdrwstr(xp,yp,mark,0.)
      endif
    1 continue
      pfPS=ipf
      end
C********************************************************************
c Plot error bars from y to y+err.
      subroutine polyerr(x,y,err,nx)
      integer nx,i
      real x(*),y(*),err(*)
      do 1 i=1,nx
         call vecw(x(i),y(i),0)
         call vecw(x(i),y(i)+err(i),1)
    1 continue
      end
C********************************************************************
c Plot error bars from y-ym*err to y+yp*err.
      subroutine polyerrs(x,y,err,nx,yp,ym)
      integer nx,i
      real x(*),y(*),err(*)
      real yp,ym
      include 'plotcom.h'
      real xn,yn
      do 1 i=1,nx
         xn=wx2nx(x(i))
         yn=wy2ny(y(i)-ym*err(i))
         call vecn(xn-.2*chrswdth,yn,0)
         call vecn(xn+.2*chrswdth,yn,1)
         yn=wy2ny(y(i)+yp*err(i))
         call vecn(xn-.2*chrswdth,yn,0)
         call vecn(xn+.2*chrswdth,yn,1)
         call vecw(x(i),y(i)-ym*err(i),0)
         call vecw(x(i),y(i)+yp*err(i),1)
    1 continue
      end
c******************************************************************
      subroutine legendline(xg,yg,nsym,string)
      real xg,yg
      integer nsym
      character*(*) string
c Draw a legendline
c At fractional position in plot box: xg, yg 
c Use symbol nsym if positive and <=256. 
c If nsym<0 use both line and symbol. If nsym=0 put only line.
c If nsym>256 use nothing but put the string at the usual place.
      include 'plotcom.h'

      real xd(2),yd(2)
      xleglen=.07
      xn=wx2nx(wxmin)*(1.-xg)+wx2nx(wxmax)*xg
      yn=wy2ny(wymin)*(1.-yg)+wy2ny(wymax)*yg
      xd(1)=xn2xw(xn)
      xd(2)=xn2xw(xn+xleglen)
      yd(1)=yn2yw(yn)
      yd(2)=yn2yw(yn)
      if(nsym.eq.0)then
         call polyline(xd,yd,2)
      elseif(nsym.eq.257)then
         call vecn(xn+xleglen,yn,0)
      elseif(nsym.eq.258)then
         call vecn(xn,yn,0)
      elseif(nsym.gt.0)then
         call polymark(xd(2),yd(2),1,nsym)
      elseif(nsym.lt.0)then
         call polymark(xd,yd,2,abs(nsym))
         call polyline(xd,yd,2)
      endif
      call drcstr(string)
c      write(*,*)xg,yg,xn,yn,xd(1),yd(1),nsym,string
      end


C********************************************************************
      subroutine polydraw(x,y,nx,drawfn)
      integer nx,i
      real x(*),y(*),xp,yp,wx2nx,wy2ny
      external drawfn
      do 1 i=1,nx
         xp=wx2nx(x(i))
         yp=wy2ny(y(i))
         call drawfn(xp,yp)
    1 continue
      end
c********************************************************************
      subroutine accircle(x,y)
      include 'plotcom.h'
      real x,y
      integer nangle
      parameter (nangle=40)
      iud=0
      do i=1,nangle
         angle=2.*3.14159*(i-1)/(nangle-1)
         xp=x+0.5*chrswdth*cos(angle)
         yp=y+0.5*chrshght*sin(angle)
         call vecn(xp,yp,iud)
         iud=1
      enddo
      end
c********************************************************************
      subroutine acx(x,y)
      include 'plotcom.h'
      real x,y
         xp=x+0.5*chrswdth
         yp=y+0.5*chrshght
         call vecn(xp,yp,0)
         xp=x-0.5*chrswdth
         yp=y-0.5*chrshght
         call vecn(xp,yp,1)
         xp=x+0.5*chrswdth
         yp=y-0.5*chrshght
         call vecn(xp,yp,0)
         xp=x-0.5*chrswdth
         yp=y+0.5*chrshght
         call vecn(xp,yp,1)
      end
c********************************************************************
      subroutine acplus(x,y)
      include 'plotcom.h'
      real x,y
         xp=x+0.5*chrswdth
         yp=y
         call vecn(xp,yp,0)
         xp=x-0.5*chrswdth
         yp=y
         call vecn(xp,yp,1)
         xp=x
         yp=y-0.5*chrshght
         call vecn(xp,yp,0)
         xp=x
         yp=y+0.5*chrshght
         call vecn(xp,yp,1)
      end
c********************************************************************
      subroutine acast(x,y)
      include 'plotcom.h'
      real x,y
      wtemp=chrswdth
      htemp=chrshght
      chrswdth=wtemp*0.707
      chrshght=htemp*0.707
      call acx(x,y)
      chrswdth=wtemp
      chrshght=htemp
      call acplus(x,y)
      end

c********************************************************************
      subroutine actrid(x,y)
      include 'plotcom.h'
      real x,y
c Inverted triangle
      call vecn(x,y-.33*sqrt(3.)*chrshght,0)
      call vecn(x+0.5*chrswdth,y+.17*sqrt(3.)*chrshght,1)
      call vecn(x-0.5*chrswdth,y+.17*sqrt(3.)*chrshght,1)
      call vecn(x,y-.33*sqrt(3.)*chrshght,1)
c      call pathfill()
      end
c********************************************************************
c General mechanism: Call acgset first then pass acgen.
      subroutine acgen(x,y)
      include 'plotcom.h'
      real x,y
c acgen data:
      integer nacg,nacgmax
      parameter (nacgmax=100)
      real xacg(nacgmax),yacg(nacgmax)
      common /acgcom/nacg,xacg,yacg,iacfill

      iud=0
      do i=1,nacg
         call vecn(x+xacg(i)*chrswdth,y+yacg(i)*chrshght,iud)
         iud=1
      enddo
      if(iacfill.ne.0)call pathfill()
      end
c********************************************************************
c Set the general accis symbol to the data in x,y, length n, scale chrs.
c if ifill ne 0 then the symbol is filled.
      subroutine acgset(x,y,n,ifill)
      real x(*),y(*)
      integer n
c acgen data:
      integer nacg,nacgmax
      parameter (nacgmax=100)
      real xacg(nacgmax),yacg(nacgmax)
      common /acgcom/nacg,xacg,yacg,iacfill
      nacg=min(n,nacgmax)
      iacfill=ifill
      do i=1,nacg
         xacg(i)=x(i)
         yacg(i)=y(i)
      enddo
      end

C********************************************************************
      subroutine poly3line(x,y,z,npts)
c Dashed line version.
      integer npts
      real x(npts),y(npts),z(npts)
      integer i
      include 'plotcom.h'
      include 'world3.h'
c Dashed line code
      real nx,ny,nz
      real vlen,dx,dy,cx,cy,plen,flen,dlen
      integer cud

c dashlen is the arc length in normalized units of the the ith line
c segment. dashdist is the starting fractional part of the ith arc.
c Segments alternate pen down, pen up.
      logical ldash
      real dashlen,dashdist
      integer MASKNO,dashmask,jmask
      parameter (MASKNO=4)
      dimension dashmask(MASKNO),dashlen(MASKNO)
      common/dashline/ldash,dashlen,dashdist,dashmask,jmask

      if(npts.le.0) return
      call vec3w(x(1),y(1),z(1),0)
      do 3 i=2,npts
         if(.not.ldash) then
            call vec3w(x(i),y(i),z(i),1)
c            call wxyz2nxyz(x(i),y(i),z(i),nx,ny,nz)
c            call trn32(nx,ny,nz,x2,y2,z2,0)
c            call optvecn(x2+xcbc2,y2+ycbc2,1)
         else
c We shall bypass vecw and go straight to normal.
            call wxyz2nxyz(x(i),y(i),z(i),nx,ny,nz)
            call trn32(nx,ny,nz,x2,y2,z2,0)
c Lengths of total vector:
            cx=x2+xcbc2
            cy=y2+ycbc2
            dx=cx- crsrx
            dy=cy- crsry
            vlen=sqrt(dx*dx+dy*dy)
c Partial length remaining:
            plen=vlen
            if(vlen.eq.0)return
c Distance to end of segment
    1       dlen=(dashlen(jmask)-dashdist)
            if(plen.gt.dlen)then
c Vector longer than this segment. Draw segment and iterate.
               dashdist=0
               plen=plen-dlen
               flen=dlen/vlen
               nx= crsrx+dx*flen
               ny= crsry+dy*flen
               cud=dashmask(jmask)
               call optvecn(nx,ny,cud)
c              call vecn(nx,ny,cud)
               jmask=mod(jmask,MASKNO)+1
               goto 1
            else
c Vector ends before segment. Draw to end of vector and quit.
               dashdist=plen+dashdist
               nx=cx
               ny=cy
               cud=dashmask(jmask)
               call optvecn(nx,ny,cud)
c              call vecn(nx,ny,cud)
            endif
         endif
    3 continue
      end
C********************************************************************
      subroutine poly3mark(x,y,z,nx,nmark)
      integer nx,nmark,i
      character*4 mark
      real x(*),y(*),z(*),xp,yp,zp
      include 'plotcom.h'
      include 'world3.h'
      ipf=pfPS
c      if(nmark.eq.3)then
c         mark='!A3'//char(0)
c  That does not quite align well. Better not to mix thinkgs up.
      if(nmark.lt.10)then
         pfPS=0
         mark=char(nmark+176)//char(0)
      elseif(nmark.eq.10)then
         mark='+'//char(0)
      elseif(nmark.eq.11)then
         mark='!AX'//char(0)
      elseif(nmark.eq.12)then
         mark='!A*'//char(0)
      elseif(nmark.eq.13)then
         mark='!A-'//char(0)
      elseif(nmark.eq.15)then
         mark='!A'//char(48)//char(0)
      else
         mark=char(nmark)//char(0)
      endif
      do 1 i=1,nx
         call wxyz2nxyz(x(i),y(i),z(i),xp,yp,zp)
         call trn32(xp,yp,zp,x2,y2,z2,0)
         cx=x2+xcbc2
         cy=y2+ycbc2
      if(nmark.eq.14)then
         call actrid(cx,cy)
      else
         call jdrwstr(cx,cy,mark,0.)
      endif
    1 continue
      pfPS=ipf
      end
c************************************************************************
      subroutine polybox(xb,y,npts)
c Plot the outline of a histogram of boxes whose (touching) boundaries are
c at xb(0:npts) and whose heights are y(npts).
      integer npts
      real y(npts),xb(0:npts)
      real xp(4),yp(4)

      yp(3)=0.
      do i=1,npts
         xp(1)=xb(i-1)
         yp(1)=yp(3)
         xp(2)=xp(1)
         yp(2)=y(i)
         xp(3)=xb(i)
         yp(3)=y(i)
         if(i.lt.npts)then
            call polyline(xp,yp,3)
         else
            xp(4)=xp(3)
            yp(4)=0.
            call polyline(xp,yp,4)
         endif
      enddo

      end
C********************************************************************
c The following are strided versions which probably should become the
c code base, with polyline and polymark being just wrappers.
C********************************************************************
      subroutine stpolyline(x,y,npts,nstx,nsty)
c Strided polyline. nstx/y are the increments of the x and y data.
c Dashed line version.
      integer npts
      real x(npts),y(npts)
      integer i
      include 'plotcom.h'
c Dashed line code
      real nx,ny
      real vlen,dx,dy,cx,cy,plen,flen,dlen
      real wx2nx, wy2ny
      integer cud

c dashlen is the arc length in normalized units of the the ith line
c segment. dashdist is the starting fractional part of the ith arc.
c Segments alternate pen down, pen up.
      logical ldash
      real dashlen,dashdist
      integer MASKNO,dashmask,jmask
      parameter (MASKNO=4)
      dimension dashmask(MASKNO),dashlen(MASKNO)
      common/dashline/ldash,dashlen,dashdist,dashmask,jmask

      if(npts.le.0) return
      call vecw(x(1),y(1),0)
      do 3 i=2,npts
         if(.not.ldash) then
            call vecw(x(1+nstx*(i-1)),y(1+nsty*(i-1)),1)
         else
c We shall bypass vecw and go straight to normal.
            nx=wx2nx(x(1+nstx*(i-1)))
            ny=wy2ny(y(1+nsty*(i-1)))
c Lengths of total vector:
            cx=nx
            cy=ny
            dx=nx- crsrx
            dy=ny- crsry
            vlen=sqrt(DX*DX+DY*DY)
c Partial length remaining:
            plen=vlen
            if(vlen.eq.0)return
c Distance to end of segment
    1       dlen=(dashlen(jmask)-dashdist)
            if(plen.gt.dlen)then
c Vector longer than this segment. Draw segment and iterate.
               dashdist=0
               plen=plen-dlen
               flen=dlen/vlen
               nx= crsrx+dx*flen
               ny= crsry+dy*flen
               cud=dashmask(jmask)
c              call optvecn(nx,ny,cud)
               call vecn(nx,ny,cud)
               jmask=mod(jmask,MASKNO)+1
               goto 1
            else
c Vector ends before segment. Draw to end of vector and quit.
               dashdist=plen+dashdist
               nx=cx
               ny=cy
               cud=dashmask(jmask)
c              call optvecn(nx,ny,cud)
               call vecn(nx,ny,cud)
            endif
         endif
    3 continue
      end
C********************************************************************
C********************************************************************
      subroutine stpolymark(x,y,nx,nstx,nsty,nmark)
c Strided polyline. nstx/y are the increments of the x and y data.
      integer nx,nmark,i,nstx,nsty
      character*4 mark
      real x(*),y(*),xp,yp,wx2nx,wy2ny
      include 'plotcom.h'
      ipf=pfPS
c      if(nmark.eq.3)then
c         mark='!A3'//char(0)
c  That does not quite align well. Better not to mix thinkgs up.
      if(nmark.lt.10)then
         pfPS=0
         mark=char(nmark+176)//char(0)
      elseif(nmark.eq.10)then
         mark='+'//char(0)
      elseif(nmark.eq.11)then
         mark='!AX'//char(0)
      elseif(nmark.eq.12)then
         mark='!A*'//char(0)
      elseif(nmark.eq.13)then
         mark='!A-'//char(0)
      elseif(nmark.eq.15)then
         mark='!A'//char(48)//char(0)
      else
         mark=char(nmark)//char(0)
      endif
      do 1 i=1,nx
         xp=wx2nx(x(1+nstx*(i-1)))
         yp=wy2ny(y(1+nsty*(i-1)))
      if(nmark.eq.14)then
         call actrid(xp,yp)
      else
         call jdrwstr(xp,yp,mark,0.)
      endif
    1 continue
      pfPS=ipf
      end
C********************************************************************
      subroutine smoothline(x,y,np,nb)
c Draw a smoothed polyline using nb to determine the range of smoothing.
c If nb is positive, boxcar average over the points before and after to
c a total extent nb.
c If nb is negative then triangular average instead.
c If nb is zero, then no smoothing is done.
c np is the number of points in the trace, which must be less than the
c locally allocated array length nl.
      integer np,nb
      real x(np),y(np)

      integer nl
      parameter (nl=10000)
      real yave(nl)

      if(np.gt.nl)then
         write(*,*)'Trace to smooth too long. Not smoothed. Split up.'
c Here we ought to put code to do repetitive smoothing and plotting
c to complete the entire trace. However, it's a bit tricky because the
c internal joins of partial traces ought to use fully smoothed parts
c of the traces, not the reduced smoothing that is the default at the
c boundaries.
         call polyline(x,y,np)
      else
         if(nb.lt.0)then
            call triangave(np,-nb,y,yave)
         else
            call boxcarave(np,nb,y,yave)
         endif
         call polyline(x,yave,np)
      endif

      end
C********************************************************************
      subroutine polygapline(x,y,n,logic)
! plot a polyline in multiple sections only where logic is true
      integer n
      real x(n),y(n)
      logical logic(n)
      logical drawing
      drawing=.true.
      id=1
      do i=1,n
      if(drawing)then
         if(.not.logic(i))then
            drawing=.false.
            call polyline(x(id),y(id),i-id)
         endif
      else
         if(logic(i))then
            drawing=.true.
            id=i
         endif
      endif
      enddo
      if(drawing)call polyline(x(id),y(id),i-id)  
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine polycolorline(x,y,npts,ipcolor)
c Version to color a polyline as you go to encode other information.
c ipcolor must run from 1-240 for gradcolor.
      integer npts
      real x(npts),y(npts)
      integer ipcolor(npts)
      include 'plotcom.h'
c Dashed line code
      real nx,ny
      real vlen,dx,dy,cx,cy,plen,flen,dlen
      real wx2nx, wy2ny
      integer cud

c dashlen is the arc length in normalized units of the the ith line
c segment. dashdist is the starting fractional part of the ith arc.
c Segments alternate pen down, pen up.
      logical ldash
      real dashlen,dashdist
      integer MASKNO,dashmask,jmask
      parameter (MASKNO=4)
      dimension dashmask(MASKNO),dashlen(MASKNO)
      common/dashline/ldash,dashlen,dashdist,dashmask,jmask

      if(npts.le.0) return
      call vecw(x(1),y(1),0)
      do 3 i=2,npts
         call gradcolor(ipcolor(i-1))
         if(.not.ldash) then
            call vecw(x(i),y(i),1)
         else
c We shall bypass vecw and go straight to normal.
            nx=wx2nx(x(i))
            ny=wy2ny(y(i))
c Lengths of total vector:
            cx=nx
            cy=ny
            dx=nx- crsrx
            dy=ny- crsry
            vlen=sqrt(DX*DX+DY*DY)
c Partial length remaining:
            plen=vlen
            if(vlen.eq.0)return
c Distance to end of segment
    1       dlen=(dashlen(jmask)-dashdist)
            if(plen.gt.dlen)then
c Vector longer than this segment. Draw segment and iterate.
               dashdist=0
               plen=plen-dlen
               flen=dlen/vlen
               nx= crsrx+dx*flen
               ny= crsry+dy*flen
               cud=dashmask(jmask)
c              call optvecn(nx,ny,cud)
               call vecn(nx,ny,cud)
               jmask=mod(jmask,MASKNO)+1
               goto 1
            else
c Vector ends before segment. Draw to end of vector and quit.
               dashdist=plen+dashdist
               nx=cx
               ny=cy
               cud=dashmask(jmask)
c              call optvecn(nx,ny,cud)
               call vecn(nx,ny,cud)
            endif
         endif
    3 continue

      end
C********************************************************************
      subroutine polycolorgapline(x,y,n,ipcolor,logic)
! plot a polycolorline in multiple sections only where logic is true
      integer n
      real x(n),y(n)
      integer ipcolor(n)
      logical logic(n)
      logical drawing
      drawing=.true.
      id=1
      do i=1,n
      if(drawing)then
         if(.not.logic(i))then
            drawing=.false.
            call polycolorline(x(id),y(id),i-id,ipcolor(id))
         endif
      else
         if(logic(i))then
            drawing=.true.
            id=i
         endif
      endif
      enddo
      if(drawing)call polycolorline(x(id),y(id),i-id,ipcolor(id))  
      end
