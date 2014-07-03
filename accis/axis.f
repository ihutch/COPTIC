c************************************************************************
c    Axis Routines:  Including GAXIS Normal units.1  23 Jan 93.
c************************************************************************
c         Automatic axis routine.
      subroutine axis()
      real first,delta
      delta=0
      call xaxis(first,delta)
      call yaxis(first,delta)
      call axbox()
      return
      end
c*********************************************************************
      subroutine axbox()
      include 'plotcom.h'
      call vecn(naxmin,naymin,0)
      call vecn(naxmax,naymin,1)
      call vecn(naxmax,naymax,1)
      call vecn(naxmin,naymax,1)
      call vecn(naxmin,naymin,1)
      return
      end
c*********************************************************************
      subroutine xaxis(first,delta)
c  Draw xaxis using the general axis routine.
      real first,delta
      include 'plotcom.h'
      real wx2nx
      integer nxfac,tog
      tog=0
      if(nrows.gt.1)then
         if( ((multype/2) .eq. 0) .and.
     $        (mod(nframe,nrows).ne.0) .and.
     $        (nxlabw.gt.0) ) tog=1  
      endif
      if(tog.eq.1)call ticlabtog()
      nxfac=0
         call gaxis(wxmin,wxmax,nxfac,first,delta, 
     $     wx2nx(wxmin),wx2nx(wxmax),naypt,naypt,.true.,lxlog)
      if(tog.eq.1)call ticlabtog()
      end
c*********************************************************************
      subroutine yaxis(first,delta)
c  Draw yaxis using the general axis routine.
      real first,delta
      include 'plotcom.h'
      real wy2ny
      integer nyfac,tog
      tog=0
      if(ncolumns.gt.1)then
         if( ((multype-2*(multype/2)) .eq. 0) .and.
     $        (nframe.gt.nrows .or. nframe.eq.0) .and.
     $        (nxlabw.gt.0) ) tog=1  
      endif
      if(tog.eq.1)call ticlabtog()
      nyfac=0
         call gaxis(wymin,wymax,nyfac,first,delta,
     $     naxpt,naxpt,wy2ny(wymin),wy2ny(wymax),.false.,lylog)
      if(tog.eq.1)call ticlabtog()
      end
c*********************************************************************
      subroutine gaxis(amin,amax,ngpow,first,delta,
     $     xgmin,xgmax,ygmin,ygmax,lpara,laxlog)
      real amin,amax,first,delta,xgmin,xgmax,ygmin,ygmax
      integer ngpow
      logical lpara,laxlog
c General axis routine.
c amin,amax           :  min and max values of axis, (without power).
c ngpow               :  power of ten for axis (label)
c first,delta         :  Set tic points. Auto if delta=0.
c xgmin/max,ygmin/max :  normal cordinates of axis ends
c lpara               :  labels parallel to axis (.true.)
c laxlog              :  log axis (true).
      real gw,dx,dy,dl,axcos,axsin,lfac,ams,tcos,tsin
      integer i
      include 'plotcom.h'
      real xdelta,x1st,xlast
      real ain,aax,x10,xm
      real xgnlog,ygnlog
c      real xgnlin,ygnlin
      integer nxfac,indx,iticnum
c Statement functions for scaling.
      xgnlog(gw)=xgmin+dx*log10(gw/ain)
      ygnlog(gw)=ygmin+dy*log10(gw/ain)
c      xgnlin(gw)=xgmin+dx*(gw-ain)
c      ygnlin(gw)=ygmin+dy*(gw-ain)
c
      iticnum=ticnum
      axcos=xgmax-xgmin
      axsin=ygmax-ygmin
      dl=sqrt(axcos*axcos+axsin*axsin)
      if(dl.eq.0) then 
        call txtmode
        write(*,*) ' GAXIS error. Axis ends coincide. '
        write(*,*)' ygmin=',ygmin,' ygmax',ygmax,' wymin',wymin,
     $       ' wymax',wymax
        write(*,*)' xgmin=',xgmin,' xgmax',xgmax,' wxmin',wxmin,
     $       ' wxmax',wxmax
        stop ' '
      endif
      axcos=axcos/dl
      axsin=axsin/dl
      if(lpara)then
         tcos=-axsin
         tsin=axcos
      else
         tcos=axsin
         tsin=-axcos
      endif
      ams=sign(1.,(amax-amin))
      xdelta=delta
      if(delta .ne. 0) x1st=first
      nxpow=ngpow
      xpow=10**nxpow
      ain=amin
      aax=amax

      if(laxlog)then
C logarithmic
         if(xdelta.eq.0)then
c Standard labelling.
            xdelta=2.
            x1st=0.
         endif
c  Log axes: Put labeled ticks at powers of ten, and ten*x1st, 
c      and unlabeled at ten*(n*xdelta) for n integer,
c      within the range of the axes.
c      Typical choices: x1st=5,xdelta=1 ; x1st=2,xdelta=2 ; x1st=0,xdelta=2.
         dy=log10(aax/ain)
         dx=(xgmax-xgmin)/dy
         dy=(ygmax-ygmin)/dy
         if(ain.gt.aax)then
            x10=ain
            ain=aax
            aax=x10
         endif
         if(ain.le.0.) stop ' GAXIS error: Log axis end negative.'
         do 1 indx=nint(log10(ain)-.49999),nint(log10(aax)+.49999)
            x10=10.**indx
            gw=x10
            if(gw.ge.(.9999*ain) .and. gw.le.(1.0001*aax))then
               call gtic(gw,2,xgnlog(gw),ygnlog(gw)
     $              ,tcos,tsin,1.,0.,lpara)
            endif
            gw=x10*x1st
            if(gw.ge.ain .and. gw.le.aax)then
               call gtic(gw,3,xgnlog(gw),ygnlog(gw)
     $              ,tcos,tsin,1.,0.,lpara)
            endif
            do 2 ixm=1,100
               xm=ixm*xdelta
               if(xm.gt.9.9999)goto 1
               gw=x10*xm
               if(gw.ge.ain .and. gw.lt.aax) then
                  call gtic(gw,0,xgnlog(gw),ygnlog(gw)
     $                 ,tcos,tsin,1.,0.,lpara)
               endif
    2       continue
    1    continue
         call vecn(xgmin,ygmin,0)
         call vecn(xgmax,ygmax,1)

      else
c Linear
         if(xdelta.eq.0)then
C Autoscale
            if(dl.lt.chrswdth*3*ticnum)then
c Axis is too short for comfort. Change ticnum
               iticnum=max(1,int(dl/(chrswdth*3)))
c               write(*,*)'Resetting ticnum',iticnum
            endif
            call fitrange(ain,aax,iticnum,
     $           nxfac,xpow,xdelta,x1st,xlast)
c Changed to -1 May 2002 for more satisfactory labeling.
            if(.not.(nxfac.gt.2.or.nxfac.lt.-1)) then
               nxfac=0
               xpow=1.
            endif
            ngpow=nxfac
            x1st=x1st/xpow
            xdelta=xdelta/xpow
            ain=ain/xpow
            aax=aax/xpow
         endif
c  Linear axes: Put ticks at first+n*delta for n integral only within axis.
         dx=(xgmax-xgmin)/(aax-ain)
         dy=(ygmax-ygmin)/(aax-ain)
         gw=x1st+ xdelta*nint(.4999+(ain-x1st)/xdelta)
c If needed increase decimal places in tic.
         ipoint=nint(log10(abs(x1st)+1.)-.49999)
     $        -nint(log10(abs(xdelta))-.49999)
c         write(*,*)'ipoint',ipoint,x1st,xdelta
         if(ipoint.gt.1)then
            inlabp=nxlabp
            nxlabp=ipoint
         endif
         if(lminor)then
c Choose the spacing of minor tics gwmt.
            xld=alog10(abs(xdelta))
            ild=nint(xld-0.4999)
            idd=nint(10.**(xld-ild))
            if(idd.le.1)then
               gwmt=0.2*10.**ild
            elseif(idd.le.3)then
               gwmt=0.5*10.**ild
            else
               gwmt=10.**ild
            endif
            if(.not.(gwmt.lt.abs(xdelta).and.gwmt.gt.0.))then
               write(*,*)'Minor tic calculation error',gwmt,xld,ild,idd
               write(*,*)xdelta,x1st
               stop
            endif
c Start earlier so as to fill in minor tics prior to xfirst:
            gwm=gw-ams*abs(xdelta)
         endif
c Iterate over major tics.
         do 4 i=1,40
            if(lminor)then
c gw is being incremented from amin to amax.
               gw1=gw+ams*abs(xdelta)
c Minor tics
 6             gwm=gwm+ams*gwmt
               if(ams*gwm.lt.ams*gw1 .and.
     $              (gwm-(aax+0.0001*(aax-ain)))*ams.le.0.)then
                  if((gwm-(ain-0.0001*(aax-ain)))*ams.ge.0.)then
                     call gtic(gwm,0,xgmin+dx*(gwm-ain),ygmin+dy*(gwm
     $                    -ain),tcos,tsin,axcos,axsin,lpara)
                  endif
                  goto 6
               endif
            endif
            if((gw-(aax+0.0001*(aax-ain)))*ams.gt.0.)goto 5
c This ought not to be necessary, I think
c           if((gw-(ain-0.0001*(aax-ain)))*ams.lt.0.)goto 5
            if(ngpow.eq.0)then
               call gtic(gw,1,xgmin+dx*(gw-ain),ygmin+dy*(gw-ain)
     $              ,tcos,tsin,axcos,axsin,lpara)
            else
               if(lpara)then
                lfac=.11*(chrswdth/0.015)*0.6/dl
                else
                lfac=.03*(chrshght/0.015)/dl
                endif
               if((gw-((1.-lfac)*aax+lfac*ain))*ams.lt.0.) then
                  call gtic(gw,1,xgmin+dx*(gw-ain),ygmin+dy*(gw-ain)
     $                 ,tcos,tsin,axcos,axsin,lpara)
               else
                  call gtic(gw,4,xgmin+dx*(gw-ain),ygmin+dy*(gw-ain),
     $                 tcos,tsin,axcos,axsin,lpara)
               endif
            endif
            gw=gw+ams*abs(xdelta)
            gwm=gw
    4    continue
    5    continue
         if(ipoint.gt.1)then
c     reset tic label places.
            nxlabp=inlabp
         endif
         if(ngpow.ne.0.and.nxlabw.gt.0)then
            call powlabel(ngpow,xgmax,ygmax,tcos,tsin,
     $                 axcos,axsin,lpara)
         endif
      endif
      call vecn(xgmin,ygmin,0)
      call vecn(xgmax,ygmax,1)
      return
      end

c*********************************************************************
      subroutine powlabel(ngpow,xg,yg,tcos,tsin,axcos,axsin,lpara)
c Put a power of ten label at xg,yg offset etc.
      real xg,yg,tcos,tsin,axcos,axsin
      integer ngpow
      logical lpara
      include 'plotcom.h'
      character*80 ticlab,label
      real cc,cs,px,py
      integer lablength

      px=xg
      py=yg
      cc=chrscos
      cs=chrssin
      call iwrite(ngpow,lablength,ticlab)
      label='!u' // ticlab(1:lablength)// '!u' // char(0)
c      ticlab='!AX!@10' //char(0)
      ticlab='!AX!@10'//label //char(0)
      if(lpara)then
         chrscos=axcos
         chrssin=axsin
         px=px+tcos*xticoff
         py=py+tsin*xticoff
      else
         chrscos=tcos
         chrssin=tsin
         px=px+tcos*(yticoff-chrswdth)
         py=py+tsin*(yticoff-chrswdth)
      endif
      if(yticoff.le.0)then
         call jdrwstr(px,py,ticlab,-0.3)
      else
         call drwstr(px,py,ticlab)
      endif
c      call drwstr(px,py,label)
      chrscos=cc
      chrssin=cs
      end

c*********************************************************************
      subroutine gtic(gw,lab,xg,yg,tcos,tsin,axcos,axsin,lpara)
      real gw,xg,yg,tcos,tsin,axcos,axsin
      integer lab
      logical lpara
c gw   value at tic
c lab  type of label
c xg,yg  tic position in normal coords.
c tcos,tsin  tic angle
c axcos,axsin  axis angle
c lpara  tic label parallel to axis.
      include 'plotcom.h'
c commons:      nxlabw   le 0 switches off the labels.
c               nxlabp   number of decimal places in labels.
c               xticlen  can reverse tics if negative (also length).
      character*80 ticlab,label
      real cc,cs
      integer nl,lablength,labp

      call vecn(xg,yg,0)
      if(lab.eq.0)then
         call vecn(xg+tcos*.6*xticlen,yg+tsin*.6*xticlen,1)
         return
      endif
      if(lab.eq.1)then
         labp=nxlabp
         call fwrite(gw,lablength,labp,ticlab)
         ticlab(lablength+1:lablength+1)=char(0)
      elseif(lab.eq.2)then
         nl=nint(-.49999+log10(gw))
         call iwrite(nl,lablength,label)
         label(lablength+1:lablength+3)= char(28) // 'u'// char(0)
         ticlab='10' // char(28) // 'u'// label(1:lablength+3)
      elseif(lab.eq.3)then
c Circumlocution to avoid mingw/wine bug.
         temp=gw/(10.**(nint(-.49999+log10(gw))))
         call fwrite(temp,lablength,0,ticlab)
      elseif(lab.eq.4)then
         call vecn(xg+tcos*xticlen,yg+tsin*xticlen,1)
         return
      endif
      call vecn(xg+tcos*xticlen,yg+tsin*xticlen,1)
c label:
      if(nxlabw.gt.0)then
      cc=chrscos
      cs=chrssin
      if(lpara)then
         chrscos=axcos
         chrssin=axsin
         call jdrwstr(xg+tcos*xticoff,yg+tsin*xticoff,ticlab,0.)
      else
         chrscos=tcos
         chrssin=tsin
         call jdrwstr(xg+tcos*yticoff,yg+tsin*yticoff,ticlab,
     $   yticoff/abs(yticoff))
      endif
      chrscos=cc
      chrssin=cs
      endif
      return
      end
c*********************************************************************
      subroutine altxaxis(xi,xa)
c  Draw alternate xaxis using the general axis routine.
c  Arguments are the factors to scale the axis ends.
c  All tic shifts must be explicit.
      real xi,xa
      real first,delta
      include 'plotcom.h'
      real wx2nx
      integer nxfac
      delta=0
      nxfac=0
         call gaxis(wxmin*xi,wxmax*xa,nxfac,first,delta, 
     $     wx2nx(wxmin),wx2nx(wxmax),naypt,naypt,.true.,lxlog)
      end
c*********************************************************************
      subroutine altyaxis(yi,ya)
c  Draw alternate yaxis using the general axis routine.
c  Arguments are the factors to scale the axis ends.
      real yi,ya
      real first,delta
      include 'plotcom.h'
      real wy2ny
      integer nyfac
      delta=0
      nyfac=0
         call gaxis(wymin*yi,wymax*ya,nyfac,first,delta,
     $     naxpt,naxpt,wy2ny(wymin),wy2ny(wymax),.false.,lylog)
      end
c*********************************************************************
      subroutine togminor()
c Toggle drawing minor tics
      include 'plotcom.h'
      lminor=.not.lminor
      end



