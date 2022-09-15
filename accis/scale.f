c****************************************************************************
        function ybyx()
c return the ratio of window height to width.
        real ybyx
        include 'plotcom.h'
        ybyx=(naymax-naymin)/(naxmax-naxmin)
        end
C********************************************************************
      subroutine pfset(isw)
c Set the plot-to-file mode. If switch=-1 ask from console.
      integer isw
      include 'plotcom.h'
c Don't switch immediately. Just set the value at the next pltinit.
      pfnextsw=isw
      goto 1
    2 write(*,*)' Plot to file? (0:no,1:hp,2:ps,3:eps)'
      read(*,*)pfnextsw
    1 if(pfnextsw.eq.-1)goto 2
      return
      end
C********************************************************************
      subroutine pfget(isw)
c Get the plot-to-file mode.
      integer isw
      include 'plotcom.h'
      isw=pfsw
      end
C********************************************************************
c   Set world to normalized scalings.
c   If min and max are both zero, leave as before.
      subroutine scalewn(wxmi,wxma,wymi,wyma,lx,ly)
      real wxmi,wxma,wymi,wyma
      logical lx,ly
      integer nerr,nemax
      include 'plotcom.h'
      data nerr/0/nemax/5/
      lxlog=lx
      lylog=ly
      if(wxmi.lt.0 .or. wxma.lt.0)lxlog=.false.
      if(wymi.lt.0 .or. wyma.lt.0)lylog=.false.
      if(wxmi.ne.0..or.wxma.ne.0.)then
         wxmin=wxmi
         wxmax=wxma
         if(wxmin.eq.wxmax)then
            if(nerr.le.nemax)write(*,*)
     $           ' SCALEWN warning: wxmin/max coincide;fixing.'
            wxmax=wxmax+1.
            nerr=nerr+1
         endif
         if(.not.lxlog)w2nx=(naxmax-naxmin)/(wxmax-wxmin)
         if(lxlog)w2nx=(naxmax-naxmin)/(log10(wxmax)-log10(wxmin))
      endif
      if(wymi.ne.0..or.wyma.ne.0.)then
         wymin=wymi
         wymax=wyma
         if(wymin.eq.wymax)then
            if(nerr.le.nemax)write(*,*)
     $           ' SCALEWN warning: wymin/max coincide;fixing.'
            wymax=wymax+1.
            nerr=nerr+1
         endif
         if(.not.lylog)w2ny=(naymax-naymin)/(wymax-wymin)
         if(lylog)w2ny=(naymax-naymin)/(log10(wymax)-log10(wymin))
      endif
      return
      end
c******************************************************************
      subroutine fitscale(xmin,xmax,ymin,ymax,lx,ly)
      real xmin,xmax,ymin,ymax
      logical lx,ly
      real xfac,xdelta,fxmin,fymin,fxmax,fymax
      integer nxfac
      if(lx)then
         fxmin=10.**(nint(log10(xmin)-0.49999))
         fxmax=10.**(nint(log10(xmax)+0.49999))
      else
         call fitrange(xmin,xmax,6,nxfac,xfac,xdelta,fxmin,fxmax)
      endif
      if(ly)then
         fymin=10.**(nint(log10(ymin)-0.49999))
         fymax=10.**(nint(log10(ymax)+0.49999))
      else
         call fitrange(ymin,ymax,6,nxfac,xfac,xdelta,fymin,fymax)
      endif
      call scalewn(fxmin,fxmax,fymin,fymax,lx,ly)
      end
c******************************************************************
c Routines for specifying colors by name.
        function idarkblue()
        idarkblue=1
        end
        function idarkgreen()
        idarkgreen=2
        end
        function iskyblue()
        iskyblue=3
        end
        function ibrickred()
        ibrickred=4
        end
        function ipurple()
        ipurple=5
        end
        function igold()
        igold=6
        end
        function igray()
        igray=7
        end
        function ilightgray()
        ilightgray=8
        end
        function iblue()
        iblue=9
        end
        function igreen()
        igreen=10
        end
        function icyan()
        icyan=11
        end
        function ired()
        ired=12
        end
        function imagenta()
        imagenta=13
        end
        function iyellow()
        iyellow=14
        end
        function iblack()
        iblack=15
        end

