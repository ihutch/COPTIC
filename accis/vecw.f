c Modifications to vecn Jan 92
C********************************************************************
      subroutine vecw(x,y,ud)
c Draw a vector in world coordinates
      real x,y
      integer ud
      real wx2nx,wy2ny
      real nx,ny
      nx=wx2nx(x)
      ny=wy2ny(y)
c      call optvecn(nx,ny,ud)
c The optvec call was the start of attempts to do hiding in projected
c 2D calls. But was not completed and gave errors. See notes.
      call vecn(nx,ny,ud)
      return
      end
C********************************************************************
c Optionally draw hiding, normalized.
      subroutine optvecn(nx,ny,ud)
      real nx,ny
      integer ud
      include 'world3.h'
      if(ihiding.ne.0)then
         call hidvecn(nx,ny,ud)
      else
         call vecn(nx,ny,ud)
      endif
      end
c*******************************************************************
c Versions for c calls
      subroutine getwx2nx(wx,xn)
      xn=wx2nx(wx)
      end
      subroutine getwy2ny(wy,yn)
      yn=wy2ny(wy)
      end
C********************************************************************
      function wx2nx(wx)
      real wx2nx,wx
      include 'plotcom.h'
      real xd
      if(lxlog)then
         if(wx.lt.0.001*wxmin .or. wx.gt.1000.*wxmax) then
            write(*,*)'ACCIS WARNING world log x value outside range:'
     $           ,wx,' plotting outside box.'
            xd=.01/w2nx
         else
            xd=log10(wx)-log10(wxmin)
         endif
      else
	 xd=wx-wxmin
      endif
      wx2nx=naxmin+xd*w2nx
      return
      end
C********************************************************************
      function wy2ny(wy)
      real wy2ny,wy
      include 'plotcom.h'
      real yd
      integer errcount
      data errcount/0/
      if(lylog)then
         if(wy.lt.0.01*wymin)then
            if(errcount.lt.5)then
               write(*,*)'ACCIS WARNING world log y value',wy
     $           ,' far below box.'
            elseif(errcount.eq.5)then
               write(*,*)'Further log value warnings suppressed'
            endif
            errcount=errcount+1
            yd=-.01/w2ny
c            write(*,*)naymin,w2ny,naymin+yd*w2ny
         elseif(wy.gt.1000.*wymax) then
            if(errcount.lt.5)then
            write(*,*)'ACCIS WARNING world log y value',wy
     $           ,' far outside box.'
            elseif(errcount.eq.5)then
               write(*,*)'Further log value warnings suppressed'
            endif
            errcount=errcount+1
c            yd=.01/w2ny
            yd=log10(1000.*wymax)
         else
            yd=log10(wy)-log10(wymin)
         endif
      else
	 yd=wy-wymin
      endif
      wy2ny=naymin+yd*w2ny
      return
      end
C********************************************************************
      function xn2xw(nx)
      real xn2xw,nx
      include 'plotcom.h'
      real xd
      xd=(nx-naxmin)/w2nx
      if(lxlog)then
         xn2xw=wxmin*(10.**xd)
      else
         xn2xw=wxmin+xd
      endif
      return
      end
C********************************************************************
      function yn2yw(ny)
      real yn2yw,ny
      include 'plotcom.h'
      real yd
      yd=(ny-naymin)/w2ny
      if(lylog)then
         yn2yw=wymin*10.**yd
      else
         yn2yw=wymin+yd
      endif
      return
      end
C********************************************************************
      subroutine trn32(x,y,z,xt,yt,zt,ifl)
      real x,y,z,xt,yt,zt
      integer ifl
c Transform point (x,y,z) to point (xt,yt) via current trans.
c If ifl.eq.1 set the transform: x,y,z: looked at, xt,yt,zt: eye.
c Projection is on to the plane through x,y,z perp to r-rt.
c Thus scaling must be done first.
c If ifl.eq.2 axonometric transform. xt=dx/dy,zt=dz/dy
c If ifl.eq.3 return the coordinates relative to center,
c          rotated to the axes in which the center to eye vector is z.
c If ifl.eq.-1 return the eye position dx,dy,dz in xt,yt,zt.
      real dx,dy,dz,rz,d,t11,t12,t13,t21,t22,t23,t31,t32,t33,dmz
      real x0,y0,z0,zsign
      save
      include 'world3.h'
      data t11,t12,t13,t21,t22,t23,t31,t32,t33
     $	 / .894427, .447214,.000000,-.182574, .365148, .912871
     $  , -.408248, .816497,-.408248/
c the next data may not be consistent with the matrix.
      data d,dx,dy,dz/10.,4.,-10.,2./
      data x0,y0,z0/0.,0.,0./
      if(ifl.eq.0)then
c Perspective, in coordinates where z-axis is the center to eye vector, is
c xt=x/(d+z) yt=y/(d+z), with d=|center-eye|. Hence we need to
c transform to these coordinates and then do this perspective scaling.
	 zt=t31*(x-x0)+t32*(y-y0)+t33*(z-z0)
	 dmz=1.+zt/d
c Prevent the perspective from amplifying too much or becoming negative.
	 if(dmz.le.0.0001)then
c	    write(*,*)' TRN32 error: point at eye'
	    return
	 endif
	 xt=(t11*(x-x0)+t12*(y-y0)+t13*(z-z0))/dmz
	 yt=(t21*(x-x0)+t22*(y-y0)+t23*(z-z0))/dmz
      elseif(ifl.eq.1)then
c Set up perspective transform.
	 dx=xt-x
	 dy=yt-y
	 dz=zt-z
         x0=x
         y0=y
         z0=z
	 d=dx*dx+dy*dy
	 rz=sqrt(d)
	 d=sqrt(d+dz*dz)
c Transformation matrix: | -c1       s1        0 |   cos,sin etc:
c Rotz  till x' perp d # | +s1c2     c1c2      s2|  c1=dy/rz s1=dx/rz
c Rotx' till z''para  d  | -s1s2    -c1s2      c2|  c2=-dz/d s2=rz/d
	 t11=-dy/rz*z3sign
	 t12=dx/rz*z3sign
	 t13=0.
	 t33=-dz/d
	 t21=t12*t33
	 t22=-t11*t33
	 t23=rz/d*z3sign
	 t31=-dx/d
	 t32=-dy/d
c	 write(*,'(3f14.6)')t11,t12,t13,t21,t22,t23,t31,t32,t33
      elseif(ifl.eq.2)then
c Axonometric setup. x=x+dxdy*y, y=z+dzdy*y, z=y, d=infinity.
	 t11=1.
	 t12=xt
	 t13=0.
	 t21=0.
	 t22=zt
	 t23=1.
	 t31=0.
	 t32=1.
	 t33=0.
	 d=1.e30
      elseif(ifl.eq.3)then
	 xt=(t11*(x-x0)+t12*(y-y0)+t13*(z-z0))
	 yt=(t21*(x-x0)+t22*(y-y0)+t23*(z-z0))
	 zt=t31*(x-x0)+t32*(y-y0)+t33*(z-z0)
      elseif(ifl.eq.-1)then
c Return current.
         xt=dx+x0
         yt=dy+y0
         zt=dz+z0
         x=x0
         y=x0
         z=z0
      endif
      end
C********************************************************************
      subroutine tn2s(px,py,sx,sy)
c Transform the coordinates. Return the screen
c transformed vector coordinates gave drwstr errors.
      real px,py
      integer sx,sy
      include 'plotcom.h'
      include 'world3.h'
      real x2,y2,z2
      external tn2shidata

      if(ihiding.lt.0)then
c calls within 1,2, or 3-planes, at position fixedn.
	 if(ihiding.eq.-1)then
	    call trn32(fixedn,px,py,x2,y2,z2,0)
	 elseif(ihiding.eq.-2)then
	    call trn32(py,fixedn,px,x2,y2,z2,0)
	 elseif(ihiding.eq.-3)then
	    call trn32(px,py,fixedn,x2,y2,z2,0)
	 elseif(ihiding.eq.-4)then
	    call trn32(fixedn,py,px,x2,y2,z2,0)
	 elseif(ihiding.eq.-5)then
	    call trn32(px,fixedn,py,x2,y2,z2,0)
	 elseif(ihiding.eq.-6)then
	    call trn32(py,px,fixedn,x2,y2,z2,0)
	 else
	    stop 'Unknown ihiding switch value'
	 endif
	 px=x2+xcbc2
	 py=y2+ycbc2
c	 x2=x2+xcbc2
c	 y2=y2+ycbc2
c	 sx=x2*scrxpix
c	 sy=scrypix-y2*n2sy
      endif
c standard call.
      sx=(px*scrxpix)
      sy=(scrypix-py*n2sy)
c         write(*,*)'scrxpix,scrypix',scrxpix,scrypix
c         write(*,*)'sx,sy,n2sy,px,py',sx,sy,n2sy,px,py
      end
C********************************************************************
c  Draw a vector in normalized coordinates.
      subroutine vecn(nx,ny,ud)
      real nx,ny
      integer ud
      include 'plotcom.h'
      real*4 prx,pry,crx,cry
      integer ret,sx,sy,ptrunc
      crx=nx
      cry=ny
      if(ltlog)then
	 prx=crsrx
	 pry=crsry
	 ret=ptrunc(prx,pry,crx,cry)
         ret2=ret/16
         ret1=ret-16*ret2
	 if(ret.ne.99)then
	    if(ret1.gt.0)then
	       call tn2s(prx,pry,sx,sy)
	       if(pfsw.ge.0)call vec(sx,sy,0)
	       if(pfsw.ne.0)call vecnp(prx,pry,0)
	    endif 
            call tn2s(crx,cry,sx,sy)
	    if(pfsw.ge.0) call vec(sx,sy,ud)
	    if(pfsw.ne.0) call vecnp(crx,cry,ud)
	    if(ret2.gt.0)then
c End point moved. Break the line here.
	       if(pfsw.ne.0)call vecnp(crx,cry,0)
	    endif
	 endif
      else
	 call tn2s(crx,cry,sx,sy)
	 if(pfsw.ge.0) call vec(sx,sy,ud)
	 if(pfsw.ne.0) call vecnp(crx,cry,ud)
      endif
      crsrx=nx
      crsry=ny
      return
      end
C********************************************************************
      subroutine truncf( x1, x2, y1, y2)
      real x1,y1,x2,y2
      include 'plotcom.h'
c   Turn on or off (if all args zero) truncation (windowing).
      trcxma=0.
      if(x1.eq.0.)then
	 if(x2.eq.0.)then
	    if(y1.eq.0.)then
	       if(y2.eq.0.)then
c Old action:	  ltlog=.false.
c New approach set truncation at the screen boundary plus20%, Not infinity.
                  ltlog=.true.
                  trcxmi=-.2
                  trcxma=1.2
                  trcymi=-.2
                  trcyma=yoverx*1.2
		  return
	       endif
	    endif
	 endif
      endif
      if(trcxma.eq.0.)then
         ltlog=.true.
         trcxmi=x1
         trcxma=x2
         trcymi=y1
         trcyma=y2
      endif
      return
      end

c**********************************************************************
c*   Truncate within the rectangle given by xma,yma,xmi,ymi in trunc.
c* Return 99 if whole vector outside, 
c* ptrunc bits 0-4 set if x1,y1 is moved; bits 5-7 set if x2,y2 moved
      function ptrunc(x1,y1,x2,y2)
      integer ptrunc
      real x1,y1,x2,y2
      include 'plotcom.h'
      real d1,d2
      integer ic
      ic=0

      d1=trcxmi-x1
      d2=trcxmi-x2
      if(d1.gt.0)then
	 if(d2.gt.0)then
	    ptrunc=99
	    return
	 endif
	 x1=trcxmi
	 y1=(-d2* y1+d1* y2)/(d1-d2)
	 ic=ic+1
      else if(d2.gt.0)then
	 x2=trcxmi
	 y2=(d2*y1-d1*y2)/(d2-d1)
         ic=ic+16
      endif

      d1=trcymi-y1
      d2=trcymi-y2
      if(d1.gt.0)then
	 if(d2.gt.0)then
	    ptrunc=99
	    return
	 endif
	 y1=trcymi
	 x1=(-d2*x1+d1*x2)/(d1-d2)
	 ic=ic+1
      else if(d2.gt.0)then
	 y2=trcymi
	 x2=(d2*x1-d1*x2)/(d2-d1)
         ic=ic+16
      endif

      d1=x1-trcxma
      d2=x2-trcxma
      if(d1.gt.0)then
	 if(d2.gt.0)then
	    ptrunc=99
	    return
	 endif
	 x1=trcxma
	 y1=(-d2*y1+d1*y2)/(d1-d2)
	 ic=ic+1
      else if(d2.gt.0)then
	 x2=trcxma
	 y2=(d2*y1-d1*y2)/(d2-d1)
         ic=ic+16
      endif

      d1=y1-trcyma
      d2=y2-trcyma
      if(d1.gt.0)then
	 if(d2.gt.0)then
	    ptrunc=99
	    return
	 endif
	 y1=trcyma
	 x1=(-d2*x1+d1*x2)/(d1-d2)
	 ic=ic+1
      else if(d2.gt.0)then
	 y2=trcyma
	 x2=(d2*x1-d1*x2)/(d2-d1)
         ic=ic+16
      endif
      ptrunc=ic
      return
      end


c***************************************************************************
	subroutine winset(wsw)
	logical wsw
      include 'plotcom.h'
      if(wsw)then
      call truncf(naxmin,naxmax,naymin,naymax)
      else
      call truncf(0.,0.,0.,0.)
      endif
      return
      end

