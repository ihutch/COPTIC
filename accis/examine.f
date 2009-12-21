c These routines are called by the X driver for interactive rotation.
c***************************************************************************
      subroutine butdown(xt,yt,zt)
      real xt,yt,zt
      real x,y,z
      real x1,y1,z1
c Get back the current position.
      call trn32(x,y,z,x1,y1,z1,-1)
c      write(*,*)'xt=',x1,y1,z1
      xt=x1
      yt=y1
      zt=z1
      end
c***************************************************************************
      subroutine butup(x2,y2,z2)
      real x2,y2,z2
      call puteye(x2,y2,z2)
      end
c***************************************************************************
      subroutine cubeupd(x2,y2,z2)
c 
      real x2,y2,z2
      real x,y,z
      data x,y,z/0.,0.,0./
      call color(0)
      call cubed(0)
      call trn32(x,y,z,x2,y2,z2,1)
      call color(15)
      call cubed(0)
      end
c***************************************************************************
      subroutine viewrot(xmoved,ymoved,x0,y0,z0,xn,yn,zn)
c return new eye position rotated around origin by xmoved,ymoved.
      real xmoved,ymoved,x0,y0,z0,xn,yn,zn
      real dist2ang
      parameter (dist2ang=.005)
      rp0=(x0*x0 + y0*y0)
      r0 = sqrt(rp0+z0*z0)
      if(r0.eq.0.) stop 'Viewrot error: zero initial length'
      tht=-dist2ang*xmoved
      xn=x0*cos(tht) - y0*sin(tht)
      yn=x0*sin(tht) + y0*cos(tht)
      zn=z0+ymoved*.1*sqrt(rp0)/r0 
      rn=sqrt(xn*xn + yn*yn + zn*zn)
      if(rn.eq.0.) stop 'Viewrot error: zero length'
      xn=r0*xn/rn
      yn=r0*yn/rn
      zn=r0*zn/rn
      end
c************************************************************************
c Convenience function for eye3d.
      function ieye3d()
      integer ieye3d
      call eye3d(ieye3d)
      end
