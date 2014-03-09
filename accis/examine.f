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
      include 'world3.h'
      z3sign=z3sign*z3signch
      z3signch=1.
      call cubeupd(x2,y2,z2)
      call puteye(x2,y2,z2)
      end
c***************************************************************************
      subroutine cubeupd(x2,y2,z2)
c Set the new transform and redraw the cube.
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
      include 'world3.h'
      z3signch=1.
      rp0=(x0*x0 + y0*y0)
      r0 = sqrt(rp0+z0*z0)
      if(r0.eq.0.) stop 'Viewrot error: zero initial length'
      tht=acos(z0/r0)
      phi=atan2(y0,x0)
      phinew=phi-z3sign*dist2ang*xmoved
      thtnew=tht-z3sign*dist2ang*ymoved
      rpnew=r0*sin(thtnew)
      if(sin(thtnew)*sin(tht).lt.0.)z3signch=-1.
      zn=r0*cos(thtnew)
      xn=rpnew*cos(phinew)
      yn=rpnew*sin(phinew)
c      write(*,'(a,9f8.4)')'xn,yn,zn,z3sign,th,ph',xn,yn,zn,z3sign,tht
c     $     ,thtnew,phinew
      end
c************************************************************************
c Convenience function for eye3d.
      function ieye3d()
      integer ieye3d
      call eye3d(ieye3d)
      end
