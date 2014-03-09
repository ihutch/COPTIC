c Test of surfdr
c Visualization of a parameterized surface.
c x and y are considered to be the parameters. z is then the dependent
c variable. However, x and y can be nonmonotonic and thus can represent
c a closed surface.

      parameter (ni=50,nj=40)
c      parameter (ni=10,nj=10)
      real x(ni,nj),y(ni,nj),z(ni,nj),work(0:ni+1,0:nj+1)
      real d(5)

c A sphere of radius 1.
      rs=1.
      do j=1,nj
         yj=-rs+2.*rs*(j-1)/(nj-1)
         rc=sqrt(rs**2-yj**2+1.e-6)
         do i=1,ni
            y(i,j)=yj
            phi=(3.141593*2.*(i-1)/(ni-1))
            x(i,j)=rc*cos(phi)
            z(i,j)=rc*sin(phi)
         enddo
      enddo

      isw=0
 1    call pltinit(0.,1.,0.,1.)
c Set the cube size.
      call setcube(.2,.1,.1,.5,.4)
c Alternative automatic call.
c      level=4
c      call surf3d(x,y,z,ni,ni,nj,level,work)
c Explicit axis control.
c     Set the perspective transform.
      call geteye(x2,y2,z2)
      call trn32(0.,0.,0.,x2,y2,z2,1)
c Set scale to maintain aspect ratios.
      call scale3(-2.,2.,-1.,1.,-1.,1.)
c Set color to zero to prevent web drawing.
      call color(0)
c      goto 10
c Draw the surface.
      call surfdr3(x,y,z,ni,ni,nj,work,isw,d)
      call color(15)
      call axproj(igetcorner())
      call ax3labels('x','y','z')
      if(ieye3d().ne.0)goto 1

c Draw using directional shading
 10   d(1)=1.
      d(2)=1.
      d(3)=1.
c 4,5 give distance limits.
      d(4)=-.8
      d(5)=.8
      call pfset(3)
      call pltinit(0.,1.,0.,1.)
c Put the pen on the drawing before setting color to make BBox correct.
      call vec3n(0.,0.,0.,0)
c Set color to zero to prevent web drawing.
      call color(0)
      call accisgradinit(22000,-40000,-40000,64000,64000,64000)
      call surfdr3(x,y,z,ni,ni,nj,work,3,d)
      call color(15)
      call axproj(igetcorner())
      call ax3labels('x','y','z')
      if(ieye3d().ne.0)goto 10
c      call pltend()
      
      end
