c***********************************************************************
      program fillgradtest
      real x(3),y(3),z(3),d(3)
      real xq(4),yq(4),zq(4),xr(4),yr(4),zr(4)
      integer ng0,ng1
      real zg0,zg1 

      data x/ 0.3,0.0,0./
      data y/ 0.3,0.0,0.8/
      data z/ 0.9,0.0,0.4/
      data d/  .0 , .0 ,1. /
      data ng0,ng1/0,230/
      data zg0,zg1/.0,1./
      data xq/0. , .8, .8,0. /
      data yq/0. ,0. , .4, .8/
      data zq/0. , .8, 0.,1. /
      data xr/0.,1.,1.,0./
      data yr/0.,0.,1.,1./
      data zr/0.,0.,1.,1./

c It is important to realize that GL filling of triangles knows only the
c three colors at the three corners. Therefore if the fill is nonlinear,
c as it is if there is saturation, then GL triangle filling over substantial
c distances will give a non-intuitive result. And in particular a quad fill
c does not come out uniform but gives an X shaped pattern.

      call pfset(3)
      call accisgradinit(0,65535,0,65535,0,65535)
      call pltinit(0.,1.,0.,1.)
c      call axis()
      call gradtri(x,y,z,d,zg0,zg1,ng0,ng1,2)
      call pltend()
      call pltinit(0.,1.,0.,1.)
      call gradquad(xq,yq,zq,d,zg0,zg1,ng0,ng1,2)
      
      call pltend()
      call pfset(3)
      call accisgradinit(-32000,-64000,0,96000,64000,128000)
      call pltinit(0.,1.,0.,1.)
      call gradquad(xr,yr,zr,d,zg0,zg1,ng0,ng1,2)
      call pltend()

      end
