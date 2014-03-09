      program arrowtest
      parameter (Li=200,imax=20,jmax=50)
      real E1(Li,Li),E2(Li,Li),x(Li),y(Li)
      real arrow(12)
      real path(3,imax)

      do i=1,imax
         x(i)=i
      enddo
      do j=1,jmax
         y(j)=j
      enddo
      do i=1,imax
         do j=1,jmax
            r=sqrt(x(i)**2+y(j)**2)
            theta=atan2(y(j),x(i))
            Emag= sin(2*theta)**2*r/(imax+(r**2/(.2*imax)))
            E1(i,j)=-Emag*sin(theta)
            E2(i,j)= Emag*cos(theta)
         enddo
      enddo
      if(.false.)then
      call pfset(3)
      call pltinit(0.,x(imax),0.,y(jmax))
      call axis()
      call arrowplot(E1,E2,imax/2.,Li,imax,jmax,x,y,1,idum,idum) 
      call pltend()

      call pltinit(0.,x(imax),0.,y(jmax))
      call axis()
      call arrowplot(E1,E2,imax/2.,Li,imax,jmax,x,y,5,4,4) 
      call pltend()

      Erange=0.
      call pltinit(0.,x(imax),0.,y(jmax))
      call axis()
      call arrowplot(E1,E2,Erange,Li,imax,jmax,x,y,5,3,6) 
      write(*,*)'Scaling Erange=',Erange
      call pltend()

      Erange=0.
      call pltinit(0.,x(imax),0.,y(jmax))
      call axis()
      call arrowplot(E1,E2,Erange,Li,imax,jmax,x,y,9,idum,idum) 
      call pltend()
      endif

c 3-D arrow plotting.
      scale=1.
      do i=1,imax
         path(1,i)=scale*.7*cos(i/3.)
         path(2,i)=scale*.7*sin(i/3.)
         path(3,i)=scale*.2*cos(i/3.)
      enddo


      do i=1,3
         arrow(i)=-.5
         arrow(i+3)=.5
      enddo
      arrow(7)=.2
c      arrow(8)=.05
c Negative means use the last two arguments.
      arrow(8)=-.1
      arrow(9)=0.2
      arrow(10)=0.1
c nangles, whether filled:
      arrow(11)=8
c Switch controlling behavior. Fill, solidshaft,
      arrow(12)=1+2
c This is pretty much a minimal 3-D plot.
      call trn32(0.,0.,0.,5.,2.5,2.,1)
      call accisgradinit(22000,-40000,-40000,64000,64000,64000)
 101  call pltinit(-1.,1.,-1.,1.)
      call setcube(.2,.2,.2,.5,.4)
      call scale3(-1.,1.,-1.,1.,-1.,1.)

c First attempt:
c      call draw3arrow(arrow)
c Improved code puts only the path into the first argument:
c      call arrow3path(arrow,2,8,0,.1,.3,.2,1)
      isw=0
      call arrow3path(path,2,6,isw,.03,.1,.05,3)
      call cubeproj(icorner)
      call axident3()
      call axproj(igetcorner())
      call eye3d(ival)
      call rotatezoom(ival)
c For subsequent calls, just redraw surface, do not recalculate.
      if(isw/4-2*(isw/8).ne.1)isw=isw+4
      if(ival.ne.0 .and. ival.ne.ichar('q'))goto 101

      end
