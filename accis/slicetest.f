      program slicetest

      integer iuds(3),ifull(3)
      integer Li
      parameter (Li=60)
      real u(Li,Li,Li),zp(Li,Li,3)
      real xn(3*Li)
      integer ixnp(4)
      integer ifixpt(3)
      data ifixpt/3*0/
      data ifull/3*Li/

      range=10.
      amp=1.
      wave=.5

      ixnp(1)=0
      do id=1,3
         iuds(id)=50
         do i=1,iuds(id)
            xn(i+ixnp(id))=range*(-1.+2.*(i-1.)/(iuds(id)-1.))
         enddo
         ixnp(id+1)=ixnp(id)+iuds(id)
      enddo

c This meaningless call corrects an apparent gfortran bug 
      write(*,'(a,$)')' '

      do k=1,iuds(3)
         z=xn(ixnp(3)+k)
         do j=1,iuds(2)
            y=xn(ixnp(2)+j)
            do i=1,iuds(1)
               x=xn(ixnp(1)+i)
               arg=wave*y
               xx=x*cos(arg)+z*sin(arg)
               yy=-x*sin(arg)+z*cos(arg)
               u(i,j,k)=1./((xx/range-.2)**2+(yy/5.)**2+.2)
c               write(*,*)x,y,z,u(i,j,k)
            enddo
         enddo
      enddo

c      iuds(3)=iuds(3)-1

      ifixpt(3)=10
      call sliceGweb(ifull,iuds,u,Li,zp,
     $        ixnp,xn,ifix,'potential:'//'!Ay!@')
      call sliceGcont(ifull,iuds,u,Li,zp,
     $        ixnp,xn,ifixpt,'potential:'//'!Ay!@')

      end
