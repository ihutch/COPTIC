      program slicetest

      integer iuds(3),ifull(3)
      integer Li
      parameter (Li=60)
      real u(Li,Li,Li),zp(Li,Li,3)
      real gradu(Li,Li,Li,3),vp(Li,Li,3,3)
      real xn(3*Li)
      integer ixnp(4)
      integer ifixpt(3)
      data ifixpt/3*0/
      data ifull/3*Li/


c      call setconlog(.true.)
c      call pfset(3)

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

c Calculate the gradient of u as a vector field to do arrow plotting.
      do k=1,iuds(3)
         do j=1,iuds(2)
            do i=1,iuds(1)
               if(i.eq.1)then
                  gradu(i,j,k,1)=(u(i+1,j,k)-u(i,j,k))/
     $                 (xn(ixnp(1)+i+1)-xn(ixnp(1)+i))
               elseif(i.eq.iuds(1))then
                  gradu(i,j,k,1)=(u(i,j,k)-u(i-1,j,k))/
     $                 (xn(ixnp(1)+i)-xn(ixnp(1)+i-1))
               else
                  gradu(i,j,k,1)=(u(i+1,j,k)-u(i-1,j,k))/
     $                 (xn(ixnp(1)+i+1)-xn(ixnp(1)+i-1))
               endif
               if(j.eq.1)then
                  gradu(i,j,k,2)=(u(i,j+1,k)-u(i,j,k))/
     $                 (xn(ixnp(2)+j+1)-xn(ixnp(2)+j))
               elseif(j.eq.iuds(2))then
                  gradu(i,j,k,2)=(u(i,j,k)-u(i,j-1,k))/
     $                 (xn(ixnp(2)+j)-xn(ixnp(2)+j-1))
               else
                  gradu(i,j,k,2)=(u(i,j+1,k)-u(i,j-1,k))/
     $                 (xn(ixnp(2)+j+1)-xn(ixnp(2)+j-1))
               endif
               if(k.eq.1)then
                  gradu(i,j,k,3)=(u(i,j,k+1)-u(i,j,k))/
     $                 (xn(ixnp(3)+k+1)-xn(ixnp(3)+k))
               elseif(k.eq.iuds(3))then
                  gradu(i,j,k,3)=(u(i,j,k)-u(i,j,k-1))/
     $                 (xn(ixnp(3)+k)-xn(ixnp(3)+k-1))
               else
                  gradu(i,j,k,3)=(u(i,j,k+1)-u(i,j,k-1))/
     $                 (xn(ixnp(3)+k+1)-xn(ixnp(3)+k-1))
               endif
            enddo
         enddo
      enddo
      

c      iuds(3)=iuds(3)-1

      ifixpt(3)=10
      ifix=0
c How to turn off pausing:
      call noeye3d(0)
      ntimes=1
      do i=1,ntimes
c And turn it on again.
         if(i.eq.ntimes)call noeye3d(9999)
c Arrowplot call, on contour plot in position 1, tell slice.
         ifix=4+16*1+64
c This is how to call with a fixed z-scale. Only the last two scale3
c arguments are relevant, adding 256 to ifix turns off z-scaling.
c         call scale3(0.,1.,0.,1.,-2.,7.)
c         ifix=4+16*1+64+256
         call sliceGweb(ifull,iuds,u,Li,zp, ixnp,xn,ifix,'potential:'/
     $        /'!Ay!@' ,gradu,vp)
      enddo
 

c      call sliceGcont(ifull,iuds,u,Li,zp,
c     $        ixnp,xn,ifixpt,'potential:'//'!Ay!@',dum,dum)

c      ifixpt(1)=-iuds(1)/2
         call sliceGcont(ifull,iuds,u,Li,zp, ixnp,xn,ifixpt,'potential:'
     $        //'!Ay!@' ,gradu,vp)

      end
