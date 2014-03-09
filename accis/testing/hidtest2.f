	program hidtest
      integer nx,ny,ud,i,j,iLx
      parameter (iLx=100,ny=99)
      real z,x,y,r,yy
      dimension z(iLx,ny),x(iLx),y(ny)
      integer nl
      parameter (nl=10)
      real cl(nl),ht
c contour storage.
      character pp(4,iLx,ny)
      parameter (ht=.02)
c      include '\fortran\graph\hid\world3.h'
	include 'world3.h'
      nx=30
      do 10 i=1,nl
	 cl(i)=ht*(-1.+i*2./nl)
   10 continue
      call pfset(-1)
      do 1 j=1,ny
	 y(j)=float(j)-ny/3. +0.5
	 yy=y(j)*y(j)
	 do 2 i=1,nx
	    x(i)=float(i)-nx/3.+0.5
	    r=sqrt(x(i)*x(i)+yy)/2.
c	    z(i,j)=ht*sin(x(i))*sin(y(j))/(x(i)*y(j))
	    z(i,j)=ht*sin(r)/r
	    ud=i-1
    2	 continue
    1 continue
      call pltinit(0.,1.,0.,1.)
c       Plot the surface. With axes. Web color 10, axis color 7.
      j=1 + 256*10 + 256*256*7
   99 call hidweb(x,y,z,iLx,nx,ny,j)
c       Draw a contour plot in perspective. Need to reset color anyway.
      call color(6)
      call axregion(-scbx3,scbx3,-scby3,scby3)
      call scalewn(x(1),x(nx),y(1),y(ny),.false.,.false.)
      call hdprset(-3,scbz3)
c       call axis   ! if desired.
      call scalewn(1.,float(nx),1.,float(ny),.false.,.false.)
c       Contour without labels, direct on mesh.
      call contourl(z,pp,iLx,nx,ny,cl,nl,r,r,16)
      call color(15)
      call pltend
      end


