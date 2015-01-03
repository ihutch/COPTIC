c Example of drawing a 3-D web and projected contour plot.
	program hidtest
	integer nx,ny,ud,i,j,iLx
	parameter (iLx=50,ny=20)
	real z,x,y,r,yy
	dimension z(iLx,ny),x(iLx),y(ny)
	integer nl
	parameter (nl=10)
	real cl(nl),ht
c       contour storage.
	character pp(4,iLx,ny)
	parameter (ht=.025)
	include 'world3.h'
c Set up data etc.
	nx=10
	do 10 i=1,nl
	   cl(i)=ht*(-1.+i*2./nl)
 10	continue
c	call pfset(-1)
	call pfset(3)
	do 1 j=1,ny
	   y(j)=-(float(j)-ny/3. +0.5)
	   yy=y(j)*y(j)
	   do 2 i=1,nx
	      x(i)=-(float(i)-nx/3.+0.5)
	      r=sqrt(x(i)*x(i)+yy)/2.
c       z(i,j)=ht*sin(x(i))*sin(y(j))/(x(i)*y(j))
	      z(i,j)=ht*sin(r)/r
	      ud=i-1
 2	   continue
 1	continue
c Start of actual plotting.
 98	call pltinit(0.,1.,0.,1.)
c       Plot the surface. With axes (1). Web color 10, axis color 7.
	j=1 + 256*10 + 256*256*7
c	write(*,*)'Starting hidweb'
	call hidweb(x,y,z,iLx,nx,ny,j)
c	write(*,*) 'Done hidweb'
c       Draw a contour plot in perspective. Need to reset color anyway.
	call color(4)
	call axregion(-scbx3,scbx3,-scby3,scby3)
	call scalewn(x(1),x(nx),y(1),y(ny),.false.,.false.)
	call hdprset(-3,scbz3)
c	write(*,*) 'Done hdprset'
c       call axis   ! if desired.
	call scalewn(1.,float(nx),1.,float(ny),.false.,.false.)
c       Contour without labels, direct on mesh.
c	write(*,*) 'Done scalewn'
c The coloring does not work properly on some displays unless initialized.
c	call accisgradinit(0,0,0,64000,64000,64000)
c Test of smooth output.
 	call contourl(z,pp,iLx,nx,ny,cl,nl,r,r,16+64)
c	call contourl(z,pp,iLx,nx,ny,cl,nl,r,r,16)
c	write(*,*) 'Done contourl'
	call color(5)
	call ax3labels('x-axis','y-axis','z-axis')
	j=1 + 256*10 + 256*256*7
	call geteye(x2,y2,z2)
c If we look from underneath, redraw the web over the contour plot.
	if(z2.lt.0.) call hidweb(x,y,z,iLx,nx,ny,j)
	call color(15)
c How to enable interactive plot rotation. Just this call instead of pltend:
	if(ieye3d().ne.0) goto 98
c	call autoplot(x,y,ny)
c	call pltend()
c	call autoplot(x,y,ny)
c	call pltend()
	end


