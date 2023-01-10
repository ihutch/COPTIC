c************************************************************************
c  Test of contouring routine.
      integer nx,ny,nlmax,nnx,nny,nl,n,i,j,snl
      parameter (nx=50,ny=50,nlmax=50)
      real z,x,y,cl
      dimension z(nx,ny),x(nx,ny),y(nx,ny),cl(nlmax)
      character ppath(4,nx,ny)
      integer consw,nxv

      write(*,'(2a)')' Enter No of contours, -ve for no labels',
     $      ', 0 for auto, and x-, y-dims <= 50.'
      read(*,*)nl,nnx,nny
      write(*,*)' Enter type of plot: 0, 1, 2'
      write(*,*)' raw-mesh, vector-rectangular, matrix-irregular mesh:'
      read(*,*)consw
      snl=0
      if(nl.lt.0) then
         snl=-1
         nl=-nl
      endif
      nxv=nnx
c Set up the arrays.
      do 100 i=1,nnx
         do 200 j=1,nny
            x(i,j)=float(i-4)
            y(i,j)=j**2*2./nny
            z(i,j)=sin((x(i,j)-.3*y(i,j))/4.)*sin(y(i,j)/6.)
  200    continue
  100 continue
      if(consw.eq.1)then
         do 400 i=1,nny
            y(i,1)=i**2*2./nny
  400    continue
      endif
      do 300 n=1,nl
         cl(n)=-1.+ 2.*n/nl
  300 continue
      write(*,*)' Contour levels,'
      write(*,'(10f7.3)')(cl(i),i=1,nl)
      write(*,'('' Mesh: nx,ny,sign'',3i4)')nnx,nny,snl
c
c Plotting calls: use built-in query.
      call pfset(-1)
c Setup call. Three types of prototype.
      if(consw.eq.0)then
c Call for raw spacing.
         call pltinit(1.,float(nnx),1.,float(nny))
         call axbox
      elseif(consw.eq.1)then
c Call for x,y vectors. Since this scale is really used, the ranges need to be
c correct.
         call pltinit(x(1,1),x(nnx,1),y(1,1),y(nny,1))
         call axis
      else
c Call for x,y arrays.
         call pltinit(x(1,1),x(nnx,1),y(1,1),y(1,nny))
         call axis()
      endif
c Contour plot
      call contourl(z,ppath,nx,nnx,nny,cl,(snl*nl),x,y,consw)
      call pltend()
      end
