c  Test of contouring routine.
      save
      integer nx,ny,nlmax
      parameter (nx=30,ny=40,nlmax=50)
      real z(nx,ny),x(nx,ny),y(nx,ny),cl(nlmax)
      character*1 ppath(nx,ny)
      integer nl,i,j,n,snl
      common/cont1stlast/c1st,clast

      write(*,'('' Enter No of contours, (<50)'')')
      read(*,*)nl
      snl=nl
      nl=abs(nl)
      do 100 i=1,nx
	 do 200 j=1,ny
	    x(i,j)=0.7*float(i-4)
	    y(i,j)=j-4
	    z(i,j)=sin((x(i,j)-.3*y(i,j))/2.)*sin(y(i,j)/4.)
 200	 continue
 100  continue
      do 300 n=1,nl
	 cl(n)=-1.+ 2.*n/nl
 300  continue
c      call pfPSset(1)
c      call glback()
      call pfset(3)
c      write(*,'(10f7.3)')(cl(i),i=1,nl)
      call multiframe(2,2,3)
c 1.Simplest call autocontours on a rectangular mesh when nl=0.
      call contrec(z,ppath,nx,ny,cl,0)
c Put axes, annotations etc on afterward if desired.
      call axis
      call axlabels('x-index','y-index')
c 2.General call. x,y not used since last arg zero.
       call pltinit(1.,float(nx),1.,float(ny))
       call contourl(z,ppath,nx,nx,ny,cl,snl,x,y,0)
c How to draw suitably scaled axes if desired:
       call scalewn(x(1,1),x(nx,1),y(1,1),y(1,ny),.false.,.false.)
       call axis
       call axbox
c       call pltend
c 3.General call. 
       call pltinit(x(1,1),x(nx,1),y(1,1),y(1,ny))
c Poor-man's gamma effect with truncated color curves.
       call accisgradinit(-40000,-40000,32000,64000,64000,64000)
c x,y matrices used since last arg 2. Coloring since +16.
c       write(*,*)'Finished gradinit'
       call contourl(z,ppath,nx,nx,ny,cl,snl,x,y,2+16)
c       write(*,*)'Finished contourl'
       call axis
       call axlabels('x','y')

c       write(*,*)c1st,clast
       call gradlegend(c1st,clast,.0,-.4,1.,-.4,.03,.true.)

       call pltend


      end

