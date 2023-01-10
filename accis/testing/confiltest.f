c  Test of contouring routine.
      save
      integer nx,ny,nlmax
      parameter (nx=30,ny=41,nlmax=50)
c      parameter (nx=8,ny=8,nlmax=50)
c      parameter  (nx=20,ny=17,nlmax=50)
      real z(nx,ny),x(nx,ny),y(nx,ny),cl(nlmax)
      character*1 ppath(nx,ny)
      integer nl,i,j,n,snl
      parameter (ipno=240)
      integer ired(ipno),igreen(ipno),iblue(ipno)

      icolsmooth=0
      nl=8
c      icolsmooth=1
c icolsmooth determines if we use gradients, and then color-step.
c      write(*,'('' Enter No of contours (<50), and colorsmooth'')')
c      read(*,*)nl,icolsmooth
      snl=nl
      nl=abs(nl)
      space=1.
      do 100 i=1,nx
         do 200 j=1,ny
            x(i,j)=0.7*float(i-4)
            y(i,j)=j-4
            z(i,j)=sin(space*(x(i,j)-.3*y(i,j))/2.)
     $           *sin(space*y(i,j)/4.)
 200     continue
 100  continue
      do 300 n=1,nl
         cl(n)=(-1.+ 2.*(n-1)/(nl-1))*1.
 300  continue
c      call pfPSset(1)
c      call pfset(3)
c 3.General call. 
       call pfset(3)
       call pltinit(x(1,1),x(nx,1),y(1,1),y(1,ny))
c Poor-man's gamma effect with truncated color curves.
c Vecx version:
c       call accisgradinit(22000,-40000,-40000,64000,64000,64000)
c Explicit version colorgradient achieving the same result.
       do i=1,ipno
          ired(i)=(22000*(ipno-i)+64000*(i-1))/(ipno-1)
          igreen(i)=(-40000*(ipno-i)+64000*(i-1))/(ipno-1)
          iblue(i)=(-40000*(ipno-i)+64000*(i-1))/(ipno-1)
       enddo
       call accisgradset(ired,igreen,iblue,ipno)
c x,y matrices used since last arg 2. Coloring since +16. +64 tricolor.
       isw=2+16
       if(icolsmooth.ne.0)isw=isw+64+256*icolsmooth
       call contourl(z,ppath,nx,nx,ny,cl,snl,x,y,isw)
c       call axis
       call axlabels('x','y')
       call pltend
       if(.false.)then
       call pltinit(x(1,1),x(nx,1),y(1,1),y(1,ny))
       npixel=ipno
       do k=1,npixel 
          if(k.le.npixel/2)then
             ired(k)=65535
             iblue(k)=65535
             igreen(k)=65535
          else
             ired(k)=55000
             iblue(k)=55000
             igreen(k)=55000
          endif
       enddo
       call accisgradset(ired,igreen,iblue,ipno)
c x,y matrices used since last arg 2. Coloring since +16. +64 tricolor.
       isw=2+16
       if(icolsmooth.ne.0)isw=isw+64+256*icolsmooth
       call contourl(z,ppath,nx,nx,ny,cl,snl,x,y,isw)
       call axis
       call axlabels('x','y')
       call pltend
       endif
c Individual contour filling.
c The switch in contourl is 
c      ifcolor=(theconsw/65536-256*(theconsw/(256*65536)))
       call pltinit(x(1,1),x(nx,1),y(1,1),y(1,ny))
       do i=4,4
          isw=2+i*65536
          call contourl(z,ppath,nx,nx,ny,cl(i),1,x,y,isw)
       enddo
       call axis
       call axlabels('x','y')
       call pltend

      stop
      end

