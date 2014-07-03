c Tektronix 4014 driver.
c Expanded to supply dummy routines for filling and rotating etc 2008.
c And for usleep, noeye3d. Routines like slicing won't actually work
c with the 4014 driver, but at least the code will run without the
c X libraries.
c********************************************************************
      integer function accis_driver()
      accis_driver=2
      end
C********************************************************************
      blockdata scrndat
      include 'plotcom.h'
c      data scrxpix,scrypix,ncolor,vmode/1024,779,15,4010/
      data scrxpix,scrypix,ncolor,vmode/4096,3120,15,4014/      
      end
C********************************************************************
c Switch to graphics mode Tek 4010/14.
      subroutine svga(scrxpix,scrypix,vmode,ncolor)
      integer scrxpix,scrypix,vmode,ncolor
      accis_nodisplay=1;
c Enter Tek mode. Modified for Xterm.
      write(*,'(1x,a)')char(27)//'[?38h'
      write(*,*)'                                 '
c      write(*,'(1x,a)')char(27)//char(12)
c      write(*,*)'                                 '
c Clear screen twice seems to do the trick for Kermit.
      write(*,'(1x,a)')char(27)//char(12)
      write(*,*)'                                 '
c Extra clear screen for non Kermit.
      write(*,'(1x,a)')char(12)
      return
      end
C********************************************************************
      subroutine txtmode
      character resp
c Write padding, so we don't lose characters on switch back.
      read(*,'(a1)')resp
      write(*,*)char(24),'                                         '
      write(*,*)char(24),'                                         '
c Fix for xterm. Did not work.
c      write(*,'(1x,a)')char(27)//'[?38l'
      end
c*********************************************************************
c 4014 vector driver.
      subroutine vec(x,y,iud)
      integer x,y,iud
      include 'plotcom.h'
      integer oylst,oyhi,oxlow,oxhi
      integer ylow,yhi,xlow,xhi,i,xlst,ylst
      integer hlab,ylab,xlab,istart
      parameter (hlab=32,ylab=96,xlab=64)
      character*6 start
c      parameter (start=char(29)//char(hlab)//char(ylab)//char(ylab)
c     $  //char(hlab)//char(xlab))
      parameter (start=' '''' @')
      character*80 outchr
      save
      data i,istart/6,1/
      data outchr(1:6)/start/
c Crunch on to the screen.
      if(x.lt.0)x=0
      if(x.gt.scrxpix)x=scrxpix
      if(y.lt.0)y=0
      if(y.gt.scrypix)y=scrypix
c Separate the coordinate 'nibbles'.
      xhi=x/128
      xlst=x-128*xhi
      xlow=xlst/4
      xlst=xlst-(xlow)*4
      xhi=xhi+hlab
      xlow=xlow+xlab
      ylst=3119-y
      yhi=ylst/128
      ylst=ylst-128*yhi
      ylow=ylst/4
      ylst=ylst-ylow*4
      yhi=yhi+hlab
      ylow=ylow+ylab
      ylst=ylab+4*ylst+xlst
      if(iud.gt.0)then
c Continuing vector. Send only the necessary parts.
         if(yhi.ne.oyhi)then
            i=i+1
            outchr(i:i)=char(yhi)
         endif
         if(ylst.ne.oylst)then
            i=i+1
            outchr(i:i)=char(ylst)
         endif
         i=i+1
         outchr(i:i)=char(ylow)
         if(xhi.ne.oxhi)then
            i=i+1
            outchr(i:i)=char(xhi)
         endif
         i=i+1
         outchr(i:i)=char(xlow)
         istart=0
      else
c Start vector.
         istart=1
      endif
      if(i.ge.74.or.istart.eq.1)then
c Finish draw and start again, if we are longer than a line or starting.
         write(*,999)outchr(1:i)
  999    format(1x,a)
         outchr(1:6)=char(29)//char(yhi)//char(ylst)
     $        //char(ylow)//char(xhi)//char(xlow)
         i=6
      endif
c Update the vector.
      oxhi=xhi
      oxlow=xlow
      oyhi=yhi
      oylst=ylst
      return
      end
c*****************************************************************
      subroutine scolor(li)
      include 'plotcom.h'
      integer li,mask
c Scale the color so that (0..15) -> (0...ncpre-1) if ncpre >16.
c But a nonzero color always makes some mark because mask>0. check this!
      if(li.eq.0)then
         mask=0
      elseif(li.eq.8)then
         mask=7
      else
         mask=mod(li,8)
      endif
c Terminate the vector so as to flush buffer
      call vecn(crsrx,crsry,0)
c Set the color.
c ANSI sequences for kermit, assumed if we are a 4014. Flush buffer first.
c      write(*,*)
c      write(*,'(1x,a1,''['',i2,''m'')')char(27),30+mask
      return
      end
c*********************************************************************

C********************************************************************
      subroutine svganodisplay(scrxpix,scrypix,vmode,ncolor)
      integer scrxpix,scrypix,vmode,ncolor
      call svga(scrxpix,scrypix,vmode,ncolor)
      end
c********************************************************************
c dummy
      subroutine vecfill()
      end
c********************************************************************
      subroutine eye3d(ival)
      ival=0
      end
c********** Use a gradient color out of 240 *************************
      subroutine acgradcolor(li)
      integer li
      integer a_gradPixno
      parameter (a_gradPixno=240)
      integer a_gradPix(a_gradPixno)
      integer a_grad_inited
      integer a_gradred(a_gradPixno)
      integer a_gradgreen(a_gradPixno)
      integer a_gradblue(a_gradPixno)
      common /a_grad/a_gradPix,a_gradred,a_gradgreen,a_gradblue
     $     ,a_grad_inited
      external a_grad_data

      if(a_grad_inited.eq.0) call accisgraddef()
      end
c***********************************************************************
      block data a_grad_data
      integer a_gradPixno
      parameter (a_gradPixno=240)
      integer a_gradPix(a_gradPixno)
      integer a_grad_inited
      integer a_gradred(a_gradPixno)
      integer a_gradgreen(a_gradPixno)
      integer a_gradblue(a_gradPixno)
      common /a_grad/a_gradPix,a_gradred,a_gradgreen,a_gradblue
     $     ,a_grad_inited
      data a_grad_inited/0/
      end
c********** Tell the current rgb color ********************************
      subroutine getrgbcolor(ipixel,red,green,blue)
      integer ipixel,red,green,blue;
      integer a_gradPixno
      parameter (a_gradPixno=240)
      integer a_gradPix(a_gradPixno)
      integer a_grad_inited
      integer a_gradred(a_gradPixno)
      integer a_gradgreen(a_gradPixno)
      integer a_gradblue(a_gradPixno)
      common /a_grad/a_gradPix,a_gradred,a_gradgreen,a_gradblue
     $     ,a_grad_inited
      if(a_grad_inited.eq.0) call accisgraddef()
      red=a_gradred(ipixel)
      green=a_gradgreen(ipixel)
      blue=a_gradblue(ipixel)
      end
c*************************************************************************
      subroutine accisgraddef()
      integer top, bot
      top=65535;
      bot=0;
      call accisgradinit(bot,bot,bot,top,top,top);
      end
c*************************************************************************
      subroutine accisgradinit(r1,g1,b1,r2,g2,b2)
      integer r1,g1,b1,r2,g2,b2
      integer i,j
c      integer ipixel,red,green,blue;
      integer a_gradPixno
      parameter (a_gradPixno=240)
      integer a_gradPix(a_gradPixno)
      integer a_grad_inited
      integer a_gradred(a_gradPixno)
      integer a_gradgreen(a_gradPixno)
      integer a_gradblue(a_gradPixno)
      common /a_grad/a_gradPix,a_gradred,a_gradgreen,a_gradblue
     $     ,a_grad_inited
      do i=0,a_gradPixno-1,1
         j=(i*r2+(a_gradPixno-1-i)*r1)/(a_gradPixno-1.)
         if(j.lt.0)then
            j=0 
         elseif(j.gt.65535) then
            j=65535
         endif
         a_gradred(i)=j
         j=(i*g2+(a_gradPixno-1-i)*g1)/(a_gradPixno-1.)
         if(j.lt.0)then
            j=0 
         elseif(j.gt.65535) then
            j=65535
         endif
         a_gradgreen(i)=j
         j=(i*b2+(a_gradPixno-1-i)*b1)/(a_gradPixno-1.)
         if(j.lt.0)then
            j=0 
         elseif(j.gt.65535) then
            j=65535
         endif
         a_gradblue(i)=j;
         a_grad_inited=2;
      enddo
      end
c*************************************************************************
      subroutine accisgradset(red,green,blue,ngcol)
c Set the gradient colors from three fortran arrays.
c Limiting the range to 0-65535, warning if the pixel number is not right.

      integer red(ngcol),green(ngcol),blue(ngcol)
      integer a_gradPixno
      parameter (a_gradPixno=240)
      integer a_gradPix(a_gradPixno)
      integer a_grad_inited
      integer a_gradred(a_gradPixno)
      integer a_gradgreen(a_gradPixno)
      integer a_gradblue(a_gradPixno)
      common /a_grad/a_gradPix,a_gradred,a_gradgreen,a_gradblue
     $     ,a_grad_inited

      do i=1,a_gradPixno
         if(i.le.ngcol)then
            a_gradred(i)=max(0,min(red(i),65535))
            a_gradgreen(i)=max(0,min(green(i),65535))
            a_gradblue(i)=max(0,min(blue(i),65535))
         else
            a_gradred(i)=0
            a_gradgreen(i)=0
            a_gradblue(i)=0
         endif
      enddo
      end
c***********************************************************************
c Dummy
      integer function igradtri(x,y,z,h,i3d)
      igradtri=0
      end
c**********************************************************************
c Dummy
      subroutine usleep(usecs)
      end
c**********************************************************************
c Dummy
      subroutine noeye3d(value)
      end
