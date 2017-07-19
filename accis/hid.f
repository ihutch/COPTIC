      subroutine hidweb(x,y,z,iLx,nx,ny,ilevel)
c Draw a 3-d web of z(x,y), dim(iLx\nx,ny), viewed from (xe,ye,ze) scaled.
c Second byte of level gives web color. Third gives axis color.
c Lowest byte: 
c   abs(level)=1 scale to fit the region, 1-d x,y vectors.
c   abs(level)=2 scale to fit region, 2-d x,y, don't hide, just wiremesh.
c   abs(level)=0 perform no scale-setting and use last perspective...
c   if bit2 (4) set, use last perspective regardless.
c   if bit3 (8) set, use last scaling regardless.
c   level.lt.0 draw no axes.
c Eye obtained from file eye.dat, or default if none exists.
      integer iLx, nx,ny,ilevel,level
      real x(*),y(*),z(iLx,ny)
      integer icorner,colw,cola
      real x2,y2,z2
      real zmin,zmax
      integer ihd
      data ihd/99/x2/0./y2/0./z2/0./
      save

      level=ilevel
      cola=level/256
      level=level-cola*256
      ipers=level/4-(level/8)*2
      iscale=level/8-(level/16)*2
      level=level-(level/4)*4
c level is lowest 8 bits.
      colw=cola-(cola/256)*256
c color of web is next 8 bits
      cola=cola/256 -(cola/65536)*256
c color of axes is next 8 bits.
      if(abs(level).ne.0 .and. ipers.eq.0)then
         call geteye(x2,y2,z2)
      endif
      if(colw.ne.0) call color(colw)
      if(abs(level).eq.1)then
         ihd=1
         ixy=1
c Set the scaling.
         xmin=x(1)
         xmax=x(nx)
         ymin=y(1)
         ymax=y(ny)
      elseif(abs(level).eq.2)then
         ihd=0
         ixy=2
c Set the scaling.       call minmax2(z,iLx,nx,ny,zmin,zmax)
         call minmax2(x,iLx,nx,ny,xmin,xmax)
         call minmax2(y,iLx,nx,ny,ymin,ymax)
      endif
c Set the top and bottom horizons.
      call hidinit(0.,1.)
c Set to non-hiding (ihd=0) or hiding (ihd=1) 3-D calls.
      call hdprset(ihd,0.)
      if(abs(level).ne.0)then
c Set the perspective transform.
         if(ipers.ne.1)call trn32(0.,0.,0.,x2,y2,z2,1)
         itics=5
         call minmax2(z,iLx,nx,ny,zmin,zmax)
         call fitrange(zmin,zmax,itics,ipow,fac10,delta,first,xlast)
         zmin=first
         zmax=xlast
         if(iscale.ne.1)call scale3(xmin,xmax,ymin,ymax,zmin,zmax)
      endif
c Draw the web
      if(ihd.eq.99)then
         write(*,*)'hidweb error. Scaling not set, but non-setting call'
         stop
      else
         if(ixy.eq.1)call webdrw(x,y,z,iLx,nx,ny,icorner)
         if(ixy.eq.2)call webdr2(x,y,z,iLx,nx,ny,icorner)
      endif
      if(cola.ne.0) call color(cola)
      if(level.ge.0)then
c Draw cube.
         call cubed(igetcubecorner())
c Draw axes.
         call axproj(igetcorner())
      endif
      end
c********************************************************************
      subroutine webdr2(x,y,z,iLx,nx,ny,icorner)
c Draw a 3-d web of z(x,y), dim(iLx\nx,ny), using current scaling.
c Return the nearest corner to eye in icorner.
c This version for x,y *arrays* rather than vectors.
      integer nx,ny,iLx
      real x(iLx,ny),y(iLx,ny),z(iLx,ny)
      integer id1,id2,ud,kx,ky,icorner
c      real d1start,d1end,d1step,z1,x2,y2,z2
      integer d1start,d1end,d1step
      real z1,x2,y2,z2
      save
c      include 'plotcom.h'
c
c Set the web order, choose the nearest corner to eye.
      do 20 id1=1,4
c This seems to be erroneous when x or y is decreasing.
         kx=1+mod(id1/2,2)*(nx-1)
         ky=1+((id1-1)/2)*(ny-1)
         xw=x(kx,ky)
         yw=y(kx,ky)
         call wxyz2nxyz(xw,yw,0.,xn,yn,zn)
         call trn32(xn,yn,zn,x2,y2,z2,0)
         if(id1.eq.1 .or. z2.lt.z1) then
            z1=z2
            icorner=id1
         endif
   20 continue
      
c Draw
      d1start=1
      d1end=nx+ny-2
      d1step=1
      if(icorner.eq.2 .or. icorner.eq.3) then
c                 Reverse the outer diagonal loop-order.
         d1start=d1end
         d1end=1
         d1step=-1
      endif
      do 1 id1=d1start,d1end,d1step
         ud=0
         do 2 id2=-min(id1,2*ny-1-id1),min(id1,2*nx-1-id1)
            kx=(id1+id2)/2 +1
            ky=(id1-id2)/2 +1
            if(icorner.eq.2 .or. icorner.eq.4) ky=ny+1-ky
c                   reverse the y-order.
            call vec3w(x(kx,ky),y(kx,ky),z(kx,ky),ud)
            ud=1
    2    continue
    1 continue
      end

c***************************************************************************
c********************************************************************
      subroutine webdrw(x,y,z,iLx,nx,ny,icorner)
c Draw a 3-d web of z(x,y), dim(iLx\nx,ny), using current scaling.
c Return the nearest corner to eye in icorner.
      integer nx,ny,iLx
      real x(nx),y(ny),z(iLx,ny)
      integer id1,id2,ud,icorner
c      real d1start,d1end,d1step,z1,x2,y2,z2
      integer d1start,d1end,d1step
      integer kx,ky
      real z1,x2,y2,z2
      save
c      include 'plotcom.h'
c
c Set the web order, choose the nearest corner to eye.
      do 20 id1=1,4
         kx=1+mod(id1/2,2)*(nx-1)
         ky=1+((id1-1)/2)*(ny-1) 
         xw=x(kx)
         yw=y(ky)
         call wxyz2nxyz(xw,yw,0.,xn,yn,zn)
         call trn32(xn,yn,zn,x2,y2,z2,0)
c        call trn32(x(kx(id1)),y(ky(id1)),0.,x2,y2,z2,0)
         if(id1.eq.1 .or. z2.lt.z1) then
            z1=z2
            icorner=id1
         endif
   20 continue
c      write(*,*) 'icorner=',icorner,z1
c Draw
      d1start=1
      d1end=nx+ny-2
      d1step=1
      if(icorner.eq.2 .or. icorner.eq.3) then
c                 Reverse the outer diagonal loop-order.
         d1start=d1end
         d1end=1
         d1step=-1
      endif
      do 1 id1=d1start,d1end,d1step
         ud=0
         do 2 id2=-min(id1,2*ny-1-id1),min(id1,2*nx-1-id1)
            kxi=(id1+id2)/2 +1
            kyi=(id1-id2)/2 +1
            if(icorner.eq.2 .or. icorner.eq.4) kyi=ny+1-kyi
c                   reverse the y-order.
            call vec3w(x(kxi),y(kyi),z(kxi,kyi),ud)
            ud=1
    2    continue
    1 continue
      end
c***************************************************************************
c Draw a vector in normalized coordinates, hiding as appropriate.
      subroutine hidvecn(x2,y2,ud)
      integer ud
      real x2,y2,x1,y1,xo,yo
c      integer ngrid
c      parameter (ngrid=1025)
c      real ytop(ngrid),ybot(ngrid)
c      common/hideln/ytop,ybot,...
      external hdcdata
      include 'hidcom.h'
      real dx,dydx,x,y
      integer signd,ix,ix1,ix2,nstate,istate
      logical lmidl
      data lmidl/.false./
      save

      if(ud.eq.0)then
         call vecn(x2,y2,0)
         x1=x2
         y1=y2
         lmidl=.false.
         return
      endif

      if((x2-x1).ne.0.)then
         xo=x1
         yo=y1
c Grid points just inside the x1 - x2 range, not outside frame.
         if(x2.ge.x1)then
            ix1=max(int(x1*(ngrid-1)),0) +1
            ix2=min(int(x2*(ngrid-1)),ngrid)
            signd=1
         else
            signd=-1
            ix1=min(int(x1*(ngrid-1)),ngrid)
            ix2=max(int(x2*(ngrid-1)),0)+1
         endif
         dx=1./(ngrid-1)
         dydx=(y2-y1)/(x2-x1)
         icount=0
         do ix=ix1,ix2,signd
            icount=icount+1
c  Do over Grid-position:
            x=dx*ix
            y=y1+(x-x1)*dydx
            nstate=0
            if(y.gt.ytop(ix))then
               nstate=1
               ytop(ix)=y
            endif
            if(y.lt.ybot(ix))then
               nstate=2
               ybot(ix)=y
            endif
c If previous call was not a moveto if(lmidl)then
c Finishing and starting segments.
            if(nstate.eq.0.and.lmidl)then
c Not drawing. If previously was, then finish segment.
               if(istate.ne.0)then
                  call vecn(xo,yo,1)
                  if(yo.gt.ytop(ix))then
                     call vecn(x,ytop(ix),1)
c                        call trihere(xo,yo,.005)
                  endif
                  if(yo.lt.ybot(ix))then
                     call vecn(x,ybot(ix),1)
c                        call trihere(xo,yo,-.005)
                  endif
c                     write(*,*)'Finish seg',ix,ix1,signd,x,y
c     $                    ,dx,dydx,xo,yo
               endif
            elseif(nstate.eq.1.and.lmidl)then
c Drawing above. If previously was not, draw start.
               if(istate.eq.0) then
                  ixms=min(max(ix-signd,1),ngrid)
                  xp=dx*ixms
                  yp=ytop(ixms)
                  call vecn(xp,yp,0)
                  call vecn(x,y,1)
               endif
            elseif(nstate.eq.2.and.lmidl)then
c Drawing below. If previously was not, draw start.
               if(istate.eq.0) then
                  ixms=min(max(ix-signd,1),ngrid)
                  xp=dx*ixms
                  yp=ybot(ixms)
                  call vecn(xp,yp,0)
                  call vecn(x,y,1)
               endif
            endif
            istate=nstate
            xo=x
            yo=y
c Set lmidl true only if we have been through this loop at least once. 
c Otherwise, starting beyond already draw xtop/bot can give spurious 
c start lines after the first vector. 
            lmidl=.true.
         enddo
c End of vector. Finish it.
         xp=x2
         yp=y2
c         call vecn(xp,yp,istate)
         call vecn(xp,yp,min(istate,1))
c         if(icount.eq.0)call trihere(xp,yp,.005)
      else
c         write(*,*)'hid xpoints coincide'
      endif
      x1=x2
      y1=y2
      end
c***************************************************************************
      subroutine trihere(xp,yp,d)
      call vecn(xp+d,yp,1)
      call vecn(xp+d,yp+d,1)
      call vecn(xp,yp,1)
      end
c***************************************************************************
      subroutine hidinit(top,bot)
      real top,bot
      include 'hidcom.h'
      integer i
      do 1 i=1,ngrid
         ytop(i)=top
         ybot(i)=bot
    1 continue
      end
c***************************************************************************
      subroutine cont3proj(z,pp,iLx,nx,ny,cl,icl,x,y,isw,f)
      real f
      real x(*),y(*)
c Project a contour plot onto a z-plane in 3-D plot, at normalized height f.
c Arguments are those of contourl.
      include 'world3.h'

      incolor=igetcolor()
      icsw=abs(isw-(isw/16)*16)
      call axregion(-scbx3,scbx3,-scby3,scby3)
      if(icsw.eq.0)then
c x, y are unused by contourl.
         call scalewn(1.,float(nx),1.,float(ny),.false.,.false.)
      elseif(icsw.eq.1)then
c x and y are really vectors.
         call scalewn(x(1),x(nx),y(1),y(ny),.false.,.false.)
      elseif(icsw.eq.2)then
c x and y are arrays.
c         write(*,*)x(1),x(nx),y(1),y(1+iLx*(ny-1)),iLx,ny
         call scalewn(x(1),x(nx),y(1),y(1+iLx*(ny-1)),.false.,.false.)
      endif
      call hdprset(-3,f*scbz3)
c       Contour without labels, direct on mesh, coloring, smooth.
      call contourl(z,pp,iLx,nx,ny,cl,icl,x,y,isw)
      call color(incolor)

c This finish is left for the caller so annotations can be added.
c      call hdprset(0,0.)
      end
