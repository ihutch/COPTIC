      subroutine hidweb(x,y,z,iLx,nx,ny,ilevel)
c Draw a 3-d web of z(x,y), dim(iLx\nx,ny), viewed from (xe,ye,ze) scaled.
c Second byte of level gives web color. Third gives axis color.
c Lowest byte: 
c   abs(level)=1 scale to fit the region, 1-d x,y vectors.
c   abs(level)=2 scale to fit region, 2-d x,y, don't hide, just wiremesh.
c   abs(level)=0 perform no scale-setting and use last perspective...
c   level.lt.0 draw no axes.
c Eye obtained from file eye.dat, or default if none exists.
      integer iLx, nx,ny,ilevel,level
      real x(*),y(*),z(iLx,ny)
      integer icorner,colw,cola
      real x2,y2,z2
      real zmin,zmax
      integer ihd
      data ihd/99/
      save

      level=ilevel
      cola=level/256
      level=level-cola*256
c level is lowest 8 bits.
      colw=cola-(cola/256)*256
c color of web is next 8 bits
      cola=cola/256 -(cola/65536)*256
c color of axes is next 8 bits.
c      write(*,*)'level=',level
      if(abs(level).ne.0)then
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
c Set the scaling.	 call minmax2(z,iLx,nx,ny,zmin,zmax)
	 call minmax2(x,iLx,nx,ny,xmin,xmax)
	 call minmax2(y,iLx,nx,ny,ymin,ymax)
      endif
c Set the top and bottom horizons.
      call hidinit(0.,1.)
c Set to non-hiding (ihd=0) or hiding (ihd=1) 3-D calls.
      call hdprset(ihd,0.)
      if(abs(level).ne.0)then
c Set the perspective transform.
         call trn32(0.,0.,0.,x2,y2,z2,1)
         itics=5
         call minmax2(z,iLx,nx,ny,zmin,zmax)
         call fitrange(zmin,zmax,itics,ipow,fac10,delta,first,xlast)
         zmin=first
         zmax=xlast
         call scale3(xmin,xmax,ymin,ymax,zmin,zmax)
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
	 icube=sign(1.,z2)*icorner
	 call cubed(icube)
c Get corner with extra switches. (Could pass this to cubed too).
         icorner=igetcorner()
c Draw axes.
	 call axproj(icorner)
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
    2	 continue
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
c	 call trn32(x(kx(id1)),y(ky(id1)),0.,x2,y2,z2,0)
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
    2	 continue
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
                  xp=dx*(ix-signd)
                  yp=ytop(ix-signd)
c                  if(ytop(ix-signd).eq.0.)
c     $                 write(*,*)ix1,ix2,ix,signd,x1,x2,x,xp
                  call vecn(xp,yp,0)
                  call vecn(x,y,1)
               endif
            elseif(nstate.eq.2.and.lmidl)then
c Drawing below. If previously was not, draw start.
               if(istate.eq.0) then
                  xp=dx*(ix-signd)
                  yp=ybot(ix-signd)
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
         call vecn(xp,yp,istate)
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
c*****************************************************************************
      subroutine surf3d(x,y,z,iLx,nx,ny,isw,work)
c Draw a 3-d surface of z(x,y), dim(iLx\nx,ny), viewed from (xe,ye,ze) scaled.
c This version designed for non-monotonic x,y, does sophisticated web order.
c work(0:iLx+1,0:ny+1) is a work array that is trashed.
c Switch: isw.lt.0 draw no axes.
c Second byte of isw gives web color. Third byte gives axis color.
c Lowest nibble levd (absolute value): 
c   =0 perform no scaling and use last perspective
c   =1 scale to fit the region, 1-d x,y
c   =2 scale to fit the region, 2-d x,y
c Second nibble (16+): isw for the surfdr3 call. bit0=1 => smooth tri.
c Eye obtained from file eye.dat, or default if none exists.
c Obsolete:
c   abs(isw)=0, 1 scale to fit the region, 1-d x,y.
c   abs(isw)=2 perform no scaling and use last perspective.
c   abs(isw)=4 scale to fit region 2-d x,y.
      integer iLx, nx,ny,isw
      real x(iLx,*),y(iLx,*),z(iLx,ny),work(0:iLx+1,0:ny+1)
      integer icorner,colw,cola
      real x2,y2,z2
      real zmin,zmax
      save

      cola=isw/256
      levd=abs(isw-(isw/16)*16)
c levd is lowest 4 bits.
      colw=abs(cola-(cola/256)*256)
c color of web is next 8 bits
      cola=abs(cola/256 -(cola/65536)*256)
c color of axes is next 8 bits.
      isw=isw/16 - (isw/256)*16
c isw is second nibble (bits 4-7).
      if(colw.gt.15)then
         call color(colw)
      else
         if(colw.gt.0) call color(colw)
      endif
c Set the top and bottom horizons.
      call hidinit(0.,1.)
      if(levd.ne.0)then
         call geteye(x2,y2,z2)
c     Set to non-hiding 3-D calls.
         call hdprset(0,0.)
c     Set the perspective transform.
         call trn32(0.,0.,0.,x2,y2,z2,1)
c     Set the scaling.
            call minmax2(z,iLx,nx,ny,zmin,zmax)
            itics=5
            call fitrange(zmin,zmax,itics,ipow,fac10,delta,first,xlast)
            zmin=first
            zmax=xlast
         if(levd.eq.2)then
            call minmax2(x,iLx,nx,ny,xmin,xmax)
            call minmax2(y,iLx,nx,ny,ymin,ymax)
         elseif(levd.eq.1)then
            call minmax(x,nx,xmin,xmax)
            call minmax(y,ny,ymin,ymax)            
         endif
         call scale3(xmin,xmax,ymin,ymax,zmin,zmax)
      endif
c Draw the surface.
      call surfdr3(x,y,z,iLx,nx,ny,work,isw,dummy)
      if(cola.ne.0) call color(cola)
      if(isw.ge.0)then
c Draw cube.
         icorner=igetcorner()
	 icube=sign(1.,z2)*icorner
c         write(*,*)icorner,x2,y2,z2
c This had a problem with extra vertical line. Ought to be fixed.
c	 call cubed(icube)
c Draw axes.
         call ax3labels('x-axis','y-axis','')
	 call axproj(icorner)
      endif
      end
c********************************************************************
c Return the nearest corner to eye in standard convention.
      function igetcorner()
c This was wrong. geteye reads the eyefile.
c      call geteye(x2,y2,z2)
      call trn32(xdum,ydum,zdum,x2,y2,z2,-1)
      if(y2.le.0.)then 
         if(x2.le.0.)then
            icorner=1
         else
            icorner=2
         endif
      else
         if(x2.le.0.)then
            icorner=4
         else
            icorner=3
         endif
      endif
c      write(*,*)'x2,y2,z2',x2,y2,z2,mod(icorner,2)
c If appropriate tell axproj to flip labels.
      if(abs(x2).gt.abs(y2).eqv.(mod(icorner,2).ne.0))
     $	   icorner=icorner+16
c If appropriate use vertical labels
c	 if(z2*z2.lt.x2*x2) icorner=icorner+32
c Trying for better results.
      xy2i=min(x2**2,y2**2)
      if(z2*z2.lt.xy2i+.2*y2**2) icorner=icorner+32
      if(z2*z2.lt.xy2i+.2*x2**2) icorner=icorner+64
      igetcorner=icorner
      end
c********************************************************************
c Attempt was made to create a webdrw that drew in the order of closest
c to furthest. However, hiding does not in the end work correctly for 
c that. So the surf3d was done instead.
c********************************************************************
      subroutine surfdr3(x,y,z,iLx,nx,ny,work,isw,d)
c Draw a 3-d surface of z(x,y), dim(iLx\nx,ny), using current scaling.
c Version using filled quadrilaterals.
c isw bit 0 set => triangular fills, else chunks. If in addition,
c isw bit 1 set => directional shading, optional argument d(5)
c    d(1-3) gives direction d(4-5) gives distance limits.
      integer nx,ny,iLx
      real x(iLx,ny),y(iLx,ny),z(iLx,ny),work(0:iLx+1,0:ny+1)
      integer isw
      real d(5)

      integer id1,id2,ud,kx,ky,icorner
      real d1start,d1end,d1step,z1,x2,y2,z2
      real dout,ddone
      parameter (dout=2.e32,ddone=1.e32)
      real xp(5),yp(5),zp(5)
      integer icx(0:4), icy(0:4), ic
      logical lchunk,ldir
      data icx/0,1,1,0,0/
      data icy/0,0,1,1,0/
      save
c
c      write(*,*)'isw',isw
      lchunk=.true.
      if(isw-isw/2 .ne.0)lchunk=.false.
      ldir=.false.
      if(isw/2-isw/4 .ne.0)ldir=.true.
      incolor=igetcolor()
c Get eye position x0,y0,z0 are unused here.
      call trn32(x0,y0,z0,x2,y2,z2,-1)
c      write(*,*)'Eye=',x2,y2,z2
c Set up the work array.
      zmin=ddone
      zmax=-ddone
      do j=1,ny-1
         do i=1,nx-1
            if(z(i,j).gt.zmax)zmax=z(i,j)
            if(z(i,j).lt.zmin)zmin=z(i,j)
            work(i,j)=0.
            do ic=0,3
               work(i,j)=work(i,j)
     $              +x(i+icx(ic),j+icy(ic))*x2
     $              +y(i+icx(ic),j+icy(ic))*y2
     $              +z(i+icx(ic),j+icy(ic))*z2

            enddo
         enddo
      enddo
      if(ldir)then
         zmin=d(4)
         zmax=d(5)
         if(.not.abs(zmax-zmin).gt.0)stop 'ERROR: surf3d d-limits'
      endif

      do k=1,(nx-1)*(ny-1)
c Search for the face not yet treated, furthest from eye
         dmin=1.e35
         do j=1,ny-1
            do i=1,nx-1
               if(work(i,j).lt.dmin)then
                  imin=i
                  jmin=j
                  dmin=work(i,j)
               endif
            enddo
         enddo
c         write(*,*)'imin,jmin,work,x,y',imin,jmin,work(imin,jmin),
c     $        x(imin,jmin),y(imin,jmin)
         do ic=0,4
            xp(1+ic)=x(imin+icx(ic),jmin+icy(ic))
            yp(1+ic)=y(imin+icx(ic),jmin+icy(ic))
            zp(1+ic)=z(imin+icx(ic),jmin+icy(ic))
         enddo
c Calculate height of centroid.
         if(lchunk)then
            zcol=0.
            do ic=0,3
               zcol=zcol+z(imin+icx(ic),jmin+icy(ic))
            enddo
            zcol=zcol/4.
            icolor=(240*(zcol-zmin))/(zmax-zmin)
            call gradcolor(icolor)
            call poly3line(xp,yp,zp,5)
            call pathfill()
         else
            if(ldir)then
               call gradquad(xp,yp,zp,d,zmin,zmax,0,239,3)
            else
               call gradquad(xp,yp,zp,zp,zmin,zmax,0,239,1)
            endif
         endif
c Here there's a problem for the very last quadrilateral. 
c All the previous ones in postscript inherit the 0 setlinewidth
c that is called from the gradquad immediately afterwards. But the 
c last one does not have it, so it has a thicker line. One option
c would be to do a stroke before 0 setlinewidth, but that makes all the 
c lines thick (as usual). 
c I'm not sure if that's what I want. Trying it; seems better.
         call color(incolor)
         if(incolor.ne.0)call poly3line(xp,yp,zp,5)
c Set face as done
         work(imin,jmin)=ddone
      enddo
      
      end







