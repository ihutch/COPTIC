c*******************************************************************
c General slicing routine with web projection.
      subroutine sliceGweb(ifull,iuds,u,nw,zp,ixnp,xn,idfix,utitle)
c Plot web-projected and/or projected contour representations 
c of quantity u on slices with
c fixed values of dimension idfix.
c The full dimensions of arrays u are
      parameter(ndims=3,nd2=2*ndims)
      integer ifull(ndims)
      real u(ifull(1),ifull(2),ifull(3))
c The used dimensions of each are
      integer iuds(ndims)
c The dimensions of square working array, zp, (nw) must be larger
c than the greatest of iuds 
      integer nw
      real zp(nw,nw)
c The node positions are given by vectors laminated into xn, whose
c starts for dimensions id are at ixnp(id)+1. 
c So the vector of positions for dimension id is 
c   ixnp(id)+1 to ixnp(id)+iuds(id) [and ixnp(id+1)-ixnp(id)>=iuds(id)]
      integer ixnp(ndims+1)
      real xn(*)
c The fixed dimension which is chosen to start, and returned after, is:
      integer idfix
c If idfix<0, then don't maintain aspect ratio.
c The plotted quantity title is
      character*(*) utitle
c Needed for perspective plot
      include 'world3.h'
c Workspace size is problematic.
      parameter (nwksp=40000)
      character*1 pp(nwksp)
c Contour levels
      real cl(30)
c Local variables:
      integer icontour,iweb,iback
      integer isw,jsw
      integer iclipping
      character*(10) cxlab,cylab
      character*(30) form1
      save nf1,nf2,nff,if1,if2,iff
      logical lfirst,laspect
      data lfirst/.true./laspect/.true./
      data iclipping/0/jsw/0/
c Tell that we are looking from the top by default.
      data ze1/1./icontour/1/iweb/1/iback/0/
      data cs/.707/sn/.707/
      save n1

      if(idfix.lt.0)then
         laspect=.false.
         idfix=abs(idfix)
      endif
      if(idfix.lt.1 .or. idfix.gt.ndims)idfix=ndims
      ips=0
      irotating=0
c Initial slice number
      call accisgradinit(64000,0,0,-64000,128000,64000)
      if(.not.lfirst)goto 19
         n1=iuds(idfix)/2
c     Plot the surface. With scaling 1. Web color 6, axis color 7.
         jsw=1 + 256*6 + 256*256*7
         iweb=1
         icontour=1
         iclipping=0
         lfirst=.false.
 20      write(*,*)' ======== Slice plotting interface:',
     $        '  arrows up/down: change slice.'
         write(*,*)
     $        ' arrows l/r: change dimension.',
     $        ' s: rescale. p: print. Drag mouse to rotate.'
         write(*,*)' i/o: move eye in/out.'
     $        ,' (jkl;) (m,./): control plotting extent in 2 axes.'
         write(*,*)
     $        ' c: contour plane position. w: toggle web. r/e: rotate.'
     $        ,' a: toggle aspect'
         write(*,*)
     $        ' t: toggle truncation.'
         write(*,*)
     $        ' d: disable interface; run continuously.',
     $        ' depress f: to interrupt running.'
 19   continue
c Start of controlled plotting loop.
 21   call pltinit(0.,1.,0.,1.)
c Set the plotting arrays for fixed dimension idfix.
      idp1=mod(idfix,3)+1
      idp2=mod(idfix+1,3)+1
      if(iclipping.eq.0)then
c Plot the full used array.
         nf1=iuds(idp1)
         nf2=iuds(idp2)
         nff=iuds(idfix)
         if1=1
         if2=1
         iff=1
         if(nf2*nf1.gt.nwksp)then
            write(*,101)nwksp,nf1,nf2
 101        format('sliceGweb error: need bigger nwksp',i6,
     $           ' smaller than',i4,' x',i4)
            return
         endif
      endif
      xdp1=xn(ixnp(idp1)+nf1)-xn(ixnp(idp1)+if1)
      xdp2=xn(ixnp(idp2)+nf2)-xn(ixnp(idp2)+if2)
c      write(*,*)'nf2,if2,xdp2',nf2,if2,xdp2
c Only works for 3-D in present implementation.
c Could be fixed to be general, I suppose.
      do i=if1,nf1
         do j=if2,nf2
            if(idfix.eq.1)then
               zp(i,j)=u(n1,i,j)
            elseif(idfix.eq.2)then
               zp(i,j)=u(j,n1,i)
            elseif(idfix.eq.3)then
               zp(i,j)=u(i,j,n1)
            endif
         enddo
      enddo

c 3D plot ranges.
      xmin=xn(ixnp(idp1)+if1)
      xmax=xn(ixnp(idp1)+nf1)
      ymin=xn(ixnp(idp2)+if2)
      ymax=xn(ixnp(idp2)+nf2)
      zmin=xn(ixnp(idfix)+iff)
      zmax=xn(ixnp(idfix)+nff)
c Web drawing. First call is needed to set scaling.
c      write(*,*)'jsw=',jsw
      if(iweb.ne.0)then
c Old buggy setting, only works for centered cube.
         if(laspect)then
            if(xdp2.gt.xdp1)then
               yc=.3
               xc=xdp1*yc/xdp2
            elseif(xdp2.lt.xdp1)then
               xc=.3
               yc=xdp2*xc/xdp1
            else
               xc=.25
               yc=.25
            endif
            zc=.2
            call setcube(xc,yc,zc,.5,.4)
         endif
c Rescale x and y (if necessary), but not z.
c         if(iclipping.ne.0)
         call scale3(xmin,xmax,ymin,ymax,wz3min,wz3max)
         call hidweb(xn(ixnp(idp1)+if1),xn(ixnp(idp2)+if2),
     $        zp(if1,if2),nw,nf1+1-if1,nf2+1-if2,jsw)
      endif
c Use this scaling until explicitly reset.
      jsw=0 + 256*6 + 256*256*7
      write(form1,'(''Dimension '',i1,'' Plane'',i4)')idfix,n1
      call drwstr(.1,.02,form1)
      call iwrite(idp1,iwidth,cxlab)
      call iwrite(idp2,iwidth,cylab)
      call ax3labels('axis-'//cxlab,'axis-'//cylab,utitle)

c Projected contouring.
      if(icontour.ne.0)then
c       Draw a contour plot in perspective. Need to reset color anyway.
         call color(4)
         call axregion(-scbx3,scbx3,-scby3,scby3)
         call scalewn(xmin,xmax,ymin,ymax,.false.,.false.)
c Calculate place of plane. 
         zplane=scbz3*(-1+(xn(ixnp(idfix)+n1)-zmin)*2./(zmax-zmin))
c accis perspective corner for axes and cube.
         icorner=igetcorner()
         if(iweb.eq.0)then
c Draw axes.
            call hdprset(0,0.)
c Ought to rescale the z-axis, but that was done in hidweb.
            call scale3(xmin,xmax,ymin,ymax,zmin,zmax)
c If we do, then we must reset jsw:
            jsw=1 + 256*6 + 256*256*7
            call axproj(icorner)
         else
c Set contour levels using the scaling of the box.
            icl=6
            do ic=1,icl
               cl(ic)=wz3min+(wz3max-wz3min)*(ic-1.)/(icl-1.)
            enddo
         endif
c Get back current eye position xe1 etc.
         call trn32(xe,ye,ze,xe1,ye1,ze1,-1)
         if(icontour.eq.1)call hdprset(-3,sign(scbz3,ze1))
         if(icontour.eq.2)call hdprset(-3,zplane)
         if(icontour.eq.3)call hdprset(-3,-sign(scbz3,ze1))
c Contour without labels, with coloring, using vector axes
         call contourl(zp(if1,if2),pp,nw,
     $        nf1+1-if1,nf2+1-if2,cl,icl,
     $        xn(ixnp(idp1)+if1),xn(ixnp(idp2)+if2),17)
         call ticlabtog()
         call axis()
         call ticlabtog()
         call axis2()
         if(iweb.ne.1)call cubed(icorner-8*(icorner/8))
      endif

      if(iweb.eq.1.and.icontour.eq.3)then
         call hidweb(xn(ixnp(idp1)+if1),xn(ixnp(idp2)+if2),
     $        zp(if1,if2),nw,nf1+1-if1,nf2+1-if2,jsw)
c This was necessary when hidweb used to change jsw.
c         jsw=0 + 256*6 + 256*256*7
      endif

      if(ips.ne.0)then
c We called for a local print of plot. Terminate and switch it off.
         call pltend()
         call pfset(0)
         ips=0
      endif

c-----------------------------------
c Double buffering switch back. Disabled here because of the compiz
c screen-writing performance issue.
      if(iback.ne.0)then 
c         call glfront()
         call accisrefresh()
      endif
c Limit framing rate to 30fps.
      call usleep(15000)
c User interface interpret key-press.
 24   call eye3d(isw)
c      write(*,*)'isw',isw
      if(isw.eq.ichar('f')) goto 24
      if(isw.eq.0) goto 23
      if(isw.eq.65364 .and. n1.gt.1) n1=n1-1
      if(isw.eq.65362 .and. n1.lt.iuds(idfix)) n1=n1+1
      if(isw.eq.ichar('q')) goto 23
      if(isw.eq.ichar('a')) laspect=.not.laspect
      if(isw.eq.ichar('d')) call noeye3d(0)
      if(isw.eq.ichar('s')) jsw=1 + 256*6 + 256*256*7
      if(isw.eq.ichar('t')) call togi3trunc()
      if(isw.eq.ichar('p'))then
         call pfset(3)
         ips=3
      endif
c Change fixed dimension, remove clipping, force scaling.
      if(isw.eq.65361)then
         idfix=mod(idfix+1,3)+1
         iclipping=0
         n1=iuds(idfix)/2
         jsw=1 + 256*6 + 256*256*7
      elseif(isw.eq.65363)then
         idfix=mod(idfix,3)+1
         iclipping=0
         n1=iuds(idfix)/2
         jsw=1 + 256*6 + 256*256*7
      endif
c Adjust clipping
      if(isw.eq.47)then
c /
         iclipping=1
         nf2=min(nf2+1,iuds(idp2))
      elseif(isw.eq.46)then
c .
         iclipping=1
         nf2=max(nf2-1,if2+1)
      elseif(isw.eq.44)then
c ,
         iclipping=1
         if2=min(if2+1,nf2-1)
      elseif(isw.eq.ichar('m'))then
c m
         iclipping=1
         if2=max(if2-1,1)
      elseif(isw.eq.ichar('l'))then
c l
         iclipping=1
         nf1=max(nf1-1,if1+1)
      elseif(isw.eq.59)then
c ;
         iclipping=1
         nf1=min(nf1+1,iuds(idp1))
      elseif(isw.eq.ichar('j'))then
         iclipping=1
         if1=max(if1-1,1)
      elseif(isw.eq.ichar('k'))then
         iclipping=1
         if1=min(if1+1,nf1-1)
      endif
      if(isw.eq.ichar('h'))goto 20
      if(isw.eq.ichar('c'))icontour=mod(icontour+1,4)
      if(isw.eq.ichar('w'))iweb=mod(iweb+1,2)
      if(isw.eq.ichar('i'))then
c Get back current eye position xe1 etc.
         call trn32(xe,ye,ze,xe1,ye1,ze1,-1)
         xe1=xe1*.9
         ye1=ye1*.9
         ze1=ze1*.9
c Move it in.
         call trn32(0.,0.,0.,xe1,ye1,ze1,1)
      elseif(isw.eq.ichar('o'))then
c Get back current eye position xe1 etc.
         call trn32(xe,ye,ze,xe1,ye1,ze1,-1)
         xe1=xe1*1.1
         ye1=ye1*1.1
         ze1=ze1*1.1
c Move it out.
         call trn32(0.,0.,0.,xe1,ye1,ze1,1)
      elseif(isw.eq.ichar('r'))then
         irotating=1
         cs=cos(.05)
         sn=sin(.05)
      elseif(isw.eq.ichar('e'))then
         irotating=1
         cs=cos(.05)
         sn=sin(-.05)
      endif
      if(irotating.gt.0)then
c         call glback()
c         iback=1
c Get back current eye position xe1 etc.
         call trn32(xe,ye,ze,xe1,ye1,ze1,-1)
c         write(*,*)'irotating',irotating,xe1,ye1,ze1,cs,sn
         xex=xe1-xe
         yex=ye1-ye
         xe1=xe+cs*xex-sn*yex
         ye1=ye+sn*xex+cs*yex
c Must tell to look at zero.
         call trn32(0.,0.,0.,xe1,ye1,ze1,1)
         irotating=irotating-1
         call trn32(xe,ye,ze,xe1,ye1,ze1,-1)
      endif
c End of user interface.
c------------------------------------
      goto 21
 23   continue
      call hdprset(0,0.)
      end
c******************************************************************
c*******************************************************************
c Three-way slice contours.
      subroutine sliceGcont(ifull,iuds,u,nw,zp,ixnp,xn,ifixpt,utitle)
c Plot projected contour representations on cuts in three dimensions 
c of quantity u on fixed values of dimensions.
c The initial intersection of the fixed planes is at ifixpt(ndims).
c The full dimensions of arrays u are
      parameter(ndims=3,nd2=2*ndims)
      integer ifull(ndims)
      real u(ifull(1),ifull(2),ifull(3))
c The used dimensions of each are
      integer iuds(ndims)
c The dimensions of square working array, zp, (nw) must be larger
c than the greatest of iuds 
      integer nw
      real zp(nw,nw,ndims)
c The node positions are given by vectors laminated into xn, whose
c starts for dimensions id are at ixnp(id)+1. 
c So the vector of positions for dimension id is 
c   ixnp(id)+1 to ixnp(id)+iuds(id) [and ixnp(id+1)-ixnp(id)>=iuds(id)]
      integer ixnp(ndims)
      real xn(*)
c The fixed dimension which is chosen to start, and returned after, is:
      integer ifixpt(ndims)
c The plotted quantity title is
      character*(*) utitle
c Needed for perspective plot
      include 'world3.h'
c Workspace size is problematic.
      parameter (nwksp=40000)
      character*1 pp(nwksp)
c Contour levels
      real cl(30)
c Make pfsw visible:
      include 'plotcom.h'
c Local variables:
      integer idire(3),idir(3)
      real cbsize(3),cmin(3),cmax(3)
      integer icontour,iweb,iback
      integer isw,jsw
      character*(10) cxlab,cylab
      character*(30) form1
      logical lfirst
      data lfirst/.true./
c Tell that we are looking from the top by default.
      data ze1/1./icontour/1/iweb/1/iback/0/
      data cs/.707/sn/.707/

      imv=1
      itri=0
      icl=0
      if(ifixpt(1).eq.0)ifixpt(1)=iuds(1)/2
      if(ifixpt(2).eq.0)ifixpt(2)=iuds(2)/2
      if(ifixpt(3).eq.0)ifixpt(3)=iuds(3)/2
      ips=0
      irotating=0
      call minmax2(u(1,1,ifixpt(3)),ifull(1),iuds(1),iuds(2),umin,umax)
      call accisgradinit(64000,0,0,-64000,128000,64000)
      if(.not.lfirst)goto 19
c Have to use goto so I can jump into the middle at 20.
c     Plot the surface. With scaling 1. Web color 6, axis color 7.
         jsw=1 + 256*6 + 256*256*7
         iweb=1
         icontour=1
         lfirst=.false.
 20      write(*,*)' ======== Slice contouring interface:',
     $        ' up/down, l/r arrows change slice.'
         write(*,*) ' d: changes l/r dimension controlled.',
     $        ' s:rescale. p:print. Drag to move view.'
         write(*,*)
     $        ' r/e: rotate. toggles: t smooth color; l labels'
 19   continue

 21   call pltinit(0.,1.,0.,1.)
      call setcube(.2,.2,.2,.5,.4)
      do id=1,3
         cmin(id)=xn(ixnp(id)+1)
         cmax(id)=xn(ixnp(id)+iuds(id))
      enddo
      call scale3(cmin(1),cmax(1),cmin(2),cmax(2),cmin(3),cmax(3))

c Set contour levels
      if(icl.eq.0)icl=10
      call fitrange(umin,umax,abs(icl),nxfac,xfac,xtic,vmin,vmax)
c      write(*,*)umin,umax,vmin,vmax,xtic
      do ic=1,abs(icl)
         cl(ic)=vmin+(ic-1)*xtic
      enddo
c      write(*,*)(cl(i),i=1,abs(icl))
c Get back current eye position xe1 etc.
      call trn32(xe,ye,ze,xe1,ye1,ze1,-1)

      idire(1)=int(sign(1.,xe1-xe))
      idire(2)=int(sign(1.,ye1-ye))
      idire(3)=int(sign(1.,ze1-ze))
      cbsize(1)=scbx3
      cbsize(2)=scby3
      cbsize(3)=scbz3

      do i=1,iuds(2)
         do j=1,iuds(3)
            zp(i,j,1)=u(ifixpt(1),i,j)
         enddo
      enddo
      do i=1,iuds(3)
         do j=1,iuds(1)
            zp(i,j,2)=u(j,ifixpt(2),i)
         enddo
      enddo
      do i=1,iuds(1)
         do j=1,iuds(2)
            zp(i,j,3)=u(i,j,ifixpt(3))
         enddo
      enddo

c This needs to be generalized to accommodate multiple levels.
c There are for each dimension three points: 1, ifixpt(id), iuds(id).
c The position of the plane is specified in 3-d normal units. 
c Each contour is called with an origin and iud in the other dims.
c Because of memory, we need to make the the origin the smaller index.
c Here we need to cycle through 4 levels of plane from back to front.

      do ilev=1,4
         if(ilev.eq.1)then
            do i=1,3
               idir(i)=-idire(i)
            enddo
         else
            do i=1,3
               idir(i)=idire(i)
            enddo
         endif
         do ifix=1,3
            id1=mod(ifix,3)+1
            id2=mod(ifix+1,3)+1
            if(ilev.eq.2)idir(id1)=-idir(id1)
            if(ilev.eq.3)idir(id2)=-idir(id2)
            call axregion(-cbsize(id1),cbsize(id1),
     $           -cbsize(id2),cbsize(id2))
            call scalewn(cmin(id1),cmax(id1),cmin(id2),cmax(id2),
     $           .false.,.false.)
            zplane=cbsize(ifix)*
     $           (-1+(xn(ixnp(ifix)+ifixpt(ifix))-cmin(ifix))*2.
     $           /(cmax(ifix)-cmin(ifix)))
            call hdprset(-ifix,zplane)
c Poor man's top lighting:            
            if(ifix.eq.3)then
               ib=10000
            elseif(ifix.eq.2)then
               ib=-5000
            else
               ib=0
            endif
            call accisgradinit(64000+ib,0+ib,0+ib,
     $           -64000+ib,128000+ib,64000+ib)
            it1=iuds(id1)*(idir(id1)+1)/2 - ifixpt(id1)*(idir(id1)-1)/2
            ib1=ifixpt(id1)*(idir(id1)+1)/2 - (idir(id1)-1)/2
            it2=iuds(id2)*(idir(id2)+1)/2 - ifixpt(id2)*(idir(id2)-1)/2
            ib2=ifixpt(id2)*(idir(id2)+1)/2 - (idir(id2)-1)/2
c Contour without labels, with coloring, using vector axes
c            write(*,*)it1,ib1,it2,ib2,ifix
            if(it1-ib1.gt.0 .and. it2-ib2.gt.0)
     $           call contourl(zp(ib1,ib2,ifix),
     $           pp,nw,it1-ib1+1,it2-ib2+1,
     $           cl,icl,
     $           xn(ixnp(id1)+ib1),xn(ixnp(id2)+ib2),17+itri*64)
            if(ilev.eq.2)idir(id1)=-idir(id1)
            if(ilev.eq.3)idir(id2)=-idir(id2)
         enddo
      enddo
      call axproj(igetcorner())
      call ax3labels('x','y','z')

      if(ips.ne.0)then
c We called for a local print of plot. Terminate and switch it off.
c Prevent pltend from querying the interface.
         pfsw=-ips
         call pltend()
         call pfset(0)
         ips=0
      endif

c-----------------------------------
c Double buffering switch back. Disabled here because of the compiz
c screen-writing performance issue.
      if(iback.ne.0)then 
c         call glfront()
         call accisrefresh()
      endif
c Limit framing rate to 30fps.
      call usleep(10000)
c User interface interpret key-press.
      call eye3d(isw)
      if(isw.eq.0) goto 23
      if(isw.eq.65364)then
         ifixpt(3)=ifixpt(3)-1
         if(ifixpt(3).lt.1)ifixpt(3)=1
      endif
      if(isw.eq.65362)then
         ifixpt(3)=ifixpt(3)+1
         if(ifixpt(3).gt.iuds(3))ifixpt(3)=iuds(3)
      endif
      if(isw.eq.ichar('q')) goto 23
      if(isw.eq.ichar('d')) imv=mod(imv,2)+1
      if(isw.eq.ichar('s')) call minmax2(u(1,1,ifixpt(3)),
     $     ifull(1),iuds(1),iuds(2),umin,umax)
      if(isw.eq.ichar('p'))then
         call pfset(3)
         ips=3
      endif
      if(isw.eq.65361)then
         ifixpt(imv)=ifixpt(imv)+1*idire(imv)
         if(ifixpt(imv).lt.1)ifixpt(imv)=1
         if(ifixpt(imv).gt.iuds(imv))ifixpt(imv)=iuds(imv)
      elseif(isw.eq.65363)then
         ifixpt(imv)=ifixpt(imv)-1*idire(imv)
         if(ifixpt(imv).lt.1)ifixpt(imv)=1
         if(ifixpt(imv).gt.iuds(imv))ifixpt(imv)=iuds(imv)
      endif
      if(isw.eq.ichar('t'))itri=mod(itri+1,2)
      if(isw.eq.ichar('l'))icl=-icl
      if(isw.eq.ichar('h'))goto 20
      if(isw.eq.ichar('c'))icontour=mod(icontour+1,4)
      if(isw.eq.ichar('w'))iweb=mod(iweb+1,2)
      if(isw.eq.ichar('i'))then
c Get back current eye position xe1 etc.
         call trn32(xe,ye,ze,xe1,ye1,ze1,-1)
         xe1=xe1*.9
         ye1=ye1*.9
         ze1=ze1*.9
c Move it in.
         call trn32(0.,0.,0.,xe1,ye1,ze1,1)
      elseif(isw.eq.ichar('o'))then
c Get back current eye position xe1 etc.
         call trn32(xe,ye,ze,xe1,ye1,ze1,-1)
         xe1=xe1*1.1
         ye1=ye1*1.1
         ze1=ze1*1.1
c Move it out.
         call trn32(0.,0.,0.,xe1,ye1,ze1,1)
      elseif(isw.eq.ichar('r'))then
         irotating=1
         cs=cos(.05)
         sn=sin(.05)
      elseif(isw.eq.ichar('e'))then
         irotating=1
         cs=cos(.05)
         sn=sin(-.05)
      endif
      if(irotating.gt.0)then
c         call glback()
c         iback=1
c Get back current eye position xe1 etc.
         call trn32(xe,ye,ze,xe1,ye1,ze1,-1)
c         write(*,*)'irotating',irotating,xe1,ye1,ze1,cs,sn
         xex=xe1-xe
         yex=ye1-ye
         xe1=xe+cs*xex-sn*yex
         ye1=ye+sn*xex+cs*yex
c Must tell to look at zero.
         call trn32(0.,0.,0.,xe1,ye1,ze1,1)
         irotating=irotating-1
      endif
c End of user interface.
c------------------------------------
      goto 21
 23   continue
      call hdprset(0,0.)
      end
c******************************************************************
