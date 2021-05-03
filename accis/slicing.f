c*******************************************************************
c General slicing routine with web projection.
      subroutine sliceGweb(ifull,iuds,u,nw,zp,ixnp,xn,idfixp,utitle
     $     ,svec,vp)
c Plot web-projected and/or projected contour representations 
c of quantity u on slices with fixed values of dimension idfix.
c The full dimensions of arrays u are
      parameter(ndims=3,nd2=2*ndims)
      integer ifull(ndims)
      real u(ifull(1),ifull(2),ifull(3))
c The used dimensions of each are
      integer iuds(ndims)
c The dimensions of square working array, zp, (nw) must be larger
c than the greatest of iuds. 
c If zp(1,1)=0 on entry, give control help.
      integer nw
!      real zp(nw,nw)
      real zp(nw*nw)
c The node positions are given by vectors laminated into xn, whose
c starts for dimensions id are at ixnp(id)+1. 
c So the vector of positions for dimension id is 
c   ixnp(id)+1 to ixnp(id)+iuds(id) [and ixnp(id+1)-ixnp(id)>=iuds(id)]
      integer ixnp(ndims+1)
      real xn(*)
c The fixed dimension which is chosen to start (extended usage below) is:
      integer idfixp
c The plotted quantity title is
      character*(*) utitle
c If abs(idfix) has bit 2 set (by adding 4), toggle plotting of svec, an
c optional vector argument on the same array of positions as u.
      real svec(ifull(1),ifull(2),ifull(3),3)
      real vp(nw,nw,2)
c These arguments ought to be present at least as dummy reals for all
c calls; otherwise the length of the utitle will not be found correctly.
c Do arrow plots of this field over contours.

c Needed for perspective plot
      include 'world3.h'
c For testing only
      include 'plotcom.h'
c Workspace size is problematic.
      parameter (nwksp=1000000,iclmax=30)
      character*1 pp(nwksp)
c Contour levels
      real cl(iclmax)
      real colorscale
c Local variables:
      integer icontour,iweb
      integer jsw
      integer iclipping,idflast,idfinlast,igradleg,icl
      integer idpa(2),idp1,idp2
      character*(30) form1
      save nf1,nf2,nff,if1,if2,iff,idfix
      logical laspect,larrow,ltellslice,logcontour
      data laspect/.true./larrow/.false./
      data ltellslice/.false./logcontour/.false./
      data iclipping/0/idflast/-9999/jsw/0/n1/0/icontour/1/iweb/1/
      data igradleg/0/
      data idfinlast/-9999/
      data colorscale/1./
      data icl/6/
c Tell that we are looking from the top by default.
      data ze1/1./
      save idp1,idp2

c Need to get more systematic with idfixin. 
c Sign negative  toggle (off) aspect ratio maintenance.
c Bits 0 (1) and 1 (2)  specify the id fixed.
c The rest become idfixf: 
c Bit 2 (4)  toggle (on) plotting of svec arrows.
c Bit 3 (8)  reinitialize.
c Bits 4,5 (16xicontour)  set the initial icontour number 0...3
c Bit 6 (64)  Toggle ltellslice 
c Bit 7 (128)  Use the first 4 values of vp on input for clipping.
c Bit 8 (256) Do no internal scaling initially.
c Bit 9 (512) Turn off contour labelling.

c Bit 10 (1024) Do surface rather than web.     Not implemented.

      ihibyte=idfixp/256
      if(idfixp/512-1024*(idfixp/1024).ne.0)then
c kluge because higher bits break something I don't understand
         iclsign=-1
c         write(*,*)'idfixin512',idfixin
         idfixin=idfixp-512
      else
         iclsign=1
         idfixin=idfixp
      endif

      idfixf=abs(idfixin)/4
      idfix=abs(idfixin)-4*(idfixf)
c Sign
      if(idfixin.lt.0)then
         laspect=.not.laspect
      endif
      if(idfix.eq.0)then
c Bits 0,1 set direction, or reinitialize and use default ndims.
         if(idflast.eq.-9999)then
            idfix=ndims
         else
            idfix=idflast
         endif
      endif
      if(idfixf-2*(idfixf/2).ne.0)then
c Bit 2 Toggle on svec 
         larrow=.not.larrow
c         write(*,*)'idfix=',idfix,'  larrow=',larrow
      endif
      idfixf=idfixf/2
      if(idfixf-2*(idfixf/2).ne.0)then 
c Bit 3=1 reinitialize:
         idfinlast=-9999
      endif
      idfixf=idfixf/2
      ips=0
      irotating=0
      if(.not.idfinlast.eq.idfixin)then
c Initialize
c Bits 4,5 set icontour
         icontour=(idfixf-4*(idfixf/4))
         idfixf=idfixf/4
c Bit 6 toggle ltellslice
         if(idfixf-2*(idfixf/2).ne.0)ltellslice=.not.ltellslice
         n1=(iuds(idfix)+1)/2
         idfixf=idfixf/2
c Bit 7 clipping
         if(idfixf-2*(idfixf/2).ne.0)then
            iclipping=1
c Set the plotting arrays for fixed dimension idfix.
            idp1=mod(idfix,3)+1
            idp2=mod(idfix+1,3)+1
            if1=max(nint(vp(1,1,1)),1)
            nf1=max(min(nint(vp(2,1,1)),iuds(idp1)),if1)
            if2=max(nint(vp(3,1,1)),1)
            nf2=max(min(nint(vp(4,1,1)),iuds(idp2)),if2)
         else
            iclipping=0
         endif
c ----------------------------------------------------------
c     Plot the surface. With scaling 1. Web color 6, axis color 7.
         jsw=1 + 256*6 + 256*256*7
         iweb=1
         write(*,*)' ======== Slice plotting interface. Hit h for help.'
         idfinlast=idfixin
         idflast=idfix
      endif
c Start of controlled plotting loop.
 21   call accisinit
c This gradient call needs to be after accisinit else nodisplay commands
c are broken.
C      call blueredgreenwhite()
      call brgwscaled(0.,colorscale)
      if(iclipping.eq.0)then
c Update the plotting arrays for fixed dimension idfix.
         idp1=mod(idfix,3)+1
         idp2=mod(idfix+1,3)+1
         idpa(1)=idp1
         idpa(2)=idp2
c Plot the full used array.
         nf1=iuds(idp1)
         nf2=iuds(idp2)
         nff=iuds(idfix)
         if1=1
         if2=1
         iff=1
         if(nf2*nf1.gt.nwksp)then
            write(*,101)nwksp,nf1,nf2
 101        format('sliceGweb error: needs bigger nwksp',i6,
     $           ' compiled smaller than',i4,' x',i4)
            write(*,'(a)')'In your accis directory execute the command'
            call iwrite(nf2*nf1,iwidth,form1)
            write(*,*)'sed -i -e "s/nwksp=1000000','[0-9]*/nwksp='
     $           ,form1(1:iwidth),'/" slicing.f'
            write(*,'(a)')'(or larger) and recompile your application.'
            return
         endif
!         if(nf1.gt.nw.or.nf2.gt.nw)then
!            write(*,102)nw,nf1,nf2
! 102        format('sliceGweb WARN: may need bigger zp dimension than'
!     $           ,i6,'^2 for',2i6)
!         endif
      endif
      xdp1=xn(ixnp(idp1)+nf1)-xn(ixnp(idp1)+if1)
      xdp2=xn(ixnp(idp2)+nf2)-xn(ixnp(idp2)+if2)
c      write(*,*)'nf2,if2,xdp2',nf2,if2,xdp2
c Only works for 3-D in present implementation. Copy to work array.
      do j=if2,nf2
         do i=if1,nf1
            if(idfix.eq.1)then
!               zp(i,j)=u(n1,i,j)   ! zp is 2-D
               zp((j-if2)*(nf1-if1+1)+(i-if1+1))=u(n1,i,j) ! zp 1-D
            elseif(idfix.eq.2)then
!               zp(i,j)=u(j,n1,i)
               zp((j-if2)*(nf1-if1+1)+(i-if1+1))=u(j,n1,i) ! zp 1-D
            elseif(idfix.eq.3)then
!               zp(i,j)=u(i,j,n1)
               zp((j-if2)*(nf1-if1+1)+(i-if1+1))=u(i,j,n1) ! zp 1-D
            endif
         enddo
      enddo
      if(larrow)then
c Set up the vector arrow plot arrays. 
         do k=1,2
            do j=if2,nf2
               do i=if1,nf1
                  if(idfix.eq.1)then
                     vp(i,j,k)=svec(n1,i,j,idpa(k))
                  elseif(idfix.eq.2)then
                     vp(i,j,k)=svec(j,n1,i,idpa(k))
                  elseif(idfix.eq.3)then
                     vp(i,j,k)=svec(i,j,n1,idpa(k))
                  endif
               enddo
            enddo
         enddo
      endif
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
!            write(*,*)'xdp1,2',xdp1,xdp2,xc,yc,zc
            call setcube(xc,yc,zc,.5,.4)
         endif
c Rescale x and y (if necessary), but not z.
c         if(iclipping.ne.0)
         call scale3(xmin,xmax,ymin,ymax,wz3min,wz3max)
         if(iweb.lt.2)then
            if(idfixin/256 -512*(idfixin/512).ne.0)then
c This call does no internal initial z-scale setting and scale3 ought to
c have been called in the external program:
               call hidweb(xn(ixnp(idp1)+if1),xn(ixnp(idp2)+if2),
!     $              zp(if1,if2),nw,nf1+1-if1,nf2+1-if2,jsw+8)
     $             zp,nf1+1-if1,nf1+1-if1,nf2+1-if2,jsw+8)  ! zp 1-D
            else
c This is the standard call that normally does internal scaling:
!               write(*,*)'nf1,if1,nf2,if2',nf1,if1,nf2,if2
               call hidweb(xn(ixnp(idp1)+if1),xn(ixnp(idp2)+if2),
!     $              zp(if1,if2),nw,nf1+1-if1,nf2+1-if2,jsw)
     $             zp,nf1+1-if1,nf1+1-if1,nf2+1-if2,jsw)  ! zp 1-D
            endif
            call color(7)
c Use this scaling until explicitly reset.
            jsw=0 + 256*6 + 256*256*7
         else
            isw=1
            if(iweb.eq.3)isw=isw+16
            call surf3d(xn(ixnp(idp1)+if1),xn(ixnp(idp2)+if2)
!     $           ,zp(if1,if2),nw,nf1+1-if1,nf2+1-if2,isw,pp)
     $           ,zp,nf1+1-if1,nf1+1-if1,nf2+1-if2,isw,pp)  ! zp 1-D
            call color(7)
            call axproj(igetcorner())
         endif
      endif

      if(ltellslice)then
         call drwstr(.1,.02
     $        ,ax3chars(idfix)(1:lentrim(ax3chars(idfix))))
         call fwrite(xn(ixnp(idfix)+n1),iwidth,2,form1)
         call drcstr(' '//form1(1:iwidth+1))
      endif
      call ax3labels(ax3chars(idp1)(1:lentrim(ax3chars(idp1)))
     $     ,ax3chars(idp2)(1:lentrim(ax3chars(idp2))),utitle)

c Projected contouring.
      if(mod(icontour,4).ne.0.and.nf1.gt.if1.and.nf2.gt.if2)then
         if(iweb.eq.0)call cubed(igetcubecorner())
c       Draw a contour plot in perspective. Need to reset color anyway.
         call axregion(-scbx3,scbx3,-scby3,scby3)
         call scalewn(xmin,xmax,ymin,ymax,.false.,.false.)
         call color(7)
c Calculate place of plane. 
         zplane=scbz3*(-1+(xn(ixnp(idfix)+n1)-zmin)*2./(zmax-zmin))
         if(iweb.eq.0)then
c Draw axes.
            call hdprset(0,0.)
c Ought to rescale the z-axis, but that was done in hidweb.
            call scale3(xmin,xmax,ymin,ymax,zmin,zmax)
c If we do, then we must reset jsw:
            jsw=1 + 256*6 + 256*256*7
            call axproj(igetcorner())
         else
c Set contour levels using the scaling of the box.
            do ic=1,icl
               cl(ic)=wz3min+(wz3max-wz3min)*(ic-1.)/(icl-1.)
            enddo
            if(logcontour)then  
               icl=max(icl,10)
               ilmax=nint(alog10(abs(wz3max))+0.5)
               do ic=2,icl-1
                  cl(icl-ic+1)=10.**(ilmax)*10**(-float((ic-1)/5))
     $                 *(1.-0.2*mod(ic+4,5))
               enddo
! Logarithmic contouring makes no sense unless >0
               if(wz3min.gt.0)cl(1)=wz3min
               cl(icl)=wz3max
            endif
!               write(*,*)icl,wz3max,wz3min,ilmax
!               write(*,*)(cl(ic),ic=1,icl)
         endif
c Get back current eye position xe1 etc.
         call trn32(xe,ye,ze,xe1,ye1,ze1,-1)
         if(mod(icontour,4).eq.1)call hdprset(-3,sign(scbz3,ze1))
         if(mod(icontour,4).eq.2)call hdprset(-3,zplane)
         if(mod(icontour,4).eq.3)call hdprset(-3,-sign(scbz3,ze1))
c Contour with coloring, using vector axes, maybe without labelling.
         iclhere=iclsign*icl
c         iclhere=icl
         icsw=17
         if(iweb.gt.1.or.icontour.ge.4)icsw=icsw-16
         call contourl(zp,pp,nf1+1-if1,
!         call contourl(zp(if1,if2),pp,nw,
     $        nf1+1-if1,nf2+1-if2,cl,iclhere,
     $        xn(ixnp(idp1)+if1),xn(ixnp(idp2)+if2),icsw)
         Erange=0.
         iasw=9
         if(igradleg.eq.1)then
            call gradlegend(cl(1),cl(icl),0.,1.1,1.,1.1,.02,.false.)
         elseif(igradleg.eq.2)then
            call gradlegend(cl(1),cl(icl),1.1,1.,1.1,0.,.02,.true.)
         endif
         call color(igray())
         if(larrow)call arrowplot(vp(if1,if2,1),vp(if1,if2,2),Erange,nw
     $        ,nf1+1-if1,nf2+1-if2,xn(ixnp(idp1)+if1),xn(ixnp(idp2)+if2)
     $        ,iasw,idum,idum) 
         call ticlabtog()
         call axis()
         call ticlabtog()
         call axis2()
      endif
      if(iweb.eq.1.and.icontour.eq.3)then
         call hidweb(xn(ixnp(idp1)+if1),xn(ixnp(idp2)+if2),
!     $        zp(if1,if2),nw,nf1+1-if1,nf2+1-if2,jsw)
     $        zp,nf1+1-if1,nf1+1-if1,nf2+1-if2,jsw)  ! zp 1-D
      endif

      if(ips.ne.0)then
c We called for a local print of plot. Terminate and switch it off.
c The problem is that pfset does not immediately turn off writing to
c the plotfile. It only sets up the switch off at the next pltinit.
c Then cubeupd calls invoked from button presses in the user interface
c are written to a closed file because pltend has closed it. Therefore
c we need an immediate turn off of the pfsw. This is fixed by calling
c pfset before pltend. On entry, pltend sees pfsw as negative, so does
c not pause. The last thing pltend does is set the pfsw from
c the pfnextsw set by pfset to zero.
         call pfset(0)
         call prtend(' ')
         write(*,*)'Terminating sliceweb ips=',ips,pfsw
         ips=0
      endif

c User interface
      call ui3d(n1,iuds,idfix,iquit,laspect,jsw,iclipping,ips,if1,if2
     $     ,nf1,nf2,idp1,idp2,icontour,iweb,ltellslice,igradleg
     $     ,colorscale,icl,logcontour)
      if(iquit.eq.0)goto 21
      call prtend(' ')
      idflast=idfix
      call hdprset(0,0.)
      end
c******************************************************************
      subroutine ui3d(n1,iuds,idfix,iquit,laspect,jsw,iclipping,ips,if1
     $     ,if2,nf1,nf2,idp1,idp2,icontour,iweb,ltellslice,igradleg
     $     ,colorscale,icl,logcontour)
c Encapsulated routine for controlling a 3-D plot.
c But many things have to be passed at present. A proper API needs
c to be designed but here's the approximate description.
c [Plane-position] n1 is controlled by up/down arrows, within the range
c 1-iuds(idfix), where iuds(3) are the used dimensional lengths, and
c idfix is the one currently being fixed (sliced).
c if1 nf1, if2 nf2 control the clipping positions  
c icontour and iweb determine the plotting of the contours and web. 
c idp1, idp2 are the the two other dimensions than idfix. 
c jsw is to do with contouring. laspect preserves aspect-ratio.
c iquit is returned as non-zero to command an end to the display.
      integer iuds(3)
      include 'world3.h'

      logical laspect,ltellslice,logcontour
c 3d display user interface.
c-----------------------------------
      iquit=0
c Limit framing rate to 30fps.
      call usleep(15000)
c User interface interpret key-press.
 24   call eye3d(isw)
c      write(*,*)'isw',isw
      if(isw.eq.ichar('f')) goto 24
      if(isw.eq.0) iquit=1
      if(isw.eq.65364 .and. n1.gt.1) n1=n1-1
      if(isw.eq.65362 .and. n1.lt.iuds(idfix)) n1=n1+1
      if(isw.eq.ichar('q')) iquit=1
      if(isw.eq.ichar('a')) laspect=.not.laspect
      if(isw.eq.ichar('x')) laspect=.false.
      if(isw.eq.ichar('z')) laspect=.false.
      if(isw.eq.ichar('v')) laspect=.false.
      if(isw.eq.ichar('d')) call noeye3d(0)
      if(isw.eq.ichar('s')) jsw=1 + 256*6 + 256*256*7
      if(isw.eq.ichar('t')) call togi3trunc()
      if(isw.eq.ichar('g')) igradleg=mod(igradleg+1,3)
      if(isw.eq.ichar('1')) read(*,*)ax3chars(1)
      if(isw.eq.ichar('2')) read(*,*)ax3chars(2)
      if(isw.eq.ichar('3')) read(*,*)ax3chars(3)
      if(isw.eq.ichar('4')) logcontour=.not.logcontour
      if(isw.eq.ichar('n')) read(*,*)icl
      if(isw.eq.ichar('p'))then
         call pfset(3)
         ips=3
      endif
      if(isw.eq.ichar('-'))then
         colorscale=colorscale*1.1
         call brgwscaled(0.,colorscale)
      elseif(isw.eq.ichar('='))then
         colorscale=colorscale/1.1
         call brgwscaled(0.,colorscale)
      endif
c Change fixed dimension, remove clipping, force scaling.
      if(isw.eq.65361)then
         idfix=mod(idfix+1,3)+1
         iclipping=0
         n1=(iuds(idfix)+1)/2
         jsw=1 + 256*6 + 256*256*7
      elseif(isw.eq.65363)then
         idfix=mod(idfix,3)+1
         iclipping=0
         n1=(iuds(idfix)+1)/2
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
         nf2=max(nf2-1,1)
         if2=min(nf2,if2)
         if(nf2.eq.if2)call linevecdoc(idfix,n1,2,nf2)
      elseif(isw.eq.44)then
c ,
         iclipping=1
         if2=min(if2+1,iuds(idp2))
         nf2=max(if2,nf2)
         if(nf2.eq.if2)call linevecdoc(idfix,n1,2,nf2)
      elseif(isw.eq.ichar('m'))then
c m
         iclipping=1
         if2=max(if2-1,1)
      elseif(isw.eq.ichar('l'))then
c l
         iclipping=1
         nf1=max(nf1-1,1)
         if1=min(nf1,if1)
         if(nf1.eq.if1)call linevecdoc(idfix,n1,1,nf1)
      elseif(isw.eq.ichar('u'))then
         ltellslice=.not.ltellslice
      elseif(isw.eq.59)then
c ;
         iclipping=1
         nf1=min(nf1+1,iuds(idp1))
c j
      elseif(isw.eq.ichar('j'))then
         iclipping=1
         if1=max(if1-1,1)
c k
      elseif(isw.eq.ichar('k'))then
         iclipping=1
         if1=min(if1+1,iuds(idp1))
         nf1=max(if1,nf1)
         if(nf1.eq.if1)call linevecdoc(idfix,n1,1,nf1)
c      elseif(isw.eq.ichar('['))then
c         iclipping=1
      else
c         write(*,*)char(isw)
      endif
      if(isw.eq.ichar('c'))icontour=mod(icontour+1,8)
      if(isw.eq.ichar('w'))iweb=mod(iweb+1,4)
      if(isw.eq.ichar('h'))then
         write(*,*)' ======== Slice plotting interface:',
     $        '  arrows up/down: change slice.'
         write(*,*)
     $        ' arrows l/r: change dimension.',
     $        ' s: rescale. p: print. Drag mouse to rotate.'
         write(*,*)' (jkl;) (m,./): control plotting extent in 2 axes.'
         write(*,*)
     $        ' c: contour plane position. w: web, solid, smooth, none'
     $        ,' a: aspect'
     $        ,' t: truncation'
         write(*,*)' u: slice-telling; g: contour gradlegend;'
     $        ,' -= scale coloring down/up'
         write(*,*)
     $        ' d: disable interface; run continuously.',
     $        ' depress f: to interrupt running.'
         write(*,*)' 1, 2, 3: Enter new label for axis 1,2,3'
     $        ,' (type in terminal window)'
         write(*,*)' n: Enter new number of contours.',
     $       ' 4: toggle logarithmic contours'
      endif
      call rotatezoom(isw)
c End of user interface.
      end

c*******************************************************************
c Three-way slice contours.
      subroutine sliceGcont(ifull,iuds,u,nw,zp,ixnp,xn,ifixpt,utitle
     $     ,svec,vp)
c Plot projected contour representations on cuts in three dimensions 
c of quantity u on fixed values of dimensions.
c The initial intersection of the fixed planes is at ifixpt(ndims).
c The full dimensions of arrays u are ifull, used iuds.
      parameter(ndims=3,nd2=2*ndims)
      integer ifull(ndims)
      real u(ifull(1),ifull(2),ifull(3))
      real svec(ifull(1),ifull(2),ifull(3),ndims)
c The used dimensions of each are
      integer iuds(ndims)
c The dimensions of square working array, zp, (nw) must be larger
c than the greatest of iuds. If zp(1,1,1)=0 on entry, give control help.
      integer nw
      real zp(nw,nw,ndims)
      real vp(nw,nw,ndims,ndims)
c The node positions are given by vectors laminated into xn, whose
c starts for dimensions id are at ixnp(id)+1. 
c So the vector of positions for dimension id is 
c   ixnp(id)+1 to ixnp(id)+iuds(id) [and ixnp(id+1)-ixnp(id)>=iuds(id)]
      integer ixnp(ndims+1)
      real xn(*)
c The plane intersection, chosen to start, and returned after, is:
      integer ifixpt(ndims)
c If ifixpt(1) is negative, then do arrowplots.
c The plotted quantity title is
      character*(*) utitle
c Needed for perspective plot
      include 'world3.h'
c Workspace size is problematic.
      parameter (nwksp=1000000)
      character*1 pp(nwksp)
c Contour levels
      real cl(30)
c Make pfsw visible:
      include 'plotcom.h'
      common/cont1stlast/c1st,clast

c Local variables:
      integer np
      parameter (np=100)
      integer narrow
      parameter (narrow=15)
      real tp(np),up(np)
      integer idire(3),idir(3)
      real cbsize(3),cmin(3),cmax(3)
      real xpl(3),xnl(3)
      integer idtx(5)
      real tx(5)
      integer icontour,iweb
      integer isw,jsw,imode,itype,ifileno
      character*(30) form1
      logical lsideplot,larrow
      data lsideplot/.false./larrow/.false./
c Tell that we are looking from the top by default.
      data ze1/1./icontour/1/iweb/1/
      data xpl/0.,0.,0./xnl/1.,0.,0./ifileno/0/
c      data jsw/1+256*6+256*256*7/
      data jsw/460289/

c      write(*,*)'ifixpt',ifixpt
      ierr=1
      imv=1
      itri=0
      icl=0
      imode=0
      itype=0
      if(ifixpt(1).lt.0)then
         larrow=.not.larrow
         ifixpt(1)=-ifixpt(1)
      endif
      if(ifixpt(1).le.0.or.ifixpt(1).gt.iuds(1))ifixpt(1)=(iuds(1)+1)/2
      if(ifixpt(2).le.0.or.ifixpt(2).gt.iuds(2))ifixpt(2)=(iuds(2)+1)/2
      if(ifixpt(3).le.0.or.ifixpt(3).gt.iuds(3))ifixpt(3)=(iuds(3)+1)/2
      ips=0
      irotating=0
c      write(*,*)'ifixpt',ifixpt
      call minmax2(u(1,1,ifixpt(3)),ifull(1),iuds(1),iuds(2),umin,umax)
      call accisgradinit(64000,0,0,-64000,128000,64000)

      call setcube(.2,.2,.2,.5,.4)
c Start of plotting loop.
 21   if(.not.lsideplot)call accisinit
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

c Store the slices in the plot arrays.      
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
      if(larrow)then
         vmax=0.
         do k=1,ndims
            do j=1,iuds(3)
               do i=1,iuds(2)
                  vh=svec(ifixpt(1),i,j,k)
                  if(abs(vh).gt.vmax)vmax=vh
                  vp(i,j,1,k)=vh
               enddo
            enddo
            do i=1,iuds(3)
               do j=1,iuds(1)
                  vh=svec(j,ifixpt(2),i,k)
                  if(abs(vh).gt.vmax)vmax=vh
                  vp(i,j,2,k)=vh
               enddo
            enddo
            do i=1,iuds(1)
               do j=1,iuds(2)
                  vh=svec(i,j,ifixpt(3),k)
                  if(abs(vh).gt.vmax)vmax=vh
                  vp(i,j,3,k)=vh
               enddo
            enddo
         enddo
      endif


c---------------------------
      if(ierr.eq.0)then
c      write(*,*)'Starting line processing'
c Line processing. Create a line of up to 5 points in order with 
c a rank index indicating how many planes it is in front of 0-3.
         do ic=1,5
            tx(ic)=0.
            idtx(ic)=0
         enddo
         tx(1)=tp(1)
         tx(2)=tp(np)
         ntx=2
c The number of plane crossings that are behind this point, which is the
c number of planes in front of which this point is. idtx
         do id=1,3
c Iterate over directions/planes. Prevent a divide by zero
            if(xnl(id).eq.0)xnl(id)=1.e-20
            tc=(xn(ixnp(id)+ifixpt(id))-xpl(id))/xnl(id)
c            write(*,*)id,xnl(id),tc,ifixpt(id),xn(ixnp(id)+ifixpt(id))
c Insert it into the sorted tx array of crossings. 
c            write(*,*)'ntx',ntx,tc,tx(1),tx(ntx)
            if(tc.gt.tx(1).and.tc.lt.tx(ntx))then
               do ix=ntx,1,-1
c Is the higher t side of this plane closer to eye?
c If so, then points at higher t are in front of this plane.
                  isn=0
                  if(xnl(id)*idire(id).gt.0.)isn=1
                  if(tc.lt.tx(ix))then
c Shift up crossing that is higher, and adjust its plane-crossing rank.
                     tx(ix+1)=tx(ix)
                     idtx(ix+1)=idtx(ix)+isn
                  else
c Put this intersection into its slot
                     if(tc.gt.tx(1))then
                        tx(ix+1)=tc
c                        idtx(ix+1)=idtx(ix+1)+isn
                        idtx(ix+1)=idtx(ix)+isn
                        ntx=ntx+1
c Mark as used up:
                        tc=tx(1)
                     endif
c Adjust the plane-crossing rank of the lower slots.
                     idtx(ix)=idtx(ix)+(1-isn)
                  endif
c                  write(*,'(a,7i4)')'ix,isn,idtx',ix,isn,idtx
               enddo
            else
c This plane perhaps behind everything. Increment ranks.
               if(xnl(id)*idire(id).gt.0.and.tc.le.tx(1) .or.
     $              xnl(id)*idire(id).lt.0.and.tc.gt.tx(ntx))then
                  do ix=1,ntx
                     idtx(ix)=idtx(ix)+1
                  enddo
               endif
            endif
         enddo
c Now tx(1:ntx) contains the ends and intersections of the line in order.
c And idtx(1:ntx) contains the number of planes behind each segment.
c tx(1) and tx(ntx) are the ends of the segment at tp=0,1.
c Draw the entire line.
         ir=-1
c         write(*,*)'idire',idire
         write(*,*)'Rankdraw with points',ntx,ir
         write(*,61)(k,idtx(k),tx(k),k=1,ntx)
 61      format(2i3,f10.4,';',2i3,f10.4,';',
     $        2i3,f10.4,';',2i3,f10.4,';',2i3,f10.4)
         call rankdraw(xpl,xnl,ntx,tx,idtx,ir)
      endif
c End of line sorting and ranking.
c---------------------------

c There are for each dimension three points: 1, ifixpt(id), iuds(id).
c The position of the plane is specified in 3-d normal units. 
c Each contour is called with an origin and iud in the other dims.
c Because of memory, we need to make the the origin the smaller index.
      ifix1=1
      ifix2=3
      if(itype.eq.1)then
         ifix1=imv
         ifix2=imv
      elseif(itype.eq.2)then
         ifix1=3
         ifix2=3
      endif
c Here we need to cycle through 4 levels of plane from back to front.
c Any plane (fixed dimension 3 of them) has four segments. 
c The furthest segment must be drawn first, then the two next most far
c in any order, then the closest. We call the distance the "level"
c The three dimensions must be interleaved with each other, drawing 
c their planes at the same level together. 
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
         do ifix=ifix1,ifix2
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
            if(it1-ib1.gt.0 .and. it2-ib2.gt.0)then
c               iclhere=iclsign*icl   ! maybe wrong here iclsign not set.
               iclhere=icl
c               write(*,*)iclsign,icl,iclhere
               call contourl(zp(ib1,ib2,ifix),
     $           pp,nw,it1-ib1+1,it2-ib2+1,
     $           cl,iclhere,
     $           xn(ixnp(id1)+ib1),xn(ixnp(id2)+ib2),17+itri*64)

               if(larrow)then
c We need a sensible arrow scale. Can't make individual plots different. 
                  Erange=narrow*1.33*vmax
                  iasw=5
c                  write(*,*)ifix,ib1,ib2,it1,it2,id1,id2
                  call color(igray())
                  call arrowplot(vp(ib1,ib2,ifix,id1)
     $                 ,vp(ib1,ib2,ifix,id2),Erange,nw
     $                 ,it1+1-ib1,it2+1-ib2
     $                 ,xn(ixnp(id1)+ib1),xn(ixnp(id2)+ib2)
     $                 ,iasw,iuds(id1)/narrow,iuds(id2)/narrow) 
               endif
            endif
            if(ilev.eq.2)idir(id1)=-idir(id1)
            if(ilev.eq.3)idir(id2)=-idir(id2)
         enddo
         call hdprset(0,0.)
         if(ilev.eq.3.and.ierr.eq.0)call rankdraw(xpl,xnl,ntx,tx,idtx,1)
      enddo
      call axproj(igetcorner())
      call ax3labels('x','y','z')
      if(ierr.eq.0)call rankdraw(xpl,xnl,ntx,tx,idtx,2)
c To do the lineout properly obscured
c 1. Draw whole before all planes. This shows it in the correct areas of 
c volumes behind any 2 (and 3 but it never shows) planes.
c 2. Draw planes at levels 1-3.
c 3. Redraw line only in places that are in front of at least 2 planes.
c 4. Draw planes at level 4 (front box).
c 5. Draw line only in places in front of all 3 planes. 

      if(ips.ne.0)then
c We called for a local print of plot. Switch it off, and terminate
c Prevent pltend from querying the interface.
         write(*,*)'Terminating sliceGweb. ips=',ips
c         pfsw=-ips
         call pfset(0)
         call prtend(' ')
         ips=0
      endif
      call gradlegend(c1st,clast,-.55,0.,-.55,1.,.03,.false.) 
c-----------------------------------
c Limit framing rate to 30fps.
      call usleep(10000)
 31   continue
      iy=0
c User interface interpret key-press.
      call eye3d(isw)
c      write(*,*)'isw=',isw
      if(isw.eq.0) goto 23
      if(isw.eq.ichar('q')) goto 23
      if(isw.eq.ichar('0'))then
         imode=0                !48
         ierr=1
      endif
      if(isw.eq.49)imode=1
      if(isw.eq.ichar('d')) imv=mod(imv,2)+1
      if(isw.eq.ichar('j'))itype=mod(itype+1,3)
      if(imode.eq.0)then
c Start of plane crossing control
         if(isw.eq.65364)then
            ifixpt(3)=ifixpt(3)-1
            if(ifixpt(3).lt.1)ifixpt(3)=1
         elseif(isw.eq.65362)then
            ifixpt(3)=ifixpt(3)+1
            if(ifixpt(3).gt.iuds(3))ifixpt(3)=iuds(3)
         elseif(isw.eq.65361)then
            ifixpt(imv)=ifixpt(imv)+1*idire(imv)
            if(ifixpt(imv).lt.1)ifixpt(imv)=1
            if(ifixpt(imv).gt.iuds(imv))ifixpt(imv)=iuds(imv)
         elseif(isw.eq.65363)then
            ifixpt(imv)=ifixpt(imv)-1*idire(imv)
            if(ifixpt(imv).lt.1)ifixpt(imv)=1
            if(ifixpt(imv).gt.iuds(imv))ifixpt(imv)=iuds(imv)
         endif
c End of plane crossing control
      elseif(imode.eq.1)then
c         write(*,*)'Plane crossing control'
         if(isw.eq.65364)then
            xpl(3)=xpl(3)-.01*(xn(ixnp(4))-xn(ixnp(3)+1))
            iy=1
         elseif(isw.eq.65362)then
            xpl(3)=xpl(3)+.01*(xn(ixnp(4))-xn(ixnp(3)+1))
            iy=1
         elseif(isw.eq.65361)then
            xpl(imv)=xpl(imv)+.01*(xn(ixnp(imv+1))-xn(ixnp(imv)+1))
            iy=1
         elseif(isw.eq.65363)then
            xpl(imv)=xpl(imv)-.01*(xn(ixnp(imv+1))-xn(ixnp(imv)+1))
            iy=1
         endif
      endif
      if(isw.eq.ichar('a'))then
         ifileno=ifileno+1
         write(form1,'(''line'',i4.4,''.dat'')')ifileno
         open(unit=10,FILE=form1,status='unknown',err=901)
         close(unit=10,status='delete',err=901)
         open(unit=10,FILE=form1,status='new',err=901)
         write(*,*)'Writing line to file ',form1
         goto 902
 901     write(10,*)'Error opening file ',form1
 902     write(10,'(a,6f10.4)')'# Line from: ',xpl
         write(10,'(a,6f10.4)')'# In direction: ',xnl
         write(10,*)np
         write(10,'(2f12.5)')(tp(k),up(k),k=1,np)
         close(10)
      endif
      if(isw.eq.ichar('s')) call minmax2(u(1,1,ifixpt(3)),
     $     ifull(1),iuds(1),iuds(2),umin,umax)
      if(isw.eq.ichar('p'))then
         call pfset(3)
         ips=3
      endif
      if(isw.eq.ichar('t'))itri=mod(itri+1,2)
      if(isw.eq.ichar('l'))icl=-icl
      if(isw.eq.ichar('h'))then
         write(*,*)' ======== Slice contouring interface:',
     $        ' up/down, <-/-> arrows change slice.'
         write(*,*) ' d: changes <-/-> dimension controlled.',
     $        ' s:rescale. p:plot-to-file.'
         write(*,*)' mode control: 1: line-rotate. 0: object-rotate.'
         write(*,*)
     $        ' toggles: t smooth color; l labels; g sideplot;'
     $        ,' j slice type.'
         write(*,*)' Drag cursor to move view.'
      endif
      if(isw.eq.ichar('c'))icontour=mod(icontour+1,4)
      if(isw.eq.ichar('w'))iweb=mod(iweb+1,2)

      if(isw.eq.ichar('g'))then
         lsideplot=.not.lsideplot
         if(lsideplot)then
c Move the 3-D plot to bottom right            
            call setcube(.1,.1,.1,.8,.2)
         else
c Move it to the center.
            call setcube(.2,.2,.2,.5,.4)
         endif
      endif
      if(imode.eq.0)then
c Rotating and Zooming of slices control interface:
         call rotatezoom(isw)
      elseif(imode.eq.1)then
c Rotating and shifting of lineout.         
         call rotateline(isw,xnl,xpl,iy)
      endif
c End of user interface.
c------------------------------------
      if(imode.eq.1)call line3out(ifull,u,ixnp,xn,xpl,xnl,np,tp,up
     $     ,ierr)
c         write(*,*)xpl,xnl
c         write(*,'(2f10.4)')(tp(k),up(k),k=1,np)

      if(lsideplot)then
c The sideplot commands.
         call axregion(.12,.5,.12,.5)
         call autoplot(tp,up,np)
         write(form1,'(a,3f7.2)')'From',xpl
         call axlabels(form1,utitle)
         write(form1,'(a,3f6.2)')'Direction',xnl
         call jdrwstr(.31,.55,form1,.0)
      endif
c Show the lineout position if it is changed.
      if(iy.ne.0)then
         call vec3w(xpl(1),xpl(2),xpl(3),0)
         call vec3w(xnl(1)+xpl(1),xnl(2)+xpl(2),xnl(3)+xpl(3),1)
         goto 31
      endif
c      write(*,*)'Returning to 21'
      goto 21
 23   continue
      call prtend(' ')
      call hdprset(0,0.)

      end
c********************************************************************
      subroutine rotateline(isw,xnl,xpl,iy)
      real xnl(3),xpl(3)
      integer irotating
c g77 does not support these intrinsics in parameters:
c      parameter (angle=3.141593/180.,cs=cos(angle),sn=sin(angle))
      parameter (angle=3.141593/180.,
     $     cs=1-angle*angle*0.5*(1.-angle*angle/12.),
     $     sn=angle*(1.-angle*angle*(1-angle*angle/20.)/6.))
      real ssn
      data irotating/0/ssn/1./
      save irotating
c      iy=0
      if(isw.eq.ichar('r'))then
c Rotate
         ssn=1.
         irotating=1
      elseif(isw.eq.ichar('e'))then
         irotating=1
         ssn=-1.
      elseif(isw.eq.ichar('y'))then
         irotating=1
         ssn=-ssn
      elseif(isw.eq.ichar('h'))then
c Help text
         write(*,*)' Line direction r/e: rotate in x/y-plane.'
         write(*,*)' z/x: elevation in z. n: enter dir. m: enter point.'
      elseif(isw.eq.ichar('z'))then
         irotating=2
         ssn=1.
      elseif(isw.eq.ichar('x'))then
         irotating=2
         ssn=-1.
      endif
      if(irotating.eq.1)then
         yn1=cs*xnl(1)-ssn*sn*xnl(2)
         yn2=ssn*sn*xnl(1)+cs*xnl(2)
         xnl(1)=yn1
         xnl(2)=yn2
         irotating=0
         iy=1
c         write(*,*)'Direction rotation',yn1,yn2
      elseif(irotating.eq.2)then
         r=sqrt(xnl(1)**2+xnl(2)**2)
         if(r.ne.0)then
            rc=xnl(1)/r
            rs=xnl(2)/r
         else
            rc=1.
            rs=0.
         endif
         z=xnl(3)
c         write(*,*)'Elevation change',r,z,ssn
         xnl(1)=(cs*r-ssn*sn*z)*rc
         xnl(2)=(cs*r-ssn*sn*z)*rs
         xnl(3)=ssn*sn*r+cs*z
         irotating=0
         iy=2
      endif
      if(isw.eq.ichar('n'))then
         write(*,'(''Enter lineout direction as 3 reals here:'',$)')
         read(*,*)xnl
         iy=3
      endif
      if(isw.eq.ichar('m'))then
         write(*,'(''Enter lineout point as 3 reals here:'',$)')
         read(*,*)xpl
         iy=4
      endif
      end
c********************************************************************
      subroutine rotatezoom(isw)
      integer irotating
      data irotating/0/
      save irotating,sn,cs
      if(isw.eq.ichar('i'))then
c Get back current eye position xe1 etc.
c Then move it in.
         call trn32(xe,ye,ze,xe1,ye1,ze1,-1)
         xe1=xe1*.9
         ye1=ye1*.9
         ze1=ze1*.9
c         call trn32(0.,0.,0.,xe1,ye1,ze1,1)
         call trn32(xe,ye,ze,xe1,ye1,ze1,1)
      elseif(isw.eq.ichar('o'))then
c Move out.
         call trn32(xe,ye,ze,xe1,ye1,ze1,-1)
         xe1=xe1/.9
         ye1=ye1/.9
         ze1=ze1/.9
c         call trn32(0.,0.,0.,xe1,ye1,ze1,1)
         call trn32(xe,ye,ze,xe1,ye1,ze1,1)
      elseif(isw.eq.ichar('r'))then
c Rotate
         irotating=1
         cs=cos(.05)
         sn=sin(.05)
      elseif(isw.eq.ichar('e'))then
         irotating=1
         cs=cos(.05)
         sn=sin(-.05)
      elseif(isw.eq.ichar('y'))then
         irotating=1
         if(cs.ne.cos(.05))then
            cs=cos(.05)
            sn=sin(.05)
         else
            sn=-sn
         endif
      elseif(isw.eq.ichar('z'))then
c Magnify the unit cube
         call getcube(cbx,cby,cbz,xcbc,ycbc)
         call setcube(cbx*1.05,cby*1.05,cbz*1.05,xcbc,ycbc)
      elseif(isw.eq.ichar('x'))then
c Shrink the unit cube
         call getcube(cbx,cby,cbz,xcbc,ycbc)
         call setcube(cbx/1.05,cby/1.05,cbz/1.05,xcbc,ycbc)
      elseif(isw.eq.ichar('v'))then
! Set view to top for 2-d plot of slice
c Get back current eye position xe1 etc.
         call trn32(xe,ye,ze,xe1,ye1,ze1,-1)
         if(abs(xe1).lt.abs(ye1))then
            call trn32(xe,ye,ze,0.,sign(0.001,ye1),10.,1)
         else
            call trn32(xe,ye,ze,sign(0.001,xe1),0.,10.,1)
         endif
      elseif(isw.eq.ichar('h'))then
c Help text
         write(*,*)' i/o: move eye in/out. r/e: rotate.'
     $        ,' z/x: zoom/shrink. v: top-view'
      endif
      if(irotating.gt.0)then
c Get back current eye position xe1 etc.
         call trn32(xe,ye,ze,xe1,ye1,ze1,-1)
         xex=xe1-xe
         yex=ye1-ye
         xe1=xe+cs*xex-sn*yex
         ye1=ye+sn*xex+cs*yex
c Must tell to look at zero. Why? Don't think so.
c         call trn32(0.,0.,0.,xe1,ye1,ze1,1)
         call trn32(xe,xe,ze,xe1,ye1,ze1,1)
c         irotating=irotating-1
         irotating=0
c Is this necessary?
      endif
      end
c******************************************************************
      integer function irotatezoom()
c Functionalized rotatezoom that calls eye3d and also sets the 
c perspective file so that calls that set the perspective based
c on the eye.dat file obey the rotation. 
      call eye3d(isw)
      call rotatezoom(isw)
      call trn32(xe,ye,ze,xe1,ye1,ze1,-1)
      call puteye(xe1,ye1,ze1)
      irotatezoom=isw
      end
c******************************************************************
c Obtain a line profile through a quantity u on n-dimensional space.
c ndims is fixed as 3 in this version.
      subroutine line3out(ifull,u,ixnp,xn,xpl,xnl,np,tp,up,ierr)
c ifull(ndims) full dimensions of u
c Node positions are given by vectors laminated into xn, whose
c starts for dimensions id are at ixnp(id)+1. 
c Line passes through point xpl with cartesian direction xnl.
c So for each dimension x=xpl+t*xnl.
c np  length of line profile array to return.
c tp, up  values of line parameter and u returned along the line.
c ierr non-zero on error.
      integer ndims
      parameter (ndims=3)
      integer ifull(ndims)
      real u(ifull(1),ifull(2),ifull(3))
      integer ixnp(ndims+1)
      real xn(*)
      real xpl(ndims)
      real xnl(ndims)
      integer np
      real tp(np),up(np)

      integer ip(ndims)
      real pf(ndims)

      ierr=0
c Calculate the end points of the line
      tmin=-1.e30
      tmax=1.e30
      do id=1,ndims
         x0=xn(ixnp(id)+1)
         x1=xn(ixnp(id+1))
c         write(*,*)x0,x1
         if(xnl(id).ne.0.)then
            t0=(x0-xpl(id))/xnl(id)
            t1=(x1-xpl(id))/xnl(id)
            if((tmin-t0)*(t0-tmax).gt.0)then
c t0 is a new bound
               if(t0.le.t1)then
                  tmin=t0
               else
                  tmax=t0
               endif
            endif
            if((tmin-t1)*(t1-tmax).gt.0)then
c t1 is a new bound
               if(t0.le.t1)then
                  tmax=t1
               else
                  tmin=t1
               endif
            endif
         endif
      enddo
c Now tmin and tmax are the end parameters of the lineout, 
c and tmax>=tmin or else there's no intersection.
      if(tmax.lt.tmin)then
c There's no intersection with the box.
         write(*,*)'Null intersection of lineout.',tmin,tmax
         ierr=1
         return
      endif
c For each point in the lineout
      do i=1,np
         tp(i)=tmin+(tmax-tmin)*(.00001+.99998*(i-1)/(np-1))
c Find the interpolated u value. First tell what box we are in.
         do id=1,ndims
            xp=xpl(id)+tp(i)*xnl(id)
            ip(id)=int(accinterp(xn(ixnp(id)+1),ixnp(id+1)-ixnp(id),
     $           xp,pf(id)))
            if(ip(id).eq.0)then
               write(*,*)'Line profile error. Outside box.',id
     $           ,ixnp(id)+1,ixnp(id+1)-ixnp(id),xp,pf(id)
               write(*,*)(xn(k),k=ixnp(id)+1,ixnp(id+1)-ixnp(id))
               ierr=2
               return
            endif
            pf(id)=pf(id)-ip(id)
         enddo
c Now we have in ip(id) and pf(id), the array position information
c Do the average over the box's vertices 
         total=0.
         do il=0,2**ndims-1
c For each box vertex use a logical offset that
c amounts to the binary representation of its sequence number, i1.
c Weight according to whether this bit is one or zero.
            thisval=u(ip(1)+ibits(il,0,1),ip(2)+ibits(il,1,1),
     $           ip(3)+ibits(il,2,1))
            i1=il
            do ik=1,ndims
               i2=i1/2
               if(i1-2*i2.ne.0)then 
                  thisval=pf(ik)*thisval
               else
                  thisval=(1.-pf(ik))*thisval
               endif
               i1=i2
            enddo
            total=total+thisval
         enddo
c Now total is the interpolated value at this point.
         up(i)=total
      enddo
c Completed the lineout.
      end


c********************************************************************
      subroutine rankdraw(xpl,xnl,ntx,tx,idtx,ir)
c Draw visibly, those segments of the line which have rank greater than
c ir.
      integer ntx,ir
      real xpl(3),xnl(3)
      real tx(ntx)
      integer idtx(ntx)

c      write(*,*)idtx
      ipen=0
      do i=1,ntx
         call vec3w(xpl(1)+xnl(1)*tx(i),xpl(2)
     $        +xnl(2)*tx(i),xpl(3)+xnl(3)*tx(i),ipen)
c         write(*,*)'drawn',i,ipen
         if(idtx(i).gt.ir)then
            ipen=1
         else
            ipen=0
         endif
      enddo

      end
c********************************************************************
c Given a monotonic function Q(x)
c on a 1-D grid x=1..nq, solve Q(x)=y for x with interpolation.
c If return is 0, then the y is outside Q's range or other error.
c The function returns the integer part of x.
      function accinterp(Q,nq,y,x)
      integer nq
      real Q(nq)
      real y,x
c
      integer iqr,iql,iqx
      real Qx,Qr,Ql

c      write(*,*)'nq=',nq
      accinterp=0
      Ql=Q(1)
      Qr=Q(nq)
      iql=1
      iqr=nq
c Circumlocution to catch y=NAN
      if(.not.((y-Ql)*(y-Qr).le.0.)) then
c Value is outside the range.
         x=0
         return
      endif
 200  if(iqr-iql.eq.1)goto 210
      iqx=(iqr+iql)/2
      Qx=Q(iqx)
c      write(*,*)y,Ql,Qx,Qr,iql,iqr
      if((Qx-y)*(Qr-y).le.0.) then
         Ql=Qx
         iql=iqx
      else
         Qr=Qx
         iqr=iqx
      endif
      goto 200
 210  continue
c Now iql and iqr, Ql and Qr bracket Q
      if(Qr-Ql.ne.0.)then
         xpart=(y-Ql)/(Qr-Ql)
         x=xpart+iql
         accinterp=iql
      else
         x=iql
         write(*,*)'****** Error!: interp coincident points'
      endif
      end
c******************************************************************
      subroutine ixnpcreate(Nx,Ny,Nz,vx,vy,vz,ixnp,xn)
c Create the ixnp and xn combined vectors from three vectors for
c different directions vx,vy,vz,Nx,Ny,Nz
      integer ixnp(4)
      real xn(*),vx(Nx),vy(Ny),vz(Nz)

c      write(*,*)'Nx,Ny,Nz',Nx,Ny,Nz
c      write(*,'(10f8.3)')vx,vy,vz
      ixnp(1)=0
      do i=1,Nx
         xn(i+ixnp(1))=vx(i)
      enddo
      ixnp(2)=ixnp(1)+Nx
      do i=1,Ny
         xn(i+ixnp(2))=vy(i)
      enddo
      ixnp(3)=ixnp(2)+Ny
      do i=1,Nz
         xn(i+ixnp(3))=vz(i)
      enddo
      ixnp(4)=ixnp(3)+Nz
      end
c******************************************************************
      subroutine linevecdoc(idfix,n1,iskip,nf)
      integer ivec(3)
      ivec(idfix)=n1
      ivec(mod(idfix+iskip-1,3)+1)=nf
      ivec(mod(idfix+mod(iskip,2),3)+1)=0
      write(*,*)'linevec',ivec
      end
c******************************************************************
c Some gradients.
c red-green gradient:
c      call accisgradinit(64000,0,0,-64000,128000,64000)
c blue purple white gradient
c      call accisgradinit(-32000,-65000,0,97000,65500,150000)
c green yellow white 
c      call accisgradinit(-32000,0,-65000,97000,150000,65500)
c red orange white
c      call accisgradinit(0,-32000,-65000,150000,97000,65500)
c*****************************************************************
      subroutine ixnpconstruct(ndims,imax,iuds,xyzmesh,ixnp,xn)
! Construct the mesh position vector from given arguments.
c Input ndims, imax, iuds, xyzmesh. Output:
c The node positions are given by vectors laminated into xn, whose
c starts for dimensions id are at ixnp(id)+1. 
c So the vector of positions for dimension id is 
c   ixnp(id)+1 to ixnp(id)+iuds(id) [and ixnp(id+1)-ixnp(id)>=iuds(id)]
      integer ndims,iuds(ndims)
      real xyzmesh(imax,ndims)
      integer ixnp(ndims+1)
      real xn(*)
      ic=0
      do i=1,ndims
         ixnp(i)=ic
         do j=1,iuds(i)
            ic=ic+1
            xn(ic)=xyzmesh(j,i)
         enddo
      enddo
      ixnp(ndims+1)=ic+1
      end
c******************************************************************
      subroutine setax3chars(a,b,c)
      character*(*) a,b,c
      include 'world3.h'
      ax3chars(1)=a
      ax3chars(2)=b
      ax3chars(3)=c
      end
c*******************************************************************
      subroutine slicesimple(ifull,iuds,u,starts,ends,nw,zp)
c Simplified call to slicegweb that creates the positional index arrays
c assuming they are uniform between starts and ends.      
      parameter(ndims=3,nd2=2*ndims)
      integer ifull(ndims)
c The used dimensions of each are
      integer iuds(ndims)
      real u(ifull(1),ifull(2),ifull(3))
      real starts(ndims),ends(ndims)
      integer nw
      real zp(nw*nw)
      integer ixnpv(ndims+1)
      parameter (nx=2000)
      real xnv(nx)
      character*20 utitle

      call simpleixnp(ndims,iuds,starts,ends,ixnpv,xnv,nx)
      idfixp=1
      utitle=''
!      write(*,*)ifull,iuds,nw,ixnpv
      call sliceGweb(ifull,iuds,u,nw,zp,ixnpv,xnv,idfixp,utitle)
!     ,svec,vp)

      end
c*********************************************************************
      subroutine simpleixnp(ndims,iuds,starts,ends,ixnpv,xnv,nx)
c Create the uniform index arrays for slicesimple
      integer iuds(ndims),ixnpv(ndims+1)
      real starts(ndims),ends(ndims),xnv(nx)
      ixnpv(1)=0
      do id=1,ndims
         ixnpv(id+1)=ixnpv(id)+iuds(id)
         do i=ixnpv(id)+1,ixnpv(id+1)
            j=i-ixnpv(id)-1
            xnv(i)=starts(id)+j*(ends(id)-starts(id))/(iuds(id)-1)
         enddo
!         write(*,*)ixnpv(id)
!         write(*,'(10f8.4)')(xnv(i),i=ixnpv(id)+1,ixnpv(id+1))
      enddo
      end
