!************************************************************************
      subroutine surfmark(Li,ni,nj,x,y,z,cij,ifix)
! Mark in 3-D the positions of those points at which cij(7,*,*) are non
! -zero, and there is a fraction<1 to the adjoining points.
      integer Li,ni,nj
      real x(ni),y(nj)
      real z(Li,nj)
      real cij(7,Li,nj)
      include 'ndimsdecl.f'
      include 'objcom.f'
      include '../accis/world3.h'

      call color(3)
      do i=1,ni
         do j=1,nj
            ipoint=int(cij(7,i,j))
            if(ipoint.ne.0)then
               call wxyz2nxyz(x(i),y(j),z(i,j),xn,yn,zn)
               do id=1,3
                  idff=mod(id-1+3-ifix,3)+1
                  do jd=1,2
                     iobj=ndata_cij*(2*(id-1)+(jd-1))+1
                     fraction=dob_cij(iobj,ipoint)
                     if(fraction.lt.1. .and. fraction.gt.0.)then
                        call vec3n(xn,yn,zn,0)
!                        write(*,*)i,j,ipoint,id,jd,fraction
                        if(idff.eq.1)then
                           xp=x(i)*(1.-fraction)+x(i+(3-2*jd))*fraction
                        else
                           xp=x(i)
                        endif
                        if(idff.eq.2)then
                           yp=y(j)*(1.-fraction)+y(j+(3-2*jd))*fraction
                        else
                           yp=y(j)
                        endif
                        zp=z(i,j)
                        if(idff.eq.3)then 
                           zp=zp+.05*(3-2*jd)*(wz3max-wz3min)
                        endif
                        call vec3w(xp,yp,zp,1)
                     endif
                  enddo
               enddo
            endif
         enddo
      enddo

      end
!*******************************************************************
      subroutine slice3web(ifull,iuds,u,cij,nw,zp,cijp,ixnp,xn,ifix,
     $           utitle,ictl)
      parameter(ndims=3,nd2=2*ndims)
! Plot web-projected and/or projected contour representations 
! of quantity u on slices with
! fixed values of dimension ifix. Mark cij boundaries.
! The full dimensions of arrays u, cij are
      integer ifull(ndims)
      real u(ifull(1),ifull(2),ifull(3))
      real cij(2*ndims+1,ifull(1),ifull(2),ifull(3))
! The used dimensions of each are
      integer iuds(ndims)
! The dimensions of square working arrays, zp, cijp, (nw) must be larger
! than the greatest of iuds 
      integer nw
      real zp(nw,nw)
      real cijp(2*ndims+1,nw,nw)
! The point positions are given by vectors laminated into xn, whose
! starts for dimensions id are at ixnp(id)+1. 
      integer ixnp(ndims)
      real xn(*)
! The fixed dimension which is chosen to start, and returned after is:
      integer ifix
! The plotted quantity title is
      character*(*) utitle
! Needed for perspective plot
      include '../accis/world3.h'
! Workspace size is problematic.
      character*1 pp(40000)
! Contour levels
      real cl(30)
! Local variables:
      integer icontour,iweb,iback
      integer isw,jsw
      character*(10) cxlab,cylab
      character*(30) form1
      logical lfirst
      data lfirst/.true./
! Tell that we are looking from the top by default.
      data ze1/1./icontour/1/iweb/1/iback/0/
      data cs/.707/sn/.707/

      if(lfirst)then
         write(*,*)' ======== Slice plotting.',
     $        ' up/down arrows change slice.'
         write(*,*) ' l/r arrows change dimension.',
     $        ' s: rescale. p: print. Drag mouse to rotate.'
         write(*,*)
     $        ' c: contour plane position. w: toggle web. r/e: rotate.'
         iweb=1
         icontour=1
         lfirst=.false.
      endif
      if(ifix.lt.1 .or. ifix.gt.ndims)ifix=ndims
      ips=0
      irotating=0
! Initial slice number
      n1=iuds(ifix)/2
!     Plot the surface. With scaling 1. Web color 6, axis color 7.
      jsw=1 + 256*6 + 256*256*7
      call accisgradinit(64000,0,0,-64000,128000,64000)
 21   call pltinit(0.,1.,0.,1.)
! Set the plotting arrays for fixed dimension ifix.
      idp1=mod(ifix,3)+1
      idp2=mod(ifix+1,3)+1
      nf1=iuds(idp1)
      call iwrite(idp1,iwidth,cxlab)
      nf2=iuds(idp2)
      call iwrite(idp2,iwidth,cylab)
! Only works for 3-D in present implementation.
! Could be fixed to be general, I suppose.
      do i=1,nf1
         do j=1,nf2
            if(ifix.eq.1)then
               zp(i,j)=u(n1,i,j)
               do k=1,nd2+1
                  cijp(k,i,j)=cij(k,n1,i,j)
               enddo
            elseif(ifix.eq.2)then
               zp(i,j)=u(j,n1,i)
               do k=1,nd2+1
                  cijp(k,i,j)=cij(k,j,n1,i)
               enddo
            elseif(ifix.eq.3)then
               zp(i,j)=u(i,j,n1)
               do k=1,nd2+1 
                  cijp(k,i,j)=cij(k,i,j,n1)
               enddo
            endif
         enddo
      enddo

! Web drawing. First call is needed to set scaling.
!      write(*,*)'jsw=',jsw
      if(iweb.eq.1)
     $     call hidweb(xn(ixnp(idp1)+1),xn(ixnp(idp2)+1),
     $        zp,nw,nf1,nf2,jsw)
! Use this scaling until explicitly reset.
      jsw=0 + 256*6 + 256*256*7
      write(form1,'(''Dimension '',i1,'' Plane'',i4)')ifix,n1
      call drwstr(.1,.02,form1)
      call ax3labels('axis-'//cxlab,'axis-'//cylab,utitle)
! Here we want to mark bounding surfaces, which are associated with those
! cijs that have non-zero pointer.
      if(ictl.eq.1 .and. iweb.eq.1) call surfmark(nw,nf1,nf2,
     $     xn(ixnp(idp1)+1),xn(ixnp(idp2)+1),
     $     zp,cijp,ifix)

! Projected contouring.
      if(icontour.ne.0)then
!       Draw a contour plot in perspective. Need to reset color anyway.
         call color(4)
         call axregion(-scbx3,scbx3,-scby3,scby3)
         xmin=xn(ixnp(idp1)+1)
         xmax=xn(ixnp(idp1)+iuds(idp1))
         ymin=xn(ixnp(idp2)+1)
         ymax=xn(ixnp(idp2)+iuds(idp2))
         zmin=xn(ixnp(ifix)+1)
         zmax=xn(ixnp(ifix)+iuds(ifix))
         call scalewn(xmin,xmax,ymin,ymax,.false.,.false.)
! Calculate place of plane. 
         zplane=scbz3*(-1+(xn(ixnp(ifix)+n1)-zmin)*2./(zmax-zmin))
! accis perspective corner for axes and cube.
         icorner=igetcorner()
         if(iweb.eq.0)then
! Draw axes.
            call hdprset(0,0.)
! Ought to rescale the z-axis, but that was done in hidweb.
            call scale3(xmin,xmax,ymin,ymax,zmin,zmax)
! If we do, then we must reset jsw:
            jsw=1 + 256*6 + 256*256*7
            call axproj(icorner)
         else
! Set contour levels using the scaling of the box.
            icl=6
            do ic=1,icl
               cl(ic)=wz3min+(wz3max-wz3min)*(ic-1.)/(icl-1.)
            enddo
         endif
! Get back current eye position xe1 etc.
         call trn32(xe,ye,ze,xe1,ye1,ze1,-1)
         if(icontour.eq.1)call hdprset(-3,sign(scbz3,ze1))
         if(icontour.eq.2)call hdprset(-3,zplane)
         if(icontour.eq.3)call hdprset(-3,-sign(scbz3,ze1))
! Contour without labels, with coloring, using vector axes
         call contourl(zp,pp,nw,nf1,nf2,cl,icl,
     $        xn(ixnp(idp1)+1),xn(ixnp(idp2)+1),17)
         call axis()
         call axis2()
         if(iweb.ne.1)call cubed(icorner-8*(icorner/8))
      endif

      if(iweb.eq.1.and.icontour.eq.3)then
         call hidweb(xn(ixnp(idp1)+1),xn(ixnp(idp2)+1),
     $        zp,nw,nf1,nf2,jsw)
! This was necessary when hidweb used to change jsw.
!         jsw=0 + 256*6 + 256*256*7
      endif

      if(ips.ne.0)then
! We called for a local print of plot. Terminate and switch it off.
         call pltend()
         call pfset(0)
         ips=0
      endif

!-----------------------------------
! Double buffering 
      if(iback.ne.0)then 
!         call glfront()
      endif
! Limit framing rate to 30fps.
      call usleep(30000)
! User interface interpret key-press.
      call eye3d(isw)
      if(isw.eq.0) goto 23
      if(isw.eq.65364 .and. n1.gt.1) n1=n1-1
      if(isw.eq.65362 .and. n1.lt.iuds(ifix)) n1=n1+1
      if(isw.eq.ichar('q')) goto 23
      if(isw.eq.ichar('s')) jsw=1 + 256*6 + 256*256*7
      if(isw.eq.ichar('p'))then
         call pfset(3)
         ips=3
      endif
      if(isw.eq.65361)then
         ifix=mod(ifix+1,3)+1
         n1=iuds(ifix)/2
      elseif(isw.eq.65363)then
         ifix=mod(ifix,3)+1
         n1=iuds(ifix)/2
      endif
      if(isw.eq.ichar('c'))icontour=mod(icontour+1,4)
      if(isw.eq.ichar('w'))iweb=mod(iweb+1,2)
      if(isw.eq.ichar('i'))then
! Get back current eye position xe1 etc.
         call trn32(xe,ye,ze,xe1,ye1,ze1,-1)
         xe1=xe1*.9
         ye1=ye1*.9
         ze1=ze1*.9
! Move it in.
         call trn32(0.,0.,0.,xe1,ye1,ze1,1)
      elseif(isw.eq.ichar('o'))then
! Get back current eye position xe1 etc.
         call trn32(xe,ye,ze,xe1,ye1,ze1,-1)
         xe1=xe1*1.1
         ye1=ye1*1.1
         ze1=ze1*1.1
! Move it out.
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
!         call glback()
         iback=1
! Get back current eye position xe1 etc.
         call trn32(xe,ye,ze,xe1,ye1,ze1,-1)
!         write(*,*)'irotating',irotating,xe1,ye1,ze1,cs,sn
         xex=xe1-xe
         yex=ye1-ye
         xe1=xe+cs*xex-sn*yex
         ye1=ye+sn*xex+cs*yex
! Must tell to look at zero.
         call trn32(0.,0.,0.,xe1,ye1,ze1,1)
         irotating=irotating-1
      endif
! End of user interface.
!------------------------------------
      goto 21
 23   continue
      call hdprset(0,0.)
      end
!******************************************************************
! Testing and examination of the intersection data.
! Assume 3-d.
! Use the knowledge that this is really the intersection with a
! sphere of center xc, and radius rc. For each used triplet,
! Find the perpendicular distance of xc from the plane, and the 
! distance from the base of the perpendicular to the point. 
! Evaluate the difference of the perpendicular length from rc,
! as a test of correct plane selection, and check base distance.
      subroutine boxintersect(ifull,iuds,cij)
      include 'ndimsdecl.f'
      parameter (mdims=3)
      integer ifull(mdims),iuds(mdims)
      real cij(ndims*2+1,ifull(1),ifull(2),ifull(3))
      include 'objcom.f'
      include 'meshcom.f'
      real xx(3),xc(3),xb(3)
      integer ijk(3)
      real fracts(2*mdims),xfr(mdims),a(mdims)
      integer ipa(mdims),ipm(mdims)
! center of sphere
      data xc/.0,.0,.0/
      data rc/.48/
      
      ermax=0.
      error=.01
      nuse=0
      nerror=0
      write(*,*)'Starting boxintersect',ndims,ifull,iuds
      do k=1,iuds(3)
         ijk(3)=k
         do j=1,iuds(2)
            ijk(2)=j
            do i=1,iuds(1)
               ijk(1)=i
               iobj=int(cij(2*ndims+1,ijk(1),ijk(2),ijk(3)))
               if(iobj.ne.0)then
!                  write(*,*)'Found object',ijk,iobj
! The point:
                  do id=1,3
                     xx(id)=xn(ixnp(id)+ijk(id))
                  enddo
! The intersection fractions (forward and backward for each dim).
                  do if=1,2*ndims
                     no=ndata_cij*(if-1)+1
!                     write(*,*)'no=',no,' iobj=',iobj
                     fracts(if)=dob_cij(no,iobj)
                  enddo
! The flags
                  iflags=idob_cij(iflag_cij,iobj)
!                  write(*,*)'xx=',xx
!                  write(*,*)' fracts=',fracts
! Treat each box with the appropriate ipmarray. 2**ndims in all.
                  do jj=1,2**ndims
                     kk=jj-1
! set the plus-minus directions for this jj. 0->+, 1->-.
! and other quantities
                     iuse=0
                     if(btest(iflags,jj-1))iuse=1
                     sumxx=0.
!                     suma=0
                     suma2=0.
                     sumxc=0.
                     do ii=1,ndims
                        ipm(ii)=(kk-2*(kk/2))
                        ipa(ii)=1-2*ipm(ii)
                        dx=xn(ixnp(ii)+ijk(ii)+ipa(ii))
     $                       -xn(ixnp(ii)+ijk(ii))
!                        if(fracts(2*ii+ipm(ii)-1).eq.1.)iuse=0
! Throw away recuts. Probably wrong but removes most errors.
!                        if(fracts(2*ii+ipm(ii)-1).eq.1.001)iuse=0
                        xfr(ii)=dx*fracts(2*ii+ipm(ii)-1)
                        a(ii)=1./xfr(ii)
                        suma2=suma2+a(ii)**2
                        sumxx=sumxx+a(ii)*xx(ii)
                        sumxc=sumxc+a(ii)*xc(ii)
                        kk=kk/2
                     enddo
                     if(iuse.eq.1)then
                        nuse=nuse+1
! Equation of plane is \sum a_i (x_i-xx_i) =1
                        amag=sqrt(suma2)
! Distance of xc from the plane is [\sum a_i(xc_i-xx_i) - 1]/|a|
                        pdist=(sumxc-sumxx-1)/amag
! Base of perp from xc to plane is xb_i = xc_i-pdist*a_i/|a|
! The distance of this from xx is the interesting thing
! It ought to be .le. sqrt(ndims)
                        xbxd=0.
                        xcxd=0.
                        do ii=1,3
                           xb(ii)=xc(ii)-(pdist/amag)*a(ii)
                           xbxd=xbxd+(xb(ii)-xx(ii))**2
                           xcxd=xcxd+(xc(ii)-xx(ii))**2
                        enddo
                        xbxd=sqrt(xbxd)
                        xcxd=sqrt(xcxd)
!                        if(xbxd.gt..05)then
                        err=abs(abs(pdist)-rc)
                        if(err.gt.ermax)ermax=err
                        if(abs(abs(pdist)-rc).gt.error)then
                           nerror=nerror+1
                           write(*,'(a,3f8.4,a,3f8.4,a,3f8.4)')
     $                          'xx=',xx,' xb=',xb
                           write(*,*)nerror,' fracts=',fracts
! Test of this case
!      do j1=1,2**ndims
!         k1=j1-1
!         do i1=1,ndims
! This could be done with bit manipulation probably more efficiently:
! direction(j)=bit(j)of(i). 0->+1, 1->-1.
!            ipa(i1)=1-2*(k1-2*(k1/2))
!            k1=k1/2
!         enddo
! Now ipa is set. Call boxedge, returning the inverse of fractions
! in fn, and the number of intersections found in npoints.
!         call boxedge(ndims,ipa,ijk,xb,npoints,1)
!         write(*,'(7i3,6f8.4)')ijk,ipa,npoints,xb,(1/xb(m),m=1,3)
!      enddo

                        endif
!                        write(*,*)iobj,ijk,ipm,pdist,xbxd,xcxd,jj
!     $                       ,(iflags/2**(m-1)-2*(iflags/2**m),m=1,8)
!     $                       ,(btest(iflags,m),m=0,7)
                     else
!                        write(*,*)iobj,ijk,ipa,'  Not used'
                     endif
                  enddo
               endif
            enddo
         enddo
      enddo
      write(*,*)'Boxintersect: Total errors',nerror,' of used',
     $     nuse,' maxerr=',ermax

      end
