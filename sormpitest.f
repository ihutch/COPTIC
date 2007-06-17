      program sormpitest
c
      include 'objcom.f'
      integer Li,ni,nj
c      parameter (Li=100,ni=40,nj=40,nk=16)
      parameter (Li=100,ni=16,nj=16,nk=20)
      integer nblks
      parameter (nblks=1)
      integer nd
      parameter (nd=ndims_sor,nd2=nd*2)
      real u(Li,Li,Li),q(Li,Li,Li),cij(nd2+1,Li,Li,Li)
      real error(Li,Li,Li)
      real zp(Li,Li),cijp(nd2+1,Li,Li)
      include 'sormesh.f'
c
      external bdyset,faddu2,cijroutine
c      real x(Li),y(Li)
      real z(Li),xp(Li)
      character*100 form1,argument
c      character*10 cxlab,cylab
      common /ctl_sor/mi_sor,xjac_sor,eps_sor,del_sor,k_sor
      logical lplot
c testing arrays
      integer iuds(nd),ifull(nd),idims(nd),ium2(nd)
      real uplot(Li,Li),uprime(Li),xprime(Li),zero(Li),uanal(Li)

      common /myidcom/myid

      include '3dcom.f'
c Common data containing the object geometric information. 
c Each object, i < 64 has: type, data(odata).
c      integer ngeomobjmax,odata,ngeomobj
c      parameter (ngeomobjmax=63,odata=16)
c      real obj_geom(odata,ngeomobjmax)
c      common /objgeomcom/ngeomobj,obj_geom

c Geometry information read in.
      call readgeom('geomtest.dat')
c First object is sphere of radius rc and potential phi.
      rc=obj_geom(5,1)
      phi=-obj_geom(10,1)/obj_geom(8,1)

      lplot=.true.
c Deal with arguments
c      if(iargc().eq.0) goto "help"
      do i=1,iargc()
         call getarg(i,argument)
         if(argument(1:3).eq.'-p ') lplot=.false.
      enddo
      ndims=nd
      idims(1)=nblks
      idims(2)=nblks
      ifull(1)=Li
      ifull(2)=Li
      iuds(1)=ni
      iuds(2)=nj
c This circumlocution quiets warnings.
      n3=3
      if(nd.eq.n3)then
         idims(n3)=nblks
         ifull(n3)=Li
         iuds(n3)=nk
      endif
      n0=nj/2
      n1=nk/2
      if(ndims.lt.3)n1=1

      iof=0
      do id=1,ndims
c Pointer to start of vector.
         ixnp(id)=iof
c Mesh size with the faces removed:
         ium2(id)=iuds(id)-2
c Mesh data:
         do i=1,iuds(id)
            xn(i+iof)=(i-1.)/(iuds(id)-1.)
         enddo
         iof=iof+iuds(id)
      enddo
      ixnp(ndims+1)=iof
c      write(*,'(10f8.4)')xn

c For possibly different jacobi convergence radius parameters.
c 7 is a pretty good value in 3-d.
      do kk=7,7,2
c This loop good for up to 3 dimensions.
         do k=1,nk
         do j=1,nj
            do i=1,ni
c initialize u, q, and cij
               u(i,j,k)=0.
c Charge rod
c               if(abs(i-ni/2).le.4 .and. abs(j-nj/2).le.4 )then
c Charge ball
               if(abs(i-ni/2).le.ni/6 .and. abs(j-nj/2).le.nj/6
     $              .and. abs(k-nk/2).le.nk/6)then
c Turned off.
c                  q(i,j,k)=1.
                  q(i,j,k)=0.
               else
                  q(i,j,k)=0.
               endif
            enddo
         enddo
         enddo

      ipoint=0
c Remove edges by starting at (2,2,2,...) and using ium2.
      call mditerate(ndims,ifull,ium2,cijroutine,
     $     cij(1,2,2,2),ipoint)
c Initialize the region flags in the object data
      call iregioninit(ndims,ifull)

c The following requires include objcom.f, which is now omitted.
c      write(*,*)'Finished mesh setup. Used No of pointers:',oi_sor
      write(*,*)'Finished mesh setup.'
c      write(*,'(a,8f8.1)')'cij(*,3,3,3)=',(cij(i,3,3,3),i=1,nd2+1)

c This jacobi radius is pretty much optimized for Poisson equation with
c fixed boundary, 2-dimensional, when kk=5. 3-d not optimized.
      xyimb=(max(ni,nj)*2.)/float(ni+nj) - 1.
      xjac_sor=1.- ((kk/5.)*4./max(10,(ni+nj)/2)**2)
     $     *(1.-0.3*xyimb)
      mi_sor=4*(ni+nj+nk)+10
      eps_sor=1.e-5
c Stop iterating immediately
c         mi_sor=3
         
c Control. Bit 1, use my sor parameters, Bit 2 use faddu.
      ictl=3
c         ictl=1
c         write(*,*)'Calling sormpi, ni,nj=',ni,nj
c The main solver call. Returns process myid.
      call sormpi(ndims,ifull,iuds,cij,u,q,bdyset,faddu2,ictl,ierr
     $     ,myid,idims)


         if(myid.eq.0)then
            if(ierr.lt.0) write(*,*)'Not converged',ierr
            write(*,*) 'mi_sor,k_sor,xjac_sor,del_sor',
     $           mi_sor,k_sor,xjac_sor,del_sor,myid
            if(ni+nj.lt.40) then
               write(form1,'(''('',i4,''f8.5)'')')nj
               write(*,*)'u='
               write(*,form1)((u(i,j,n1),j=1,nj),i=1,ni)
            endif
         endif
      enddo
         
      if(myid.eq.0.and.lplot)then
c         call cijplot(ndims,ifull,iuds,cij)
      endif

      if(myid.eq.0)then
c intersection examination
c         call boxintersect(ndims,ifull,iuds,cij)

      endif


      if(myid.eq.0)then
c-------------------------------------------------------------------
c Write some data out for cross checking
c Delete the file first to help with nfs problems.
         open(21,file='smt.round',status='unknown')
         close(21,status='delete')
         open(21,file='smt.round',status='unknown')
         write(21,*)ni,nj,nk/2
         write(21,'(6e12.4)')((u(i,j,nk/2),i=1,ni/2),j=1,nj/2)
         close(21)
c This file has full accuracy.
         open(20,file='smt.out',status='unknown')
         close(20,status='delete')
         open(20,file='smt.out',status='unknown')
         write(20,*)ni,nj,nk/2
         write(20,*)((u(i,j,nk/2),i=1,ni/2),j=1,nj/2)
         close(20)
c-------------------------------------------------------------------
c Do some analytic checking of the case with a fixed potential sphere
c inside a logarithmic derivative boundary condition. 1/r solution.
c         rc=.2
         errmax=0.
         errvar=0.
         count=0
         do i=2,ni-1
            do j=2,nj-1
               do k=2,nk-1
                  count=count+1.
                  r=0.
                  xr=xn(i)
                  yr=xn(iuds(1)+j)
                  zr=xn(iuds(1)+iuds(2)+k)
                  r=sqrt((xr-.5)**2+(yr-.5)**2+(zr-.5)**2)
                  if(u(i,j,k).gt.1.e-4 .and.
     $                 r.ge.rc)then
                     e=u(i,j,k)-phi*rc/r
                     error(i,j,k)=e
                     errvar=errvar+e**2
                     if(abs(e).gt.abs(errmax))errmax=e
c                  if(abs(e).gt.1.e-3) 
c                  write(*,'(3i3,10f8.4)')i,j,k,
c     $                 xr,yr,zr,r,u(i,j,k),2.*rc/r,e
                  else
                     error(i,j,k)=0.
                  endif
               enddo
            enddo
         enddo
         errvar=errvar/count
       write(*,*)'Max error=',errmax,' Standard Deviation=',sqrt(errvar)
c      iform=7
c      uscale=10000000.
c      call udisplay(ndims,u,ifull,iuds,iform,uscale)
c-------------------------------------------------------------------
c Start of plotting section.
         if(lplot .and. abs(errmax).lt..1) then
            call yautoplot(u(1,n0,n1),ni)
            do i=1,ni
               xr=xn(i)
               yr=xn(iuds(1)+n0)
               zr=xn(iuds(1)+iuds(2)+n1)
               r=sqrt((xr-.5)**2+(yr-.5)**2+(zr-.5)**2)
               z(i)=phi*rc/r
               xp(i)=i
            enddo
            call polyline(xp,z,ni)
            call pltend()
         endif
         if(lplot)then
c Plotting slices.
            ifix=3
            call slice3web(ifull,iuds,u,cij,Li,zp,cijp,ixnp,xn,ifix,
     $           'potential:'//'!Ay!@')
c-------------------------------------------------------------------
c Start of gradient testing. Do a contour plot of u in a fixed plane
c Then for an array of points finer than the original array, do
c arrow plot showing gradient in this plane.
            idf=3
            id1=mod(idf,3)+1
            id2=mod(idf+1,3)+1
            ifixed=6
            if(.true.)then
            do i=1,iuds(id1)
               do j=1,iuds(id2)
                  ium2(idf)=ifixed
                  ium2(id1)=i
                  ium2(id2)=j
                  uplot(i,j)=u(ium2(1),ium2(2),ium2(3))
               enddo
               write(*,'(10f8.4)')(uplot(i,j),j=1,iuds(id2))
            enddo
            write(*,*)idf,id1,id2,ium2(1),ium2(2),ium2(3)
            call autocolcont(uplot,Li,iuds(id1),iuds(id2))
c            call contourl(uplot,cworka,Li,,iym,zclv,icl,x,y,icsw)
            call color(15)
            call drwstr(.3,.3,'hello')
            size=0.02
            angle=.5
            call charsize(size,0.5*size)
            call charangl(180.*angle/3.141593)
            call jdrwstr(wx2nx(10.),wy2ny(10.),
     $           char(ichar('_')+128)//char(0),0.)
            call charsize(0.,0.)
            call charangl(0.)
            call pltend()
            endif
c

c Dimension along which we are interpolating: idf
            ium2(id1)=iuds(id1)/2
            ium2(id2)=iuds(id2)/2
            write(*,*)'Starting uprime calculation'
            do i=1,Li
               zero(i)=0.
c Offset to start of idf position array.
               ioff=ixnp(idf)
c xn is the position array for each dimension arranged linearly.
c An array of size Li crossing the idf coordinate range
               xprime(i)=xn(ioff+2)+
     $              (xn(ioff+iuds(idf)-1)-xn(ioff+2))*(i-.99)/(Li-.98)

c Packaging the interpolation appears to be working. 
c Gives identical output.
               icinc=nd2+1
               iuinc=1
               do k=1,idf-1
                  icinc=icinc*ifull(k)
                  iuinc=iuinc*ifull(k)
               enddo
               ium2(idf)=1
               call gradinterpcorrect(
     $              cij(nd2+1,ium2(1),ium2(2),ium2(3))
     $              ,u(ium2(1),ium2(2),ium2(3))
     $              ,idf,icinc,iuinc,
     $              xprime(i),uprime(i),ix,xm)
c Analytic comparison.
               uanal(i)=-phi*rc/((xprime(i)-0.5)*abs(xprime(i)-0.5))
               write(*,*)i,idf,ioff,ix,xm,uprime(i)
            enddo
            write(*,*)'Ended uprime calculation'
            call autoplot(xprime,uprime,Li)
c            call polymark(xprime,uprime,Li,1)
            call polymark(xn(ioff+1),zero,iuds(idf),10)
            call winset(.true.)
            call dashset(1)
            call polyline(xprime,uanal,Li)
c Crude differencing to show the improvement in our interpolation.
            do kk=2,nk-2
               x0=xn(ioff+kk)
               x1=xn(ioff+kk+1)
               ium2(idf)=kk
               u0=u(ium2(1),ium2(2),ium2(3))
               ium2(idf)=kk+1
               u1=u(ium2(1),ium2(2),ium2(3))
               dudx=(u1-u0)/(x1-x0)
               xq=(x1+x0)*0.5
c               write(*,*)kk,x0,x1,xq
               call polymark(xq,dudx,1,4)
            enddo
            call pltend()
c-------------------------------------------------------------------
         endif
c End of plotting.
      endif
c-------------------------------------------------------------------
      call MPI_FINALIZE(ierr)
      end
c**********************************************************************
c**********************************************************************
c      subroutine bdyset0(ndims,ifull,iuds,cij,u,q)
c     Null version
c      end
c**********************************************************************
      subroutine bdyset(ndims,ifull,iuds,cij,u,q)
      integer ndims,ifull(ndims),iuds(ndims)
      real cij(*),u(*),q(*)
      external bdy3slope
c set the derivative to zero on boundaries 3.
      ipoint=0
      call mditerate(ndims,ifull,iuds,bdy3slope,u,ipoint)

c set the second derivative to zero on max j
c      do i=2,ni-1
c         u(i,nj)=relax*(2.*u(i,nj-1)-u(i,nj-2)) +(1.-relax)*u(i,nj)
c      enddo
      end
c**********************************************************************
      real function faddu0(u,fprime)
      faddu0=0.
      slope=1000.
      faddu=(u)*slope
      fprime=slope
      end
c**********************************************************************
      real function faddu(u,fprime)
      real u,fprime
      real*8 slope,expu,temp
      faddu1=0
      slope=1000.
      temp=u
      expu=exp(temp)
      faddu=(expu-1.D0)*slope
      fprime=faddu+slope
      end
c**********************************************************************
      real function faddu2(u,fprime)
      faddu2=0.
      faddu=0.
      fprime=0.
      end

c************************************************************************
      subroutine bdy3slope(inc,ipoint,indi,ndims,iused,u)
c Version of bdyroutine that sets derivative=0 on 3-boundary.
      integer ipoint,inc
      integer indi(ndims),iused(ndims)
      real u(*)

      parameter (mdims=10)
c Structure vector needed for finding adjacent u values.
      integer iLs(mdims+1)
      common /iLscom/iLs

c Algorithm: if on a boundary face of dimension >1, steps of 1 (dim 1).
c Otherwise steps of iused(1)-1 or 1 on first or last (of dim 1).
      inc=1
      do n=ndims,2,-1
         if(indi(n).eq.0)then
c On boundary face 0 of dimension n>1. Break.
c This is where we put boundary setting for n>1
            u(ipoint+1)=0.
            if(n.eq.3)then
c Second derivative is zero:
c               u(ipoint+1)=2.*u(ipoint+1+iLs(n))-u(ipoint+1+2*iLs(n)
c First derivative is zero:
               u(ipoint+1)=u(ipoint+1+iLs(n))
            endif
            goto 101
         elseif(indi(n).eq.iused(n)-1)then
c On boundary face iused(n) of dimension n>1. Break.
            u(ipoint+1)=0.
            if(n.eq.3) u(ipoint+1)=u(ipoint+1-iLs(n))
            goto 101
         endif
      enddo
c     We are not on any higher boundary.
c This is where the boundary setting is done for n=1
      u(ipoint+1)=0.
      if(indi(n).eq.0)then
         inc=iused(1)-1
      elseif(indi(n).eq.iused(n)-1)then
         inc=1
      else
         write(*,*)'BDY Error. We should not be here',
     $        n,ipoint,indi
         stop
      endif
 101  continue
c      write(*,*)'indi,inc,iused,ipoint',indi,inc,iused,ipoint

      end
c************************************************************************
      subroutine surfmark(Li,ni,nj,x,y,z,cij,ifix)
c Mark in 3-D the positions of those points at which cij(7,*,*) are non
c -zero, and there is a fraction<1 to the adjoining points.
      integer Li,ni,nj
      real x(ni),y(nj)
      real z(Li,nj)
      real cij(7,Li,nj)
      include 'objcom.f'

      call color(3)
      do i=1,ni
         do j=1,nj
            ipoint=cij(7,i,j)
            if(ipoint.ne.0)then
               call wxyz2nxyz(x(i),y(j),z(i,j),xn,yn,zn)
               do id=1,3
                  idff=mod(id-1+3-ifix,3)+1
                  do jd=1,2
                     iobj=ndata_sor*(2*(id-1)+(jd-1))+1
                     fraction=dob_sor(iobj,ipoint)
                     if(fraction.lt.1. .and. fraction.gt.0.)then
c                        write(*,*)i,j,ipoint,id,jd,fraction
                        ip=i
                        jp=j
                        zp=z(i,j)
                        if(idff.eq.1)ip=ip+(3-2*jd)
                        if(idff.eq.2)jp=jp+(3-2*jd)
                        if(idff.eq.3)zp=zp+.05*(3-2*jd)
                        call vec3n(xn,yn,zn,0)
                        call vec3w(x(i)*(1.-fraction)+x(ip)*fraction,
     $                       y(j)*(1.-fraction)+y(jp)*fraction,
     $                       z(i,j)*(1.-fraction)+zp*fraction,1)
                     endif
                  enddo
               enddo
            endif
         enddo
      enddo

      end
c*******************************************************************
      subroutine slice3web(ifull,iuds,u,cij,nw,zp,cijp,ixnp,xn,ifix,
     $           utitle)
      parameter(ndims=3,nd2=2*ndims)
c Plot web-projected and/or projected contour representations 
c of potential u on slices with
c fixed values of dimension ifix. Mark cij boundaries.
c The full dimensions of arrays u, cij are
      integer ifull(ndims)
      real u(ifull(1),ifull(2),ifull(3))
      real cij(2*ndims+1,ifull(1),ifull(2),ifull(3))
c The used dimensions of each are
      integer iuds(ndims)
c The dimensions of square working arrays, zp, cijp, (nw) must be larger
c than the greatest of iuds 
      integer nw
      real zp(nw,nw)
      real cijp(2*ndims+1,nw,nw)
c The point positions are given by vectors laminated into xn, whose
c starts for dimensions id are at ixnp(id)+1. 
      integer ixnp(ndims)
      real xn(*)
c The fixed dimension which is chosen to start, and returned after is:
      integer ifix
c The plotted quantity title is
      character*(*) utitle
c Needed for perspective plot
      include '/home/hutch/accis/world3.h'
c Workspace size is problematic.
      character*1 pp(40000)
c Contour levels
      real cl(30)
c Local variables:
      character*(10) cxlab,cylab
      character*(30) form1

      write(*,*)' Slice plotting.',
     $           ' up/down arrows change slice.'
      write(*,*) ' left/right arrows change dimension.',
     $           ' s: rescale. p: print. Drag mouse to rotate.'
      write(*,*) ' c: toggle contour plane. w: toggle web. r: rotate'
      iweb=1
      icontour=1
      if(ifix.lt.1 .or. ifix.gt.ndims)ifix=ndims
      ips=0
      irotating=0
c Initial slice number
      n1=iuds(ifix)/2
c     Plot the surface. With scaling 1. Web color 6, axis color 7.
      jsw=1 + 256*6 + 256*256*7
 21   call pltinit(0.,1.,0.,1.)
c Set the plotting arrays for fixed dimension ifix.
      idp1=mod(ifix,3)+1
      idp2=mod(ifix+1,3)+1
      nf1=iuds(idp1)
      call iwrite(idp1,iwidth,cxlab)
      nf2=iuds(idp2)
      call iwrite(idp2,iwidth,cylab)
c Only works for 3-D in present implementation.
c Could be fixed to be general, I suppose.
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

c Web drawing
      if(iweb.eq.1)call hidweb(xn(ixnp(idp1)+1),xn(ixnp(idp2)+1),
     $        zp,nw,nf1,nf2,jsw)
c Use this scaling until explicitly reset.
      jsw=0 + 256*6 + 256*256*7
      write(form1,'(''Dimension '',i1,'' Plane'',i4)')ifix,n1
      call drwstr(.1,.02,form1)
      call ax3labels('axis-'//cxlab,'axis-'//cylab,utitle)
c Here we want to mark bounding surfaces, which are associated with those
c cijs that have non-zero pointer.
      if(iweb.eq.1) call surfmark(nw,nf1,nf2,
     $     xn(ixnp(idp1)+1),xn(ixnp(idp2)+1),
     $     zp,cijp,ifix)

c Projected contouring.
      if(icontour.ne.0)then
c       Draw a contour plot in perspective. Need to reset color anyway.
         call color(4)
         call accisgradinit(64000,0,0,-64000,128000,64000)
         call axregion(-scbx3,scbx3,-scby3,scby3)
         xmin=xn(ixnp(idp1)+1)
         xmax=xn(ixnp(idp1)+iuds(idp1))
         ymin=xn(ixnp(idp2)+1)
         ymax=xn(ixnp(idp2)+iuds(idp2))
         zmin=xn(ixnp(ifix)+1)
         zmax=xn(ixnp(ifix)+iuds(ifix))
         call scalewn(xmin,xmax,ymin,ymax,.false.,.false.)
c Calculate place of plane. 
         zplane=scbz3*(-1+(xn(ixnp(ifix)+n1)-zmin)*2./(zmax-zmin))
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
         if(icontour.eq.1)call hdprset(-3,scbz3)
         if(icontour.eq.2)call hdprset(-3,zplane)
         call scalewn(1.,float(nf1),1.,float(nf2),.false.,.false.)
c Contour without labels, with coloring, direct on mesh.
         call contourl(zp,pp,nw,nf1,nf2,cl,icl,
     $        xn(ixnp(idp1)+1),xn(ixnp(idp2)+1),16)
         call cubed(icorner-8*(icorner/8))
      endif

      if(ips.ne.0)then
c We called for a local print of plot. Terminate and switch it off.
         call pltend()
         call pfset(0)
         ips=0
      endif

c User interface:
      if(irotating.gt.0)then
c Get back current eye position xe1 etc.
         call trn32(xe,ye,ze,xe1,ye1,ze1,-1)
         cs=cos(.1)
         sn=sin(.1)
         write(*,*)'irotating',irotating,xe1,ye1,ze1,cs,sn
         xex=xe1-xe
         yex=ye1-ye
         xe1=xe+cs*xex-sn*yex
         ye1=ye+sn*xex+cs*yex
c Rotation has a bug in it here. ze1 somehow gets changed.
         write(*,*)'setting',xe1,ye1,ze1
         call trn32(xe,ye,ze,xe1,ye1,ze1,1)
         irotating=irotating-1
         goto 21
      endif
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
      if(isw.eq.ichar('c'))icontour=mod(icontour+1,3)
      if(isw.eq.ichar('w'))iweb=mod(iweb+1,2)
      if(isw.eq.ichar('i'))then
c Get back current eye position xe1 etc.
         call trn32(xe,ye,ze,xe1,ye1,ze1,-1)
         xe1=xe1*.9
         ye1=ye1*.9
         ze1=ze1*.9
c Move it in.
         call trn32(xe,ye,ze,xe1,ye1,ze1,1)
      endif
      if(isw.eq.ichar('o'))then
c Get back current eye position xe1 etc.
         call trn32(xe,ye,ze,xe1,ye1,ze1,-1)
         xe1=xe1*1.1
         ye1=ye1*1.1
         ze1=ze1*1.1
c Move it out.
         call trn32(xe,ye,ze,xe1,ye1,ze1,1)
      endif
      if(isw.eq.ichar('r')) irotating=10
c End of user interface.
      goto 21
 23   continue
      end
c******************************************************************
c Testing and examination of the intersection data.
c Assume 3-d.
c Use the knowledge that this is really the intersection with a
c sphere of center xc, and radius rc. For each used triplet,
c Find the perpendicular distance of xc from the plane, and the 
c distance from the base of the perpendicular to the point. 
c Evaluate the difference of the perpendicular length from rc,
c as a test of correct plane selection, and check base distance.
      subroutine boxintersect(ndims,ifull,iuds,cij)
      integer ndims
      parameter (mdims=3)
      integer ifull(mdims),iuds(mdims)
      real cij(ndims*2+1,ifull(1),ifull(2),ifull(3))
      include 'objcom.f'
      include 'sormesh.f'
      real xx(3),xc(3),xb(3)
      integer ijk(3)
      real fracts(2*mdims),xfr(mdims),a(mdims)
      integer ipa(mdims),ipm(mdims)
c center of sphere
      data xc/.5,.5,.5/
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
               iobj=cij(2*ndims+1,ijk(1),ijk(2),ijk(3))
               if(iobj.ne.0)then
c                  write(*,*)'Found object',ijk,iobj
c The point:
                  do id=1,3
                     xx(id)=xn(ixnp(id)+ijk(id))
                  enddo
c The intersection fractions (forward and backward for each dim).
                  do if=1,2*ndims
                     no=ndata_sor*(if-1)+1
c                     write(*,*)'no=',no,' iobj=',iobj
                     fracts(if)=dob_sor(no,iobj)
                  enddo
c The flags
                  iflags=idob_sor(iflag_sor,iobj)
c                  write(*,*)'xx=',xx
c                  write(*,*)' fracts=',fracts
c Treat each box with the appropriate ipmarray. 2**ndims in all.
                  do jj=1,2**ndims
                     kk=jj-1
c set the plus-minus directions for this jj. 0->+, 1->-.
c and other quantities
                     iuse=0
                     if(btest(iflags,jj-1))iuse=1
                     sumxx=0.
                     suma=0
                     suma2=0.
                     sumxc=0.
                     do ii=1,ndims
                        ipm(ii)=(kk-2*(kk/2))
                        ipa(ii)=1-2*ipm(ii)
                        dx=xn(ixnp(ii)+ijk(ii)+ipa(ii))
     $                       -xn(ixnp(ii)+ijk(ii))
c                        if(fracts(2*ii+ipm(ii)-1).eq.1.)iuse=0
c Throw away recuts. Probably wrong but removes most errors.
c                        if(fracts(2*ii+ipm(ii)-1).eq.1.001)iuse=0
                        xfr(ii)=dx*fracts(2*ii+ipm(ii)-1)
                        a(ii)=1./xfr(ii)
                        suma2=suma2+a(ii)**2
                        sumxx=sumxx+a(ii)*xx(ii)
                        sumxc=sumxc+a(ii)*xc(ii)
                        kk=kk/2
                     enddo
                     if(iuse.eq.1)then
                        nuse=nuse+1
c Equation of plane is \sum a_i (x_i-xx_i) =1
                        amag=sqrt(suma2)
c Distance of xc from the plane is [\sum a_i(xc_i-xx_i) - 1]/|a|
                        pdist=(sumxc-sumxx-1)/amag
c Base of perp from xc to plane is xb_i = xc_i-pdist*a_i/|a|
c The distance of this from xx is the interesting thing
c It ought to be .le. sqrt(ndims)
                        xbxd=0.
                        xcxd=0.
                        do ii=1,3
                           xb(ii)=xc(ii)-(pdist/amag)*a(ii)
                           xbxd=xbxd+(xb(ii)-xx(ii))**2
                           xcxd=xcxd+(xc(ii)-xx(ii))**2
                        enddo
                        xbxd=sqrt(xbxd)
                        xcxd=sqrt(xcxd)
c                        if(xbxd.gt..05)then
                        err=abs(abs(pdist)-rc)
                        if(err.gt.ermax)ermax=err
                        if(abs(abs(pdist)-rc).gt.error)then
                           nerror=nerror+1
                           write(*,'(a,3f8.4,a,3f8.4,a,3f8.4)')
     $                          'xx=',xx,' xb=',xb
                           write(*,*)nerror,' fracts=',fracts
c Test of this case
c      do j1=1,2**ndims
c         k1=j1-1
c         do i1=1,ndims
c This could be done with bit manipulation probably more efficiently:
c direction(j)=bit(j)of(i). 0->+1, 1->-1.
c            ipa(i1)=1-2*(k1-2*(k1/2))
c            k1=k1/2
c         enddo
c Now ipa is set. Call boxedge, returning the inverse of fractions
c in fn, and the number of intersections found in npoints.
c         call boxedge(ndims,ipa,ijk,xb,npoints,1)
c         write(*,'(7i3,6f8.4)')ijk,ipa,npoints,xb,(1/xb(m),m=1,3)
c      enddo

                        endif
c                        write(*,*)iobj,ijk,ipm,pdist,xbxd,xcxd,jj
c     $                       ,(iflags/2**(m-1)-2*(iflags/2**m),m=1,8)
c     $                       ,(btest(iflags,m),m=0,7)
                     else
c                        write(*,*)iobj,ijk,ipa,'  Not used'
                     endif
                  enddo
               endif
            enddo
         enddo
      enddo
      write(*,*)'Boxintersect: Total errors',nerror,' of used',
     $     nuse,' maxerr=',ermax

      end
