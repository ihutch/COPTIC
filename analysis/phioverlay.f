      program phioverlay
c Overlay contours obtained using xfig2trace on a contour plot of 
c potential. 

      include 'examdecl.f'
      parameter (ntheta=4,nr=38)
c      real thetadist(nr,ntheta),thetaval(ntheta)
      real rvalneg(0:nr),rval(0:nr)
c      integer ithetacount(nr,ntheta)

      parameter (nzmax=500)
      real rzdist(0:nr,nzmax),irzcount(nr,nzmax)
c      real rzmval(nzmax),rzmpos(nzmax),zpos(nzmax)
      character*1 ppath(0:nr,nzmax)
      character*20 string

      real cl(100)
      real xl(2),yl(2)
c xfig2trace parameters.
      parameter (np=200,nt=50)
      real xtraces(np,nt), ytraces(np,nt)
      integer nl(nt)
      real xyia(4)
      character*(100)overfile
c 
      call examargs(rp,iover,overfile)
      ied=1
      call array3read(phifilename,ifull,iuds,ied,u,ierr)
      if(ierr.eq.1)stop
      write(*,'(a,3i4,$)')'On grid',iuds
      write(*,'(a,2f7.2,a,2f7.2,a,2f7.2)')
     $     (',',xn(ixnp(kk)+1),xn(ixnp(kk+1)),kk=1,3)

      call accisgradinit(64000,0,0,-64000,128000,64000)
      x0=(xn(ixnp(1)+1)+xn(ixnp(2)))/2.
      y0=(xn(ixnp(2)+1)+xn(ixnp(3)))/2.
      rx=max(abs(xn(ixnp(1)+1)-x0),abs(xn(ixnp(2))-x0))
      ry=max(abs(xn(ixnp(2)+1)-x0),abs(xn(ixnp(3))-y0))
      rs=sqrt(rx**2+ry**2)
      call rzbin(nr,rzdist,irzcount,rval,rvalneg,phimax)
c Now rzdist(nr,nz) is the rz-distribution of the (average) potential.
c Contour it.
      call pfset(3)
      nc=30
      call fitrange(-phimax,phimax,nc,ipow,fac10,delta,first,xlast)
c         write(*,*)ipow,fac10,delta,first,xlast
      do k=1,nc
         cl(k)=(first+k*delta)
      enddo
      nc=abs(2.*phimax)/abs(delta)
c         write(*,*)' Phimax', phimax,' fac10',fac10
c         write(*,*)' Contours',nc,(cl(k),k=1,nc)
      iconsw=1+64
c Change z1 and z2, and nrp if you want to crop the region plotted.
      z1=1
      z2=iuds(3)
      nz=z2+1-z1
      nrp=nr-10
c      call pltinaspect(-rval(nrp),rval(nrp),
      call pltinit(0.,rval(nrp),
     $     xn(ixnp(3)+z1),xn(ixnp(3)+z2))
      call color(15)
      call contourl(rzdist(0,z1),ppath,nr+1,nrp+1,nz,
     $     cl,nc,rval,xn(ixnp(3)+z1),iconsw)
c      call contourl(rzdist(0,z1),ppath,nr+1,nrp+1,nz,
c     $     cl,nc,rvalneg,xn(ixnp(3)+z1),iconsw)
      call color(15)
c      call axis()
c      call axlabels('r','z')
c         write(*,*)'ipow',ipow
      call fwrite(delta,iwd,max(-ipow,0)+1,string)
      call boxtitle('!Af!@-contours ('//string(1:iwd)
     $     //'T!de!d spaced)')
      call gradlegend(-phimax,phimax,-.35,0.,-.35,1.,-.02,.false.) 
c Indicate rectangle limits.
      yl(1)=xn(ixnp(3)+z1)
      yl(2)=xn(ixnp(3)+z2)
      call dashset(2)
      xl(1)=rx
      xl(2)=rx
      call polyline(xl,yl,2)
      xl(1)=-rx
      xl(2)=-rx
      call polyline(xl,yl,2)
      if(ry.ne.rx)then
         xl(1)=ry
         xl(2)=ry
         call polyline(xl,yl,2)
         xl(1)=-ry
         xl(2)=-ry
c         call polyline(xl,yl,2)
      endif
      call dashset(0)
      xl(1)=0.
      xl(2)=0.
      call polyline(xl,yl,2)
      call fwrite(vd,iwd,2,string)
c         call jdrwstr(.02,.31,'v!dd!d='//string(1:iwd),1.)
      call fwrite(debyelen,iwd,1,string)
c         call jdrwstr(.02,.25,'!Al!@!dDe!d='//string(1:iwd),1.)
      call fwrite(phip,iwd,2,string)
c         call jdrwstr(.02,.28,'!Af!@!dp!d='//string(1:iwd),1.)
      call scalewn(0.,rval(nrp)/debyelen,
     $     xn(ixnp(3)+z1)/debyelen,xn(ixnp(3)+z2)/debyelen
     $     ,.false.,.false.)
      call axis()
      call axlabels('r/!Al!@!dDe!d','z/!Al!@!dDe!d')
c-------------------------------------------------------------------
      if(iover.gt.0)call xfig2trace(np,nt,xtraces,ytraces,il,nl,xyia
     $     ,overfile)
      write(*,*)'Return from xfig2',il,(nl(k),k=1,il),xyia
      if(il.gt.0)then
         do k=1,il
c            do j=1,nl(k)
c               xtraces(j,k)=xtraces(j,k)*debyelen
c               ytraces(j,k)=ytraces(j,k)*debyelen
c            enddo
            write(*,*)k,(xtraces(j,k),j=1,nl(k))
            write(*,*)k,(ytraces(j,k),j=1,nl(k))
            call winset(.true.)
            call dashset(2)
            call color(15)
            call polyline(xtraces(1,k),ytraces(1,k),nl(k))
         enddo
      endif
c------------------------------------------------------------------

      call pltend()

      call exit(0)
 101  write(*,*)'Error writing output'
      end

c*************************************************************
      subroutine examargs(rp,iover,overfile)
      character*100 overfile
      include 'examdecl.f'

      ifull(1)=na_i
      ifull(2)=na_j
      ifull(3)=na_k

c Defaults and silence warnings.
      phifilename=' '
      zp(1,1,1)=0.

c Deal with arguments
      if(iargc().eq.0) goto 201
      do i=1,iargc()
         call getarg(i,argument)
         if(argument(1:1).eq.'-')then
         if(argument(1:13).eq.'--objfilename')
     $        read(argument(14:),'(a)',err=201)objfilename
         if(argument(1:2).eq.'-f')
     $        read(argument(3:),'(a)',err=201)phifilename
         if(argument(1:2).eq.'-r')
     $        read(argument(3:),'(f8.4)',err=201)rp
         if(argument(1:2).eq.'-h')goto 203
         if(argument(1:2).eq.'-?')goto 203
         if(argument(1:2).eq.'-o')then
            read(argument(3:),*,end=14,err=14)overfile
            iover=1
 14         continue
         endif
         else
            read(argument(1:),'(a)',err=201)phifilename
         endif
         
      enddo
      goto 202
c------------------------------------------------------------
c Help text
 201  continue
      write(*,*)'=====Error reading command line argument'
 203  continue
 301  format(a,i5)
      write(*,301)'Usage: phiexamine [switches] <phifile>'
      write(*,301)' --objfile<filename>  set name of object data file.'
     $     //' [ccpicgeom.dat'
      write(*,301)' -h -?   Print usage.'
      call exit(0)
 202  continue
      if(lentrim(partfilename).lt.5)goto 203
      end
c*****************************************************************
      subroutine rzbin(nr,rzdist,irzcount,rval,rvalneg,phimax)

      include 'examdecl.f'
      real rzdist(0:nr,iuds(3)),irzcount(nr,iuds(3))
      real rvalneg(0:nr),rval(0:nr)

c r-z binning and plotting
      nz= ixnp(3+1)-ixnp(3)
c Set up r-mesh.
      do j=1,nr
         do i=1,nz
            rzdist(j,i)=0.
            irzcount(j,i)=0
         enddo
         rval(j)=(rs)*(j-0.5)/float(nr)
         rvalneg(j)=-rval(j)
c            write(*,*)j,rval(j)
      enddo

c Bin values
      rmin=1.e10
      phimax=0.
      do k=1,iuds(3)
         z=xn(ixnp(3)+k)
         do j=1,iuds(2)
            do i=1,iuds(1)
               x=xn(ixnp(1)+i)
               y=xn(ixnp(2)+j)
               r=sqrt(x**2+y**2)
               if(r.lt.rmin)rmin=r
               ir=nint(0.5+nr*(r+.00001)/(rs+.00002))
               if(ir.le.nr .and. ir.ge.1)then
                  irzcount(ir,k)=irzcount(ir,k)+1
                  rzdist(ir,k)=rzdist(ir,k)+u(i,j,k)
               endif
               if(u(i,j,k).gt.phimax)phimax=u(i,j,k)                  
            enddo
         enddo
      enddo

c         write(*,*)((j,rval(j),irzcount(j,k),j=1,nr),k=1,1)
      do j=1,nr
         do i=1,iuds(3)
            if(irzcount(j,i).ne.0)then
               rzdist(j,i)=rzdist(j,i)/irzcount(j,i)
            endif
         enddo
      enddo
c Fix up zero values if necessary.
      do j=1,nr
         do i=1,iuds(3)
            if(irzcount(j,i).eq.0)then
               jp=j+1
               if(jp.gt.nr)jp=j
               jm=j-1
               if(jm.lt.1)jm=1
               if(irzcount(jm,i).ne.0)then
                  rzdist(j,i)=rzdist(jm,i)
               elseif(irzcount(jp,i).ne.0)then
                  rzdist(j,i)=rzdist(jp,i)
               else
                  write(*,*)'Failed setting rzdist',j,i,jm,jp
     $                 ,' Too coarse a mesh? rmin=',rmin
c     $                    ,rval(jm),rval(jp)
               endif
            endif
         enddo
      enddo
c Fill in along the axis.
      rval(0)=0.
      rvalneg(0)=0.
      do k=1,iuds(3)
         rzdist(0,k)=rzdist(1,k)
      enddo

      end
c*****************************************************************
      include 'xfig2trace.f'
