      program phiexamine

      include 'examdecl.f'

      parameter (ntheta=4,nr=38)
      real thetadist(nr,ntheta),thetaval(ntheta),rval(0:nr)
      real rvalneg(0:nr)
      integer ithetacount(nr,ntheta)
      real rzdist(0:nr,na_k),irzcount(nr,na_k)
      real rzmval(na_k),rzmpos(na_k),zpos(na_k)
      character*1 ppath(0:nr,na_k)
      character*20 string

      real oneoverr(100),huckel(100),ro(100),cl(100)
      real xl(2),yl(2)
c 
      call examargs(rp)
         
c      write(*,*)(u(16,16,k),k=1,36) 
      ied=1
      call array3read(phifilename,ifull,iuds,ied,u,ierr)
      if(ierr.eq.1)stop

      write(*,'(a,3i4,$)')'On grid',iuds
      write(*,'(a,2f7.2,a,2f7.2,a,2f7.2)')
     $     (',',xn(ixnp(kk)+1),xn(ixnp(kk+1)),kk=1,3)

c      write(*,*)(u(16,16,k),k=1,36) 
      ifix=1
c      call noeye3d(0)
      call sliceGweb(ifull,iuds,u,na_m,zp,
     $     ixnp,xn,ifix,'potential:'//'!Af!@')
c      call sliceGcont(ifull,iuds,u,na_m,zp,
c     $     ixnp,xn,ifix,'potential:'//'!Af!@')

      iplot=2
c Default spherical r plot.
      do k=1,ndims_mesh-1
         if(xn(ixnp(k+1)).ne.xn(ixnp(k+2)))then 
            if(k.eq.1)then
               write(*,*)'Unequal x-y mesh dimensions. No radial plot.'
               iplot=0
               goto 1
            else
c               write(*,*)'Unequal y-z mesh dimensions. r-z plot.'
               iplot=1
            endif
         endif
      enddo
 1    continue
      
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if(iplot.eq.2)then
c r-theta spherical plots 
c For boxes
         rs=sqrt(3.)*rs
         do j=1,nr
            do i=1,ntheta
               thetadist(j,i)=0.
               ithetacount(j,i)=0
               if(i.eq.1)thetaval(i)=-1.+ (i-0.5)/float(ntheta)
            enddo
            rval(j)=1.+(rs-1.)*(j-0.5)/float(nr)
c         write(*,*)j,rval(j)
         enddo

c plot potential versus radius.
c         write(*,*)rs
         call pltinit(0.,rs,u(iuds(1)/2,iuds(2)/2,iuds(3)/2),0.)
         call axis()
         call axlabels('radius','potential')
         call charsize(.005,.005)
         phimin=0.
         rmin=1.e10
         redge=0.
         do k=1,iuds(3)
            do j=1,iuds(2)
               do i=1,iuds(1)
                  x=xn(ixnp(1)+i)
                  y=xn(ixnp(2)+j)
                  z=xn(ixnp(3)+k)
                  r=sqrt(x**2+y**2+z**2)
                  call polymark(r,u(i,j,k),1,10)
                  if(r.gt.rs .and. u(i,j,k).ne.0)then
c                  write(*,'(4f12.6,3i3)')x,y,z,u(i,j,k),i,j,k
                  endif
                  if(u(i,j,k).le.phimin.and.r.le.rmin)then
                     rmin=r
                     if(abs(u(i,j,k)-phimin).gt.1.e-5
     $                    .or.r.gt.redge)redge=r
                     phimin=u(i,j,k)
                  endif
               enddo
            enddo
         enddo
c If redge is nearly 1, guess it is exactly 1.
         if(abs(redge-1.).lt..1)redge=1.
         do i=1,100
c            ro(i)=1.+(rs-1.)*i/100
            ro(i)=rmin+(rs-rmin)*(i-1.)/(100-1.)
            oneoverr(i)=phimin*redge/ro(i)
            huckel(i)=phimin*(redge/ro(i))*
     $           exp(-max(0.,(ro(i)-redge))/debyelen)
         enddo
         write(*,*)rmin,redge,phimin,debyelen
c Average together.
         denmin=0.
         do k=1,iuds(3)
            do j=1,iuds(2)
               do i=1,iuds(1)
                  x=xn(ixnp(1)+i)
                  y=xn(ixnp(2)+j)
                  z=xn(ixnp(3)+k)
                  r=sqrt(x**2+y**2+z**2)
                  c=z/r
                  itheta=nint(0.5+ntheta*(c+1.000001)/2.00001)
                  ir=nint(0.5+nr*(r-.999999)/(rs-1.00001))
                  if(ir.le.nr .and. ir.ge.1)then
                     ithetacount(ir,itheta)=ithetacount(ir,itheta)+1
                     thetadist(ir,itheta)=thetadist(ir,itheta)+u(i,j,k)
                  endif
                  if(u(i,j,k).lt.denmin)denmin=u(i,j,k)
               enddo
            enddo
         enddo
         do j=1,nr
            do i=1,ntheta
               if(ithetacount(j,i).ne.0)then
                  thetadist(j,i)=thetadist(j,i)/ithetacount(j,i)
               endif
            enddo
         enddo
         do j=1,nr
            do i=1,ntheta
               if(ithetacount(j,i).eq.0)then
                  ip=mod(i,ntheta)+1
                  im=mod(i-2,ntheta)+1
                  if(ithetacount(j,ip).ne.0)then
                     thetadist(j,i)=thetadist(j,ip)
                  elseif(ithetacount(j,im).ne.0)then
                     thetadist(j,i)=thetadist(j,im)
                  else
                     write(*,*)'Failed setting',j,i,im,ip,ntheta
                     write(*,'(20i4)')(ithetacount(j,ii),ii=1,ntheta)
                  endif
               endif
            enddo
         enddo

c Plot binned data:
         do i=1,ntheta
            call color(mod(i,16))
            call polyline(rval(1),thetadist(1,i),nr)
         enddo

c Uncomment for written output.
c         write(*,*)'r, phi distribution (ir,itheta)'
c         write(*,*)nr, ntheta
c         do j=1,nr
c            write(*,'(10f8.4)')rval(j),(thetadist(j,i),i=1,ntheta)
c         enddo

         call charsize(0.,0.)
         call color(2)
         call dashset(2)
         call winset(.true.)
         call polyline(ro,oneoverr,100)
         call legendline(.5,.1,0,'Coulomb Potential')
         call color(3)
         call dashset(3)
         call polyline(ro,huckel,100)
         call legendline(.5,.15,0,'Yukawa Potential')
         call dashset(0)
         call pltend()
c%%%%%%%%%%%%%%% End of spherical r-theta plotting %%%%%%%%%%%
      elseif(iplot.eq.1)then
c r-z binning and plotting
         x0=(xn(ixnp(1)+1)+xn(ixnp(2)))/2.
         rx=max(abs(xn(ixnp(1)+1)-x0),abs(xn(ixnp(2))-x0))
         y0=(xn(ixnp(2)+1)+xn(ixnp(3)))/2.
         ry=max(abs(xn(ixnp(2)+1)-x0),abs(xn(ixnp(3))-y0))
         rs=sqrt(rx**2+ry**2)

         call cyl3bin(nr,rzdist,irzcount,rval,rvalneg,phimax,u)

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
         call pltinaspect(-rval(nr),rval(nr),
     $        xn(ixnp(3)+1),xn(ixnp(3)+iuds(3)))
         call color(15)
         call contourl(rzdist,ppath,nr+1,nr+1,iuds(3),
     $        cl,nc,rval,xn(ixnp(3)+1),iconsw)
         call contourl(rzdist,ppath,nr+1,nr+1,iuds(3),
     $        cl,nc,rvalneg,xn(ixnp(3)+1),iconsw)
         call color(15)
         call axis()
         call axlabels('r','z')
c         write(*,*)'ipow',ipow
         call fwrite(delta,iwd,max(-ipow,0)+1,string)
         call boxtitle('!Af!@-contours ('//string(1:iwd)
     $        //'T!de!d spaced)')
         call gradlegend(-phimax,phimax,-.6,0.,-.6,1.,-.04,.false.) 
c Indicate rectangle limits.
         yl(1)=xn(ixnp(3)+1)
         yl(2)=xn(ixnp(3)+iuds(3))
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
            call polyline(xl,yl,2)
         endif
         call dashset(0)
         xl(1)=0.
         xl(2)=0.
         call polyline(xl,yl,2)
c         call fwrite(vd,iwd,2,string)
c         call jdrwstr(.02,.31,'v!dd!d='//string(1:iwd),1.)
c         call fwrite(debyelen,iwd,1,string)
c         call jdrwstr(.02,.25,'!Al!@!dDe!d='//string(1:iwd),1.)
c         call fwrite(phip,iwd,2,string)
c         call jdrwstr(.02,.28,'!Af!@!dp!d='//string(1:iwd),1.)
         write(*,'(''vd='',f7.3,'' debyelen='',f7.3,''  phip='',f7.3)')
     $        vd,debyelen,phip
         call pltend()
      endif

c Now we want to find the r-position of the peak potential at each
c z. 

      write(*,*)'Not finding potential peaks, but could do so ...'
      call exit(1)

      iz0=0
      do i=1,iuds(3)
         rzmval(i)=0.
         rzmpos(i)=0.
         jrm=0
         zpos(i)=xn(ixnp(3)+i)/debyelen
         if(xn(ixnp(3)+i).gt.0. .and. iz0.eq.0)iz0=i
         do j=1,nr
            if(rzmval(i).lt.rzdist(j,i))then
               jrm=j
               rzmval(i)=rzdist(j,i)
               rzmpos(i)=rval(j)/debyelen
            endif
         enddo
      enddo
      izt=iuds(3)+1-iz0
c      call autoplot(rzmpos(iz0),xn(ixnp(3)+iz0),izt)
      call pltinaspect(0.,rval(nr)/debyelen,0.,xn(ixnp(3)+iuds(3))
     $     /debyelen)
      call axis()
      call polyline(rzmpos(iz0),zpos(iz0),izt)
      call axlabels('peak-!Af!@ radius [/!Al!@!dDe!d]'
     $     ,'z/!Al!@!dDe!d')
      call pltend
      open(22,file='peakfile.dat',status='unknown',err=101)
      close(22,status='delete')
      open(22,file='peakfile.dat',status='new',err=101)
      write(22,'(a,f5.2)')'legend: M=',vd
      write(22,*)izt
      do i=0,izt-1
c         write(22,*)xn(ixnp(3)+iz0+i),rzmpos(iz0+i)
c         write(22,*)rzmpos(iz0+i),xn(ixnp(3)+iz0+i)
         write(22,*)rzmpos(iz0+i),zpos(iz0+i)
      enddo
      close(22)

      call exit(0)
 101  write(*,*)'Error writing output'
      end

c*************************************************************
      subroutine examargs(rp)
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
c Bin rectangular mesh values into radial (cylindrical) 
c mesh values and average in angle.
      subroutine cyl3bin(nr,rzdist,irzcount,rval,rvalneg,phimax,up)
      implicit none
      integer nr
      include 'examdecl.f'
      real rzdist(0:nr,na_k),irzcount(nr,na_k)
      real rval(0:nr),rvalneg(0:nr)
      real phimax
c Passed quantity whose mesh is declared in examdecl. E.g. potential.
c It is defined for iuds.
      real up(ifull(1),ifull(2),ifull(3))
c Local variables:
      real rmin,x,y,z,r
      integer i,j,k,ir,nz,jp,jm

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
                     rzdist(ir,k)=rzdist(ir,k)+up(i,j,k)
                  endif
                  if(up(i,j,k).gt.phimax)phimax=up(i,j,k)                  
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
     $                    ,' Too coarse a mesh? rmin=',rmin
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
