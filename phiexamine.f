      program phiexamine

      include 'examdecl.f'

      parameter (ntheta=4,nr=38)
      real thetadist(nr,ntheta),thetaval(ntheta),rval(0:nr)
      real rvalneg(0:nr)
      integer ithetacount(nr,ntheta)

      parameter (nzmax=500)
      real rzdist(0:nr,nzmax),irzcount(nr,nzmax)
      character*1 ppath(0:nr,nzmax)
      character*20 string

      real oneoverr(100),ro(100),cl(100)
c 
c silence warnings:
      fluxfilename=' '

      call examargs
      
      call array3read(phifilename,ifull,iuds,u,ierr)
      if(ierr.eq.1)stop

      write(*,'(a,3i4,$)')'On grid',iuds
      write(*,*)(',',xn(ixnp(kk)+1),xn(ixnp(kk+1)),kk=1,3)

      ifix=1
      call sliceGweb(ifull,iuds,u,Li,zp,
     $     ixnp,xn,ifix,'potential:'//'!Ay!@')

      iplot=2
c Default spherical r plot.
      do k=1,ndims_mesh-1
         if(xn(ixnp(k+1)).ne.xn(ixnp(k+2)))then 
            if(k.eq.1)then
               write(*,*)'Unequal x-y mesh dimensions. No radial plot.'
               iplot=0
               goto 1
            else
               write(*,*)'Unequal y-z mesh dimensions. r-z plot.'
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
         write(*,*)rs
         call pltinit(0.,rs,u(iuds(1)/2,iuds(2)/2,iuds(3)/2),0.)
         call axis()
         call axlabels('radius','potential')
         call charsize(.001,.001)
         phimin=0.
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
                  if(u(i,j,k).lt.phimin)phimin=u(i,j,k)
               enddo
            enddo
         enddo
         do i=1,100
            ro(i)=1.+(rs-1.)*i/100
            oneoverr(i)=phimin/ro(i)
         enddo
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
                     write(*,*)'Failed setting',j,i,im,ip
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

         write(*,*)'r, phi distribution (ir,itheta)'
         write(*,*)nr, ntheta
         do j=1,nr
            write(*,'(10f8.4)')rval(j),(thetadist(j,i),i=1,ntheta)
         enddo

         call charsize(0.,0.)
         call color(2)
         call polyline(ro,oneoverr,100)
         call legendline(.5,.1,0,'Coulomb Potential')
         call pltend()
c%%%%%%%%%%%%%%% End of spherical r-theta plotting %%%%%%%%%%%
      elseif(iplot.eq.1)then
c r-z binning and plotting
         x0=(xn(ixnp(1)+1)+xn(ixnp(2)))/2.
         rx=max(abs(xn(ixnp(1)+1)-x0),abs(xn(ixnp(2))-x0))
         y0=(xn(ixnp(2)+1)+xn(ixnp(3)))/2.
         ry=max(abs(xn(ixnp(2)+1)-x0),abs(xn(ixnp(3))-y0))
         rs=sqrt(rx**2+ry**2)
         nz= ixnp(3+1)-ixnp(3)
         write(*,*)x0,rx,y0,ry,rs,nz
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
         phimax=0.
         do k=1,iuds(3)
            z=xn(ixnp(3)+k)
            do j=1,iuds(2)
               do i=1,iuds(1)
                  x=xn(ixnp(1)+i)
                  y=xn(ixnp(2)+j)
                  r=sqrt(x**2+y**2)
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
                     write(*,*)'Failed setting',j,i,jm,jp
                     write(*,'(20i4)')(irzcount(j,ii),ii=1,iuds(3))
                  endif
               endif
            enddo
         enddo
c Fill in along the axis.
            rval(0)=0.
            rvalneg(0)=0.
         do k=1,iuds(3)
 6          rzdist(0,k)=rzdist(1,k)
         enddo

c Now rzdist(nr,nz) is the rz-distribution of the (average) potential.
c Contour it.
         nc=30
         call fitrange(-phimax,phimax,nc,ipow,fac10,delta,first,xlast)
         write(*,*)ipow,fac10,delta,first,xlast
         do k=1,nc
            cl(k)=(first+k*delta)
         enddo
         nc=abs(2.*phimax)/abs(delta)
         write(*,*)' Phimax', phimax,' fac10',fac10
         write(*,*)' Contours',nc,(cl(k),k=1,nc)
         iconsw=1+64
         call pltinaspect(-rval(nr),rval(nr),
     $        xn(ixnp(3)+1),xn(ixnp(3)+iuds(3)))
         call contourl(rzdist,ppath,nr+1,nr+1,iuds(3),
     $        cl,nc,rval,xn(ixnp(3)+1),iconsw)
         call contourl(rzdist,ppath,nr+1,nr+1,iuds(3),
     $        cl,nc,rvalneg,xn(ixnp(3)+1),iconsw)
         call color(15)
         call axis()
         call axlabels('r','z')
         call fwrite(delta,iwd,2,string)
         call boxtitle('!Af!@-contours ('//string(1:iwd)
     $        //'T!de!d spaced)')
         call gradlegend(-phimax,phimax,-.6,0.,-.6,1.,-.04,.false.) 
         call pltend()

      endif


      end


c*************************************************************
      subroutine examargs()
      include 'examdecl.f'

      do i=1,ndims
         ifull(i)=Li
      enddo

c Defaults and silence warnings.
      phifilename=' '
      fluxfilename=' '
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
