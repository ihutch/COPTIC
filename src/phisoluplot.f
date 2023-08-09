      subroutine vaccheck(ifull,iuds,cij,u,thetain,nth,rs,ltestplot)
      integer ifull(*),iuds(*)
      real cij(*),u(*),thetain
      logical ltestplot
      include 'ndimsdecl.f'
      include '3dcom.f'
! Silence warning:
      r=rs
!      write(*,*)r,obj_geom(otype,2),ltestplot
! Do some analytic checking of the case with a fixed potential sphere
! inside an object 2 sphere with its radius and BC.
      if(obj_geom(otype,2).ne.257)return
! Also write out some data for checking.
      rc=obj_geom(oradius,1)
      phip=-obj_geom(oabc+2,1)/obj_geom(oabc,1)
! Calculate the object-2 outer boundary correction to phiinf
      a=obj_geom(oabc,2)
      b=obj_geom(oabc+1,2)
      c=obj_geom(oabc+2,2)
      r=obj_geom(oradius,2)
      if(r.eq.0)return
      phiinf=-(a*phip*rc/r -b*phip*rc/r**2+c)
     $     /(-a*rc/r + a + b*rc/r**2)
!         write(*,*)'Vacuum phiinfty=',phiinf,'  rabc=',r,a,b,c,oabc
      write(*,'(a,f7.4,a,f8.4,a,f7.4,a,f8.4)')
     $     ' Vacuum solution: rc=',rc,' phip=',phip,
     $     ' inside r=',r,' phi_infty=',phiinf
      call spherecheck(ifull,iuds,u,phip,rc,phiinf)
      if(ltestplot)then
! Plot some of the initial-solver data.
         call solu3plot(ifull,iuds,u,cij,
     $        phip,phiinf,rc,thetain,nth)
         write(*,*)'Return from solu3plot.'
         stop
      endif
      end

!********************************************************************
      subroutine circleplot(xc,yc,r)
      integer nangle
      parameter (nangle=100,pi=3.1415927)

      ip=0
      do i=1,nangle
         t=2.*pi*(i-1)/(nangle-1)
         xp=xc+r*cos(t)
         yp=yc+r*sin(t)
         call vecw(xp,yp,ip)
         if(i.eq.1)ip=1
      enddo
      end
!********************************************************************
      subroutine spherecheck(ifull,iuds,u,phi,rc,phiinf)
! Do some analytic checking of the case with a fixed potential sphere
! inside a logarithmic derivative boundary condition. 1/r solution.
      include 'ndimsdecl.f'
      include 'meshcom.f'
      integer iuds(ndimsmax),ifull(ndimsmax)
      real u(ifull(1),ifull(2),ifull(3))
      errmax=0.
      errvar=0.
      count=0
      do i=2,iuds(1)-1
         do j=2,iuds(2)-1
            do k=2,iuds(3)-1
               xr=xn(i)
               yr=xn(iuds(1)+j)
               zr=xn(iuds(1)+iuds(2)+k)
               r=sqrt(xr**2+yr**2+zr**2)
               if(abs(u(i,j,k)).gt.1.e-4 .and.
     $              r.ge.rc)then
                  count=count+1.
                  e=u(i,j,k)-(phiinf+(phi-phiinf)*rc/r)
!                     error(i,j,k)=e
                  errvar=errvar+e**2
                  if(abs(e).gt.abs(errmax))then
                     errmax=e
                     iem=i
                     jem=j
                     kem=k
                  endif
               endif
            enddo
         enddo
      enddo
!      write(*,*)errvar,count,errmax,iuds,u(3,3,3)
      errvar=errvar/count
      write(*,'(a,f8.5,a,f10.6)')' Maximum nodal vacuum phi error='
     $     ,errmax,'   Standard Deviation=',sqrt(errvar)
!      stop
      if(sqrt(errvar).gt.0.4)then
         write(*,*)
     $     '********** Standard Deviation too large *************'
      endif
! Rarely needed printout of u:
!      iform=7
!      uscale=10000000.
!      call udisplay(ndims,u,ifull,iuds,iform,uscale)
! Write some data out for cross checking
! Delete the file first to help with nfs problems.
         open(21,file='smt.round',status='unknown')
         close(21,status='delete')
         open(21,file='smt.round',status='unknown')
         write(21,*)iuds(1),iuds(2),iuds(3)/2
         write(21,'(6e12.4)')
     $        ((u(i,j,iuds(3)/2),i=1,iuds(1)/2),j=1,iuds(2)/2)
         close(21)
! This file has full accuracy.
         open(20,file='smt.out',status='unknown')
         close(20,status='delete')
         open(20,file='smt.out',status='unknown')
         write(20,*)iuds(1),iuds(2),iuds(3)/2
         write(20,*)
     $        ((u(i,j,iuds(3)/2),i=1,iuds(1)/2),j=1,iuds(2)/2)
         close(20)
      end
!***************************************************************
! Master plotting routine.
      subroutine solu3plot(ifull,iuds,u,cij,
     $     phi,phiinf,rc,thetain,nth)
      real thetain
      include 'ndimsdecl.f'
      integer ifull(ndimsmax),iuds(ndimsmax)
      integer iLs(ndims+1)
      real u(*),cij(*)
      iLs(1)=1
      do i=1,3
         iLs(i+1)=iLs(i)*ifull(i)
      enddo
      call fixedline(ifull,iuds,u,phi,phiinf,rc)
!      call slicesolu(ifull,iuds,u,cij)
      call gradradial(ifull,iLs,u,cij,phi,phiinf,rc,thetain,nth)
!      call phiradial(ifull,iLs,cij,phi,phiinf,rc,thetain,nth,rs)
      call eremesh(ifull,iLs,iuds,u,cij,phi,phiinf,rc)
      end
!***************************************************************
! Packaged version of plotting.
      subroutine slicesolu(ifull,iuds,u,cij)

      include 'ndimsdecl.f'
      parameter (nd2=ndimsmax*2)
      integer ifull(ndimsmax),iuds(ndimsmax)
      real cij(*),u(*)
      integer ifmax
      parameter (ifmax=300,Li=ifmax)
      real cijp(2*ndims+1,ifmax,ifmax)
      real zp(ifmax,ifmax)
      save cijp,zp ! Suppress warnings about size.
      include 'meshcom.f'
      include '3dcom.f'

! Plotting slices.
      ifix=3
      call slice3web(ifull,iuds,u,cij,Li,zp,cijp,ixnp,xn,ifix,
     $     'potential:'//'!Ay!@',1)
      end
!*******************************************************************
! Packaged version of plotting.
      subroutine gradradial(ifull,iLs,u,cij,phi,phiinf,rc,thetain,nth)
!     $     ,rs)
      include 'ndimsdecl.f'
      parameter (nd2=ndimsmax*2)
      integer ifull(ndimsmax)
!     integer iuds(ndims)
      real cij(2*ndimsmax+1,ifull(1),ifull(2),ifull(3))
      real u(ifull(1),ifull(2),ifull(3))
      integer ifmax
      parameter (ifmax=500,Li=ifmax)
      real uanal(ifmax)
      real rfield(ifmax),tfield(ifmax)
      real rprime(ifmax),xprime(ndimsmax,ifmax)
      real xfrac(ndims),xff(ndims),upnd(ndimsmax,ifmax)
      real upsimple(ndimsmax,ifmax)
      integer itemp(ndimsmax)
      real rsimple(ifmax),region(ifmax)
      integer iLs(ndims+1)

      character*40 form1
      include 'meshcom.f'
      include '3dcom.f'      
!-------------------------------------------------------------------
! Field plotting for various different lines:
! Spherical angles in 3-D
      graderr=0.
      tgraderr=0.
      do iti=1,nth
         theta=thetain*iti
         varphi=0.
         ct=cos(theta)
         st=sin(theta)
         cp=cos(varphi)
         sp=sin(varphi)
!         write(*,*)'rs=',rs
         do i=1,Li
!            zero(i)=0.
!            rp=rs*i/(Li-1)
!            rp=.99999*rs*i/Li
            rp=.99999*5.*i/Li
            rprime(i)=rp
! Coordinates relative to center of first object (sphere).
            xprime(1,i)=rp*st*cp + obj_geom(ocenter,1)  
            xprime(2,i)=rp*st*sp + obj_geom(ocenter+1,1)  
            xprime(3,i)=rp*ct + obj_geom(ocenter+2,1)  
            
            iregion=insidemask(ndims,xprime(1,i))
! Calculate fractional mesh positions of this point, always positive.
! Thus the origin of the box is below point in all dimensions.
            do id=1,ndims
! Offset to start of idf position array.
               ioff=ixnp(id)
! xn is the position array for each dimension arranged linearly.
! Find the index of xprime in the array xn:
               ix=interp(xn(ioff+1),ixnp(id+1)-ioff,xprime(id,i),xm)
!                  ix=xm
               if(xm.ge.(ixnp(id+1)-ixnp(id)))then
                  write(*,*)'****** ccpicplot xm too big',id,xm
                  xm=ixnp(id+1)-ixnp(id)-1.e-3
               endif
               xff(id)=xm
               xfrac(id)=xm-ix
               itemp(id)=ix
! attempt to fix getsimple.
               if(itemp(id).le.1)then
                  itemp(id)=itemp(id)+1
                  xfrac(id)=xfrac(id)-1.
               elseif(itemp(id).ge.ixnp(id+1)-ixnp(id)-1)then
                  itemp(id)=itemp(id)-1
                  xfrac(id)=xfrac(id)+1.
               endif
            enddo
            
! Get the ndims field components at this point.
!               write(*,*)'xfrac',(xfrac(kk),kk=1,ndims)
!     $              ,(xprime(kk,i),kk=1,ndims)
            do idf=1,ndims
               ioff=ixnp(idf)
! New approach mirrors padvnc.
               call getfield(
     $              cij(nd2+1,1,1,1)
     $              ,u,iLs,xn(ixnp(idf)+1),idf
     $              ,xff,iregion,upnd(idf,i))

               call getsimple3field(
     $              ndims,u(itemp(1),itemp(2),itemp(3))
     $              ,iLs,xn(ioff+itemp(idf)),idf
     $              ,xfrac,upsimple(idf,i))

               if(.not.abs(upnd(idf,i)).lt.1.e6)then
                  write(*,*)'Field BIG/weird!',upnd(idf,i),i,idf,xff
               endif

            enddo
            if(ifmax.lt.40)write(*,'(3f8.4,i3,a,3f8.4)')xff,
     $           iregion,' Field',(upnd(k,i),k=1,3)

! Radial component of field
            rfield(i)=
     $           upnd(1,i)*st*cp +
     $           upnd(2,i)*st*sp +
     $           upnd(3,i)*ct
            rsimple(i)=
     $           upsimple(1,i)*st*cp +
     $           upsimple(2,i)*st*sp +
     $           upsimple(3,i)*ct
! Tangential component (magnitude) of field
            tfield(i)=-sqrt(max(0.,
     $           upnd(1,i)**2+upnd(2,i)**2+upnd(3,i)**2
     $           -rfield(i)**2))

! Analytic comparison.
            uanal(i)=(phi-phiinf)*rc/(rprime(i)**2)
!               write(*,'(''i,rprime,rfield,uanal(i)'',i4,4f10.5)')
!     $              i,rprime(i),rfield(i),uanal(i)
!     $              ,rsimple(i)
! Region of point
            region(i)=insideall(ndims,xprime(1,i))
! iregion is the masked region.
            if(iregion.eq.2)then 
               gerr=abs(rfield(i)-uanal(i))
               if(gerr.gt.graderr)then 
                  graderr=gerr
                  rgraderr=rp
                  tgerr=theta
               endif
               terr=abs(tfield(i))
               if(terr.gt.tgraderr)then
                  tgraderr=terr
                  rtgraderr=rp
                  ttgraderr=theta
               endif
            endif
            region(i)=-0.1*region(i)
         enddo
!         write(*,*)'Ended uprime calculation'
!         write(*,'(10f8.4)')rfield
         call dashset(0)
         if(iti.eq.1)then
            call autoplot(rprime,rfield(1),Li)
            call axis2()
            call axlabels('r','E!d(r)!d=-grad!Af!@')
            call legendline(.5,.1,0,'getfield')
            call winset(.true.)
            call dashset(4)
            call color(ired())
            call polyline(rprime,uanal,Li)
            call legendline(.5,.05,0,'analytic')
            call winset(.false.)
            call color(iblue())
            call dashset(3)
            call polyline(rprime,rsimple(1),Li)
            call legendline(.5,.15,0,'simplefield')            
            call color(13)
            call dashset(0)
            call legendline(.5,.25,0,'region/10')
            call color(idarkgreen())
            call dashset(2)
            call legendline(.5,.2,0,'tangential')
         else
            call polyline(rprime,rfield(1),Li)
         endif
         call dashset(2)
         call color(idarkgreen())
         call polyline(rprime,tfield(1),Li)
         call color(13)
         call dashset(0)
         call polyline(rprime,region,Li)
         call color(15)
         call winset(.false.)
         form1='!Aq!@='
         call fwrite(180*theta/3.1415926,iwdth,1,form1(7:))
         call jdrwstr(.01,.05+.03*iti,form1,1.)
      enddo
      call dashset(0)
      write(*,*)'Along the fine-stepped diagnostic line(s):'
      write(*,*)'Max radial field error=    ',graderr,'  at r,th='
     $     ,rgraderr,tgerr
      write(*,*)'Max tangential field error=',tgraderr,'  at r,th='
     $     ,rtgraderr ,ttgraderr
      call pltend()
      end
!-------------------------------------------------------------------
!-------------------------------------------------------------------
!***************************************************************
! Packaged version of potential plotting.
      subroutine phiradial(ifull,iLs,cij,
     $     phi,phiinf,rc,thetain,nth,rs)
      include 'ndimsdecl.f'
      parameter (nd2=ndimsmax*2)
      integer ifull(ndimsmax)
      real cij(*)
      integer ifmax
      parameter (ifmax=1000,Li=ifmax)
      real uanal(ifmax)
      real rfield(ifmax),tfield(ifmax)
      real rprime(ifmax),xprime(ndimsmax,ifmax)
      real xfrac(ndims),xff(ndims)
!      real upsimple(ndims,ifmax),upnd(ndims,ifmax)
      integer itemp(ndims)
      real rsimple(ifmax),region(ifmax)
      integer iLs(ndims+1)

      character*40 form1
      include 'meshcom.f'
      include '3dcom.f'
! Potential plotting for different lines:
! Spherical angles in 3-D
      phierr=0
      do iti=1,nth
         theta=thetain*iti
         varphi=0.
         ct=cos(theta)
         st=sin(theta)
         cp=cos(varphi)
         sp=sin(varphi)
!         write(*,*)'Starting x calculation'
         do i=1,Li
!            zero(i)=0.
!            rp=rs*i/(Li-1)
            rp=0.99999*rs*i/Li
            rprime(i)=rp
! Coordinates relative to center of first object (sphere).
            xprime(1,i)=rp*st*cp + obj_geom(ocenter,1)  
            xprime(2,i)=rp*st*sp + obj_geom(ocenter+1,1)  
            xprime(3,i)=rp*ct + obj_geom(ocenter+2,1)  
            
            iregion=insidemask(ndims,xprime(1,i))
! Calculate fractional mesh positions of this point, always positive.
! Thus the origin of the box is below point in all dimensions.
            do id=1,ndims
! Offset to start of idf position array.
               ioff=ixnp(id)
! xn is the position array for each dimension arranged linearly.
! Find the index of xprime in the array xn:
               ix=interp(xn(ioff+1),ixnp(id+1)-ioff,xprime(id,i),xm)
               if(ix.eq.0.or.ix.ge.ixnp(id+1)-ioff)then
                  write(*,*)'ccpicplot interp outside range'
     $                 ,xprime(id,i)
     $                 ,ioff,(xn(ioff+k),k=1,ixnp(id+1)-ioff),' end '
               endif
!                  ix=xm
               xff(id)=xm
               xfrac(id)=xm-ix
               itemp(id)=ix
            enddo

!               write(*,*)'xfrac',(xfrac(kk),kk=1,ndims)
!     $              ,(xprime(kk,i),kk=1,ndims)            

! Analytic comparison.
            uanal(i)=phiinf+(phi-phiinf)*rc/rprime(i)
            if(uanal(i).lt.phi)uanal(i)=phi
!               write(*,'(''i,rprime,rfield,uanal(i)'',i4,4f10.5)')
!     $              i,rprime(i),rfield(i),uanal(i)
!     $              ,rsimple(i)

            rfield(i)=getpotential(u,cij,iLs,xff,iregion,2)
            rsimple(i)=getpotential(u,cij,iLs,xff,iregion,1)
            tfield(i)=getpotential(u,cij,iLs,xff,iregion,3)
! Region of point
            region(i)=insideall(ndims,xprime(1,i))
!            write(*,*)'region=',region(i)
            if(iregion.eq.2)then 
               perr=abs(rfield(i)-uanal(i))
               if(perr.gt.phierr)then 
                  phierr=perr
                  rphierr=rp
                  tphierr=theta
               endif
            endif
            region(i)=-0.1*region(i)
         enddo
         call dashset(0)
         if(iti.eq.1)then
            call pltinit(0.,rprime(Li),phi*1.02,0.)
            call axis()
            call axis2()
            call axlabels('r','!Af!@')
            call legendline(.5,.1,0,'getpotential')
            call winset(.true.)
            call polyline(rprime,rfield(1),Li)
            call dashset(4)
            call color(ired())
            call polyline(rprime,uanal,Li)
            call legendline(.5,.05,0,'analytic')
            call winset(.false.)
            call color(iblue())
            call dashset(3)
            call polyline(rprime,rsimple(1),Li)
            call legendline(.5,.15,0,'simple')            
            call color(13)
            call dashset(0)
            call legendline(.5,.25,0,'region/10')
            call color(idarkgreen())
            call dashset(2)
            call legendline(.5,.2,0,'nearest')
         else
            call polyline(rprime,rfield(1),Li)
         endif
         call dashset(2)
         call color(idarkgreen())
         call polyline(rprime,tfield(1),Li)
         call color(13)
         call dashset(0)
         call polyline(rprime,region,Li)
         call color(15)
         call winset(.false.)
         form1='!Aq!@='
         call fwrite(180*theta/3.1415926,iwdth,1,form1(7:))
         call jdrwstr(.01,.05+.03*iti,form1,1.)
      enddo
      call dashset(0)
      write(*,*)'Max potential error=',phierr,'  at',rphierr,tphierr
      call pltend()
      end
!********************************************************************
      subroutine fixedline(ifull,iuds,u,phi,phiinf,rc)
      integer iuds(3),ifull(3)
      real u(ifull(1),ifull(2),ifull(3))
      include 'ndimsdecl.f'
      include 'meshcom.f'
      integer ifmax
      parameter (ifmax=1000,Li=ifmax)
      real z(ifmax),xp(ifmax)
      write(*,*)'In solu3plot',ifull,iuds,rc
      n0=max(iuds(2)/2,1)
      n1=max(iuds(3)/2,1)
      idf=3
      do i=1,iuds(1)
         xr=xn(i)
         yr=xn(iuds(1)+n0)
         zr=xn(iuds(1)+iuds(2)+n1)
         r=sqrt((xr-.0)**2+(yr-.0)**2+(zr-.0)**2)
         z(i)=phiinf+(phi-phiinf)*rc/r
         xp(i)=xr
!         write(*,*)i,xr
      enddo
      call autoplot(xp,u(1,n0,n1),iuds(1))
      call polymark(xp,u(1,n0,n1),iuds(1),1)
      call color(ibrickred())
      call winset(.true.)
      call polyline(xp,z,iuds(1))
      call winset(.false.)
      call color(15)
      call axlabels('x mesh (at half mesh in y, z)',
     $     'u(x) and !Af!@(r!dc!d)r!dc!d/r')
      call pltend()
      end
!********************************************************************
!***************************************************************
! Packaged version of potential error finding and plotting.
      subroutine eremesh(ifull,iLs,iuds,u,cij,phi,phiinf,rc)
      include 'ndimsdecl.f'
      parameter (nd2=ndimsmax*2)
      integer ifull(ndimsmax),iuds(ndimsmax)
      real cij(*)
      real u(*)
      integer ifmax
      parameter (ifmax=100,Li=ifmax)
      real xprime(ndimsmax),xff(ndims)
      real ere(ifmax,ifmax),xcont(ifmax),ycont(ifmax)
      character cwork(ifmax,ifmax)
      real zclv(30)
      character*6 xtit,ytit

      integer iLs(ndims+1)
      include 'meshcom.f'
!      include '3dcom.f'

! Silence warnings:,thetain,nth,rs)
!      x0=13.8
!      y0=12.8
!      x1=17.8
!      y1=16.8
      x0=12.*iuds(1)/32
      y0=12.*iuds(2)/32
      z0=12.*iuds(3)/32
      x1=21.*iuds(1)/32
      y1=21.*iuds(2)/32
      z1=21.*iuds(3)/32

      idf=2
      id1=mod(idf,3)+1
      write(xtit,'(''axis-'',i1)')id1
      id2=mod(idf+1,3)+1
      write(ytit,'(''axis-'',i1)')id2
      ide=id1

      k1=1
      k2=ifmax
      errvtot=0.
 1    errvar=0.
      count=0.
      errmax=0.
      kerrmax=k1
      do k=k1,k2
         f3=(k-0.999)/(ifmax-0.998)
         xff(idf)=(1.-f3)*z0 + f3*z1
!         write(*,*)'f3,xff(idf)',f3,xff(idf),z0,z1,errmax
         ixff=int(xff(idf))
         ff=xff(idf)-ixff
         xprime(idf)=xn(ixnp(idf)+ixff)*(1.-ff)
     $        +xn(ixnp(idf)+ixff+1)*ff
      do j=1,ifmax
         f2=(j-0.999)/(ifmax-0.998)
         xff(id2)=(1.-f2)*y0 + f2*y1
         ixff=int(xff(id2))
         ff=xff(id2)-ixff
         xprime(id2)=xn(ixnp(id2)+ixff)*(1.-ff)
     $        +xn(ixnp(id2)+ixff+1)*ff
         ycont(j)=xprime(id2)
         do i=1,ifmax
            f1=(i-0.999)/(ifmax-0.998)
            xff(id1)=(1.-f1)*x0 + f1*x1
            ixff=int(xff(id1))
            ff=xff(id1)-ixff
            xprime(id1)=xn(ixnp(id1)+ixff)*(1.-ff)
     $           +xn(ixnp(id1)+ixff+1)*ff
            xcont(i)=xprime(id1)
            iregion=insidemask(ndims,xprime)
!            write(*,*),i,j,iregion,xprime

            call getfield(
!     $           cij(nd2+1,1,1,1)
     $           cij(nd2+1)
     $           ,u,iLs,xn(ixnp(ide)+1),ide
     $           ,xff,iregion,ere(i,j))

! Analytic calculation of field. E= phi _r_ /r^3
! Assuming sphere is centered on origin. 
            r2=(xprime(idf)**2+xprime(id1)**2+xprime(id2)**2)
            if(ide.eq.id1)then
               anal=xprime(id1)*(phi-phiinf)*rc/(r2*sqrt(r2))
            elseif(ide.eq.id2)then
               anal=xprime(id2)*(phi-phiinf)*rc/(r2*sqrt(r2))
            else
               anal=xprime(idf)*(phi-phiinf)*rc/(r2*sqrt(r2))
            endif
            if(iregion.ne.2)anal=0.
            ere(i,j)=ere(i,j)-anal
            err=abs(ere(i,j))
            errvar=errvar+err*err
            count=count+1
            if(err.gt.errmax)then
               errmax=err
               errs=ere(i,j)
               f3errmax=f3
               kerrmax=k
            endif
         enddo
      enddo

      enddo
      errvar=errvar/(count-1)
      if(k1.ne.k2)then
         k1=kerrmax
         k2=k1
         errvtot=errvar
         errmtot=errmax
         goto 1
      endif
      icl=20
      call fitrange(0.,errmtot,icl/2,ipow,fac10,delta,first,flast)
      icsw=1
!      write(*,*)errmtot,delta,first,delta
      do i=1,icl
         zclv(i)=float(i-icl/2)*delta
      enddo
      write(*,'(a,f8.4)')'Contour spacing=',delta
      call pltinit(xcont(1),xcont(ifmax),ycont(1),ycont(ifmax))
      call contourl(ere,cwork,ifmax,ifmax,ifmax,zclv,icl,
     $           xcont,ycont,icsw) 
      call axis()
      call axlabels(xtit,ytit)
      call winset(.true.)
      do j=1,iuds(id2)
         do i=1,iuds(id1)
            xp=xn(ixnp(id1)+i)
            yp=xn(ixnp(id2)+j)
            call polymark(xp,yp,1,1)
         enddo
      enddo
      write(*,'(2f8.4,a,f8.4,a,f8.4,a,2f8.4)')xff(idf),xprime(idf)
     $     ,' Max Field Error=',errs,' dir ',ide
     $     ,', sd=',sqrt(errvar),sqrt(errvtot)
!     $     ,' Contours',zclv(2)-zclv(1)
      icode=ieye3d()
      if(icode.eq.ichar('u'))then
         k1=k1+1
         k2=k1
!         write(*,*)xff(idf),xprime(idf)
         goto 1 
      elseif(icode.eq.ichar('d'))then
         k1=k1-1
         k2=k1
!         xff(idf)=xff(idf)-.1
!         write(*,*)xff(idf),xprime(idf)
         goto 1
      endif

      end
!*******************************************************************
