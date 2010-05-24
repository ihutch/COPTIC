      subroutine vaccheck(ifull,iuds,cij,u,thetain,nth,rs,ltestplot)
      logical ltestplot
      include '3dcom.f'
c Do some analytic checking of the case with a fixed potential sphere
c inside an object 2 sphere with its radius and BC.
c Also write out some data for checking.
      rc=obj_geom(oradius,1)
      phip=-obj_geom(oabc+2,1)/obj_geom(oabc,1)
c Calculate the object-2 outer boundary correction to phiinf
      a=obj_geom(oabc,2)
      b=obj_geom(oabc+1,2)
      c=obj_geom(oabc+2,2)
      r=obj_geom(oradius,2)
      if(r.eq.0)return
      phiinf=-(a*phip*rc/r -b*phip*rc/r**2+c)
     $     /(-a*rc/r + a + b*rc/r**2)
c         write(*,*)'Vacuum phiinfty=',phiinf,'  rabc=',r,a,b,c
      write(*,'(a,f7.4,a,f8.4,a,f7.4,a,f8.4)')
     $     ' Vacuum solution: rc=',rc,' phip=',phip,
     $     ' inside r=',r,' phi_infty=',phiinf
      call spherecheck(ifull,iuds,u,phip,rc,phiinf)
      if(ltestplot)then
c Plot some of the initial-solver data.
         call solu3plot(ifull,iuds,u,cij,
     $        phip,phiinf,rc,thetain,nth,rs)
         write(*,*)'Return from solu3plot.'
         stop
      endif
      end

c********************************************************************
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
c********************************************************************
      subroutine spherecheck(ifull,iuds,u,phi,rc,phiinf)
c Do some analytic checking of the case with a fixed potential sphere
c inside a logarithmic derivative boundary condition. 1/r solution.
      include 'meshcom.f'
      integer iuds(3),ifull(3)
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
c                     error(i,j,k)=e
                  errvar=errvar+e**2
                  if(abs(e).gt.abs(errmax))errmax=e
c                  else
c                     error(i,j,k)=0.
               endif
            enddo
         enddo
      enddo
      errvar=errvar/count
      write(*,'(a,f8.5,a,f10.6)')' Maximum vacuum phi error=',errmax,
     $     '   Standard Deviation=',sqrt(errvar)
      if(sqrt(errvar).gt.0.1)write(*,*)
     $     '********** Standard Deviation too large *************'
c Rarely needed printout of u:
c      iform=7
c      uscale=10000000.
c      call udisplay(ndims,u,ifull,iuds,iform,uscale)
c Write some data out for cross checking
c Delete the file first to help with nfs problems.
         open(21,file='smt.round',status='unknown')
         close(21,status='delete')
         open(21,file='smt.round',status='unknown')
         write(21,*)iuds(1),iuds(2),iuds(3)/2
         write(21,'(6e12.4)')
     $        ((u(i,j,iuds(3)/2),i=1,iuds(1)/2),j=1,iuds(2)/2)
         close(21)
c This file has full accuracy.
         open(20,file='smt.out',status='unknown')
         close(20,status='delete')
         open(20,file='smt.out',status='unknown')
         write(20,*)iuds(1),iuds(2),iuds(3)/2
         write(20,*)
     $        ((u(i,j,iuds(3)/2),i=1,iuds(1)/2),j=1,iuds(2)/2)
         close(20)
      end
c***************************************************************
c Packaged version of plotting.
      subroutine solu3plot(ifull,iuds,u,cij,
     $     phi,phiinf,rc,thetain,nth,rs)
      parameter (ndims=3,nd2=ndims*2)
      integer ifull(ndims),iuds(ndims)
      real cij(2*ndims+1,ifull(1),ifull(2),ifull(3))
      real u(ifull(1),ifull(2),ifull(3))
      integer ifmax
      parameter (ifmax=1000,Li=ifmax)
      real cijp(2*ndims+1,ifmax,ifmax)
      real zp(ifmax,ifmax)
      real z(ifmax),xp(ifmax)
c      real uplot(ifmax,ifmax)
c      real zero(ifmax)
      real uanal(ifmax)
      real rfield(ifmax),tfield(ifmax)
      real rprime(ifmax),xprime(ndims,ifmax)
      real xfrac(ndims),xff(ndims),upnd(ndims,ifmax)
      real upsimple(ndims,ifmax)
      integer itemp(ndims)
      real rsimple(ifmax),region(ifmax)
      parameter (mdims=10)
      integer iLs(mdims+1)
      common /iLscom/iLs

      character*40 form1
      include 'meshcom.f'
      include '3dcom.f'

c Start of plotting section.

c      write(*,*)'In solu3plot',ifull,iuds,rc,thetain,nth
      n0=iuds(2)/2
      n1=iuds(3)/2
      idf=3
c      id1=mod(idf,3)+1
c      id2=mod(idf+1,3)+1
c      ifixed=iuds(3)/2
      do i=1,iuds(1)
         xr=xn(i)
         yr=xn(iuds(1)+n0)
         zr=xn(iuds(1)+iuds(2)+n1)
         r=sqrt((xr-.0)**2+(yr-.0)**2+(zr-.0)**2)
         z(i)=phiinf+(phi-phiinf)*rc/r
         xp(i)=xr
c         write(*,*)i,xr
      enddo
      if(.false.)then
      call autoplot(xp,u(1,n0,n1),iuds(1))
      call polymark(xp,u(1,n0,n1),iuds(1),1)
      call color(ibrickred())
      call polyline(xp,z,iuds(1))
      call color(15)
      call axlabels('x mesh (at half mesh in y, z)',
     $     'u(x) and !Af!@(r!dc!d)r!dc!d/r')
      call pltend()
c Plotting slices.
      ifix=3
      call slice3web(ifull,iuds,u,cij,Li,zp,cijp,ixnp,xn,ifix,
     $     'potential:'//'!Ay!@',1)
      endif
c-------------------------------------------------------------------
c Field plotting for various different lines:
c Spherical angles in 3-D
      do iti=1,nth
         theta=thetain*iti
         varphi=0.
         ct=cos(theta)
         st=sin(theta)
         cp=cos(varphi)
         sp=sin(varphi)
         graderr=0
c         write(*,*)'Starting uprime calculation'
         do i=1,Li
c            zero(i)=0.
c            rp=rs*i/(Li-1)
            rp=.99999*rs*i/Li
            rprime(i)=rp
c Coordinates relative to center of first object (sphere).
            xprime(1,i)=rp*st*cp + obj_geom(ocenter,1)  
            xprime(2,i)=rp*st*sp + obj_geom(ocenter+1,1)  
            xprime(3,i)=rp*ct + obj_geom(ocenter+2,1)  
            
            iregion=insidemask(ndims,xprime(1,i))
c Calculate fractional mesh positions of this point, always positive.
c Thus the origin of the box is below point in all dimensions.
            do id=1,ndims
c Offset to start of idf position array.
               ioff=ixnp(id)
c xn is the position array for each dimension arranged linearly.
c Find the index of xprime in the array xn:
               ix=interp(xn(ioff+1),ixnp(id+1)-ioff,xprime(id,i),xm)
c                  ix=xm
               if(xm.ge.(ixnp(id+1)-ixnp(id)))then
                  write(*,*)'****** ccpicplot xm too big',id,xm
                  xm=ixnp(id+1)-ixnp(id)-1.e-3
               endif
               xff(id)=xm
               xfrac(id)=xm-ix
               itemp(id)=ix
c attempt to fix getsimple.
               if(itemp(id).le.1)then
                  itemp(id)=itemp(id)+1
                  xfrac(id)=xfrac(id)-1.
               elseif(itemp(id).ge.ixnp(id+1)-ixnp(id)-1)then
                  itemp(id)=itemp(id)-1
                  xfrac(id)=xfrac(id)+1.
               endif
            enddo
            
c Get the ndims field components at this point.
c               write(*,*)'xfrac',(xfrac(kk),kk=1,ndims)
c     $              ,(xprime(kk,i),kk=1,ndims)
            do idf=1,ndims
               ioff=ixnp(idf)
c New approach mirrors padvnc.
               call getfield(
     $              ndims
     $              ,cij(nd2+1,1,1,1)
     $              ,u,iLs,xn(ixnp(idf)+1),idf
     $              ,xff,iregion,upnd(idf,i))

               call getsimple3field(
     $              ndims,u(itemp(1),itemp(2),itemp(3))
     $              ,iLs,xn(ioff+itemp(idf)),idf
     $              ,xfrac,upsimple(idf,i))

               if(.not.abs(upnd(idf,i)).lt.1.e6)then
                  write(*,*)'Field weird!',upnd(idf,i),i,idf,xff
               endif

            enddo


c Radial component of field
            rfield(i)=
     $           upnd(1,i)*st*cp +
     $           upnd(2,i)*st*sp +
     $           upnd(3,i)*ct
            rsimple(i)=
     $           upsimple(1,i)*st*cp +
     $           upsimple(2,i)*st*sp +
     $           upsimple(3,i)*ct
c Tangential component (magnitude) of field
            tfield(i)=-sqrt(max(0.,
     $           upnd(1,i)**2+upnd(2,i)**2+upnd(3,i)**2
     $           -rfield(i)**2))

c Analytic comparison.
            uanal(i)=(phi-phiinf)*rc/(rprime(i)**2)
c               write(*,'(''i,rprime,rfield,uanal(i)'',i4,4f10.5)')
c     $              i,rprime(i),rfield(i),uanal(i)
c     $              ,rsimple(i)
c Region of point
            region(i)=insideall(ndims,xprime(1,i))
            if(region(i).eq.2)then 
               gerr=abs(rfield(i)-uanal(i))
               if(gerr.gt.graderr)graderr=gerr
            endif
            region(i)=-0.1*region(i)
         enddo
c         write(*,*)'Ended uprime calculation'
c         write(*,'(10f8.4)')rfield
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
      call pltend()
c-------------------------------------------------------------------
c-------------------------------------------------------------------
c Potential plotting for different lines:
c Spherical angles in 3-D
      do iti=1,nth
         theta=thetain*iti
         varphi=0.
         ct=cos(theta)
         st=sin(theta)
         cp=cos(varphi)
         sp=sin(varphi)
         phierr=0
c         write(*,*)'Starting x calculation'
         do i=1,Li
c            zero(i)=0.
c            rp=rs*i/(Li-1)
            rp=0.99999*rs*i/Li
            rprime(i)=rp
c Coordinates relative to center of first object (sphere).
            xprime(1,i)=rp*st*cp + obj_geom(ocenter,1)  
            xprime(2,i)=rp*st*sp + obj_geom(ocenter+1,1)  
            xprime(3,i)=rp*ct + obj_geom(ocenter+2,1)  
            
            iregion=insidemask(ndims,xprime(1,i))
c Calculate fractional mesh positions of this point, always positive.
c Thus the origin of the box is below point in all dimensions.
            do id=1,ndims
c Offset to start of idf position array.
               ioff=ixnp(id)
c xn is the position array for each dimension arranged linearly.
c Find the index of xprime in the array xn:
               ix=interp(xn(ioff+1),ixnp(id+1)-ioff,xprime(id,i),xm)
               if(ix.eq.0.or.ix.ge.ixnp(id+1)-ioff)then
                  write(*,*)'ccpicplot interp outside range'
     $                 ,xprime(id,i)
     $                 ,ioff,(xn(ioff+k),k=1,ixnp(id+1)-ioff),' end '
               endif
c                  ix=xm
               xff(id)=xm
               xfrac(id)=xm-ix
               itemp(id)=ix
            enddo

c               write(*,*)'xfrac',(xfrac(kk),kk=1,ndims)
c     $              ,(xprime(kk,i),kk=1,ndims)            

c Analytic comparison.
            uanal(i)=phiinf+(phi-phiinf)*rc/rprime(i)
            if(uanal(i).lt.phi)uanal(i)=phi
c               write(*,'(''i,rprime,rfield,uanal(i)'',i4,4f10.5)')
c     $              i,rprime(i),rfield(i),uanal(i)
c     $              ,rsimple(i)

            rfield(i)=getpotential(u,cij,iLs,xff,iregion,2)
c Confusing tests.
c            rfcomp=potentialatpoint(x,u,cij,iLs)
c            if(rfcomp.ne.rfield(i))write(*,*)
c     $           'getpotential vs potential at point difference:'
c     $           ,rfield(i),rfcomp
            rsimple(i)=getpotential(u,cij,iLs,xff,iregion,1)
            tfield(i)=getpotential(u,cij,iLs,xff,iregion,3)
c Region of point
            region(i)=insideall(ndims,xprime(1,i))
            if(region(i).eq.2)then 
               perr=abs(rfield(i)-uanal(i))
               if(perr.gt.phierr)phierr=perr
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
      call pltend()
c-------------------------------------------------------------------
      write(*,*)'Max potential error=',phierr
      write(*,*)'Max field error=    ',graderr
      end
c********************************************************************
c Packaged version of plotting.
      subroutine orbit3plot(ifull,iuds,u,phi,rc,rs)
      parameter (ndims=3)
      integer ifull(ndims),iuds(ndims),itemp(ndims)
      real u(ifull(1),ifull(2),ifull(3))
      integer ifmax
      parameter (ifmax=100,Li=ifmax)
      real uplot(ifmax,ifmax)
      character cwork(ifmax,ifmax)

      include 'meshcom.f'
      include 'partcom.f'
c Silence warnings with spurious assigment.
      zclv=phi
c Calculate some stuff for contour plot.
      idf=3
      id1=mod(idf,3)+1
      id2=mod(idf+1,3)+1
      ifixed=iuds(idf)/2
      if(.true.)then
         do i=1,iuds(id1)
            do j=1,iuds(id2)
               itemp(idf)=ifixed
               itemp(id1)=i
               itemp(id2)=j
               uplot(i,j)=u(itemp(1),itemp(2),itemp(3))
            enddo
c               write(*,'(10f8.4)')(uplot(i,j),j=1,iuds(id2))
         enddo
      endif
      call dashset(0)
      nf1=iuds(id1)
      nf2=iuds(id2)
      call pltinit(-rs,rs,-rs,rs)
c Contour without labels, with coloring, using vector coordinates.
      zclv=20.
      icl=0
      call contourl(uplot,cwork,Li,nf1,nf2,zclv,icl,
     $        xn(ixnp(id1)+1),xn(ixnp(id2)+1),17)
      call color(15)
      call axis()
      call color(13)
      call circleplot(0.,0.,rs)
      call circleplot(0.,0.,rc)
      do kk=1,norbits
         call color(kk)
         call polyline(xorbit(1,kk),yorbit(1,kk),iorbitlen(kk))
         call polymark(xorbit(1,kk),yorbit(1,kk),iorbitlen(kk),3)
      enddo
      call pltend()
c      write(*,*)'Returning from orbit3plot.'
      end
