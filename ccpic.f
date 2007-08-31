      program ccpic
c
      include 'objcom.f'
c Storage array spatial count size
      integer Li,ni,nj,nk
c      parameter (Li=100,ni=40,nj=40,nk=16)
      parameter (Li=100,ni=16,nj=16,nk=20)
c      parameter (Li=100,ni=32,nj=32,nk=40)
      real u(Li,Li,Li),q(Li,Li,Li),cij(2*ndims_sor+1,Li,Li,Li)
      real psum(Li,Li,Li)
c Used dimensions, Full dimensions. Used dims-2
      integer iuds(ndims_sor),ifull(ndims_sor),ium2(ndims_sor)
c Mesh spacing description structure
      include 'meshcom.f'
c Processor cartesian geometry
      integer nblksi,nblksj,nblksk
      parameter (nblksi=1,nblksj=1,nblksk=1)
      integer idims(ndims_sor)
c MPI information.
      common /myidcom/myid,nprocs
c Common data containing the BC-object geometric information
      include '3dcom.f'
c Structure vector needed for finding adjacent u values.
c Set by mditerate.
      parameter (mdims=10)
      integer iLs(mdims+1)
      common /iLscom/iLs
c Particle common data
      include 'partcom.f'
c Plasma common data
      include 'plascom.f'
      integer ndims,nd2
      parameter (ndims=ndims_sor,nd2=ndims*2)

c Diagnostics
c      real error(Li,Li,Li)
      real zp(Li,Li),cijp(2*ndims_sor+1,Li,Li)
c
      external bdyset,faddu,cijroutine,psumtoq
      real z(Li),xp(Li)
      character*100 form1,argument
      common /ctl_sor/mi_sor,xjac_sor,eps_sor,del_sor,k_sor
      logical lplot,l1plot

c testing arrays
      parameter (ntests=1)
      real uplot(Li,Li),zero(Li),uanal(Li)
      real rfield(Li,ntests),tfield(Li,ntests)
      real rprime(Li),xprime(ndims,Li)
      real xfrac(ndims),upnd(ndims,Li),upsimple(ndims,Li),itemp(ndims)
c      real xcenter(ndims),upregion(Li),uprime(Li),xnd(ndims)
      integer id1,id2,idf
      real rsimple(Li,ntests)

c Initialize ids to silence spurious warnings.
      data id1,id2,idf/0,0,0/
c-------------------------------------------------------------

c Geometry information read in.
      call readgeom('ccpicgeom.dat')
c First object is sphere of radius rc and potential phi.
      rc=obj_geom(5,1)
      phi=-obj_geom(10,1)/obj_geom(8,1)
      write(*,*)'rc=',rc,'  phi=',phi
c Second object is bounding sphere of radius rs. 
c But use a tad more for the mesh size
      rs=obj_geom(5,2)*1.00001

c Default values, may be overridden by arguments.
c Set mesh and mpi parameters.
      idims(1)=nblksi
      idims(2)=nblksj
      idims(3)=nblksk
      ifull(1)=Li
      ifull(2)=Li
      ifull(3)=Li
      iuds(1)=ni
      iuds(2)=nj
      iuds(3)=nk

      thetain=.1
      nth=1
      lplot=.true.
      l1plot=.false.

c Places where plotting occurs.
      n0=nj/2
      n1=nk/2
      if(ndims.lt.3)n1=1
c--------------------------------------------------------------
c Deal with arguments
c      if(iargc().eq.0) goto "help"
      do i=1,iargc()
         call getarg(i,argument)
         if(argument(1:3).eq.'-p ') lplot=.false.
         if(argument(1:3).eq.'-p1') l1plot=.true.
         if(argument(1:2).eq.'-t')read(argument(3:),*)thetain
         if(argument(1:2).eq.'-n')read(argument(3:),*)nth
           
      enddo
c---------------------------------------------------------------
c Construct the mesh vector(s) and ium2
      iof=0
      do id=1,ndims
c Pointer to start of vector.
         ixnp(id)=iof
c Mesh data. Make it extend from -rs to +rs
         do i=1,iuds(id)
            xn(i+iof)=rs*(2*(i-1.)/(iuds(id)-1.) - 1.)
         enddo
         iof=iof+iuds(id)
c Mesh size with the faces removed:
         ium2(id)=iuds(id)-2
      enddo
      ixnp(ndims+1)=iof

c----------------------------------------------------------------
c INITIALIZATIONS
c Initialize cij:
      ipoint=0
c Remove edges by starting at (2,2,2,...) and using ium2.
      call mditerate(ndims,ifull,ium2,cijroutine,
     $     cij(1,2,2,2),ipoint)
c Initialize the region flags in the object data
      call iregioninit(ndims,ifull)
c Text graphic of slice through cij
c      write(*,*)'iregion:'
      write(form1,'(''('',i2,''i1)'')')iuds(1)
      write(*,form1)((ireg3(j,iuds(2)/2,k,ifull,cij),j=1,iuds(1)),
     $           k=1,iuds(3))
c The following requires include objcom.f
      write(*,*)'Finished mesh setup. Used No of pointers:',oi_sor
c      write(*,*)'Finished mesh setup.'
c      write(*,'(a,8f8.1)')'cij(*,3,3,3)=',(cij(i,3,3,3),i=1,nd2+1)
      if(lplot)call cijplot(ndims,ifull,iuds,cij,rs)

c Initialize charge
      call zero3array(q,iLs,ni,nj,nk)

c Set the SOR parameters.
c For possibly different jacobi convergence radius parameters.
c 7 is a pretty good value in 3-d.
      kk=7
c This jacobi radius is pretty much optimized for Poisson equation with
c fixed boundary, 2-dimensional, when kk=5. 3-d not optimized.
      xyimb=(max(ni,nj)*2.)/float(ni+nj) - 1.
      xjac_sor=1.- ((kk/5.)*4./max(10,(ni+nj)/2)**2)
     $     *(1.-0.3*xyimb)
      mi_sor=4*(ni+nj+nk)+20
      eps_sor=1.e-5

c Demonstration of initialization call that returns myid if I need it.
c Control bit 4 pure initialization, no communication. But don't let
c it tell us to use internal sor params (1).
         ictl=5
c Returns process myid.
      call sormpi(ndims,ifull,iuds,cij,u,q,bdyset,faddu,ictl,ierr
     $     ,myid,idims)
c      write(*,*)'sormpi Initialization returns',myid
c End of initializations
c---------------------------------------------------------------         
c Control. Bit 1, use my sor parameters, Bit 2 use faddu.
c      ictl=3
      ictl=1
c         write(*,*)'Calling sormpi, ni,nj=',ni,nj
c The main solver call.
      call sormpi(ndims,ifull,iuds,cij,u,q,bdyset,faddu,ictl,ierr
     $     ,myid,idims)

      if(myid.eq.0)then
c------------------------------------------------------------------
         if(ierr.lt.0)then
            write(*,*)'Not converged.'
     $        ,ierr
            write(*,*) 'mi_sor,k_sor,xjac_sor,del_sor',
     $           mi_sor,k_sor,xjac_sor,del_sor,myid
c         if(ni+nj.lt.40) then
c            write(form1,'(''('',i4,''f8.5)'')')nj
c            write(*,*)'u='
c            write(*,form1)((u(i,j,n1),j=1,nj),i=1,ni)
c         endif
         endif
c-------------------------------------------------------------------
c-------------------------------------------------------------------
c Do some analytic checking of the case with a fixed potential sphere
c inside a logarithmic derivative boundary condition. 1/r solution.
c Also write out some data for checking.
      call spherecheck(ifull,iuds,u,phi,rc)
      call solu3plot(ifull,iuds,u,cij,phi,rc,thetain,nth,iregion
     $     ,rs)
c Replaced with packaged call.
      if(.false.)then
c-------------------------------------------------------------------
c Start of plotting section.
c Calculate some stuff for autocolorcontour.
       idf=3
       id1=mod(idf,3)+1
       id2=mod(idf+1,3)+1
       ifixed=nk/2
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
c         if(lplot .and. abs(errmax).lt..1) then
         if(lplot) then
            call yautoplot(u(1,n0,n1),ni)
            do i=1,ni
               xr=xn(i)
               yr=xn(iuds(1)+n0)
               zr=xn(iuds(1)+iuds(2)+n1)
               r=sqrt((xr-.0)**2+(yr-.0)**2+(zr-.0)**2)
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
         endif
         if(l1plot)then
c Different lines:
c Spherical angles in 3-D
            do iti=1,nth
            itest=1
            theta=thetain*iti
            varphi=0.
            ct=cos(theta)
            st=sin(theta)
            cp=cos(varphi)
            sp=sin(varphi)
            write(*,*)'Starting uprime calculation'
            do i=1,Li
               zero(i)=0.
               rp=rs*(i)/(Li-1)
               rprime(i)=rp
c Coordinates relative to center of first object (sphere).
               xprime(1,i)=rp*st*cp + obj_geom(ocenter,1)  
               xprime(2,i)=rp*st*sp + obj_geom(ocenter+1,1)  
               xprime(3,i)=rp*ct + obj_geom(ocenter+2,1)  
               
               iregion=insideall(ndims,xprime(1,i))
c Calculate fractional mesh positions of this point, always positive.
c Thus the origin of the box is below point in all dimensions.
               do id=1,ndims
c Offset to start of idf position array.
                  ioff=ixnp(id)
c xn is the position array for each dimension arranged linearly.
c Find the index of xprime in the array xn:
                  ix=interp(xn(ioff+1),ixnp(id+1)-ioff,xprime(id,i),xm)
c                  ix=xm
                  xfrac(id)=xm-ix
                  itemp(id)=ix
               enddo

c Get the ndims field components at this point.
c               write(*,*)'xfrac',(xfrac(kk),kk=1,ndims)
c     $              ,(xprime(kk,i),kk=1,ndims)
               do idf=1,ndims
                  ioff=ixnp(idf)
                  call getfield(
     $                 ndims
     $                 ,cij(nd2+1,itemp(1),itemp(2),itemp(3))
     $                 ,u(itemp(1),itemp(2),itemp(3))
     $                 ,iLs
     $                 ,xn(ioff+itemp(idf)),idf
     $                 ,xfrac,iregion,upnd(idf,i))
                  call getsimple3field(
     $                 ndims,u(itemp(1),itemp(2),itemp(3))
     $                 ,iLs,xn(ioff+itemp(idf)),idf
     $                 ,xfrac,upsimple(idf,i))
               enddo

c Radial component of field
               rfield(i,itest)=
     $              upnd(1,i)*st*cp +
     $              upnd(2,i)*st*sp +
     $              upnd(3,i)*ct
               rsimple(i,itest)=
     $              upsimple(1,i)*st*cp +
     $              upsimple(2,i)*st*sp +
     $              upsimple(3,i)*ct
c Tangential component (magnitude) of field
               tfield(i,itest)=-sqrt(max(0.,
     $              upnd(1,i)**2+upnd(2,i)**2+upnd(3,i)**2
     $              -rfield(i,itest)**2))

c Analytic comparison.
               uanal(i)=phi*rc/(rprime(i)**2)
c               write(*,'(''i,rprime,rfield,uanal(i)'',i4,4f10.5)')
c     $              i,rprime(i),rfield(i,itest),uanal(i)
c     $              ,rsimple(i,itest)

            enddo
            write(*,*)'Ended uprime calculation'
            call dashset(0)
            if(iti.eq.1)then
               call autoplot(rprime,rfield(1,itest),Li)
               call winset(.true.)
               call dashset(1)
               call color(ired())
               call polyline(rprime,uanal,Li)
               call winset(.false.)
c            call polyline(xprime,upregion,Li)
               call color(iblue())
               call dashset(3)
               call polyline(rprime,rsimple(1,itest),Li)
            else
               call polyline(rprime,rfield(1,itest),Li)
            endif
            call dashset(2)
            call color(idarkgreen())
            call polyline(rprime,tfield(1,itest),Li)
c            do itest=2,ntests
c               call polyline(rprime,rfield(1,itest),Li)
c               call dashset(2)
c            enddo
            call color(15)
            call winset(.false.)
            form1='!Aq!@='
            call fwrite(180*theta/3.1415926,iwdth,1,form1(7:))
            call jdrwstr(.01,.1,form1,1.)
            enddo
            call dashset(0)
            call pltend()
c-------------------------------------------------------------------
         endif
c End of plotting.
      endif
      endif
c------------------------------------------------------------------
c Orbit testing:
      n_part=1
      if_part(1)=1
c We are going to populate this region; so identify it.
      x_part(1,1)=.9*rs
      x_part(2,1)=0.
      x_part(3,1)=0.
      iregion_part=insideall(ndims,x_part(1,1))
      call srand(myid)
      write(*,*)'Calling pinit'
      call pinit()
      write(*,*)'Return from pinit'
      norbits=1
      dt=.025 
      x_part(1,1)=2.
      x_part(2,1)=0.
c      x_part(1,1)=0.
c      x_part(2,1)=-.4
      x_part(3,1)=0.
      x_part(4,1)=2*dt
      x_part(5,1)=1.1
      x_part(6,1)=0.
c Already set.
c      do id=1,ndims
c         ium2(id)=iuds(id)-2
c      enddo
      rhoinf=1.     !for now.
      do j=1,30
         write(*,'(i4,i4''  x_p='',6f10.5)')j,ierr, (x_part(k,1),k=1,6)
         call zero3array(psum,iLs,ni,nj,nk)
         call chargetomesh(psum,iLs,diags)
c Convert psums to charge, q. Remember external psumtoq!
         call mditerarg(ndims,ifull,ium2,psumtoq,
     $        0,psum(2,2,2),q(2,2,2),rhoinf)
c Some diagnostics.
c         write(*,*)'Psum:'
c         call diag3array(psum,iLs,ni,nj,nk)
c         write(*,*)'q:'
c         call diag3array(q,iLs,ni,nj,nk)
         call sormpi(ndims,ifull,iuds,cij,u,q,bdyset,faddu,ictl,ierr
     $     ,myid,idims)
c         write(*,*)'Sormpi iterations:',ierr
         if(lplot)
     $        call slice3web(ifull,iuds,u,cij,Li,zp,cijp,ixnp,xn,ifix,
     $           'potential:'//'!Ay!@')
c
         call padvnc(ndims,cij,u,iLs)
      enddo
c      write(*,*)iorbitlen(1),(xorbit(k,1),k=1,10)

      call dashset(0)
      call autocolcont(uplot,Li,iuds(id1),iuds(id2))
      call scalewn(-rs,rs,-rs,rs,.false.,.false.)
      call ticset(-.01,-.01,-.03,-.02,0,0,0,0)
      call color(15)
      call axis()
      call ticset(.0 ,.0 ,.0,.0,0,0,0,0)
      call color(12)
      call circleplot(0.,0.,rs)
      call circleplot(0.,0.,rc)
      call color(5)
      call polyline(xorbit,yorbit,iorbitlen(1))
c      call autoplot(xorbit,yorbit,iorbitlen(1))
      call polymark(xorbit,yorbit,iorbitlen(1),1)
      call pltend()

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
c     L(u) + f(u) = q(x,y,...), 
c     where L is a second order elliptical differential operator 
c     represented by a difference stencil of specified coefficients,
c     f is some additional function, and q is the "charge density".
c f is exp(u) here for Boltzmann electrons and densities normalized
c to unity at infinity.
      real function faddu(u,fprime)
      real u,fprime
      fprime=exp(u)
      faddu=fprime
      end
c************************************************************************
      subroutine bdy3slope(inc,ipoint,indi,ndims,iused,u)
c Version of bdyroutine that sets derivative=0 on 3-boundary.
      integer ipoint,inc
      integer indi(ndims),iused(ndims)
      real u(*)

      parameter (mdims=10)
c Structure vector needed for finding adjacent u values.
c Can't be passed here because of mditerate argument conventions.
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
c*********************************************************************
      subroutine zero3array(array,iLs,ni,nj,nk)
      real array(*)
      integer iLs(4)
      do k=1,nk
         do j=1,nj
            do i=1,ni
               index=(i+iLs(2)*(j-1)+iLs(3)*(k-1))
               array(index)=0.
            enddo
         enddo
      enddo
      end
c*********************************************************************      
      subroutine diag3array(array,iLs,ni,nj,nk)
      real array(*)
      integer iLs(4)
      do k=1,nk
         do j=1,nj
            do i=1,ni
               index=(i+iLs(2)*(j-1)+iLs(3)*(k-1))
               x=array(index)
               if(x.ne.0)
     $              write(*,'(a,4i7,f14.5)')
     $              'diag3array at  ',index,i,j,k,x
            enddo
         enddo
      enddo
      end

c********************************************************************
c Default plasma data.
      block data plascomset
      include 'plascom.f'
      data debylen,Ti,vd,rs/1.,1.,0.,5./
c      data debylen,Ti,vd,rs/1.,1.,0.,.45/

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
      subroutine spherecheck(ifull,iuds,u,phi,rc)
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
                  e=u(i,j,k)-phi*rc/r
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
      write(*,*)'Max phi error=',errmax,
     $     ' Standard Deviation=',sqrt(errvar)
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
      subroutine solu3plot(ifull,iuds,u,cij,phi,rc,thetain,nth,iregion
     $     ,rs)
      parameter (ndims=3,nd2=ndims*2)
      integer ifull(ndims),iuds(ndims)
      real cij(2*ndims+1,ifull(1),ifull(2),ifull(3))
      real u(ifull(1),ifull(2),ifull(3))
      integer ifmax
      parameter (ifmax=100,Li=ifmax)
      real cijp(2*ndims+1,ifmax,ifmax)
      real zp(ifmax,ifmax)
      real z(ifmax),xp(ifmax)
c      real uplot(ifmax,ifmax)
      real zero(ifmax),uanal(ifmax)
      real rfield(ifmax),tfield(ifmax)
      real rprime(ifmax),xprime(ndims,ifmax)
      real xfrac(ndims),upnd(ndims,ifmax)
      real upsimple(ndims,ifmax),itemp(ndims)
      real rsimple(ifmax)
      parameter (mdims=10)
      integer iLs(mdims+1)
      common /iLscom/iLs

      character*40 form1
      include 'meshcom.f'
      include '3dcom.f'

c Start of plotting section.
c Calculate some stuff for autocolorcontour.
c Places where plotting occurs.

      write(*,*)'In solu3plot',ifull,iuds,rc,thetain,nth
      n0=iuds(2)/2
      n1=iuds(3)/2
      idf=3
      id1=mod(idf,3)+1
      id2=mod(idf+1,3)+1
      ifixed=iuds(3)/2
      call yautoplot(u(1,n0,n1),iuds(1))
      do i=1,iuds(1)
         xr=xn(i)
         yr=xn(iuds(1)+n0)
         zr=xn(iuds(1)+iuds(2)+n1)
         r=sqrt((xr-.0)**2+(yr-.0)**2+(zr-.0)**2)
         z(i)=phi*rc/r
         xp(i)=i
      enddo
      call polyline(xp,z,iuds(1))
      call pltend()
c Plotting slices.
      ifix=3
      call slice3web(ifull,iuds,u,cij,Li,zp,cijp,ixnp,xn,ifix,
     $     'potential:'//'!Ay!@')
c-------------------------------------------------------------------
c Different lines:
c Spherical angles in 3-D
      do iti=1,nth
         theta=thetain*iti
         varphi=0.
         ct=cos(theta)
         st=sin(theta)
         cp=cos(varphi)
         sp=sin(varphi)
         write(*,*)'Starting uprime calculation'
         do i=1,Li
            zero(i)=0.
            rp=rs*(i)/(Li-1)
            rprime(i)=rp
c Coordinates relative to center of first object (sphere).
            xprime(1,i)=rp*st*cp + obj_geom(ocenter,1)  
            xprime(2,i)=rp*st*sp + obj_geom(ocenter+1,1)  
            xprime(3,i)=rp*ct + obj_geom(ocenter+2,1)  
            
            iregion=insideall(ndims,xprime(1,i))
c Calculate fractional mesh positions of this point, always positive.
c Thus the origin of the box is below point in all dimensions.
            do id=1,ndims
c Offset to start of idf position array.
               ioff=ixnp(id)
c xn is the position array for each dimension arranged linearly.
c Find the index of xprime in the array xn:
               ix=interp(xn(ioff+1),ixnp(id+1)-ioff,xprime(id,i),xm)
c                  ix=xm
               xfrac(id)=xm-ix
               itemp(id)=ix
            enddo
            
c Get the ndims field components at this point.
c               write(*,*)'xfrac',(xfrac(kk),kk=1,ndims)
c     $              ,(xprime(kk,i),kk=1,ndims)
            do idf=1,ndims
               ioff=ixnp(idf)
               call getfield(
     $              ndims
     $              ,cij(nd2+1,itemp(1),itemp(2),itemp(3))
     $              ,u(itemp(1),itemp(2),itemp(3))
     $              ,iLs
     $              ,xn(ioff+itemp(idf)),idf
     $              ,xfrac,iregion,upnd(idf,i))
               call getsimple3field(
     $              ndims,u(itemp(1),itemp(2),itemp(3))
     $              ,iLs,xn(ioff+itemp(idf)),idf
     $              ,xfrac,upsimple(idf,i))
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
            uanal(i)=phi*rc/(rprime(i)**2)
c               write(*,'(''i,rprime,rfield,uanal(i)'',i4,4f10.5)')
c     $              i,rprime(i),rfield(i),uanal(i)
c     $              ,rsimple(i)

         enddo
         write(*,*)'Ended uprime calculation'
         call dashset(0)
         if(iti.eq.1)then
            call autoplot(rprime,rfield(1),Li)
            call winset(.true.)
            call dashset(1)
            call color(ired())
            call polyline(rprime,uanal,Li)
            call winset(.false.)
c            call polyline(xprime,upregion,Li)
            call color(iblue())
            call dashset(3)
            call polyline(rprime,rsimple(1),Li)
         else
            call polyline(rprime,rfield(1),Li)
         endif
         call dashset(2)
         call color(idarkgreen())
         call polyline(rprime,tfield(1),Li)
         call color(15)
         call winset(.false.)
         form1='!Aq!@='
         call fwrite(180*theta/3.1415926,iwdth,1,form1(7:))
         call jdrwstr(.01,.1,form1,1.)
      enddo
      call dashset(0)
      call pltend()
c-------------------------------------------------------------------
      end
