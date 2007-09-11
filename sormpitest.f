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
      include 'meshcom.f'
c
      external bdyset,faddu2,cijroutine
c      real x(Li),y(Li)
      real z(Li),xp(Li)
      character*100 form1,argument
c      character*10 cxlab,cylab
      common /ctl_sor/mi_sor,xjac_sor,eps_sor,del_sor,k_sor
      logical lplot
c testing arrays
      parameter (ntests=3)
      integer iuds(nd),ifull(nd),idims(nd),ium2(nd)
      real uplot(Li,Li),uprime(Li),xprime(Li),zero(Li),uanal(Li)
      real upregion(Li),tfield(Li,ntests)
      real xnd(nd),xfrac(nd)

      common /myidcom/myid,nprocs

      include '3dcom.f'
c Common data containing the object geometric information. 
c Each object, i < 64 has: type, data(odata).
c      integer ngeomobjmax,odata,ngeomobj
c      parameter (ngeomobjmax=63,odata=16)
c      real obj_geom(odata,ngeomobjmax)
c      common /objgeomcom/ngeomobj,obj_geom

c Structure vector needed for finding adjacent u values.
c Set by mditerate.
      parameter (mdims=10)
      integer iLs(mdims+1)
      common /iLscom/iLs
c Similar set here:
      integer iLsc(mdims+1)

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


      iLsc(1)=nd2+1
      do n=1,ndims
         iLsc(n+1)=iLsc(n)*ifull(n)
      enddo

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

      write(*,*)'iregion:'
      write(form1,'(''('',i2,''i1)'')')iuds(2)
c      write(*,form1)((ireg3(iuds(1)/2,j,k,ifull,cij),j=1,iuds(2)),
c     $           k=1,iuds(3))
      write(*,form1)((ireg3(j,iuds(2)/2,k,ifull,cij),j=1,iuds(1)),
     $           k=1,iuds(3))


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


         if(myid.eq.1)then
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
c         call cijplot(ndims,ifull,iuds,cij,.5)
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
c               write(*,'(10f8.4)')(uplot(i,j),j=1,iuds(id2))
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
c Different lines:
            ium2(id1)=iuds(id1)/2 -3
            xnd(id1)=xn(ixnp(id1)+ium2(id1))
            xnd(id2)=xn(ixnp(id2)+ium2(id2))
            xfrac(id1)=0.
            xfrac(id2)=0.
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
c               icinc=nd2+1
c               iuinc=1
c               do k=1,idf-1
c                  icinc=icinc*ifull(k)
c                  iuinc=iuinc*ifull(k)
c               enddo
c Alternative:
               iuinc=iLs(idf)
               icinc=iLsc(idf)

               ium2(idf)=1
               call gradinterpcorrect(
     $              cij(nd2+1,ium2(1),ium2(2),ium2(3))
     $              ,u(ium2(1),ium2(2),ium2(3))
     $              ,idf,icinc,iuinc,
     $              xprime(i),uprime(i),ix,xm)

c Testing of gradinterpregion
               xnd(idf)=xprime(i)
               iregion=insideall(ndims,xnd)
               call gradinterpregion(
     $              cij(nd2+1,ium2(1),ium2(2),ium2(3))
     $              ,u(ium2(1),ium2(2),ium2(3))
     $              ,idf,icinc,iuinc,
     $              xprime(i),upregion(i),iregion,ix,xm)

c Testing of gradlocalregion
c Offset to start of idf position array.
               ioff=ixnp(idf)
c xn is the position array for each dimension arranged linearly.
c Find the index of xprime in the array xn:
               ix=interp(xn(ioff+1),ixnp(idf+1)-ioff,xprime(i),xm)
               ix=nint(xm)
               xm=xm-ix
               ium2(idf)=ix
               xf=xm
               xfrac(idf)=xf

              call gradlocalregion(
     $              cij(nd2+1,ium2(1),ium2(2),ium2(3))
     $              ,u(ium2(1),ium2(2),ium2(3))
     $              ,idf,icinc,iuinc,xn(ioff+ix)
     $              ,xf,uplocal,iregion,ixout,xm)
 
               call getfield(
     $              ndims
     $              ,cij(nd2+1,ium2(1),ium2(2),ium2(3))
     $              ,u(ium2(1),ium2(2),ium2(3))
c     $              ,iLsc
     $              ,iLs
     $              ,xn(ioff+ix),idf
     $              ,xfrac,iregion,field)

               write(*,*)iregion,ix,xm,xprime(i),uprime(i),uplocal
     $              ,field
c     $              ,upregion(i)

               do itest=1,ntests
                  xfrac(id1)=(itest-1.)/(ntests)
                  call getfield(
     $                 ndims
     $                 ,cij(nd2+1,ium2(1),ium2(2),ium2(3))
     $                 ,u(ium2(1),ium2(2),ium2(3))
c     $                 ,iLsc
     $                 ,iLs
     $                 ,xn(ioff+ix),idf
     $                 ,xfrac,iregion,tfield(i,itest))
               enddo
               write(*,*)uplocal,(tfield(i,itest),itest=1,ntests)

c Analytic comparison.
               uanal(i)=-phi*rc/((xprime(i)-0.5)*abs(xprime(i)-0.5))
c               write(*,*)i,idf,ioff,ix,xm,uprime(i)
            enddo
            write(*,*)'Ended uprime calculation'
            call autoplot(xprime,uprime,Li)
c            call polymark(xprime,uprime,Li,1)
            call polymark(xn(ioff+1),zero,iuds(idf),10)
            call winset(.true.)
            call dashset(1)
            call polyline(xprime,uanal,Li)
            call dashset(2)
            call color(ired())
            call polyline(xprime,upregion,Li)
            call color(iblue())
            do itest=1,ntests
               call polyline(xprime,tfield(1,itest),Li)
            enddo
            call color(15)
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
      real function faddu2(u,fprime)
      faddu2=0.
      faddu=0.
      fprime=0.
      end

