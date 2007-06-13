      program sormpitest
c
      integer Li,ni,nj
      parameter (Li=100,ni=40,nj=40,nk=16)
      integer nd
      parameter (nd=3,nd2=nd*2)
      real u(Li,Li,Li),q(Li,Li,Li),cij(nd2+1,Li,Li,Li)


c object-data storage. Guess at needed size. 
c      integer nobj
c      parameter (nobj=nd2+2)
c      real obj(nobj,2*Li*Li)
c      common objcom/obj/
      include 'objcom.f'
      include 'sormesh.f'
c
      external bdyset,faddu2,cijroutine
      real x(Li),y(Li)
      character*100 form1,argument
      common /sorctl/isor_mi,sor_jac,sor_eps,sor_del,isor_k
      logical lplot
c testing arrays
      integer iuds(nd),ifull(nd),idims(nd),ium2(nd)

      common /myidcom/myid

      lplot=.true.
      call getarg(1,argument)
      if(argument(1:2).eq.'-p') lplot=.false.

      ndims=nd

      dx2=1./ni**2
      dy2=1./nj**2
      dz2=1./nk**2

      idims(1)=2
      idims(2)=2
      ifull(1)=Li
      ifull(2)=Li
      iuds(1)=ni
      iuds(2)=nj
c This circumlocution quiets warnings.
      n3=3
      if(nd.eq.n3)then
         idims(n3)=2
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
            xn(i+iof)=(i-1.)/iuds(id)
         enddo
         iof=iof+iuds(id)
      enddo
c      write(*,'(10f8.4)')xn

c For possibly different jacobi convergence radius parameters.
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
c               q(i,j,k)=1.
c               do n=1,ndims
c                  cij(2*n-1,i,j,k)=iuds(n)**2
c                  cij(2*n,i,j,k)=iuds(n)**2
c               enddo
c No objects for now.
c               cij(2*ndims+1,i,j,k)=0.
            enddo
         enddo
         enddo

      ipoint=0
c      call mditerate(ndims,ifull,iuds,cijroutine,cij,ipoint)
c Remove edges by starting at (2,2,2,...) and using ium2.
      call mditerate(ndims,ifull,ium2,cijroutine,
     $     cij(1,2,2,2),ipoint)

c Example of object data pure fabrication.
c Turned off.
      nnodes=0
      do nn=1,nnodes
         i=ni/2
         j=nj/2
         k=nk/2
         cij(nd2+1,i,j,k)=nn
         do no=1,nobj
            obj(no,nn)=1
         enddo
         obj(1,nn)=.9
         obj(3,nn)=.4
         obj(5,nn)=.5
         obj(7,nn)=.8
         obj(9,nn)=.7
         obj(11,nn)=.8
      enddo


c This jacobi radius is pretty much optimized for Poisson equation with
c fixed boundary, 2-dimensional, when kk=5. 3-d not optimized.
      xyimb=(max(ni,nj)*2.)/float(ni+nj) - 1.
      sor_jac=1.- ((kk/5.)*4./max(10,(ni+nj)/2)**2)
     $     *(1.-0.3*xyimb)
      isor_mi=4*(ni+nj+nk)+10
      sor_eps=1.e-5
c Stop iterating immediately
c         isor_mi=3
         
c Control. Bit 1, use my sor parameters, Bit 2 use faddu.
      ictl=3
c         ictl=1
c         write(*,*)'Calling sormpi, ni,nj=',ni,nj
c The main solver call. Returns process myid.
      call sormpi(ndims,ifull,iuds,cij,u,q,bdyset,faddu2,ictl,ierr
     $     ,myid,idims,nobj,obj)


         if(myid.eq.0)then
            if(ierr.lt.0) write(*,*)'Not converged',ierr
            write(*,*) 'isor_mi,isor_k,sor_jac,sor_del',
     $           isor_mi,isor_k,sor_jac,sor_del,myid
            if(ni+nj.lt.40) then
               write(form1,'(''('',i4,''f10.7)'')')nj
               write(*,*)'u='
               write(*,form1)((u(i,j,n1),j=1,nj),i=1,ni)
            endif
         endif
      enddo

         
      if(myid.eq.0.and.lplot)then
         call cijplot(ndims,ifull,iuds,cij)
      endif

c      write(*,*)'Returned myid=',myid

c      iform=7
c      uscale=10000000.
c      call udisplay(ndims,u,ifull,iuds,iform,uscale)
c no plotting for now.
c         myid=1
c-------------------------------------------------------------------
c Start of plotting section.
      if(myid.eq.0)then
         if(lplot) then
            call yautoplot(u(1,n0,n1),ni)
            call pltend()

c     Web surface plot.
            do i=1,ni
               x(i)=i
            enddo
            do j=1,nj
               y(j)=j
            enddo
 22         continue
 21         call pltinit(0.,1.,0.,1.)
c     Plot the surface. With axes (1). Web color 6, axis color 7.
            jsw=1 + 256*6 + 256*256*7
c     write(*,*)'Starting hidweb'
            call hidweb(x,y,u(1,1,n1),Li,ni,nj,jsw)
            write(form1,'(''Level z='',i4)')n1
            call drwstr(.1,.1,form1)
            if(ieye3d().ne.0) goto 21
            n1=n1+1
            if(n1.lt.nk .and. ndims.gt.2) goto 22
         endif
c-------------------------------------------------------------------
c Write some data out for cross checking
c Delete the file first to help with nfs problems.
         open(21,file='smt.round',status='unknown')
         close(21,status='delete')
         open(21,file='smt.round',status='unknown')
         write(21,*)ni,nj,nk/2
         write(21,'(6e12.4)')((u(i,j,nk/2),i=1,ni/2),j=1,nj/2)
         close(21)
c Delete the file first to help with nfs problems.
         open(20,file='smt.out',status='unknown')
         close(20,status='delete')
         open(20,file='smt.out',status='unknown')
         write(20,*)ni,nj,nk/2
         write(20,*)((u(i,j,nk/2),i=1,ni/2),j=1,nj/2)
         close(20)
      endif
      call MPI_FINALIZE(ierr)
      end
c**********************************************************************
c**********************************************************************
      subroutine bdyset0(ndims,ifull,iuds,cij,u,q)
c     Null version
      end
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
