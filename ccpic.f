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
c Don't use the mditerate common. It might not be right.
      integer iLs(ndims_sor+1)
c Particle common data
      include 'partcom.f'
c Plasma common data
      include 'plascom.f'
      integer ndims,nd2
      parameter (ndims=ndims_sor,nd2=ndims*2)

      external bdyset,faddu,cijroutine,cijedge,psumtoq
      character*100 form1,argument
      common /ctl_sor/mi_sor,xjac_sor,eps_sor,del_sor,k_sor
      logical lplot,lcijplot,lsliceplot,lorbitplot

c Diagnostics
c      real error(Li,Li,Li)
      real zp(Li,Li),cijp(2*ndims_sor+1,Li,Li)
c testing arrays
      real xir(ndims)
c Set up the structure vector.
      data iLs/1,Li,(Li*Li),(Li*Li*Li)/
c Mesh and mpi parameter defaults:
      data idims/nblksi,nblksj,nblksk/
      data ifull/Li,Li,Li/
      data iuds/ni,nj,nk/
c Point which lies in the plasma region:
      data xir/2.,2.,2./
c Data for plotting etc.
      data lplot,lcijplot,lsliceplot,lorbitplot/
     $     .false.,.false.,.false.,.false./
      data thetain,nth/.1,1/

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

c--------------------------------------------------------------
c Deal with arguments
c      if(iargc().eq.0) goto "help"
      do i=1,iargc()
         call getarg(i,argument)
         if(argument(1:3).eq.'-p ') lplot=.true.
         if(argument(1:3).eq.'-pc') lcijplot=.true.
         if(argument(1:3).eq.'-ps') lsliceplot=.true.
         if(argument(1:3).eq.'-po') lorbitplot=.true.
         if(argument(1:2).eq.'-t')read(argument(3:),*)thetain
         if(argument(1:2).eq.'-n')read(argument(3:),*)nth
           
      enddo
c---------------------------------------------------------------
c Construct the mesh vector(s) and ium2
      call meshconstruct(ndims,iuds,ium2,rs)
c----------------------------------------------------------------
c INITIALIZATIONS
      write(*,*)'Initializing the stencil data cij'
c Initialize cij:
c New scheme. Pass whole cij but give offset of one in each dimension.
      ipoint=iLs(1)+iLs(2)+iLs(3)
c      write(*,*)ipoint,iLs
      call mditerate(ndims,ifull,ium2,cijroutine,
     $     cij(1,1,1,1),ipoint)

c Set an object pointer for all the edges so their regions get
c set by the following call
      ipoint=0
      call mditerate(ndims,ifull,iuds,cijedge,cij,ipoint)
      
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
      if(lcijplot)call cijplot(ndims,ifull,iuds,cij,rs)

c Initialize charge
      call zero3array(q,iLs,ni,nj,nk)
c---------------------------------------------------------------         
c Control. Bit 1, use my sor parameters (not here). Bit 2 use faddu.
c Control bit 4 pure initialization, no communication (to get myid).
c      ictl=2
      ictl=0
c         write(*,*)'Calling sormpi, ni,nj=',ni,nj
c An initial solver call.
      call sormpi(ndims,ifull,iuds,cij,u,q,bdyset,faddu,ictl,ierr
     $     ,myid,idims)
      write(*,*)ierr

      if(myid.eq.0)then
c------------------------------------------------------------------
c         if(ierr.lt.0)then
c            write(*,*)'Not converged.',ierr
c            write(*,*) 'mi_sor,k_sor,xjac_sor,del_sor',
c     $           mi_sor,k_sor,xjac_sor,del_sor,myid
c         endif
c-------------------------------------------------------------------
c-------------------------------------------------------------------
c Do some analytic checking of the case with a fixed potential sphere
c inside a logarithmic derivative boundary condition. 1/r solution.
c Also write out some data for checking.
         if(lplot)then
            call spherecheck(ifull,iuds,u,phi,rc)
c Plot some of the initial-solver data.
            call solu3plot(ifull,iuds,u,cij,phi,rc,thetain,nth
     $           ,rs)
         endif
c End of plotting.
      endif
c------------------------------------------------------------------
c We are going to populate the region in which xir lies.
      iregion_part=insideall(ndims,xir)
c With a specified number of particles.
      n_part=2
      call srand(myid)
      write(*,*)'Initializing',n_part,' particles'
      call pinit()
      write(*,*)'Return from pinit'
c A special orbit.
c Pinit resets x_part. So set it for the special first particle.
      x_part(1,1)=2.
      x_part(2,1)=0.
      x_part(3,1)=0.
      dt=.1
c Prior half-step radial velocity
      x_part(4,1)=0.5*dt*(abs(phi)/x_part(1,1)**2)
c Tangential velocity of circular orbit at r=4.
      x_part(5,1)=sqrt(abs(phi/x_part(1,1))-x_part(4,1)**2)
      x_part(6,1)=0.
      nsteps=10.
      write(*,*)'dt=',dt,' dtheta=',dt*x_part(5,1)/x_part(1,1),
     $     ' steps=',nsteps,' total theta/2pi='
     $     ,nsteps*dt*x_part(5,1)/x_part(1,1)/2./3.1415927

c Follow the first norbits particles.
      norbits=5

c Turn off self field.
      rhoinf=1.e20     !for now.
      do j=1,nsteps
         write(*,'(i4,i4''  x_p='',7f9.5)')j,ierr, (x_part(k,1),k=1,6)
     $        ,sqrt(x_part(1,1)**2+x_part(2,1)**2)
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
         if(lsliceplot)
     $        call slice3web(ifull,iuds,u,cij,Li,zp,cijp,ixnp,xn,ifix,
     $           'potential:'//'!Ay!@')
c
         call padvnc(ndims,cij,u,iLs)
      enddo
c      write(*,*)iorbitlen(1),(xorbit(k,1),k=1,10)
      
      if(lorbitplot)call orbit3plot(ifull,iuds,u,phi,rc,rs)
c-------------------------------------------------------------------
      call MPI_FINALIZE(ierr)
      end
c**********************************************************************
c********************************************************************
c Default plasma data.
      block data plascomset
      include 'plascom.f'
      data debylen,Ti,vd,rs/1.,1.,0.,5./
c      data debylen,Ti,vd,rs/1.,1.,0.,.45/

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
c*******************************************************************
      subroutine meshconstruct(ndims,iuds,ium2,rs)
      integer iuds(ndims),ium2(ndims)
      include 'meshcom.f'

      iof=0
      do id=1,ndims
c Pointer to start of vector.
         ixnp(id)=iof
c Mesh data. Make it extend from -rs to +rs
         do i=1,iuds(id)
            xn(i+iof)=rs*(2*(i-1.)/(iuds(id)-1.) - 1.)
            if(.false.)then
c Two region Non-uniform assuming iuds even.
            iu2=iuds(id)/2
c Uniform -1 to 1:
            xi=(i-iu2-0.5)/(iu2-0.5)
c Half of mesh at each spacing
            xmid=0.5
c First half of mesh extends to this fraction of the total rs.
            xs=0.25
            if(abs(xi).lt.xmid)then
               xn(i+iof)=rs*xs*xi/xmid
            else
               xn(i+iof)=rs*sign(((abs(xi)-xmid)+xs*(1.-abs(xi))),xi)
     $              /(1-xmid)
            endif
            endif
         enddo
         iof=iof+iuds(id)
c Mesh size with the faces removed:
         ium2(id)=iuds(id)-2
      enddo
      ixnp(ndims+1)=iof

      end
