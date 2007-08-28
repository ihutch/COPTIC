c Solve an elliptical problem in ndims dimensions
c by simultaneous overrelaxation
c     L(u) + f(u) = q(x,y,...), 
c     where L is a second order elliptical differential operator 
c     represented by a difference stencil of specified coefficients,
c     f is some additional function, and q is the "charge density".
c Ian Hutchinson, January 2006, Dec 2006
c
      subroutine sormpi(ndims,ifull,iuds,cij,u,q,bdyset,faddu,ictl,ierr
     $     ,mpiid,idims)

c The used number of dimensions. But this must be equal to ndims 
c parameter used for bbdydecl dimensions.
      integer ndims
c     The number of dimensions declared in bbdydecl.f
      parameter (ndimsdecl=3)
c     ifull full dimensions, iuds used dimensions
      integer ifull(ndims)
c     iuds(ndims) declared in bbdydecls
c     cij coefficients of the finite difference scheme, in the order
c         east,west,north,south ... [regarding i,j, as x,y].
c         real cij(nd2+1,Li,nj)
c     The last cij value of the first index is a pointer to object info
c     and if it is zero (null) then no object code affects this point.
c     Objects cut the mesh between nodes and e.g. specify u there.
c     They can thus define objects embedded in the mesh.
c
c     In this mpi routine we must use linear addressing for cij,u,q
c     so that the pointers can be used for each block.
      real cij(*)
c     u  potential to be solved for (initialized on entry).
c        its boundaries are at 1,ni; 1,nj; 
c        which are normally set on entry or by bdyset,
c        but but boundary values are never changed inside sormpi.
      real u(*)
c     q  "charge density" input
      real q(*)
c User supplied functions which should be declared external in the 
c calling routine.
c     bdyset subroutine that evaluates the boundary conditions and
c               deposits into the edge values of the potential.
c               bdyset(Li,ni,nj,cij,u,q)
c               If non-constant coefficients are required, then
c               bdyset could be used to adjust their values.
c     faddu(u,fprime)  
c               real function that returns the additional component
c               f, and as parameter fprime=f'.
      external faddu
c     ictl  integer control switches
c           bit 1: user-defined iteration parameters (default no)
c           bit 2: use faddu (default no)
c           bit 3: just initialize bbdy (return mpiid).
c     ierr  Returns Positive: number of iterations to convergence.
c           Negative: Not-converged maximum iterations.
c           Zero: something bad.
c     mpiid Returns my MPI process number for this mpi version.
      integer ictl,ierr
c     integer idims(*)  number of blocks, declared in bbdydecl.f


c Other things that we might want control over include the maximum
c number of iterations, the Jacobi radius, and the convergence size.
c k_sor is the sor iteration index, for diagnostics.
      common /ctl_sor/mi_sor,xjac_sor,eps_sor,del_sor,k_sor

      logical lconverged
      logical laddu
c      logical lmpisplit

c Declarations fo MPI calls
      include 'mpif.h'

c bbdydecl declares most things for bbdy, using parameter ndimsdecl.
      include 'bbdydecl.f'

      real delta,umin,umax
c      data lmpisplit/.false./
c This saves data between calls so we can use separate initialization
c But at present that causes seg faults!
      save

c-------------------------------------------------------------------
      if(ndims.ne.ndimsdecl)then
         write(*,*)'Wrong number of dimensions in sormpi call',
     $        ndims,ndimsdecl
         stop
      endif
c      write(*,*)'In sormpi',ictl,ndims,ifull,iuds,idims,ierr
c      return
      do nd=1,ndims
         lperiod(nd)=.false.
      enddo
c--------------------------------------------------------------------
c Control Functions:
      ictlh=ictl
c First bit of ictlh indicates if sorctl preset or defaults.
      if(mod(ictlh,2).eq.0)then
         xyimb=(max(iuds(1),iuds(2))*2.)/float(iuds(1)+iuds(2)) - 1.
         xjac_sor=1.- (4./max(10,(iuds(1)+iuds(2))/2)**2)
     $        *(1.-0.3*xyimb)
         mi_sor=2.*(iuds(1)+iuds(2))+10
         eps_sor=1.e-5
      endif
c Second bit of ictlh indicates if there's additional term.
      ictlh=ictlh/2
      if(mod(ictlh,2).ne.0)then
         laddu=.true.
      else
         laddu=.false.
      endif
c Third bit of ictlh indicates we are just initializing.
      ictlh=ictlh/2
      if(mod(ictlh,2).ne.0)then
c (Re)Initialize the block communications:
c Set boundary conditions (and conceivably update cij).
         k_sor=-2
         call bbdy(iLs,ifull,iuds,u,k_sor,ndims,idims,lperiod,
     $        icoords,iLcoords,myside,myorig,
     $        icommcart,mycartid,myid)
         mpiid=myid
c         if(mpiid.eq.0)write(*,*)'bbdy (re)initialized'
         return
      endif
c End of control functions.
c      write(*,*) 'Entering bbdy second',ndims,lperiod,iLs,
c     $        icoords,iLcoords,myside,myorig,icommcart,mycartid,myid
c------------------------------------------------------------------
      ierr=0
      omega=1.
      oaddu=0.
c This makes cases with strong faddu effects converge better.
      underrelax=1.2
c Main iteration      
      do k_sor=1,mi_sor
         delta=0.
         umin=1.e30
         umax=0.
         relax=(omega+oaddu)/(1.+underrelax*oaddu) 
         oaddu=0.
c Set boundary conditions (and conceivably update cij).
         call bdyset(ndims,ifull,iuds,cij,u,q)
c Do block boundary communications, returns block info icoords...myid.
         call bbdy(iLs,ifull,iuds,u,k_sor,ndims,idims,lperiod,
     $        icoords,iLcoords,myside,myorig,
     $        icommcart,mycartid,myid)
c If this is found to be an unused node, jump to barrier.
         if(mycartid.eq.-1)goto 999
c Do a relaxation.
c            write(*,*) 'Calling sorrelaxgen',k_sor,myorig,ndims,myside
            call sorrelaxgen(k_sor,ndims,iLs,myside,
     $           cij(1+(2*ndims+1)*(myorig-1)),
     $           u(myorig),
     $           q(myorig),
     $           laddu,faddu,oaddu,relax,delta,umin,umax)
c            write(*,*) 'Return from sorrelaxgen0',k_sor,delta,umin,umax
c Test convergence
         call testifconverged(eps_sor,delta,umin,umax,
     $        lconverged,icommcart)
c         if(myid.eq.0)
c     $        write(*,*)k_sor,delta,umin,umax,lconverged,relax
         if(lconverged.and.k_sor.ge.2)goto 11
c Chebychev acceleration:
         if(k_sor.eq.1)then
            omega=1./(1.-0.5*xjac_sor**2)
         else
            omega=1./(1.-0.25*xjac_sor**2*omega)
         endif
      enddo
c We finished the loop, implies we did not converge.
      k_sor=-mi_sor
c-------------------------------------------------------------------
 11   continue
c Do the final mpi_gather [or allgather if all processes need
c the result].
      nk=-1
      call bbdy(iLs,ifull,iuds,u,nk,ndims,idims,lperiod,
     $        icoords,iLcoords,myside,myorig,
     $        icommcart,mycartid,myid)
c Boundary conditions need to be reset based on the gathered result.
c But that's not sufficient when there's a relaxation so be careful!
      call bdyset(ndims,ifull,iuds,cij,u,q)
      del_sor=delta
      ierr=k_sor
 999  continue
c Indirection is needed here because otherwise the finalize call
c seems to cause the return to fail. Probably unnecessary.
      mpiid=myid
c mpi version needs gracious synchronization when some processes
c are unused by the iteration.
      call MPI_BARRIER(MPI_COMM_WORLD,ierrmpi)
c unfortunately this is too early, and causes an exit.
c      if(lmpisplit)call MPI_FINALIZE(ierrmpi)
      end
c**********************************************************************
c***********************************************************************
c The challenge here is to ensure that all processes decide to end
c at the same time. If not then a process will hang waiting for message.
c So we have to universalize the convergence test. All block must be
c converged, but the total spread depends on multiple blocks.
      subroutine testifconverged(eps,delta,umin,umax,lconverged,
     $     icommcart)
      include 'mpif.h'
      logical lconverged
      real convgd(3)
      convgd(1)=abs(delta)
      convgd(2)=-umin
      convgd(3)=umax
c Here we need to allreduce the data, selecting the maximum values,
c doing it in place.
c      write(*,*)'convgd,icommcart',convgd,icommcart
      call MPI_ALLREDUCE(MPI_IN_PLACE,convgd,3,MPI_REAL,
     $     MPI_MAX,icommcart,ierr)
      if(convgd(1).lt.eps*(convgd(2)+convgd(3))) then
         lconverged=.true.
      else
         lconverged=.false.
      endif
      delta=sign(convgd(1),delta)
      end
c***********************************************************************
c Cut here and throw the rest away for the basic routines
c***********************************************************************
