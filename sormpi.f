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
c         real cij(nd2+1,ifull(1),ifull(2),...)
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
c        but boundary values are never changed inside sormpi.
      real u(*)
c     q  "charge density" input
      real q(*)
c User supplied functions which should be declared external in the 
c calling routine.
c     bdyset subroutine that evaluates the boundary conditions and
c               deposits into the edge values of the potential.
c               bdyset(ndims,ifull,iuds,cij,u,q)
c               If non-constant coefficients are required, then
c               bdyset could be used to adjust their values.
c     faddu(u,fprime,index)  
c               real function that returns the additional component
c               f, and as parameter fprime=df/du.
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

c bbdydecl declares most things for bbdy, using parameter ndimsdecl.
      include 'bbdydecl.f'

      real delta,umin,umax
c This saves data between calls so we can use separate initialization
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
c This jacobi radius is pretty much optimized for Poisson equation with
c fixed boundary, 3-D, but is supposed to be general dimensional.
         maxlen=0
         sumlen=0
         do k=1,ndims
            sumlen=sumlen+iuds(k)
            if(iuds(k).gt.maxlen)maxlen=iuds(k)
         enddo
         xyimb=ndims*maxlen/sumlen - 1.
         xjac_sor=1.- (5./max(10.,sumlen/ndims)**2)*(1.-0.3*xyimb)
         mi_sor=int(2.*sumlen+20)
         eps_sor=1.e-5
      endif
c Second bit of ictlh indicates if there's additional term.
      ictlh=ictlh/2
      if(mod(ictlh,2).ne.0)then
         laddu=.true.
      else
         laddu=.false.
      endif
c      write(*,*)'laddu=',laddu
c Third bit of ictlh indicates we are just initializing.
      ictlh=ictlh/2
      if(mod(ictlh,2).ne.0)then
c (Re)Initialize the block communications:
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
c Experiments:
c      xjac_sor=1.-.6*(1.-xjac_sor)
c      write(*,*)'xjac_sor=',xjac_sor
c      do k_sor=1,mi_sor*2
      do k_sor=1,mi_sor
         delta=0.
         umin=1.e30
         umax=0.
         relax=(omega+oaddu)/(1.+underrelax*oaddu) 
         oaddu=0.
c Set boundary conditions (and conceivably update cij).
c Only needed every other step, and gives identical results.
         if(mod(k_sor,2).eq.1)call bdyset(ndims,ifull,iuds,cij,u,q)
c Do block boundary communications, returns block info icoords...myid.
         call bbdy(iLs,ifull,iuds,u,k_sor,ndims,idims,lperiod,
     $        icoords,iLcoords,myside,myorig,
     $        icommcart,mycartid,myid)
c If this is found to be an unused node, jump to barrier.
         if(mycartid.eq.-1)goto 999
c Do a relaxation.
c         write(*,*)'At sorrelax',myside,myorig,icoords,iLs
c            if(k_sor.le.2)
c         write(*,*) 'Calling sorrelaxgen',delta,oaddu,relax
         call sorrelaxgen(k_sor,ndims,iLs,myside,
     $        cij(1+(2*ndims+1)*(myorig-1)),u(myorig),q(myorig),myorig,
     $        laddu,faddu,oaddu,
     $        relax,delta,umin,umax)

c          call checkdelta(delta,deltaold)
c          if(k_sor.le.2)
c            write(*,'(''sorrelaxgen'',i4,4f13.5)')k_sor,delta,oaddu
c     $           ,umin,umax
c Test convergence
         call testifconverged(eps_sor,delta,umin,umax,
     $        lconverged,icommcart)
c         if(myid.eq.0)
c     $        write(*,*)k_sor,delta,umin,umax,lconverged
c     $        ,relax,omega
         if(lconverged.and.k_sor.ge.2)goto 11

         if(k_sor.eq.1)then
c Chebychev acceleration doesn't work very well at start.
c Better not to start with such a big value of omega.
c It tends to go unstable. But for now, I'm leaving the original
c for checking purposes.
            omega=1./(1.-0.5*xjac_sor**2)
c            omega=1./(1.-0.45*xjac_sor**2)
         else
            omega2=1./(1.-0.25*xjac_sor**2*omega)
            omega=omega2
         endif
      enddo
c We finished the loop, implies we did not converge.
      k_sor=-mi_sor
      ierr=-1
c      write(*,*)'Finished sorrelax loop unconverged',k_sor,delta
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
      call bdyset(ndims,ifull,iuds,cij,u,q)
      del_sor=delta
      ierr=k_sor
 999  continue
c Indirection is needed here because otherwise the finalize call
c seems to cause the return to fail. Probably unnecessary.
      mpiid=myid
      end
c**********************************************************************
c***********************************************************************
c The challenge here is to ensure that all processes decide to end
c at the same time. If not then a process will hang waiting for message.
c So we have to universalize the convergence test. All blocks must be
c converged, but the total spread depends on multiple blocks.
      subroutine testifconverged(eps,delta,umin,umax,lconverged,
     $     icommcart)
      implicit none
      real eps,delta,umin,umax
      integer ierr,icommcart
      logical lconverged
      real convgd(3)
      convgd(1)=abs(delta)
      convgd(2)=-umin
      convgd(3)=umax
c Here we need to allreduce the data, selecting the maximum values,
      call mpiconvgreduce(convgd,icommcart,ierr)
c........
      if(convgd(1).lt.eps*(convgd(2)+convgd(3))) then
         lconverged=.true.
      else
         lconverged=.false.
      endif
      delta=sign(convgd(1),delta)
      end
c**********************************************************************
c This is of dubious sense.
      subroutine avederiv(x,xval,xd,ns)
c return value and change of x averaged over the last ns (<ntot) steps.
c If ns=1, return the value of x now, and the change since last call.
c Be careful with the first call which might have rubbish in it.
c To initialize to xval and xd, call twice with ns=1.
      data xave,xdav,xlast/0.,0.,0./
      save xave,xdav,xlast
      
      if(ns.eq.0)stop 'avederiv error ns=0'
      xave=(xave*(ns-1) + x)/ns
      xval=xave
      xdn=x-xlast
      xdav=(xdav*(ns-1) + xdn)/ns
      xd=xdav
      xlast=x
      end
c***********************************************************************
