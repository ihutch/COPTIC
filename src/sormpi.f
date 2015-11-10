c Solve an elliptical problem in ndims dimensions
c by successive overrelaxation
c     L(u) + f(u) = q(x,y,...), 
c     where L is a second order elliptical differential operator 
c     represented by a difference stencil of specified coefficients,
c     f is some additional function, and q is the "charge density".
c Ian Hutchinson, 2006-2014
c
      subroutine sormpi(ndims,ifull,iuds,cij,u,q,bdyshare,bdyset,faddu
     $     ,ictl,ierr,mpiid,idims)
c The used number of dimensions. But this must be equal to ndims 
c parameter used for bbdydecl dimensions.
      integer ndims
c     The number of dimensions declared in bbdydecl.f
      parameter (ndimsbbdy=3)
c     ifull full dimensions, iuds used dimensions
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
c     bdyshare subroutine sets boundary conditions. It has the same
c              argument list as bbdy, declared in bbdydecl.f. 
c              If it returns idone(1)=0, fall back to the following.
c     bdyset subroutine that evaluates the boundary conditions and
c               deposits into the edge values of the potential.
c               bdyset(ndims,ifull,iuds,cij,u,q)
c     faddu(u,fprime,index)  
c               real function that returns the additional component
c               f, and as parameter fprime=df/du.
      external faddu,bdyshare
c     ictl  integer control switches
c           bit 1: user-defined iteration parameters (default no)
c           bit 2: use faddu (default no)
c           bit 3: just initialize bbdy (return mpiid).
c           bits ib=4 on: if 1 then dimension ib-3 is periodic.
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
      logical ldebugs

c bbdydecl declares most things for bbdy, using parameter ndimsbbdy.
      include 'bbdydecl.f'
c iLs has been removed from bbdydecl.f
      integer iLs(ndimsbbdy+1)

c Scratch arrays for bdyshare communications
      integer idone(ndimsbbdy)
      integer zeros(ndimsbbdy)
      integer ones(ndimsbbdy)

      real delta,umin,umax
      data ldebugs/.false./
c This saves data between calls so we can use separate initialization
      save
c-------------------------------------------------------------------
      if(ndims.gt.ndimsbbdy)then
         write(*,*)'Too many dimensions in sormpi call',
     $        ndims,ndimsbbdy
         stop
      endif
      if(ldebugs)write(*,*)'In sormpi',ictl,ndims,ifull,iuds,idims,ierr
c      return
c--------------------------------------------------------------------
c Decide if any of the dimensions is periodic
      ictlh=ictl
      ictlh=ictlh/8
c Initialize and set accordingly.
      do nd=1,ndims
         idone(nd)=0
         zeros(nd)=0
         ones(nd)=1
         if(mod(ictlh,2).eq.0)then
            lperiod(nd)=.false.
         else
c            write(*,*)'Dimension',nd,'  periodic.'
            lperiod(nd)=.true.
         endif
         ictlh=ictlh/2
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
c Plain Poisson convergence with electrons is often slower.
         mi_sor=mi_sor*2
         laddu=.false.
      endif
c      write(*,*)'laddu=',laddu
c Third bit of ictlh indicates we are just initializing.
      ictlh=ictlh/2
      if(mod(ictlh,2).ne.0)then
c bbdy() is the only place the parallel nature of this routine affects
c A serial version is obtained by linking with nonmpibbdy.f, rather than
c mpibbdy.f so as to replace the mpi calls with dummies. Those files
c also contain a few other routines called here.
c (Re)Initialize the block communications:
         k_sor=-2
         call bbdy(iLs,ifull,iuds,u,k_sor,ndims,idims,
     $        icoords,iLcoords,myside,myorig,
     $        icommcart,mycartid,myid,lperiod)
         mpiid=myid
         if(mpiid.eq.0)write(*,*)'bbdy (re)initialized'
         return
      endif
c End of control functions.
c------------------------------------------------------------------
      ierr=0
      omega=1.
      oaddu=0.
c This makes cases with strong faddu effects converge better.
      underrelax=1.2
c------------------------------------------------------------------
c Main iteration      
      do k_sor=1,mi_sor
         delta=0.
         umin=1.e30
         umax=0.
         relax=(omega+oaddu)/(1.+underrelax*oaddu) 
         oaddu=0.
c------------------------------------------------------------------
c Do block boundary communications, returns block info icoords...myid.
         call bbdy(iLs,ifull,iuds,u,k_sor,ndims,idims
     $        ,icoords,iLcoords,myside,myorig,
     $        icommcart,mycartid,myid,lperiod)
c If this is found to be an unused node, jump to barrier.
         if(mycartid.eq.-1)goto 999
c If we are running without MPI then set periodic conditions explicitly.
         if(icommcart.eq.0)idone(2)=1
c------------------------------------------------------------------
c Set boundary conditions (and conceivably update cij).
c Only needed every other step, and gives identical results.
         if(mod(k_sor,2).eq.1)then
c The parallelized boundary setting routine
            idone(1)=0
            call bdyshare(ifull,iuds,u,idone,ndims,idims,
     $        icoords,iLcoords,myside,myorig,
     $        icommcart,mycartid,myid,lperiod)
c If this did not succeed. Fall back to global setting.
c            write(*,*)idone(1)
            if(idone(1).eq.0)call bdyset(ndims,ifull,iuds,cij,u,q)
         endif
c-------------------------------------------------------------------
c Do a relaxation.
c         write(*,*)'At sorrelax',myside,myorig,icoords,iLs
c            if(k_sor.le.2)
c         write(*,*) 'Calling sorrelaxgen',delta,oaddu,relax
         call sorrelaxgen(k_sor,ndims,iLs,myside,
     $        cij(1+(2*ndims+1)*(myorig-1)),
     $        cij(1+(2*ndims+1)*(myorig-1)),u(myorig),q(myorig),myorig,
     $        laddu,faddu,oaddu,
     $        relax,delta,umin,umax)

c          call checkdelta(delta,deltaold)
         if(ldebugs.and.k_sor.le.2) write(*,
     $        '(''sorrelaxgen returned'',i4,4f13.5)')
     $        k_sor,delta,oaddu ,umin,umax
c Test convergence
         call testifconverged(eps_sor,delta,umin,umax,
     $        lconverged,icommcart)
c         if(myid.eq.0) write(*,'(i5,f10.6,2f8.4,l3,2f8.4)')k_sor,delta
c     $        ,umin,umax,lconverged ,relax,omega
         if(lconverged.and.k_sor.ge.2)goto 11

         if(k_sor.eq.1)then
c Chebychev acceleration doesn't work very well at start.
c Better not to start with such a big value of omega.
c Otherwise it tends to go unstable. 
            omega=1./(1.-0.45*xjac_sor**2)
c One might leave the original for checking purposes:
c            omega=1./(1.-0.5*xjac_sor**2)
         else
            omega2=1./(1.-0.25*xjac_sor**2*omega)
            omega=omega2
         endif
      enddo
c------------------------------------------------------------------
c We finished the loop, implies we did not converge.
      k_sor=-mi_sor
      ierr=-1
c      write(*,*)'Finished sorrelax loop unconverged',k_sor,delta
c-------------------------------------------------------------------
 11   continue
c Do the final mpi_gather [or allgather if all processes need
c the result].
      nk=-1
      call bbdy(iLs,ifull,iuds,u,nk,ndims,idims,
     $        icoords,iLcoords,myside,myorig,
     $        icommcart,mycartid,myid,lperiod)
c Boundary conditions need to be fully set based on the gathered result
c Call the parallelized boundary setting routine but lie to it that 
c it is the only process. Also insist on explicit setting.
      idone(2)=1
      idone(1)=0
c Every process must do this.
      call bdyshare(ifull,iuds,u,idone,ndims,ones,
     $        zeros,iLcoords,iuds,ones,
     $        icommcart,mycartid,myid,lperiod)
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
      real cscale
c This determines a minimum convergence scale for uniform cases.
      parameter (cscale=0.01)
      convgd(1)=abs(delta)
      convgd(2)=-umin
      convgd(3)=umax
c Here we need to allreduce the data, selecting the maximum values,
      call mpiconvgreduce(convgd,icommcart,ierr)
c........
      if(convgd(1).le.eps*max(convgd(2)+convgd(3),cscale)) then
         lconverged=.true.
      else
         lconverged=.false.
      endif
      delta=sign(convgd(1),delta)
      end
c**********************************************************************

