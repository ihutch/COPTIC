! Solve an elliptical problem in ndims dimensions
! by successive overrelaxation
!     L(u) + f(u) = q(x,y,...), 
!     where L is a second order elliptical differential operator 
!     represented by a difference stencil of specified coefficients,
!     f is some additional function, and q is the "charge density".
! Ian Hutchinson, 2006-2014
!
      subroutine sormpi(ndims,ifull,iuds,cij,u,q,bdyshare,bdyset,faddu
     $     ,ictl,ierr,mpiid,idims)
! The used number of dimensions. But this must be equal to ndims 
! parameter used for bbdydecl dimensions.
      integer ndims
!     The number of dimensions declared in bbdydecl.f
      parameter (ndimsbbdy=3)
!     ifull full dimensions, iuds used dimensions
!     cij coefficients of the finite difference scheme, in the order
!         east,west,north,south ... [regarding i,j, as x,y].
!         real cij(nd2+1,ifull(1),ifull(2),...)
!     The last cij value of the first index is a pointer to object info
!     and if it is zero (null) then no object code affects this point.
!     Objects cut the mesh between nodes and e.g. specify u there.
!     They can thus define objects embedded in the mesh.
!
!     In this mpi routine we must use linear addressing for cij,u,q
!     so that the pointers can be used for each block.
      real cij(*)
!     u  potential to be solved for (initialized on entry).
!        its boundaries are at 1,ni; 1,nj; 
!        which are normally set on entry or by bdyset,
!        but boundary values are never changed inside sormpi.
      real u(*)
!     q  "charge density" input
      real q(*)
! User supplied functions which should be declared external in the 
! calling routine.
!     bdyshare subroutine sets boundary conditions. It has the same
!              argument list as bbdy, declared in bbdydecl.f. 
!              If it returns idone(1)=0, fall back to the following.
!     bdyset subroutine that evaluates the boundary conditions and
!               deposits into the edge values of the potential.
!               bdyset(ndims,ifull,iuds,cij,u,q)
!     faddu(u,fprime,index)  
!               real function that returns the additional component
!               f, and as parameter fprime=df/du.
      external faddu,bdyshare
!     ictl  integer control switches
!           bit 1: user-defined iteration parameters (default no)
!           bit 2: use faddu (default no)
!           bit 3: just initialize bbdy (return mpiid).
!           bits ib=4 on: if 1 then dimension ib-3 is periodic.
!     ierr  Returns Positive: number of iterations to convergence.
!           Negative: Not-converged maximum iterations.
!           Zero: something bad.
!     mpiid Returns my MPI process number for this mpi version.
      integer ictl,ierr
!     integer idims(*)  number of blocks, declared in bbdydecl.f


! Other things that we might want control over include the maximum
! number of iterations, the Jacobi radius, and the convergence size.
! k_sor is the sor iteration index, for diagnostics.
      common /ctl_sor/mi_sor,xjac_sor,eps_sor,del_sor,k_sor

      logical lconverged
      logical laddu
      logical ldebugs

! bbdydecl declares most things for bbdy, using parameter ndimsbbdy.
      include 'bbdydecl.f'
! iLs has been removed from bbdydecl.f
      integer iLs(ndimsbbdy+1)

! Scratch arrays for bdyshare communications
      integer idone(ndimsbbdy)
      integer zeros(ndimsbbdy)
      integer ones(ndimsbbdy)

      real delta,umin,umax
      data ldebugs/.false./
! This saves data between calls so we can use separate initialization
      save

!-------------------------------------------------------------------
      if(ndims.gt.ndimsbbdy)then
         write(*,*)'Too many dimensions in sormpi call',
     $        ndims,ndimsbbdy
         stop
      endif
      if(ldebugs)write(*,*)'In sormpi',ictl,ndims,ifull,iuds,idims,ierr
!      return
!--------------------------------------------------------------------
! Decide if any of the dimensions is periodic
      ictlh=ictl
      ictlh=ictlh/8
! Initialize and set accordingly.
      do nd=1,ndims
         idone(nd)=0
         zeros(nd)=0
         ones(nd)=1
         if(mod(ictlh,2).eq.0)then
            lperiod(nd)=.false.
         else
!            write(*,*)'Dimension',nd,'  periodic.'
            lperiod(nd)=.true.
         endif
         ictlh=ictlh/2
      enddo
!--------------------------------------------------------------------
! Control Functions:
      ictlh=ictl
! First bit of ictlh indicates if sorctl preset or defaults.
      if(mod(ictlh,2).eq.0)then
! This jacobi radius is pretty much optimized for Poisson equation with
! fixed boundary, 3-D, but is supposed to be general dimensional.
         maxlen=0
         sumlen=0
         do k=1,ndims
            sumlen=sumlen+iuds(k)
            if(iuds(k).gt.maxlen)maxlen=iuds(k)
         enddo
         xyimb=ndims*maxlen/sumlen - 1.
         xjac_sor=1.- (5./max(10.,sumlen/ndims)**2)*(1.-0.3*xyimb)
! More accurate calculation of the Jacobi radius is probably better.
         if(.true.)then
            call getxjac(ndims,ifull,iuds,cij,xjac)
!            write(*,'(a,f10.6,a,f10.6)')
!     $           'Old xjac',xjac_sor,' New xjac',xjac
            xjac_sor=xjac
         endif
         mi_sor=int(2.*sumlen+20)
         eps_sor=1.e-5
      endif
! Second bit of ictlh indicates if there's additional term.
      ictlh=ictlh/2
      if(mod(ictlh,2).ne.0)then
         laddu=.true.
      else
! Plain Poisson convergence with electrons is often slower.
         mi_sor=mi_sor*2
         laddu=.false.
      endif
!      write(*,*)'laddu=',laddu
! Third bit of ictlh indicates we are just initializing.
      ictlh=ictlh/2
      if(mod(ictlh,2).ne.0)then
! bbdy() is the only place the parallel nature of this routine affects
! A serial version is obtained by linking with nonmpibbdy.f, rather than
! mpibbdy.f so as to replace the mpi calls with dummies. Those files
! also contain a few other routines called here.
! (Re)Initialize the block communications:
         k_sor=-2
         call bbdy(iLs,ifull,iuds,u,k_sor,ndims,idims,
     $        icoords,iLcoords,myside,myorig,
     $        icommcart,mycartid,myid,lperiod)
         mpiid=myid
         if(mpiid.eq.0)write(*,*)'bbdy (re)initialized'
         return
      endif
! End of control functions.
!------------------------------------------------------------------
      ierr=0
      omega=1.
      oaddu=0.
! This makes cases with strong faddu effects converge better.
      underrelax=1.2
!------------------------------------------------------------------
! Main iteration      
      do k_sor=1,mi_sor
         delta=0.
         umin=1.e30
         umax=0.
         relax=(omega+oaddu)/(1.+underrelax*oaddu) 
         oaddu=0.
!------------------------------------------------------------------
! Do block boundary communications, returns block info icoords...myid.
         call bbdy(iLs,ifull,iuds,u,k_sor,ndims,idims
     $        ,icoords,iLcoords,myside,myorig,
     $        icommcart,mycartid,myid,lperiod)
! If this is found to be an unused node, jump to barrier.
         if(mycartid.eq.-1)goto 999
! If we are running without MPI then set periodic conditions explicitly.
         if(icommcart.eq.0)idone(2)=1
!------------------------------------------------------------------
! Set boundary conditions (and conceivably update cij).
! Only needed every other step, and gives identical results.
         if(mod(k_sor,2).eq.1)then
! The parallelized boundary setting routine
            idone(1)=0
            call bdyshare(ifull,iuds,u,idone,ndims,idims,
     $        icoords,iLcoords,myside,myorig,
     $        icommcart,mycartid,myid,lperiod)
! If this did not succeed. Fall back to global setting.
!            write(*,*)idone(1)
            if(idone(1).eq.0)call bdyset(ndims,ifull,iuds,cij,u,q)
         endif
!-------------------------------------------------------------------
! Do a relaxation.
!         write(*,*)'At sorrelax',myside,myorig,icoords,iLs
!            if(k_sor.le.2)
!         write(*,*) 'Calling sorrelaxgen',delta,oaddu,relax
         call sorrelaxgen(k_sor,ndims,iLs,myside,
     $        cij(1+(2*ndims+1)*(myorig-1)),
     $        cij(1+(2*ndims+1)*(myorig-1)),u(myorig),q(myorig),myorig,
     $        laddu,faddu,oaddu,
     $        relax,delta,umin,umax)

!          call checkdelta(delta,deltaold)
         if(ldebugs.and.k_sor.le.2) write(*,
     $        '(''sorrelaxgen returned'',i4,4f13.5)')
     $        k_sor,delta,oaddu ,umin,umax
! Test convergence
         call testifconverged(eps_sor,delta,umin,umax,
     $        lconverged,icommcart)
!         if(myid.eq.0) write(*,'(i5,f10.6,2f8.4,l3,2f8.4)')k_sor,delta
!     $        ,umin,umax,lconverged ,relax,omega
         if(lconverged.and.k_sor.ge.2)goto 11

         if(k_sor.eq.1)then
! Chebychev acceleration doesn't work very well at start.
! Better not to start with such a big value in omega.
! Otherwise it tends to go unstable. 
            omega=1./(1.-0.45*xjac_sor**2)
! One might leave the original for checking purposes:
!            omega=1./(1.-0.5*xjac_sor**2)
         else
! The value 0.25 is too close to instability. So
            omega2=1./(1.-0.245*xjac_sor**2*omega)
            omega=omega2
         endif
      enddo
!------------------------------------------------------------------
! We finished the loop, implies we did not converge.
      k_sor=-mi_sor
      ierr=-1
!      write(*,*)'Finished sorrelax loop unconverged',k_sor,delta
!-------------------------------------------------------------------
 11   continue
! Do the final mpi_gather [or allgather if all processes need
! the result].
      nk=-1
      call bbdy(iLs,ifull,iuds,u,nk,ndims,idims,
     $        icoords,iLcoords,myside,myorig,
     $        icommcart,mycartid,myid,lperiod)
! Boundary conditions need to be fully set based on the gathered result
! Call the parallelized boundary setting routine but lie to it that 
! it is the only process. Also insist on explicit setting.
      idone(2)=1
      idone(1)=0
! Every process must do this.
      call bdyshare(ifull,iuds,u,idone,ndims,ones,
     $        zeros,iLcoords,iuds,ones,
     $        icommcart,mycartid,myid,lperiod)
      del_sor=delta
      ierr=k_sor
 999  continue
! Indirection is needed here because otherwise the finalize call
! seems to cause the return to fail. Probably unnecessary.
      mpiid=myid
! If we need to broadcast more widely, do it.
      call meshbroadcast(u,ndims,ifull,iuds,iLs,0,icommcart)
      end      
!****************************************************************
!***********************************************************************
! The challenge here is to ensure that all processes decide to end
! at the same time. If not then a process will hang waiting for message.
! So we have to universalize the convergence test. All blocks must be
! converged, but the total spread depends on multiple blocks.
      subroutine testifconverged(eps,delta,umin,umax,lconverged,
     $     icommcart)
      implicit none
      real eps,delta,umin,umax
      integer ierr,icommcart
      logical lconverged
      real convgd(3)
      real cscale
! This determines a minimum convergence scale for uniform cases.
      parameter (cscale=0.01)
      convgd(1)=abs(delta)
      convgd(2)=-umin
      convgd(3)=umax
! Here we need to allreduce the data, selecting the maximum values,
      call mpiconvgreduce(convgd,icommcart,ierr)
!........
      if(convgd(1).le.eps*max(convgd(2)+convgd(3),cscale)) then
         lconverged=.true.
      else
         lconverged=.false.
      endif
      delta=sign(convgd(1),delta)
      end
!**********************************************************************
      subroutine  getxjac(ndims,ifull,iuds,cij,xjac)
! Get the Jacobi convergence radius assuming uniform (anisotropic) mesh
! M squared = [sum_i(1/dx_i^2)]/[sum_i(1/(dx_i^2*N_i^2))]
! basing the 1/dx_i^2 on the cij for each direction. 
      integer ndims,ifull(ndims),iuds(ndims)
      real cij(*)
!      real cij(7,160,129,3)
      parameter (ndimsmax=3)
      real ctyp(ndimsmax)

      if(ndims.gt.ndimsmax) stop 'getxjac error ndims exceeds ndimsmax'

      istep=2*ndims+1
      ipointer=0
      iref=2
! Point to cij position 1,iref,iref,iref
      do i=1,ndims
         ipointer=ipointer+(iref-1)*istep
         istep=istep*ifull(i)
      enddo
! Characterize the mesh spacings in different dimensions by
! ctyp is 1/Dx^2 in each direction
      dxsum=0.
      dxNsum=0.
!      write(*,*)(cij(ipointer+i),i=1,2*ndims)
      do i=1,ndims
         ctyp(i)=cij(ipointer+2*i)
!         write(*,*)'i, ipointer, ctyp, iuds',i,ipointer,ctyp(i),iuds(i)
         dxsum=dxsum+ctyp(i)
         dxNsum=dxNsum+ctyp(i)/iuds(i)**2         
      enddo
! Jacobi radius is 1- (\pi^2/2) dxNsum/dxsum 
      xjac=1.-3.1415926**2/2. *dxNsum/dxsum
      end
