!! Optional test program main
!      call ffttest
!      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      subroutine fftsolve(M,N,inout,xL,yL)
! Solve the poisson equation by fourier transforms for a periodic 2D domain
! of lengths xL, yL. The matrix is complex*8 inout(M,N) entirely filled.
! \nabla^2\phi = -rho/epsilon0
! On input it contains rho/epsilon0. On output it contains \phi.
      include 'fftw3.f'
      integer M,N
      double complex inout(M,N)
      integer*8 planforward,planbackward
      real k2
      integer idebug,ifirst
      logical lreuse
      data idebug/0/ifirst/0/lreuse/.true./
      save
      if(idebug.gt.1)write(*,*)'Input'
      if(idebug.gt.1)write(*,'(10f8.4)')real(inout)

      if(.not.lreuse.or.ifirst.eq.0)then
         call dfftw_plan_dft_2d(planforward, M,N, inout,inout,
     &        FFTW_FORWARD, FFTW_ESTIMATE)
         call dfftw_plan_dft_2d(planbackward, M,N, inout,inout,
     &        FFTW_BACKWARD, FFTW_ESTIMATE)
         ifirst=ifirst+1
      endif

      call dfftw_execute_dft(planforward, inout, inout)
      if(idebug.gt.1)write(*,*)'Transformed'
      if(idebug.gt.1)write(*,'(10f8.4)')inout

      scale=M*N*(2*3.1415926)**2
      do j=1,N
         ky=j-1
         if(ky.gt.N/2)ky=ky-N
!         write(*,*)ky
         do i=1,M
            kx=i-1
            if(kx.gt.M/2)kx=kx-M
! This divides by -(k^2*M*N) to invert the Laplacian. 
! k_j=2\pi (j-1)/(N*DY) etc. N*DY= Y-length. And M*N to scale.
            k2=scale*((kx/xL)**2+(ky/yL)**2)
            if(k2.eq.0.)then
! Set the average charge to zero. Otherwise nonphysical results. 
               inout(i,j)=0.
            else
               inout(i,j)=-inout(i,j)/k2
            endif
         enddo
      enddo
      call dfftw_execute_dft(planbackward, inout, inout)

      if(idebug.gt.1)write(*,*)'Transformed back, in='
      if(idebug.gt.1)write(*,'(10f8.4)')real(inout)

! If we reuse the plans we can't destroy on return. 
! However, reusing plan for big arrays saves only a few percent of time.
      if(.not.lreuse)then
         call dfftw_destroy_plan(planforward)
         call dfftw_destroy_plan(planbackward)
      endif

      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine fftphisolve2d(ifull,iuds,u,q,cscratch,xL,yL,ierr)
! Copy to/from complex scratch and solve poisson's equation.
! On input the charge is in q(2:iuds(1)-1,2:iuds(2)-1) guard cells not used
! On output the potential is in u(1:iuds(1),1:iuds(1)) including guards.
! That is why the arguments in the call must be q(1,1,2) u(1,1,1)
      integer ifull(2),iuds(2)
      real u(ifull(1),ifull(2),3),q(ifull(1),ifull(2))
      real xL,yL
      double complex cscratch(iuds(1)-2,iuds(2)-2)
      integer M,N

      M=iuds(1)-2
      N=iuds(2)-2
      do j=1,N       ! Copy q without the guard cells to scratch.
         do i=1,M
            cscratch(i,j)=-q(i+1,j+1)
         enddo
      enddo
      call fftsolve(M,N,cscratch,xL,yL) ! Solve only on true cells.
      do j=1,N+2
         do i=1,M+2   ! Fill back u including periodic guard cells.
            u(i,j,2)=real(cscratch(mod(i-2+M,M)+1,mod(j-2+N,N)+1),4)
            u(i,j,1)=u(i,j,2)  ! Also need z-guards.
            u(i,j,3)=u(i,j,2)
         enddo
      enddo
      ierr=0
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine ffttest
      integer M,N
      parameter (M=500,N=300)      ! Mesh counts
      double complex inout(M,N)  ! Charge density
      real worka(M,N),workb(M,N)
      real xL,yL   ! Lengths of domain
      real xn,yn   ! Wavelength numbers
      data xL/1/yL/1/xn/3/yn/2/


      sk=2.*3.1415926*sqrt((xn/xL)**2+(yn/yL)**2)

!      do ic=1,100   ! Time for 100 500x300 calls 1.8s including planning.
      do j=1,N
         do i=1,M
            worka(i,j)=cos(2*3.1415926*xn*(i-1.)/M)
     $                *cos(2*3.1415926*yn*(j-1.)/N)
         enddo
      enddo
      inout=worka

      call fftsolve(M,N,inout,xL,yL)
      workb=real(inout,4)
      call minmax2(workb,M,M,N,zmin,zmax)
      write(*,*)'phimax should nearly equal 1/k**2',zmax,1/sk**2
!      enddo

      call multiframe(2,1,3)
      call autocolcont(worka,M,M,N)
      call autocolcont(workb,M,M,N)
      call pltend

      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
