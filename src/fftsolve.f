!! Optional test program main
!      call ffttest
!      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      subroutine fftsolve3d(M,N,O,inout,xL,yL,zL)
! Solve the poisson equation by fourier transforms for a periodic 2D domain
! of lengths xL, yL. The matrix is complex*8 inout(M,N) entirely filled.
! \nabla^2\phi = -rho/epsilon0
! On input it contains rho/epsilon0. On output it contains \phi.
      include 'fftw3.f'
      integer M,N,O
      double complex inout(M,N,O)
      integer*8 planforward,planbackward
      real k2
      integer idebug,ifirst
      logical lreuse
      data idebug/0/ifirst/0/lreuse/.true./
      save
      if(idebug.gt.1)write(*,*)'Input'
      if(idebug.gt.1)write(*,'(10f8.4)')real(inout)

      if(.not.lreuse.or.ifirst.eq.0)then
         call dfftw_plan_dft_3d(planforward, M,N,O, inout,inout,
     &        FFTW_FORWARD, FFTW_ESTIMATE)
         call dfftw_plan_dft_3d(planbackward, M,N,O, inout,inout,
     &        FFTW_BACKWARD, FFTW_ESTIMATE)
         ifirst=ifirst+1
      endif

      call dfftw_execute_dft(planforward, inout, inout)
      if(idebug.gt.1)write(*,*)'Transformed'
      if(idebug.gt.1)write(*,'(10f8.4)')inout

      scale=M*N*O
      do k=1,O
         kz=k-1
         if(kz.gt.O/2)kz=kz-O
         skz=2.*sin(kz*3.1415926/O)
         do j=1,N
            ky=j-1
            if(ky.gt.N/2)ky=ky-N
            sky=2.*sin(ky*3.1415926/N)
            do i=1,M
               kx=i-1
               if(kx.gt.M/2)kx=kx-M
! The inefficiency of recalculating sine costs about 20% in speed.
               skx=2.*sin(kx*3.1415926/M)
               k2=scale*((skx*M/xL)**2+(sky*N/yL)**2+(skz*O/zL)**2)
               if(k2.eq.0.)then
! Set the average charge to zero. Otherwise nonphysical results. 
                  inout(i,j,k)=0.
               else
                  inout(i,j,k)=-inout(i,j,k)/k2
               endif
            enddo
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
      subroutine fftphisolve3d(ifull,iuds,u,q,cscratch,xL,yL,zL,ierr)
! Copy to/from complex scratch and solve poisson's equation.
! On input the charge is in q(2:iuds(1)-1,2:iuds(2)-1) guard cells not used
! On output the potential is in u(1:iuds(1),1:iuds(1)) including guards.
      integer ifull(3),iuds(3)
      real u(ifull(1),ifull(2),ifull(3)),q(ifull(1),ifull(2),ifull(3))
      real xL,yL,zL
      double complex cscratch(iuds(1)-2,iuds(2)-2,iuds(3)-2)
      integer M,N,O

      M=iuds(1)-2
      N=iuds(2)-2
      O=iuds(3)-2
      do k=1,O 
         do j=1,N               ! Copy q without the guard cells to scratch.
            do i=1,M
               cscratch(i,j,k)=-q(i+1,j+1,k+1)
            enddo
         enddo
      enddo
      call fftsolve3d(M,N,O,cscratch,xL,yL,zL) ! Solve only on true cells.
      do k=1,O+2
         do j=1,N+2
            do i=1,M+2          ! Fill back u including periodic guard cells.
               u(i,j,k)=real(cscratch(mod(i-2+M,M)+1,mod(j-2+N,N)+1
     $              ,mod(k-2+O,O)+1),4)
            enddo
         enddo
      enddo
      ierr=0
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This testing routine shows there's no real speed advantage in a 2d call
! relative to a 3d call with dimension 1 in the final. 
      subroutine ffttest
      integer M,N,O
      parameter (M=500,N=300,O=1)      ! Mesh counts
      double complex inout(M,N)  ! Charge density
      double complex inout3(M,N,O)  ! Charge density
      equivalence (inout,inout3)
      real worka(M,N),workb(M,N)
      real xL,yL,zL   ! Lengths of domain
      real xn,yn,zn   ! Wavelength numbers
      integer n2cycles,n3cycles
      data xL/1/yL/1/zL/1/xn/3/yn/2/zn/0/
      data n2cycles/0/n3cycles/1000/

      sk=2.*3.1415926*sqrt((xn/xL)**2+(yn/yL)**2)


! Now the same using 3d routines
      do ic=1,n3cycles ! Time for 1000 500x300 calls 10.6s. The same.
         ! With the corrected sinc expression 12.5s.
      do j=1,N
         do i=1,M
            worka(i,j)=cos(2*3.1415926*xn*(i-1.)/M)
     $                *cos(2*3.1415926*yn*(j-1.)/N)
         enddo
      enddo
      inout=worka

      call fftsolve3d(M,N,O,inout3,xL,yL,zL)
      workb=real(inout,4)
      call minmax2(workb,M,M,N,zmin,zmax)
      write(*,*)'3d phimax should nearly equal 1/k**2',zmax,1/sk**2
      enddo

      if(n3cycles.eq.1)then
      call multiframe(2,1,3)
      call autocolcont(worka,M,M,N)
      call autocolcont(workb,M,M,N)
      call pltend
      endif

!      write(*,*)'If the plots and the above are the same we are good.'

      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
