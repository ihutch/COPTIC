!! Optional test program main
!      call ffttest
!      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      subroutine fftsolve3d(M,N,O,inout,xL,yL,zL)
! Solve the poisson equation by fourier transforms for a periodic 3D domain
! of lengths xL, yL, zL. The matrix is complex*8 inout(M,N,O) entirely filled.
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

!      call ffttrid(ifull,iused,u,q,cscratch,xL,yL,zL)! Or
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
      integer M,N,O,MP2,NP2,OP2
      parameter (M=200,N=300,O=1)      ! Mesh counts
      parameter (MP2=M+2,NP2=N+2,OP2=O+2)
      double complex inout(M,N)  ! Charge density
      double complex inout3(M,N,O)  ! Charge density
      double complex cscratch(N,O,M) ! Reordered!      
      real uinout(MP2,NP2,OP2), qin(MP2,NP2,OP2)
      equivalence (inout,inout3), (inout3,cscratch)
      real worka(M,N),workb(M,N),workc(M,N)
      real xL,yL,zL   ! Lengths of domain
      real xn,yn,zn   ! Wavelength numbers
      integer n2cycles,n3cycles
      integer ifull(3),iused(3)
      data ifull/MP2,NP2,OP2/iused/MP2,NP2,OP2/
      data xL/1/yL/1/zL/1/xn/4/yn/1/zn/0/
      data n2cycles/0/n3cycles/1000/

      n3cycles=1
      sk=2.*3.1415926*sqrt((xn/xL)**2+(yn/yL)**2)


! Now the same using 3d routines
      do ic=1,n3cycles ! Time for 1000 500x300 calls 10.6s. The same.
         ! With the corrected sinc expression 12.5s.
      do j=1,N
         do i=1,M
            worka(i,j)=sin(2*3.1415926*xn*(i-1/2.)/M)
     $                *cos(2*3.1415926*yn*(j-1/2.)/N)
         enddo
      enddo
! 3D version
! Older test direct call. 
!         inout=-worka
!         call fftsolve3d(M,N,O,inout3,xL,yL,zL)
!         workb=real(inout,4)
! Test full external call.
      qin=0.
      qin(2:M+1,2:N+1,2)=worka
      call fftphisolve3d(ifull,iused,uinout,qin,cscratch, xL,yL,zL,ier)
      workb=uinout(2:M+1,2:N+1,2)
! 2D transform version
      call ffttrid(ifull,iused,uinout,qin,inout,xL,yL,zL)
      workc=uinout(2:M+1,2:N+1,2)
      call minmax2(workb,M,M,N,zmin,zmax)
      call minmax2(workc,M,M,N,cmin,cmax)
      write(*,'(a,g12.6,a,g12.6)')'3D FFT phimax ',zmax
     $     ,' should nearly equal 1/k**2 ',1/sk**2
      write(*,'(a,g12.6,a,g12.6)')'2D-Tri phimax ',cmax
     $     ,' should nearly equal 1/k**2 ',1/sk**2
      enddo
      call minmax2(workb-workc,M,M,N,dmin,dmax)
      write(*,'(a,g14.6,f10.6)')'max difference and frac',dmax,dmax/zmax

      if(n3cycles.eq.1)then
      call multiframe(3,1,3)
      call autocolcont(workb,M,M,N)
      call autocolcont(workc,M,M,N)
      call autocolcont(workb-workc,M,M,N)
      call pltend
      endif


!      write(*,*)'If the plots and the above are the same we are good.'

      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine ffttrid(ifull,iuds,u,q,cscratch,xL,yL,zL)
! Solve the poisson equation by fourier transforms for periodic 2D directions
! y, z, but by tridiagonal elimination and back substitution in x direction.
! In x-direction domain is full; in y,z, it is short 1/2 mesh.
! Lengths are xL, yL, zL. Work matrix complex*8 cscratch(M+2,N+2,O+2).
! \nabla^2\phi = -rho/epsilon0; 
! q is rho/epsilon0 on entry, u is phi on exit with guard values.
      include 'fftw3.f'
      include 'griddecl.f'  ! Maybe use to allocate complex diagonals.
      include 'ndimsdecl.f'
      include 'partcom.f'   ! for dt
!      parameter (na_m=1000)  ! Instead doing this for now.
      integer ifull(3),iuds(3)
      real u(ifull(1),ifull(2),ifull(3)),q(ifull(1),ifull(2),ifull(3))
      double complex cscratch(0:iuds(2)-1,0:iuds(3)-1,0:iuds(1)+3) ! Reordered!
      integer M,N,O
      complex C(na_m),D(na_m),E(na_m),B(na_m)  ! cgtsl diagonals
      integer*8 planforward,planbackward
      real k2
      integer idebug,ifirst
      logical lreuse,loutward
      data idebug/0/ifirst/0/lreuse/.true./
      data loutward/.true./ ! Outward propagating else log slope match
      save

      loutward=.false.
      M=iuds(1)-2
      N=iuds(2)-2
      O=iuds(3)-2
! Create plans
! The FFTW_UNALIGNED flag makes it safe to use other arrays with same plan.
      if(.not.lreuse.or.ifirst.eq.0)then
         call dfftw_plan_dft_2d(planforward, N,O, cscratch,cscratch,
     &        FFTW_FORWARD, FFTW_ESTIMATE+FFTW_UNALIGNED)
         call dfftw_plan_dft_2d(planbackward, N,O, cscratch,cscratch,
     &        FFTW_BACKWARD, FFTW_ESTIMATE+FFTW_UNALIGNED)
         ifirst=ifirst+1
      endif
! Initialize to zero the additional storage planes
      if(ifirst.eq.0)cscratch(0:N+1,0:O+1,M+2:M+5)=0.

! For each x, 2D-FFT q giving a 3D array in which x
! is still position (not k).
      do i=1,M
         do k=1,O 
            do j=1,N   ! Copy reordered q without the guard cells to scratch.
               cscratch(j,k,i)=-q(i+1,j+1,k+1)
            enddo
         enddo
         call dfftw_execute_dft(planforward, cscratch(1,1,i) ,
     $        cscratch(1,1,i))
      enddo

      scale=N*O
      cscratch(1:N,1:O,1:M)=cscratch(1:N,1:O,1:M)/scale
! For each ky,kz, calculate the effective (piecewise linear) k2 (perp) and 
! solve in the x-direction by tridiagonal elimination. 
      dx=xL/M
      ci=(1/dx)**2
      ri=0. ! Make edge potential zero.
      ri=1. ! Logarithmic slope =-1/dx
      ri=5. ! Use as upper limit of log slope inverse.
      do k=1,O
         kz=k-1
         if(kz.gt.O/2)kz=kz-O
         skz=2.*sin(kz*3.1415926/O)
         do j=1,N
            ky=j-1
            if(ky.gt.N/2)ky=ky-N
            sky=2.*sin(ky*3.1415926/N)
            k2=((sky*N/yL)**2+(skz*O/zL)**2)
            di=-2*ci-k2
            G=0.
            if(k2.gt.0)then 
               G=min(5.,dt/(dx*sqrt(k2)))
               r=min(1/(dx*sqrt(k2)),ri)
            else
               r=ri ! set k2=0 ratio
            endif
            H=(1.-G)/(1.+G)
            do i=1,M
! Create the initial RHS for this ky,kz, which is single-complex 
! Forward Fourier transform of -q
               B(i)= cmplx(cscratch(j,k,i))
! Create the diagonals and subdiagonals
               C(i)= ci
               D(i)= di
               E(i)= ci
            enddo
! Fix the Boundary corners. Non propagating version:
            if(loutward)then
               stable=1.1 ! Stabilizing hack.
! Outward propagating boundary conditions:
               D(1)=stable*((-2.-H)*ci-k2)
               D(M)=stable*((-2.-H)*ci-k2)
! Now we need to fix B at ends by adding FT-u terms from prior timestep
! Requires the extra cscratch storage.
               B(M)=B(M)+cmplx(H*cscratch(j,k,M+5)+cscratch(j,k,M+4))*ci
               B(1)=B(1)+cmplx(H*cscratch(j,k,M+2)+cscratch(j,k,M+3))*ci
            else
! Setting no B modification is consistent if the D value is adjusted
! but otherwise is assuming effectively phiG=0.  

! phi + r dphi/dx=0 requires (pN+pG)/2=-r(pG-pN)/dx (r is positive).
! i.e. pG=(2r/dx-1)/(2r/dx+1)pN. So the resulting last row equation is
! pNm-2pN+(2r/dx-1)/(2r/dx+1)pN=B*dx^2 i.e. with k2 correction: 
!               D(M)=(-2.-k2/ci+(2.*r/dx-1)/(2.*r/dx+1))*ci
!               D(1)=(-2.-k2/ci+(2.*r/dx-1)/(2.*r/dx+1))*ci
! Alternative with uncentered value. pG=pN*r/(1+r), so last row is
! PNm-2pN+pN*r/(1+r) with k2 correction
               D(M)=(-2.-k2/ci+r/(1.+r))*ci
               D(M)=(-2.-k2/ci+r/(1.+r))*ci
            endif
! Install in the temporary Guard FT u-values equal to -q at edge
            cscratch(j,k,0)=cscratch(j,k,1)
            cscratch(j,k,M+1)=cscratch(j,k,M)
! Call linpack complex tridiagonal solver cgtsl
            call CGTSL(M,C,D,E,B,INFO)
            if(INFO.ne.0) stop 'CGTSL Zero diagonal trap'
! Install solution into in-mesh cscratch
            cscratch(j,k,1:M)= B(1:M)
            if(loutward)then
! Turn the temporary guard values into the guard FT u-values
! consistent with the poisson equation at mesh edge.
               cscratch(j,k,0)=cscratch(j,k,0)/ci
     $              +(2.+k2/ci)*cscratch(j,k,1)-cscratch(j,k,2)
               cscratch(j,k,M+1)=cscratch(j,k,M+1)/ci
     $              +(2.+k2/ci)*cscratch(j,k,M)-cscratch(j,k,M-1)
            else 
! Set guard r=-phi/phi'dx= -phiG/(phiG-phiM) [uncentered]
! So phiG*(r+1)=phiM*r, i.e. 
               cscratch(j,k,M+1)=cscratch(j,k,M)*r/(1.+r)
               cscratch(j,k,0)=cscratch(j,k,1)*r/(1.+r)
            endif
! Save the 0,1,M,M+1 planes of the fourier transformed solution
! (guard and edge) for use in the next timestep.
            cscratch(j,k,M+2)=
     $           cscratch(mod(j-2+N,N)+1,mod(k-2+O,O)+1,0)
            cscratch(j,k,M+3)=
     $           cscratch(mod(j-2+N,N)+1,mod(k-2+O,O)+1,1)
            cscratch(j,k,M+4)=
     $           cscratch(mod(j-2+N,N)+1,mod(k-2+O,O)+1,M)
            cscratch(j,k,M+5)=
     $           cscratch(mod(j-2+N,N)+1,mod(k-2+O,O)+1,M+1)
         enddo
      enddo

! Then do a 2D backward transform at each x-position to arrive at a
! 3D spatial array when copied into u.
      do i=0,M+1
         call dfftw_execute_dft(planbackward, cscratch(1,1,i) ,
     $        cscratch(1,1,i))
         do k=0,O+1             ! Copy cscratch to u with periodic guard cells.
            do j=0,N+1  
               u(i+1,j+1,k+1)=
     $              real(cscratch(mod(j-1+N,N)+1,mod(k-1+O,O)+1,i),4)
            enddo
         enddo
      enddo
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE CGTSL (N, C, D, E, B, INFO)
C***BEGIN PROLOGUE  CGTSL
C***PURPOSE  Solve a tridiagonal linear system.
C***LIBRARY   SLATEC (LINPACK)
C***CATEGORY  D2C2A
C***TYPE      COMPLEX (SGTSL-S, DGTSL-D, CGTSL-C)
C***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, SOLVE, TRIDIAGONAL
C***AUTHOR  Dongarra, J., (ANL)
C***DESCRIPTION
C
C     CGTSL given a general tridiagonal matrix and a right hand
C     side will find the solution.
C
C     On Entry
C
C        N       INTEGER
C                is the order of the tridiagonal matrix.
C
C        C       COMPLEX(N)
C                is the subdiagonal of the tridiagonal matrix.
C                C(2) through C(N) should contain the subdiagonal.
C                On output C is destroyed.
C
C        D       COMPLEX(N)
C                is the diagonal of the tridiagonal matrix.
C                On output D is destroyed.
C
C        E       COMPLEX(N)
C                is the superdiagonal of the tridiagonal matrix.
C                E(1) through E(N-1) should contain the superdiagonal.
C                On output E is destroyed.
C
C        B       COMPLEX(N)
C                is the right hand side vector.
C
C     On Return
C
C        B       is the solution vector.
C
C        INFO    INTEGER
C                = 0 normal value.
C                = K if the K-th element of the diagonal becomes
C                    exactly zero.  The subroutine returns when
C                    this is detected.
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   780814  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  CGTSL
      INTEGER N,INFO
      COMPLEX C(*),D(*),E(*),B(*)
C
      INTEGER K,KB,KP1,NM1,NM2
      COMPLEX T
      COMPLEX ZDUM
      REAL CABS1
      CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
C***FIRST EXECUTABLE STATEMENT  CGTSL
         INFO = 0
         C(1) = D(1)
         NM1 = N - 1
         IF (NM1 .LT. 1) GO TO 40
            D(1) = E(1)
            E(1) = (0.0E0,0.0E0)
            E(N) = (0.0E0,0.0E0)
C
            DO 30 K = 1, NM1
               KP1 = K + 1
C
C              FIND THE LARGEST OF THE TWO ROWS
C
               IF (CABS1(C(KP1)) .LT. CABS1(C(K))) GO TO 10
C
C                 INTERCHANGE ROW
C
                  T = C(KP1)
                  C(KP1) = C(K)
                  C(K) = T
                  T = D(KP1)
                  D(KP1) = D(K)
                  D(K) = T
                  T = E(KP1)
                  E(KP1) = E(K)
                  E(K) = T
                  T = B(KP1)
                  B(KP1) = B(K)
                  B(K) = T
   10          CONTINUE
C
C              ZERO ELEMENTS
C
               IF (CABS1(C(K)) .NE. 0.0E0) GO TO 20
                  INFO = K
                  GO TO 100
   20          CONTINUE
               T = -C(KP1)/C(K)
               C(KP1) = D(KP1) + T*D(K)
               D(KP1) = E(KP1) + T*E(K)
               E(KP1) = (0.0E0,0.0E0)
               B(KP1) = B(KP1) + T*B(K)
   30       CONTINUE
   40    CONTINUE
         IF (CABS1(C(N)) .NE. 0.0E0) GO TO 50
            INFO = N
         GO TO 90
   50    CONTINUE
C
C           BACK SOLVE
C
            NM2 = N - 2
            B(N) = B(N)/C(N)
            IF (N .EQ. 1) GO TO 80
               B(NM1) = (B(NM1) - D(NM1)*B(N))/C(NM1)
               IF (NM2 .LT. 1) GO TO 70
               DO 60 KB = 1, NM2
                  K = NM2 - KB + 1
                  B(K) = (B(K) - D(K)*B(K+1) - E(K)*B(K+2))/C(K)
   60          CONTINUE
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
  100 CONTINUE
C
      RETURN
      END
