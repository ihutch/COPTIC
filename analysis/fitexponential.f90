!Optional main
!call testfitexponential
!end
!*******************************************************************
subroutine fitexponential(n,z,A,B,ierr)
! To a complex trace z of length n, use FFT to fit a complex exponential form
! to find A and B such that z(j)~=A*exp(i B*(j-1)).
! Uses fftpack routines, or fftw3 routines.
  integer :: n,ierr
  complex :: z(n),A,B

  integer :: error,inc,k
  integer :: ml1(1)
  real :: j(n)                 ! Automatic arrays.
  complex :: zlocal(n),zinv(n),wk(n) ! Automatic arrays.
  real, allocatable :: WSAVE(:),WORK(:)  !fftpack Allocatable arrays
  complex :: F,G,P,Q,R,S,T,D
  integer :: idebug=1

  j=(/(k,k=0,n-1)/)
  inc=1
! Fourier transform z locally
  zlocal=z
!FFTW
!  call sfftw_plan_dft_1d(iplan,n,zlocal,zlocal,-1,64) ! forward, estimate
!  call sfftw_execute_dft(iplan,zlocal,zlocal)
!  call sfftw_destroy_plan(iplan)
!  zlocal=zlocal/n
!fftpack
  lenc=2*n
  lenwrk=2*n
  lensav=2*n+30
  allocate(wsave(lensav),stat=error)
  allocate(work(lenc),stat=error)
  if (error.ne.0) write(*,*)'****** fitexponential Allocation error',error
  call CFFT1I (n, WSAVE, LENSAV, error)
  call CFFT1F (n,inc,zlocal,n,WSAVE,LENSAV,WORK,LENWRK,error)  
  if(error.ne.0)write(*,*)'***** fitexponential CFFT error',error
!!!!!!!!!!!!!!!!!!!!!!!!!!!! Find the peak and the range nearby
! Find the maximum amplitude frequency position
  ml1=maxloc(abs(zlocal))
  kr=ml1(1)
  zm=abs(zlocal(kr))
  zfac=.6        ! Probably it would be good to autotune this value
  do kw=0,n/4    ! Find kw as the half-width zfac-height of the zlocal peak.
     if(kr-kw.le.1.or.kr+kw.ge.n)exit
     if(abs(zlocal(kr+kw)).lt.zm*zfac.and.abs(zlocal(kr-kw)).lt.zm*zfac)exit
  enddo
! Near this position 1/zlocal ~ F*j+G 
  zinv=1./zlocal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
! The fit minimizes \sum |F*j+G-zinv(j)|^2 which requires 
! P*F+Q*G-R=0 and Q*F+S*G-T=0, where
! P=Sum(j^2), Q=Sum(j), R=Sum(j*zinv(j)), S=Sum(1), T=Sum(zinv(j))
! whose solution is F=-(Q*T-S*R)/D, G=(P*T-Q*R)/D, D=(P*S-Q**2). 
  P=sum(j(kr-kw:kr+kw)**2)
  Q=sum(j(kr-kw:kr+kw))
  R=sum(j(kr-kw:kr+kw)*zinv(kr-kw:kr+kw))
  S=1+kw*2
  T=sum(zinv(kr-kw:kr+kw))
  D=P*S-Q**2
  if(D.eq.0)then
     write(*,*)'fitexponential D=0 Error. kw=',kw
     ierr=1
     return
  endif
  F=-(Q*T-S*R)/D
  G= (P*T-Q*R)/D
  ! Now we have a fit to zinv=F*j+G, where zinv=1/zlocal=1/FFT(z)
  wk=F*j+G
  ! One can show that (to a good approximation)
  B=-2.*3.1415926*G/F/n
  A=-2.*3.1415926*complex(0.,1.)/(exp(complex(0.,1.)*B*n)-1)/F
  ! However, there is ambiguity 2pi in the value of Br, and probably one
  ! should insist that |B|<=pi. So if Br>pi Br=Br-2pi
  if(real(B).gt.3.1415926)B=B-2.*3.1415926
  if(idebug.ne.0)then
     write(*,*)'kw=',kw
     write(*,*)'F,G',F,G
     call autoplot(j,abs(zlocal),n)
     call axlabels('j','|FFT|')
     call color(4)
     call polyline(j(kr-kw:kr+kw),abs(zlocal(kr-kw:kr+kw)),1+2*kw)
     call color(13)
     call polyline(j(kr-kw:kr+kw),abs(1/wk(kr-kw:kr+kw)),1+2*kw)
     call color(15)
     call pltend
!     call autoplot(j,real(zinv),n)
!     call axlabels('j','1/FFT')
!     call polyline(j,imag(zinv),n)
     !  call autoinit(j(kr-kw:kr+kw),real(zinv(kr-kw:kr+kw)),1+2*kw)
     !  call axis
!     call color(12)
!     call polyline(j(kr-kw:kr+kw),real(zinv(kr-kw:kr+kw)),1+2*kw)
!     call polyline(j(kr-kw:kr+kw),imag(zinv(kr-kw:kr+kw)),1+2*kw)
!     call polyline(j(kr-kw:kr+kw),real(wk(kr-kw:kr+kw)),1+2*kw)
!     call polyline(j(kr-kw:kr+kw),imag(wk(kr-kw:kr+kw)),1+2*kw)
!     call pltend
  endif
  ierr=0
end subroutine fitexponential
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Optional main program as a subroutine.
subroutine testfitexponential
integer, parameter :: nz=200
complex :: z(nz),A,B,Bin,zout(nz),AR,BR,zcreal(nz)
real :: zreal(nz),zrealout(nz)

Bin=complex(50./nz,3./nz)
do j=1,nz
   z(j)=exp(complex(0.,1.)*Bin*(j-1)) + .1*(rand()-0.5+complex(0.,rand()-0.5))
   zreal(j)=2*real(z(j))
enddo
zcreal=zreal
call fitexponential(nz,z,A,B,ierr)
write(*,*)'Bin=',Bin
write(*,*)'A=',A,' B=',B
call fitexponential(nz,zcreal,AR,BR,ierr)
do j=1,nz
   zout(j)=A*exp(complex(0.,1.)*B*(j-1))
   zrealout(j)=real(AR*exp(complex(0.,1.)*BR*(j-1)))
enddo
call yautoplot(real(z),nz)
call axlabels('index','real(z)')
call color(5)
call ypolyline(real(zout),nz)
call color(7)
call ypolyline(zrealout,nz)
call color(15)
call pltend
end
