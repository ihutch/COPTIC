      integer ntmax,nt,nxmax,nx,nmodes
      parameter (ntmax=2000,nxmax=1500,nmodes=11,npmax=20)
      complex phimodes(ntmax,nxmax,nmodes)
      real rphimodes(ntmax,nxmax,nmodes),cphimodes(ntmax,nxmax,nmodes)
      integer iudsphi(3),ifullphi(3),ixnps(4)
      integer ml1(1),mloc(2),it(npmax)
      real xns(nxmax)
      real time(ntmax),xn(nxmax),omega(ntmax)
      parameter (nex=1)
      real zp(ntmax,ntmax),yp(ntmax,nmodes),yp2(ntmax,nmodes-nex)
      real ap(ntmax,nmodes),rap(ntmax,nmodes)
      complex cap(ntmax,nmodes)
      real xmodenums(nmodes)
      character*20 string


      INTEGER NMAX, INC, LENSAV, LENWRK, IER
      PARAMETER (NMAX=2000,INC=1,LENSAV=(2*NMAX+30),LENWRK=2*NMAX)
      REAL       WSAVE(LENSAV), WORK(LENWRK)


      data xmodenums/0,1,2,3,4,5,6,7,8,9,10/
      data ifullphi/ntmax,nxmax,nmodes/

      call unsavemodes(ntmax,nt,nxmax,nx,nmodes,nmr,phimodes,time,xn
     $     ,iyperr)

! Use this for rounded over nonlinear cases:
      nt=.9*nt

      if(iyperr.ne.0) Stop 'No savedmodes.dat file'

! Testing verification of correct recall.
      do i=1,4
         write(*,'(8g10.2)')(phimodes(i,j,2),j=1,4)
      enddo

      clip=5.
      iclipped=nint(clip/(xn(2)-xn(1))) ! Go clip Debye lengths.
      ix1=max(nx/2+1-iclipped,1)
      ix2=min(nx/2+iclipped,nx)
! Slicing call needs:
      iudsphi(1)=nt
      iudsphi(2)=ix2-ix1+1
      iudsphi(3)=nmodes
      call ixnpcreate(iudsphi(1),iudsphi(2),iudsphi(3)
     $     ,time,xn(ix1),xmodenums,ixnps,xns)
      call setax3chars('time','x','mode')
      rphimodes=0.
      rphimodes(1:nt,:,:)=real(phimodes(1:nt,:,:))
!      cphimodes=abs(phimodes)
      cphimodes=abs(rphimodes)
      yp=0.
      cap=0.
      ap=0.
      do n=ix1,ix2
         yp=yp+abs(xn(n))*cphimodes(:,n,:) ! Makes sense only for absolute.
         ap=ap+xn(n)*rphimodes(:,n,:)
         cap=cap+xn(n)*phimodes(:,n,:)
      enddo
! yp is now the absolute amplitude and ap the amplitude of real part.

! Find the maxima of the oscillations starting at the overall maximum
!      yp2=yp(:,2:nmodes)
      yp2=yp(:,1+nex:nmodes)
      mloc=maxloc(yp2)
      write(*,*)'mloc',mloc
      mode=mloc(2)  ! The mode gt 0 that grows biggest.
      call autoplot(time,yp(1,mode+nex),nt)
      call axlabels('time','amplitude')
      call iwrite(mode+nex-1,iwidth,string)
      call legendline(.1,.8,0,'mode '//string)
      call dashset(4)
      call polyline(time,yp(1,mode+nex-1),nt)
      call iwrite(mode+nex-2,iwidth,string)
      call legendline(.1,.75,0,'mode '//string)
      call dashset(5)
      call polyline(time,yp(1,mode+nex+1),nt)
      call iwrite(mode+nex,iwidth,string)
      call legendline(.1,.85,0,'mode '//string)
      it(1)=mloc(1)
      istep=nint(12./(time(2)-time(1)))
      istep=1
      write(*,*)'istep,dt=',istep,time(2)-time(1)
      do i=1,npmax-1 ! The peak number counting back from max
         do j=4,50   ! Step difference back from current peak.
           write(*,'(a,3i4,2i5,2f8.4)')'i,j=',i,j,it(i),
     $           it(i)-j,istep,yp2(it(i),mode),yp2(it(i)-j,mode)
           if(yp2(it(i)-j,mode).gt.yp2(it(i)-j+istep,mode))exit
         enddo    ! Found adjacent minimum.
         if(.not.yp2(it(i)-j,mode)/yp2(it(i),mode).lt..8)exit ! Too shallow
         ml1=maxloc(yp2(1:it(i)-j,mode))
         it(i+1)=ml1(1)         
         if(it(i+1).lt.50)exit
! Now it(i+1) is the position of the next peak down. 
      enddo
      write(*,*)'mode',mode
      write(*,*)'imin it   time      dt   maximum    gamma   1/gamma'
     $     ,'   omega'
      do j=1,i
! exponential growth rate gamma=log(y(n+1)/y(n))/dt
         dt=-time(it(min(j+1,i)))+time(it(j))
         growthrate=-alog(yp2(it(min(j+1,i)),mode)/yp2(it(j),mode))/dt
         write(*,'(2i4,2f8.1,4f9.4)')j,it(j),time(it(j)),dt,yp2(it(j)
     $        ,mode),growthrate,1/growthrate,3.14159/dt
         call polymark(time(it(j)),yp2(it(j),mode),1,1)
      enddo
      if(i-1.gt.3)then
         period=2.*(time(it(2))-time(it(i-1)))/((i-1)-(2))
         gamma=-alog(yp2(it(3),mode)/yp2(it(i-1),mode))
     $        /(time(it(i-1))-time(it(3)))
         write(*,*)'period=',period,' gamma=',gamma,1./gamma
      endif
      call pltend
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(.false.)then   
      write(*,*)'rphimodes(i,j,2)'
      do i=1,4
         write(*,'(6g10.2)')(rphimodes(i,j,2),j=1,4)
      enddo
      call sliceGweb(ifullphi,iudsphi,rphimodes(1,ix1,1),ntmax,zp,
     $              ixnps,xns,3+64,'Amplitude  ',dum,dum)   
!      call autoplot(time,cphimodes(1,nx/2+2,mode),nt)
!      call polyline(time,rphimodes(1,nx/2+2,mode),nt)
!      call pltend
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(.true.)then  ! Not debugged properly
! Use fft to find frequency of different modes
      do nm=1,nmodes
         call CFFT1I (nt, WSAVE, LENSAV, IER)
         call CFFT1F (nt,INC,cap(:,nm),nt,WSAVE,LENSAV,WORK,LENWRK,IER)
      enddo
c B is a straight sum, F is sum divided by N.
c      call CFFT1B (nt,INC,cap,nt,WSAVE,LENSAV,WORK,LENWRK,IER)
      if(ier.ne.0)write(*,*)'CFFT1F error',ier

!      do nm=1,nmodes
!         call realtocomplexfft(nt,ap(:,nm),cap(:,nm))
!      enddo
!      write(*,'(11f7.3)')((real(cap(i,j)),j=1,11),i=1,nt)
      rap=abs(cap)
! Now rap is abs(FT) of the time dependent modes. 
! Figure out the frequency range/array. 
      dt=time(2)-time(1)                          !=T/N
      domega=2.*3.1415926/(time(nt)-time(1))      !=2\pi/T
      ntop=nt/6
      do i=1,ntop
         omega(i)=i*domega
      enddo
 
      do nm=mode+nex,mode+nex
      write(string,'(a,f3.0)')'mode=',xmodenums(nm)
         call dashset(0)
         call autoplot(time,ap(1,nm),nt)
         call axlabels('time','amplitude')
         call boxtitle(string)
         call pltend
      ioff=1
      ml1=maxloc(rap(1:ntop,nm))
      write(*,*)'i=',ml1(1),' omega=',omega(ml1(1)),
     $     ' period=',2*3.14159/omega(ml1(1))
      call autoinit(omega(ioff),rap(ioff,nm),ntop)
      call axis
      call axlabels('omega','abs(mode)')
      call boxtitle(string)
      call dashset(0)
!      call smoothline(omega(ioff),rap(ioff,nm),ntop-ioff,3)
!      call color(2)
      call polyline(omega(ioff),rap(ioff,nm),ntop-ioff)
      call polymark(omega(ml1(1)),rap(ml1(1),nm),1,1)
      call pltend
      call autoplot(omega,imag(cap(ml1(1),nm)/cap(:,nm)),ntop)
      call polymark(omega,imag(cap(ml1(1),nm)/cap(:,nm)),ntop,1)
      call color(2)
      call polyline(omega,real(cap(ml1(1),nm)/cap(:,nm)),ntop)
      call pltend

      enddo

      endif


      end

include 'modesaving.f'

!******************************************************************
      subroutine realtocomplexfft(n,rarray,cfft)
! Calculate the complex fft cfft of a real array rarray of length N.
! Library fftpack5.1 is used. Maximum length of array nmax parameter.
      integer n
      real rarray(n)
      complex cfft(n)

      INTEGER NMAX, INC, LENC, LENSAV, LENWRK, IER
      PARAMETER (NMAX=2000,INC=1,LENSAV=(2*NMAX+30),LENWRK=2*NMAX)
      REAL       WSAVE(LENSAV), WORK(LENWRK)

      if(n.gt.nmax)then
         write(*,*)'realtocomplexfft array length too great',n
         return
      endif
      LENC=n
      do j=1,n
         cfft(j)=rarray(j)
      enddo
      call CFFT1I (n, WSAVE, LENSAV, IER)
      call CFFT1F (n,INC,cfft,LENC,WSAVE,LENSAV,WORK,LENWRK,IER)
c B is a straight sum, F is sum divided by N.
c      call CFFT1B (n,INC,cfft,LENC,WSAVE,LENSAV,WORK,LENWRK,IER)
      if(ier.ne.0)then
         write(*,*)'CFFT1F error',ier
      endif
c$$$ fftpack5.1 Documentation:
c$$$ Output Arguments
c$$$ 
c$$$ C       For index J*INC+1 where J=0,...,N-1 (that is, for the Jth 
c$$$         element of the sequence),
c$$$ 
c$$$            C(J*INC+1) = 
c$$$            N-1
c$$$            SUM C(K*INC+1)*EXP(-I*J*K*2*PI/N)
c$$$            K=0
c Actually by experiment, forward transform is divided by N, backward not.
c$$$         where I=SQRT(-1).
c$$$ IER     =  0 successful exit
c$$$         =  1 input parameter LENC   not big enough
c$$$         =  2 input parameter LENSAV not big enough
c$$$         =  3 input parameter LENWRK not big enough
c$$$         = 20 input error returned by lower level routine
      end
