      integer ntmax,nt,nxmax,nx,nmodes
      parameter (ntmax=2000,nxmax=1024,nmodes=11,npmax=20)
      complex phimodes(ntmax,nxmax,nmodes)
      real rphimodes(ntmax,nxmax,nmodes),cphimodes(ntmax,nxmax,nmodes)
      integer iudsphi(3),ifullphi(3),ixnps(4)
      integer ml1(1),mloc(2),it(npmax)
      real xns(nxmax)
      real time(ntmax),xn(nxmax)
      parameter (nex=1)
      real zp(ntmax,ntmax),yp(ntmax,nmodes),yp2(ntmax,nmodes-nex)
      real xmodenums(nmodes)
      character*20 string
      data xmodenums/0,1,2,3,4,5,6,7,8,9,10/
      data ifullphi/ntmax,nxmax,nmodes/

      call unsavemodes(ntmax,nt,nxmax,nx,nmodes,nmr,phimodes,time,xn
     $     ,iyperr)

      if(iyperr.ne.0) Stop 'No savedmodes.dat file'

! Testing verification of correct recall.
      do i=1,4
         write(*,'(8g10.2)')(phimodes(i,j,2),j=1,4)
      enddo

      clip=3.   
      iclipped=nint(clip/(xn(2)-xn(1))) ! Go clip Debye lengths.
      ix1=max(nx/2+1-iclipped,1)
      ix2=min(nx/2+iclipped,nx)
! Slicing call
      iudsphi(1)=nt
      iudsphi(2)=ix2-ix1+1
      iudsphi(3)=nmodes
      call ixnpcreate(iudsphi(1),iudsphi(2),iudsphi(3)
     $     ,time,xn(ix1),xmodenums,ixnps,xns)
      call setax3chars('time','x','mode')
         
      rphimodes=real(phimodes)
!      cphimodes=sqrt(real(phimodes)**2+imag(phimodes)**2)
      cphimodes=abs(phimodes)
      yp=0.
      do n=ix1,ix2
         yp=yp+cphimodes(:,n,:)
      enddo
! Find the maxima of the oscillations starting at the overall maximum
!      yp2=yp(:,2:nmodes)
      yp2=yp(:,1+nex:nmodes)
      mloc=maxloc(yp2)
      write(*,*)mloc
      mode=mloc(2)  ! The mode gt 0 that grows biggest.
      call autoplot(time,yp2(1,mode),nt)
      call axlabels('time','amplitude')
      call iwrite(mode+nex-1,iwidth,string)
      call legendline(.1,.8,0,'mode '//string)
      call dashset(4)
      call polyline(time,yp(1,mode),nt)
      call iwrite(mode+nex-2,iwidth,string)
      call legendline(.1,.75,0,'mode '//string)
      call dashset(5)
      call polyline(time,yp(1,mode+2),nt)
      call iwrite(mode+nex,iwidth,string)
      call legendline(.1,.85,0,'mode '//string)
      it(1)=mloc(1)
      istep=3
      do i=1,npmax-1
         do j=4,50
!           write(*,*)'i,j',i,j,it(i)-j,yp2(it(i),mode),yp2(it(i)-j
!     $           ,mode)
            if(yp2(it(i)-j,mode).gt.yp2(it(i)-j+istep,mode))exit
         enddo    ! Found adjacent minimum.
         if(.not.yp2(it(i)-j,mode)/yp2(it(i),mode).lt..8)exit ! Too shallow
         ml1=maxloc(yp2(1:it(i)-j,mode))
         it(i+1)=ml1(1)         
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

      if(.true.)then   

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


      end

include 'modesaving.f'
