c Main program to get (not plot) phasespace data from file[s] on command line.
c At the same time get the potential u, and calculate energy \int E^2dx
      implicit none
      include '../src/ndimsdecl.f'
      include '../src/phasecom.f'
      include '../src/partcom.f'
      include '../accis/plotcom.h'
      character*100 phasefilename
      real x(npsbuf),u(npsbuf)!,uave(npsbuf)
      integer j,k,nsmooth,knmax,idk,nxmax,ixstep,nxua,icl,icsw
      parameter (nxmax=512)
      real energy(npsbuf),time(npsbuf),xlen,dlnen(npsbuf),smooth(npsbuf)
      real uarray(npsbuf,nxmax),xuarray(nxmax),worka(npsbuf,nxmax)
      real umin,umax,zclv(2)
      real t,slope,enl,enkn,enu,tknl,tknmax,tknu
      real phirange,phirangeinit
      parameter (phirangeinit=0.5)
      integer i,ii,n,Np,Nave,Nastep,idone,irun
      character*12 nlabel(2)
      logical ldebug
      data nlabel/' !Bn!di!d!@',' !Bn!de!d!@'/
      data idone/0/irun/0/
      ldebug=.false.
      nsmooth=160
      Nave=1
      Nastep=1
      phirange=phirangeinit
      j=0
      call pfset(3)
      do i=1,iargc()
         call getarg(i,phasefilename)
         if(ldebug)write(*,*)phasefilename(1:50)
! Deal with possible switches
         if(phasefilename(1:2).eq.'-N')then
c Set the starting number of plot file writing to be N
            read(phasefilename(3:),*,err=2,end=2)Np
            pfilno=Np
            write(*,*)'Plot file number offset:',Np
            goto 1
 2          write(*,*)'Garbled -N flag (needs number):',phasefilename
            goto 1
         endif
         if(phasefilename(1:2).eq.'-A')then
            read(phasefilename(3:),*,err=3,end=3)Nave
            write(*,*)'Averaging potential over steps:',Nave
            goto 1
 3          write(*,*)'Garbled -a flag (needs number):',phasefilename
            goto 1
         elseif(phasefilename(1:3).eq.'-ns')then
            read(phasefilename(4:),*,err=6,end=6)nsmooth
            write(*,*)'dln(energy)/dt smoothing over range',nsmooth
            goto 1
 6          write(*,*)'Garbled -ns flag (needs number):',phasefilename
            goto 1
         elseif(phasefilename(1:2).eq.'-d')then
            ldebug=.not.ldebug
         elseif(phasefilename(1:2).eq.'-h')then
            goto 4
         elseif(phasefilename(1:2).eq.'-q')then
            irun=1
            call pfset(-3)
         else
! Not a known switch. Interpret as a file name, and read it.
            j=j+1
            if(j.gt.npsbuf)then
               write(*,*)'File number exceeds npsbuf',npsbuf
               goto 5
            endif
            n=npsbuf
            call phaseread(phasefilename,n,x,u,t)
! Store in possibly limited x-range x,u versus time arrays.
            ixstep=1+(n-1)/nxmax
            k=0
            do ii=1,n,ixstep
               k=k+1
               uarray(j,k)=u(ii)
               if(j.eq.1)xuarray(k)=x(ii)
            enddo
            nxua=k
! Find total average energy density at this time.
            energy(j)=0.
            xlen=abs(x(n)-x(1))
            do ii=2,n
               energy(j)=energy(j)+(u(ii)-u(ii-1))**2/(x(ii)-x(ii-1))
            enddo
            energy(j)=energy(j)/xlen  ! Thus energy density
            time(j)=t
         endif
 1       continue
      enddo
      write(*,'(10f8.2)')(xuarray(k),k=1,nxua) ! Check xuarray.
      if(j.eq.0)goto 4
 5    continue  
! Calculate the time derivative of ln(energy)
      do i=2,j  
         dlnen(i)=(alog(energy(i))-alog(energy(i-1)))
     $        /(time(i)-time(i-1))
      enddo
! Smooth dlnen over a range of nsmooth+1 points. (Should be adjustable.)
! And find the peak smoothed value. 
      knmax=1
      do i=1+nsmooth/2,j-nsmooth/2  ! Triangular smoothing
         smooth(i)=0.
         do k=i-nsmooth/2,i+nsmooth/2
            smooth(i)=smooth(i)+dlnen(k)*(1-abs(k-i)/(nsmooth/2.))
         enddo
         smooth(i)=smooth(i)/(nsmooth/2.)
         if(smooth(i).gt.smooth(knmax))knmax=i
      enddo
! Now we have the steepest slope of the smoothed dlnen.
      slope=smooth(knmax)
! Get the ends of the straight line with this slope to draw
      idk=nint(.6*knmax)
      tknmax=time(knmax)
      tknu=time(knmax+idk)
      tknl=time(knmax-idk)
      enkn=energy(knmax)
      enl=enkn*exp(-slope*(tknmax-tknl))
      enu=enkn*exp(-slope*(tknmax-tknu))
      write(*,'(a,f10.5a,f7.2)')'Peak growth rate:',slope
     $ ,'   Growth time:',1/slope
      if(ldebug)then
         write(*,*)'knmax=',knmax,' idk=',idk,' slope=',slope
         write(*,*)'tknl,tknu',tknl,tknu
         write(*,*)'enl ,enu ',enl,enu
         call autoplot(time,smooth,j) ! Check the smoothed log slope.
         call axlabels('time','smoothed dln<E!u2!u>/dt')
         call pltend
      endif
      call multiframe(2,1,1)
! Do the growth plot
      call lautoplot(time,energy,j,.false.,.true.)
      call axlabels('time','field energy density <E!u2!u>')
      call axis2
! Overplot the slope line
      call color(4)
      call polyline([tknl,tknu],[enl,enu],2)
      call color(15)
!      call pltend
! Put contour plot of uarray(t,x) here.
      call pltinit(time(1),time(j),xuarray(1),xuarray(nxua))
      call minmax2(uarray,npsbuf,j,nxua,umin,umax)
      zclv(1)=umin/2
      zclv(2)=umax/2
      icl=2
      icsw=1+16+32
!      call blueredgreenwhite()
      call contourl(uarray,worka,npsbuf,j,nxua,zclv,icl,
     $    time,xuarray,icsw)
      call color(2)
      call axis
      call axis2
      call axlabels('time','x')
      call gradlegend(zclv(1),zclv(2),1.04,0.,1.04,.9,.02,.true.)
      call legendline(1.03,.95,258,'!Af!@')
      call pltend
      call multiframe(0,0,0)
      call exit
 4    continue
      write(*,*)'Obtain the evolution of the field energy from',
     $     ' one dimensional pps files'
      write(*,*)'Usage:  EnergyPhase [Options] file1 [file2 ....]'
      write(*,*)'Options: -A<nnn> average-number',' -N starting-number'
      write(*,*)' -ns<nnn> smoothing range -q no screen plots or halts.'
      write(*,*)' -d toggle debugging'
      end
 
