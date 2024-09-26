C MAIN9999 program to get phasespace data from file[s] on command line.
c At the same time get the potential u, and calculate energy \int E^2dx
c Then contour plot phi(x,t), f(v,t) with energy.
      implicit none
      include '../src/ndimsdecl.f'
      include '../src/phasecom.f'
      include '../src/partcom.f'
      include '../src/fvcom.f'
      include '../accis/plotcom.h'
      character*100 phasefilename
      real x(npsbuf),u(npsbuf)!,uave(npsbuf)
      integer j,k,nsmooth,knmax,idk,nxmax,ixstep,nxua,icl,icsw,jx,iv
      parameter (nxmax=512)
      real energy(npsbuf),time(npsbuf),xlen,dlnen(npsbuf),smooth(npsbuf)
      real uarray(npsbuf,nxmax),xuarray(nxmax),worka(npsbuf,nxmax)
      real faveofv(npsbuf,npsv),varrayoft(npsbuf,npsv),tarr(npsbuf,npsv)
      real umin,umax,zclv(2),pmin,pmax,emin,emax
!      real tbar,tot,v2bar,vbar
      integer ispecies,iwidth
      real t,slope,enl,enkn,enu,tknl,tknmax,tknu
      real phirange,phirangeinit
      real vrange,vrange1
      parameter (phirangeinit=0.5)
      integer i,ii,n,Np,Nave,Nastep,idone,irun,lentrim
      character*12 nlabel(2)
      character*20 string
      character*100 title
      logical ldebug,lrgb
      real wx2nx,wy2ny
      data nlabel/' !Bn!di!d!@',' !Bn!de!d!@'/
      data idone/0/irun/0/
      ldebug=.false.
      nsmooth=160
      Nave=1
      Nastep=1
      phirange=phirangeinit
      j=0
      ispecies=1
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
            goto 1
         elseif(phasefilename(1:2).eq.'-t')then
            read(phasefilename(3:),*)title
            goto 1
         elseif(phasefilename(1:2).eq.'-c')then
            lrgb=.true.
            goto 1
         else
! Not a known switch. Interpret as a file name, and read it.
            j=j+1
            if(j.gt.npsbuf)then
               write(*,*)'File number exceeds npsbuf',npsbuf
               goto 5
            endif
            n=npsbuf
            call phaseread(phasefilename,n,x,u,t)
            if(n.le.0)then
               write(*,*)'Read failed for ', phasefilename
               goto 1
            endif
! Store in possibly limited x-range x,u versus time arrays.
            ixstep=1+(n-1)/nxmax
            k=0
            do ii=1,n,ixstep
               k=k+1
               uarray(j,k)=u(ii)
               if(j.eq.1)xuarray(k)=x(ii)
            enddo
!            if(j.eq.1)write(*,'(10f8.3)')x(1:n)
            nxua=k
! Find total average energy density at this time.
            energy(j)=0.
!            write(*,*)'n,npsbuf',n,npsbuf
            xlen=abs(x(n)-x(1))
            do ii=2,n
               energy(j)=energy(j)+(u(ii)-u(ii-1))**2/(x(ii)-x(ii-1))
            enddo
            energy(j)=energy(j)/xlen  ! Thus energy density
            time(j)=t
         endif
! Integrate wrt x to get fave as a function of t and v.
         do iv=1,npsv
            faveofv(j,iv)=0.
            varrayoft(j,iv)=psv(iv,ispecies)
            tarr(j,iv)=time(j)
            do jx=1,npsx
               faveofv(j,iv)=faveofv(j,iv)+psfxv(jx,iv,ispecies)
            enddo
         enddo
! Rescale?
         if(j.eq.1)vrange1=varrayoft(j,npsv)-varrayoft(j,1)
         vrange=varrayoft(j,npsv)-varrayoft(j,1)
         faveofv(j,:)=(faveofv(j,:)*vrange1/vrange)
!         faveofv(j,:)=alog10(faveofv(j,:)*vrange1/vrange+1000.)
 1       continue
      enddo
! End of command line arguments 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      write(*,'(10f8.3)')time(1:j)
!      write(*,'(10f8.2)')(xuarray(k),k=1,nxua) ! Check xuarray.
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
      if(nsmooth.gt.j)nsmooth=j-1
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
      write(*,'(aa,f10.5a,f7.2)')title(1:lentrim(title))
     $     ,'  Peak growth rate:',slope,'   Growth time:',1/slope
      if(ldebug)then
         write(*,*)'knmax=',knmax,' idk=',idk,' slope=',slope
         write(*,*)'tknl,tknu',tknl,tknu
         write(*,*)'enl ,enu ',enl,enu
         call autoplot(time,smooth,j) ! Check the smoothed log slope.
         call axlabels('time','smoothed dln<E!u2!u>/dt')
         call pltend
      endif


      call multiframe(3,1,1)
! Do the growth plot
      call minmax(energy,j,emin,emax)
      call dcharsize(0.02,0.02)
      call pltinit(time(1),time(j),emin,emax)
      emin=10.**(nint(log10(emin)-0.49999))
      emax=10.**(nint(log10(emax)+0.49999))      
      call scalewn(time(1),time(j),emin,emax,.false.,.true.)
      call axis
      call polyline(time,energy,j)
      call axlabels('time','field energy density <E!u2!u>')
      call axis2
      call legendline(.1,1.1,258,title(1:lentrim(title)))
      if(emax/emin.gt.20.)then
! Overplot the slope line
         call winset(.true.)
         call color(4)
         call polyline([tknl,tknu],[enl,enu],2)
         call fwrite(1/slope,iwidth,1,string)
         call legendline(.05,.92,258,'Growth time '//string(1:iwidth))
         call color(15)
      endif
! Put contour plot of uarray(t,x) here.
      call pltinit(time(1),time(j),xuarray(1),xuarray(nxua))
      call minmax2(uarray,npsbuf,j,nxua,umin,umax)
      zclv(1)=umin/2
      zclv(2)=umax/2
! Contour uarray^1/3 to enhance sensitivity near zero. 
      zclv(2)=umax**.3333
      zclv(1)=-abs(umin)**.3333
      icl=2
      icsw=1+16+32
      if(lrgb)call blueredgreenwhite()
!      call contourl(uarray,worka,npsbuf,j
      call contourl(sign(abs(uarray)**.3333,uarray),worka,npsbuf,j
     $     ,nxua,zclv,icl,time,xuarray,icsw)
      call color(2)
      call axis
      call axis2
      call axlabels('time','x')
      call gradlegend(zclv(1),zclv(2),1.04,0.06,1.04,.9,.02,.true.)
      call legendline(1.02,.97,258,'!Af!@!u1/3!u')
      call color(10)
      call jdrwstr(wx2nx(time(j)),wy2ny(0.)-.01,'-c!ds0!d',-1.5)
      call winset(.true.)
      call polyline([time(j)-xuarray(nxua)*sqrt(1836.),time(j)]
     $     ,[xuarray(nxua),0.],2)
      call color(6)
      call polyline([time(j)-xuarray(nxua),time(j)]
     $     ,[xuarray(nxua),0.],2)
      call winset(.false.)
      call jdrwstr(wx2nx(time(j)-xuarray(nxua)),wy2ny(xuarray(nxua))+.02
     $     ,'-v!dte!d',-0.5)

!      call pltend
! Plot contours of the evolution of the ispecies distribution function.
      call blueredgreenwhite()
      call pltinit(time(1),time(j),psv(1,ispecies),psv(npsv,ispecies))
      call minmax2(faveofv,npsbuf,j,npsv,pmin,pmax)
      write(*,*)'pmin,pmax',pmin,pmax
      zclv(1)=pmin
      zclv(2)=pmax
      icl=2
      icsw=2+16+32
      call contourl(faveofv,worka,npsbuf,j,npsv,zclv,icl,
     $    tarr,varrayoft,icsw)
      call color(2)
      call axis
      call axis2
      call axlabels('time','v')
      call gradlegend(zclv(1),zclv(2),1.04,0.,1.04,.75,.02,.true.)
      call legendline(1.01,.92,258,'f!de!d(v)dv')
      call pltend
!      call multiframe(0,0,0)


      call exit
 4    continue
      write(*,*)'Obtain the evolution of the field energy from',
     $     ' one dimensional pps files'
      write(*,*)'and plot t,x contours of phi, t,v contours of f' 
      write(*,*)'Usage:  EnergyPhase [Options] file1 [file2 ....]'
      write(*,*)'Options: -A<nnn> average-number',' -N starting-number'
      write(*,*)' -ns<nnn> smoothing range. -q no screen plot or halt.'
      write(*,*)' -d toggle debugging. -t<Title> set title. -c color'
      end
 
