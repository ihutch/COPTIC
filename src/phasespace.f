c Code for phasespace accumulation, reading, writing, plotting.
c**********************************************************************
      block data phaseblockdata
      include 'ndimsdecl.f'
      include 'phasecom.f'
      data psfmax/nspeciesmax*0./ipsftri/0/
      end
c**********************************************************************
      subroutine pszero(ispecies)
      integer ispecies
      include 'ndimsdecl.f'
      include 'phasecom.f'
      do i=1,npsx
         do j=1,npsv
            psfxv(i,j,ispecies)=0.
         enddo
      enddo
      end
c**********************************************************************
      subroutine psaccum(ispecies,id)
      implicit none
c Accumulate all particles of species ispecies
c into phase-space x,v for dimension id
c Summed over the other dimensions. 
c The bins are uniform from psvmin to psvmax and xmeshstart to end.
      include 'ndimsdecl.f'
      include 'phasecom.f'
      include 'meshcom.f'
      include 'partcom.f'
      include 'plascom.f'
      include 'fvcom.f'
      include 'myidcom.f'
      integer ispecies,id
      integer i,ierr,ixbin,ivbin,isp
      real v,x,vs,vt,sv
      logical linitedps(nspeciesmax)
      real vplim(2,nspeciesmax)
      data linitedps/nspeciesmax*.false./
      save vplim,linitedps

!      write(*,*)ispecies,linitedps(ispecies)
      if(.not.linitedps(ispecies))then
c Ensure the limits etc of the phase space array are set.
         psxmin=xmeshstart(id)
         psxmax=xmeshend(id)
         vplim(:,ispecies)=-1.e5
         if(nc(ispecies).ne.0)then
            isp=ispecies
            vs=max(maxval(vsc(1:nc(isp),isp))
     $           ,maxval(vsc(1:nc(isp),isp)))
            vt=maxval(vtc(1:nc(isp),isp))
            psvmax(isp)=sqrt(abs(eoverms(isp)))*(3*vt+vs)
            vs=min(minval(vsc(1:nc(isp),isp))
     $           ,minval(vsc(1:nc(isp),isp)))
            psvmin(isp)=sqrt(abs(eoverms(isp)))*(-3*vt+vs)
         elseif(.false.)then ! problematic.
! Assumed unshifted maxwellian when not using -fp argument
            psvmax(ispecies)=3.*sqrt(abs(eoverms(ispecies))
     $           *Ts(ispecies))
            psvmin=-psvmax
         endif
      endif
!      write(*,*)'psaccum',ispecies,psvmin(1:2),psvmin(1:2)


 1    continue
      if(.not.linitedps(ispecies))then
c The centers of the bins in phase space (redundancy negligible).
         do i=1,npsx
            psx(i)=psxmin+(i-0.5)*(psxmax-psxmin)/npsx
         enddo
         do i=1,npsv
            psv(i,ispecies)=psvmin(ispecies)+(i-0.5)*(psvmax(ispecies)
     $           -psvmin(ispecies))/npsv
         enddo
         call fvinfincalc(ispecies)
         linitedps(ispecies)=.true.
      endif

      sv=1.5
      call pszero(ispecies)
c Accumulate
      do i=iicparta(ispecies),iocparta(ispecies)
         if(x_part(iflag,i).ne.0)then
            v=x_part(id+ndims,i)
! Auto-upscaling information
            if(v.gt. sv*vplim(1,ispecies))vplim(1,ispecies)= v/sv
            if(v.lt.-sv*vplim(2,ispecies))vplim(2,ispecies)=-v/sv
            x=x_part(id,i)-0.5*x_part(idtp,i)*v
            ixbin=ceiling(.99999*(x-psxmin)/(psxmax-psxmin)*float(npsx)
     $           +1)
            ivbin=floor(.99999*(v-psvmin(ispecies))/(psvmax(ispecies)
     $           -psvmin(ispecies))*float(npsv)+1)
c Wrap periodically the x-position bins in case of exit.
            if(ixbin.lt.1)ixbin=npsx+ixbin
            if(ixbin.gt.npsx)ixbin=ixbin-npsx
            if(ivbin.gt.0 .and. ivbin.le.npsv)then
c Not for velocity beyond the vrange or x beyond x-range.
               psfxv(ixbin,ivbin,ispecies)=psfxv(ixbin,ivbin,ispecies)+1
            endif
         endif
      enddo
c All reduce to sum the distributions from all processes.
      call mpiallreducesum(psfxv(1,1,ispecies),npsx*npsv,ierr)
! Upscaling 
      call mpiallreducemax(vplim(1,ispecies),2,ierr)
      if(vplim(1,ispecies).gt.psvmax(ispecies))then
         if(myid.eq.0)write(*,'(a,i2,2f10.5)')' Rescale vmax',ispecies
     $        ,vplim(1,ispecies),psvmax(ispecies)
         psvmax(ispecies)=vplim(1,ispecies)+(psvmax(ispecies)
     $        -psvmin(ispecies))*.1
         linitedps(ispecies)=.false.
         goto 1
      endif
      if(vplim(2,ispecies).gt.-psvmin(ispecies))then
         if(myid.eq.0)write(*,'(a,i2,2f10.5)')' Rescale vmin',ispecies,
     $        -vplim(2,ispecies),psvmin(ispecies)
         psvmin(ispecies)=-vplim(2,ispecies)-(psvmax(ispecies)
     $        -psvmin(ispecies))*.1
         linitedps(ispecies)=.false.
         goto 1
      endif
!      write(*,*)'psaccumfinal',ispecies,psvmax(1:2),psvmin(1:2)
!     $     ,linitedps
      end
c***********************************************************************
      subroutine psnaccum(ispecies,id)
      implicit none
c Accumulate the densities of ispecies into phasespace x-bins psn
      include 'ndimsdecl.f'
      include 'phasecom.f'
      include 'meshcom.f'
      include 'partcom.f'
      include 'plascom.f'
      integer ispecies,id
      integer i,ixbin,ierr
      real x
      psn(:,ispecies)=0. ! zero density array. Then check that psx params set:
      if(psxmin.ne.xmeshstart(id))Stop 'psnaccum called before psaccum'
      do i=iicparta(ispecies),iocparta(ispecies)
         x=x_part(id,i)
         ixbin=int(.99999*(x-psxmin)/(psxmax-psxmin)*float(npsx)+1)
         psn(ixbin,ispecies)=psn(ixbin,ispecies)+1.
         psvave(ixbin,ispecies)=
     $        psvave(ixbin,ispecies)+x_part(id+ndims,i)
      enddo
      call mpiallreducesum(psn(1,ispecies),npsx,ierr)
      call mpiallreducesum(psvave(1,ispecies),npsx,ierr)
      end
c***********************************************************************
      subroutine phasewrite(phasefilename,nu,x,u,t)
c Write file with phasespace data plus u(x) if length nu != 0.
      character*(*) phasefilename
      integer nu
      real x(nu),u(nu),t
      include 'ndimsdecl.f'
      include 'phasecom.f'
      include 'partcom.f' ! For nspecies
      include 'fvcom.f'   ! For nc

      open(12,file=phasefilename,status='unknown',form='unformatted',err
     $     =101)
      write(12)npsx,npsv,nu,nspecies
      write(12)psfxv(:,:,1:nspecies),psvmax,psvmin,psxmax,psxmin,psx,psv
      write(12)psn(:,1:nspecies) ! Save the density as fn of x
      if(nu.ne.0)write(12)(x(i),u(i),i=1,nu),t
      do isp=1,nspecies
         if(nc(isp).ne.0)then
            write(12)nc(isp)
            write(12)(finfofv(i,isp),i=1,npsv)
         endif
      enddo
! Addition 3 Oct 2024 of average velocity as a fn of position.
      write(12)psvave(:,1:nspecies)
      close(12)
      return
 101  write(*,*)'Error opening file:',
     $     phasefilename(1:lentrim(phasefilename))
      end
c***********************************************************************
      subroutine phaseread(phasefilename,nu,x,u,t)
      character*(*) phasefilename
      integer nu
      real x(nu),u(nu)
      include 'ndimsdecl.f'
      include 'phasecom.f'
      include 'partcom.f'
      include 'fvcom.f'   ! For nc
      ipsversion=0
      nuin=nu
      open(13,file=phasefilename,status='old',form='unformatted',err
     $     =101)
      read(13)npsxf,npsvf,nu,nspecies
      if(npsx.ne.npsxf.or.npsv.ne.npsvf)then
         write(*,'(a,2i4,a,2i4)')
     $        'Incorrect phase-space array length allocation npsx,npsv',
     $        npsx,npsv,' file requires',npsxf,npsvf
         stop
      else
         read(13,err=102,end=102)psfxv(:,:,1:nspecies),psvmax,psvmin
     $        ,psxmax,psxmin,psx,psv
         read(13,err=102,end=102)psn(:,1:nspecies)
         if(nu.gt.nuin)stop 'phaseread file nu length too great'
         if(nu.ne.0)read(13)(x(i),u(i),i=1,nu),t
         do isp=1,nspecies
            read(13,end=103)nc(isp)
            if(nc(isp).ne.0)read(13)(finfofv(i,isp),i=1,npsv)
         enddo
! Incremental version expansion.
         read(13,err=104,end=104)psvave(:,1:nspecies)
         ipsversion=1
 103  continue
      endif
      close(13)
      return
 101  write(*,*)'Phase read ERROR opening file: ',
     $     phasefilename(1:lentrim(phasefilename))
      nu=0
      return
 102  write(*,*)'Fatal Read Error for file ',phasefilename
      psn(1,1)=0.
      close(13)
      stop
 104  continue
      write(*,'(a,i3,2a)')'PS-Version',ipsversion,' No psvave in file '
     $     ,phasefilename(1:lentrim(phasefilename))
      psvave=0.
      close(13)
      end
c***********************************************************************
      subroutine phaseplot(ispecies)
      include 'ndimsdecl.f'
      include 'phasecom.f'
      character*30 string
      real faveofv(npsv)
      integer ifcolor(npsv)
      logical logspec

!      write(*,*)'phaseplot psvmin,psvmax',ispecies,psvmin(ispecies)
!     $     ,psvmax(ispecies)
      call pltinit(psxmin,psxmax,psvmin(ispecies),psvmax(ispecies))
      call color(15)
      write(string,'('' f!d'',i1,''!d'')')ispecies
      call jdrwstr(wx2nx(psxmax),wy2ny(psv(40,ispecies))
     $     ,string(1:lentrim(string)),1.)
      call blueredgreenwhite()
      write(string,'(''v!d'',i1,''!d'')')ispecies
      call axlabels('x',string(1:lentrim(string)))
! Set adjust f-scale if necessary
      call minmax2(psfxv(1,1,ispecies),npsx,npsx,npsv,pmin,pmax)
      if(pmax.lt.psfmax(ispecies)*0.9)then
         psfmax(ispecies)=pmax
      elseif(pmax.gt.psfmax(ispecies)*1.03)then
         psfmax(ispecies)=pmax
      endif
      logspec=.false.
      if(ispecies.eq.ilogspec)logspec=.true.
      call logphasecont(ispecies,logspec)
      call color(15)
      if(ipsversion.ge.1)call polyline(psx,psvave(1,ispecies),npsx)
      call polyline(psxmax+.3*(psxmax-psxmin)*finfofv(:,ispecies)
     $        /finfmax(ispecies),psv(1,ispecies),npsv)
! Integrate wrt x to get fave as a function of v.
      if(lsideplot)then
      fapeak=0.
      do i=1,npsv
         faveofv(i)=0.
         do j=1,npsx
            faveofv(i)=faveofv(i)+psfxv(j,i,ispecies)
         enddo
         if(faveofv(i).gt.fapeak)fapeak=faveofv(i)
      enddo
      vbar=0.
      tot=0.
      do i=1,npsv
         faveofv(i)=faveofv(i)/fapeak
         vbar=vbar+psv(i,ispecies)*faveofv(i)
         v2bar=v2bar+psv(i,ispecies)**2*faveofv(i)
         tot=tot+faveofv(i)
      enddo 
      vbar=vbar/tot
      v2bar=v2bar/tot
      tbar=v2bar-vbar**2
! And plot it.
      call color(4)
      call polyline(psxmax+.3*(psxmax-psxmin)*faveofv,psv(1,ispecies)
     $     ,npsv)
!      write(*,*)'vbar',vbar
      call polyline([psxmax,psxmax*1.04],[vbar,vbar],2)
      call drcstr('!pv!q!o-!o')
      call color(15)
      endif
      if(.false.)then
! The hard part is scaling the colors (cfac) to what is being used in 
! the histogram plot. The histogram just plots the number of particles
! whereas the side plot is the f(v), normalized to unity integral. 
! Not yet figured out.
         ifcolor=nint(finfofv(:,ispecies)*cfac)
         call polycolorline(psxmax+.3*psxmax*finfofv(:,ispecies),psv(1
     $        ,ispecies),npsv,ifcolor)
      endif
c If needed, do pltend externally.
      end
c**********************************************************************
      subroutine psfmaxset(fac,ispecies)
c Set the value of psfmax to fac times the maximum in the current array.
      real fac
      include 'ndimsdecl.f'
      include 'phasecom.f'
      call minmax2(psfxv(1,1,ispecies),npsx,npsx,npsv,pmin,pmax)
      psfmax=pmax*fac
      end
c**********************************************************************
      subroutine psftri
c Toggle the phase contour triangular gradients
      include 'ndimsdecl.f'
      include 'phasecom.f'
      if(ipsftri.eq.0)then
         ipsftri=64
      else
         ipsftri=0
      endif
      end
c*********************************************************************
      subroutine fvinfincalc(ispecies)
      include 'ndimsdecl.f'
      include 'fvcom.f'
      include 'phasecom.f'
      include 'plascom.f'
c      write(*,'(a,i2,9f10.4)')'species,dcc,vsc,vtc',ispecies ,(dcc(i
c     $     ,ispecies),vsc(i,ispecies),vtc(i,ispecies),i=1,nc(ispecies))
      finfmax(ispecies)=0.
      do i=1,npsv
         fi=0.
         vinf=psv(i,ispecies)/sqrt(abs(eoverms(ispecies)))
         do j=1,nc(ispecies)
            fi=fi+dcc(j,ispecies)*exp(-0.5*(vinf-vsc(j,ispecies))**2
     $           /vtc(j,ispecies)**2)/(vtc(j,ispecies)*sqrt(2.
     $           *3.1415926))
            if(fi.gt.finfmax(ispecies))finfmax(ispecies)=fi
         enddo
         finfofv(i,ispecies)=fi
      enddo
      end
c*********************************************************************
      subroutine logphasecont(ispecies,llog)
      integer ispecies
      logical llog
      include 'ndimsdecl.f'
      include 'phasecom.f'
      real cworka(npsx,npsv),zclv(2)
      icl=2
      icsw=1+16+32
      if(.not.llog)then
c Set extrema of coloring range from psfmax.
         zclv(1)=0.
         zclv(2)=psfmax(ispecies)
         call contourl(psfxv(1,1,ispecies),cworka,npsx,npsx,npsv,zclv
     $        ,icl,psx,psv(1,ispecies),icsw)
         call color(ilightgray())
      else
         zclv(1)=alog10(psfmax(ispecies)/1.e3)
         zclv(2)=alog10(psfmax(ispecies))
         call contourl(alog10(psfxv(:,:,ispecies)+.5),cworka,npsx,npsx
     $        ,npsv,zclv,icl,psx,psv(1,ispecies),icsw)
         call color(ilightgray())
         call legendline(.17,.94,258,'log!d10!d(!Bf!@)')
      endif
      call gradlegend(zclv(1),zclv(2),.3,.94,.7,.94,.05,.true.)
      call color(15)
      call axis()
      call axis2

      end
