c Code for phasespace accumulation, reading, writing, plotting.
c**********************************************************************
      block data phaseblockdata
      include 'ndimsdecl.f'
      include 'phasecom.f'
      data psfmax/0./ipsftri/0/
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
      integer ispecies,id
      integer i,ierr,ixbin,ivbin,isp
      real v,x,vs,vt
      logical linitedps(nspeciesmax)
      data linitedps/nspeciesmax*.false./

      if(.not.linitedps(ispecies))then
c Ensure the limits etc of the phase space array are set.
         psxmin=xmeshstart(id)
         psxmax=xmeshend(id)
         if(nc(ispecies).ne.0)then
            isp=ispecies
            vs=max(maxval(vsc(1:nc(isp),isp))
     $           ,maxval(vsc(1:nc(isp),isp)))
            vt=maxval(vtc(1:nc(isp),isp))
            psvmax(isp)=sqrt(abs(eoverms(isp)))*(3*vt+vs)
            vs=min(minval(vsc(1:nc(isp),isp))
     $           ,minval(vsc(1:nc(isp),isp)))
            psvmin(isp)=sqrt(abs(eoverms(isp)))*(-3*vt+vs)
         else
            psvmax(ispecies)=3.*sqrt(abs(eoverms(ispecies))
     $           *Ts(ispecies))
         endif
c The centers of the bins in phase space (redundancy negligible).
         do i=1,npsx
            psx(i)=psxmin+(i-0.5)*(psxmax-psxmin)/npsx
         enddo
         do i=1,npsv
            psv(i,ispecies)=psvmin(ispecies)+(i-0.5)*(psvmax(ispecies)
     $           -psvmin(ispecies))/npsv
         enddo
         call fvinfincalc(ispecies)
      endif

      call pszero(ispecies)
c Accumulate
      do i=iicparta(ispecies),iocparta(ispecies)
         if(x_part(iflag,i).ne.0)then
            v=x_part(id+ndims,i)
            x=x_part(id,i)-0.5*x_part(idtp,i)*v
            ixbin=int(.99999*(x-psxmin)/(psxmax-psxmin)*float(npsx)+1)
            ivbin=int(.99999*(v-psvmin(ispecies))/(psvmax(ispecies)
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
      enddo
      call mpiallreducesum(psn(1,ispecies),npsx,ierr)
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
 103  continue
      endif
      close(13)
      return
 101  write(*,*)'Phase read ERROR opening file:',
     $     phasefilename(1:lentrim(phasefilename))
      nu=0
      return
 102  write(*,*)'Reading Error for file',phasefilename
      psn(1,1)=0.
      close(13)
      stop
      end
c***********************************************************************
      subroutine phaseplot(ispecies)
      include 'ndimsdecl.f'
      include 'phasecom.f'
      character*30 string
      real cworka(npsx,npsv),zclv(2)
      integer ifcolor(npsv)

      call pltinit(psxmin,psxmax,psvmin(ispecies),psvmax(ispecies))
      call color(15)
      write(string,'('' f!d'',i1,''!d'')')ispecies
      call jdrwstr(wx2nx(psxmax),wy2ny(psv(40,ispecies))
     $     ,string(1:lentrim(string)),1.)
      call blueredgreenwhite()
      write(string,'(''v!d'',i1,''!d'')')ispecies
      call axlabels('x',string(1:lentrim(string)))
c If unset, set psfmax for less than full range. But better set earlier.
      call minmax2(psfxv(1,1,ispecies),npsx,npsx,npsv,pmin,pmax)
      if(pmax.gt.psfmax*1.1)call psfmaxset(1.1,ispecies)
c Set extrema of coloring range from psfmax.
      zclv(1)=0.
      zclv(2)=psfmax
      icl=2
      icsw=1+16+32+ipsftri
c Using triangular gradients +64 gives too large ps output.
      call contourl(psfxv(1,1,ispecies),cworka,npsx,npsx,npsv,zclv,icl
     $     ,psx,psv(1,ispecies),icsw) 
      call axis()
      call axis2
      call color(ilightgray())
      call gradlegend(zclv(1),zclv(2),.3,.94,.7,.94,.05,.true.)
      call color(15)
      call polyline(psxmax+.3*psxmax*finfofv(:,ispecies),psv(1,ispecies)
     $     ,npsv)
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
      do i=1,npsv
         fi=0.
         vinf=psv(i,ispecies)/sqrt(abs(eoverms(ispecies)))
         do j=1,nc(ispecies)
            fi=fi+dcc(j,ispecies)*exp(-0.5*(vinf-vsc(j,ispecies))**2
     $           /vtc(j,ispecies)**2)/(vtc(j,ispecies)*sqrt(2.
     $           *3.1415926))
         enddo
         finfofv(i,ispecies)=fi
      enddo
      end
