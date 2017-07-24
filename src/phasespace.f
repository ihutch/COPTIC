c Code for phasespace accumulation, reading, writing, plotting.
c**********************************************************************
      block data phaseblockdata
      include 'phasecom.f'
      data psfmax/0./
      end
c**********************************************************************
      subroutine pszero
      include 'phasecom.f'
      do i=1,npsx
         do j=1,npsv
            psfxv(i,j)=0.
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
      include 'phasecom.f'
      include 'ndimsdecl.f'
      include 'meshcom.f'
      include 'partcom.f'
      integer ispecies,id
      integer i,ierr,ixbin,ivbin
      real v,x

      call pszero()
c Ensure the limits etc of the phase space array are set.
      psxmin=xmeshstart(id)
      psxmax=xmeshend(id)
      psvmax=3.
      psvmin=-3.
c The centers of the bins in phase space (redundancy negligible).
      do i=1,npsx
         psx(i)=psxmin+(i-0.5)*(psxmax-psxmin)/npsx
      enddo
      do i=1,npsv
         psv(i)=psvmin+(i-0.5)*(psvmax-psvmin)/npsv
      enddo

c Accumulate
      do i=iicparta(ispecies),iocparta(ispecies)
         x=x_part(id,i)
         v=x_part(id+ndims,i)
         ixbin=int(.99999*(x-psxmin)/(psxmax-psxmin)*float(npsx)+1)
         ivbin=int(.99999*(v-psvmin)/(psvmax-psvmin)*float(npsv)+1)
         if(ivbin.gt.0 .and. ivbin.le.npsv
     $        .and.x_part(iflag,i).ne.0)then
c Not for velocity beyond the vrange
            psfxv(ixbin,ivbin)=psfxv(ixbin,ivbin)+1.
         endif
      enddo
c All reduce to sum the distributions from all processes.
      call mpiallreducesum(psfxv,npsx*npsv,ierr)
      end
c***********************************************************************
      subroutine phasewrite(phasefilename,nu,x,u,t)
c Write file with phasespace data plus u(x) if length nu != 0.
      character*(*) phasefilename
      integer nu
      real x(nu),u(nu),t
      include 'phasecom.f'

      open(12,file=phasefilename,status='unknown',form='unformatted',err
     $     =101)
      write(12)npsx,npsv,nu
      write(12)psfxv,psvmax,psvmin,psxmax,psxmin,psx,psv
      if(nu.ne.0)write(12)(x(i),u(i),i=1,nu),t
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
      include 'phasecom.f'

      nuin=nu
      open(13,file=phasefilename,status='old',form='unformatted',err
     $     =101)
      read(13)npsxf,npsvf,nu
      if(npsx.ne.npsxf.or.npsv.ne.npsvf)then
         write(*,'(a,2i4,a,2i4)')
     $        'Incorrect phase-space array length allocation npsx,npsv',
     $        npsx,npsv,' file requires',npsxf,npsvf
         stop
      else
         read(13)psfxv,psvmax,psvmin,psxmax,psxmin,psx,psv
         if(nu.gt.nuin)stop 'phaseread file nu length too great'
         if(nu.ne.0)read(13)(x(i),u(i),i=1,nu),t
      endif
      close(13)
      return
 101  write(*,*)'Error opening file:',
     $     phasefilename(1:lentrim(phasefilename))
      end
c***********************************************************************
      subroutine phaseplot
      include 'phasecom.f'
      real cworka(npsx,npsv),zclv(2)

      call pltinit(psxmin,psxmax,psvmin,psvmax)
      call blueredgreenwhite()
      call axis()
      call axlabels('x','v')
c If unset, set psfmax for less than full range. But better set earlier.
      if(psfmax.eq.0.)call psfmaxset(1.5)
c Set extrema of coloring range from psfmax.
      zclv(1)=0.
      zclv(2)=psfmax
      icl=2
      icsw=1+16+32
c Using triangular gradients gives too large a ps file. +64 not.
      call contourl(psfxv,cworka,npsx,npsx,npsv,zclv,icl,psx,psv,icsw) 
      call gradlegend(zclv(1),zclv(2),.3,1.15,.7,1.15,.003,.true.)
c If needed, do pltend externally.
      end
c**********************************************************************
      subroutine psfmaxset(fac)
c Set the value of psfmax to fac times the maximum in the current array.
      real fac
      include 'phasecom.f'
      call minmax2(psfxv,npsx,npsx,npsv,pmin,pmax)
      psfmax=pmax*fac
      end
