c Initialize the flux data, determining what we are saving and where.
c Routine depends upon knowledge of what the actual problem is.
c So this is just an example, giving the effective template.
      subroutine fluxdatainit()
      include '3dcom.f'
c-----------------------------------------------
c Initialize here to avoid giant block data program.
      nf_step=0
      do i=1,nf_quant
         do j=1,nf_obj
            nf_posno(i,j)=0
            do k=1-nf_posdim,nf_maxsteps
               nf_address(i,j,k)=0
            enddo
         enddo
      enddo
      do i=1,nf_datasize
         ff_data(i)=0.
      enddo
c
c Example: save two fluxes (i=2) for only one object (j=1).
c------------------------------------------------
c Initialize a non-zero mapping from 3dobjects to flux surfaces.
      nf_map(1)=1
c Arbitrary choice of a parameter.
      nfluxes=50
c Set the correct numbers of quantities and objects and the posnos
c for all these quantities/objects.
      mf_quant=2
      mf_obj=1
c There are nfluxes positions for each quantity.
      nf_posno(1,1)=nfluxes
      nf_posno(2,1)=nfluxes
c-------------------------------------------------
c Now we create the addressing arrays etc.
      call nfaddressinit()
c After which, nf_address(i,j,k) points to the start of data for 
c quantity i, object j, step k. 
c So we can pass nf_data(nf_address(i,j,k)) as a vector start.
c-------------------------------------------------
c The k=1-nf_posdim to k=0
c slots exist for us to put descriptive information, such
c as the angle-values that provide positions to correspond to the fluxes.
c For example putting angle data into k=0:
      do ip=1,nf_posno(1,1)
         c=-1.+(ip-0.5)*2./nfluxes
         ff_data(nf_address(nf_flux,1,0)+ip-1)=c
      enddo

      end
c******************************************************************
      subroutine nfaddressinit()
      include '3dcom.f'
c General iteration given correct settings of nf_posno. Don't change!
c Zero nums to silence incorrect warnings.
      numdata=0
      numobj=0
      nf_address(1,1,0)=1
      do k=1-nf_posdim,nf_maxsteps
         if(k.gt.1-nf_posdim)
     $        nf_address(1,1,k)=nf_address(1,1,k-1)+numobj
         numobj=0
         do j=1,mf_obj
            if(j.gt.1)nf_address(1,j,k)=nf_address(1,j-1,k)+numdata
            numdata=0
            do i=1,mf_quant
               if(i.gt.1)nf_address(i,j,k)=
     $              nf_address(i-1,j,k)+nf_posno(i-1,j)
               numdata=numdata+nf_posno(i,j)
            enddo
            numobj=numobj+numdata
         enddo
      enddo
c Check if we might overrun the datasize.
      if(nf_address(mf_quant,mf_obj,nf_maxsteps)+numobj
     $     .gt.nf_datasize)then
         write(*,*)'DANGER: data from',mf_quant,mf_obj,nf_maxsteps,
     $        ' would exceed nf_datasize',nf_datasize
      else
         write(*,'(a,i10,a,i1,i2,i5,a,i10)')'Maximum nf_address:',
     $     nf_address(mf_quant,mf_obj,nf_maxsteps)+numobj,
     $        ' (',mf_quant,mf_obj,nf_maxsteps,
     $        ') is safely below nf_datasize:',nf_datasize
      endif
      end
c******************************************************************
      subroutine tallyexit(i,idiffreg)
c Document the exit of this particle just happened. 
c Assign the exit to a specific object, and bin on object.      
c On entry
c        i is particle number, 
c        idiffreg is the difference between its region and active.
c
      include 'partcom.f'
      include '3dcom.f'

      idiff=idiffreg
c Determine the object crossed. Maximum obj number if multiple.
      iobj=0
 1    if(idiff.eq.0) goto 2
      iobj=iobj+1
      idiff=idiff/2
      goto 1
 2    continue
      call objsect(i,iobj,ierr)
      if(ierr.ne.0)then
         write(*,*)'Tallyexit error',ierr,i,iobj
      endif
      
      end
c******************************************************************
c****************************************************************
      subroutine objsect(j,iobj,ierr)
c Find the intersection of the last step of particle j with  
c object iobj, and update the positioned-fluxes accordingly.
c
c Currently implemented only for object which is
c ndims-dimensional spheroid of semi-radii rc(ndims), center xc.
c
      include '3dcom.f'
      include 'partcom.f'

      real xc(npdim),rc(npdim),x1(npdim),x2(npdim)

      ierr=0

c Do nothing for untracked object
      if(nf_map(iobj).eq.0)return

      itype=obj_geom(1,iobj)
c Use only bottom 8 bits:
      itype=itype-256*(itype/256)

      if(itype.eq.1)then
c Sphere intersection.
         A=0.
         B=0.
         C=-1.
         D=-1.
c x1 and x2 are the coordinates in system in which sphere 
c has center 0 and radius 1.
         do i=1,npdim
            xc(i)=obj_geom(ocenter+i-1,iobj)
            rc(i)=obj_geom(oradius+i-1,iobj)
            x1(i)=(x_part(i,j)-dt*x_part(i+3,j)-xc(i))/rc(i)
            x2(i)=(x_part(i,j)-xc(i))/rc(i)
            A=A+(x2(i)-x1(i))**2
            B=B+x1(i)*(x2(i)-x1(i))
            C=C+x1(i)**2
            D=D+x2(i)**2
         enddo
c This condition tests for a sphere crossing.
         if(D.ne.0. .and. D*C.le.0.)then
            if(B.ge.0. .and. A*C.le.0.) then
               fraction=(-B+sqrt(B*B-A*C))/A
c            elseif(B.lt.0. .and. A*C.ge.0.)then
            else
               fraction=(-B-sqrt(B*B-A*C))/A
            endif
c That should exhaust the possibilities.
c Get the actual position and bin it.
c
c Here we need code that decides which of the nf_posno for this object
c to update corresponding to this crossing, and then update it. 
c Example: bin by cos(theta)=x12(3) for nf_posno(nf_flux,iobj) bins. 
c Flux only.
            infobj=nf_map(iobj)
            z12=(1.-fraction)*x1(3)+fraction*x2(3)
            ibin=nf_posno(nf_flux,infobj)*(0.999999*z12+1.)*0.5
     $           + nf_address(nf_flux,infobj,nf_step)
            ff_data(ibin)=ff_data(ibin)+1
         else
c Did not intersect!
            ierr=1
         endif
      else
c Unknown object type.
         ierr=99
      endif

      end
c*********************************************************************
      subroutine timeave(nu,u,uave,ictl)
c Average a quantity u(nu) over steps with a certain decay number
c into uave.
c ictl controls the actions as follows:
c       0    do the averaging.
c       1    do nothing except set nstep=1.
c       2    do nothing except set nstave=nu
c       3    do both 1 and 2. 
c       99   do nothing except increment nstep.
c The 99 call should be done at the end of all usage of this routine
c for the present step.
      real u(nu),uave(nu)

      integer nstep,nstave
      data nstep/1/nstave/20/
c Normal call.
      if(ictl.eq.0)then
         do i=1,nu
            uave(i)=(uave(i)*(nstep-1)+u(i))/nstep
         enddo
         return
      endif
      
      if(ictl.eq.1 .or. ictl.eq.3)then
         nstep=1
      endif
      if(ictl.eq.2 .or. ictl.eq.3)then
         nstave=nu
      endif
      if(ictl.ge.99)then
         if(nstep.le.nstave) nstep=nstep+1
      endif

      end
c***********************************************************************
      subroutine fluxdiag()
      include '3dcom.f'
      
      sum=0
      do i=1,nf_posno(1,1)
         sum=sum+ff_data(nf_address(1,1,nf_step)+i-1)
      enddo
      write(*,*)'Total flux number',sum
      write(*,'(10f7.1)')(ff_data(nf_address(1,1,nf_step)+i-1),
     $     i=1,nf_posno(1,1))

      end
c***********************************************************************
c Averaging the flux data over all positions.
c The positions might be described by more than one dimension, but
c that is irrelevant to the averaging, though not to plotting.
      subroutine fluxave()
      include '3dcom.f'
      parameter (nfluxmax=200)
      real flux(nfluxmax),angle(nfluxmax)

      do i=1,nf_posno(1,1)
         flux(i)=0.
      enddo
      sum=0
      do i=1,nf_posno(1,1)
         do is=1,nf_step
            flux(i)=flux(i)+ff_data(nf_address(1,1,is)+i-1)
            sum=sum+ff_data(nf_address(1,1,is)+i-1)
         enddo
         flux(i)=flux(i)/nf_step
      enddo

c From here on is non-general and is mostly for testing.
      do i=1,nf_posno(1,1)
c Here's the assumption that k=0 is angle information, and all different
c we could make this more general by binning everything with the same
c angle together.
         angle(i)=ff_data(nf_address(1,1,0)+i-1)
      enddo
      write(*,*) 'Average Flux'
      write(*,'(10f8.4)')(flux(i),i=1,nf_posno(1,1))

      call autoplot(angle,flux,nf_posno(1,1))
      call axlabels('angle cosine','average counts')
c      do i=1,nf_step
c         call polyline(angle,ff_data(nf_address(1,1,i)),nf_posno(1,1))
c      enddo
      call pltend()

      end
c*******************************************************************
      subroutine outputflux(name)
c File name:
      character*(*) name
c Common data containing the BC-object geometric information
      include '3dcom.f'
c Particle common data
      include 'partcom.f'
c Plasma common data
      include 'plascom.f'
c Zero the name first. Very Important!
c Construct a filename that contains many parameters
c Using the routines in strings_names.f
      name='flux_'
      call nameappendexp(name,'T',Ti,1)
      call nameappendint(name,'v',nint(100*vd),3)
      if(rs.ge.100)then
         call nameappendint(name,'R',ifix(rs/10.),2)
      else
         call nameappendint(name,'r',ifix(rs),2)
      endif
      call nameappendint(name,'P',ifix(abs(phip)),2)
      call nameappendexp(name,'L',debyelen,1)
      np=nbcat(name,'.dat')
      write(*,*)name

      open(22,file=name,status='unknown',err=101)
      close(22,status='delete')
      open(22,file=name,status='new',form='unformatted',err=101)
      write(22)nf_step,mf_quant,mf_obj
c      write(*,*)nf_step,mf_quant,mf_obj
      write(22)((nf_posno(i,j),i=1,mf_quant),j=1,mf_obj)
      write(22)(((nf_address(i,j,k),i=1,mf_quant),j=1,mf_obj),
     $     k=1-nf_posdim,nf_step+1)
      write(22)(ff_data(i),i=1,nf_address(1,1,nf_step+1)-1)

      close(22)

      write(*,*)'Wrote flux data to ',name(1:lentrim(name))
      return

 101  continue
      write(*,*)'Error opening file:',name
      close(22,status='delete')

      end

c*****************************************************************
      subroutine readfluxfile(name)
      character*(*) name
      include '3dcom.f'

      open(23,file=name,status='old',form='unformatted',err=101)
      read(23)nf_step,mf_quant,mf_obj
      read(23)((nf_posno(i,j),i=1,mf_quant),j=1,mf_obj)
      read(23)(((nf_address(i,j,k),i=1,mf_quant),j=1,mf_obj),
     $     k=1-nf_posdim,nf_step+1)
      read(23)(ff_data(i),i=1,nf_address(1,1,nf_step+1)-1)
      close(23)

      write(*,*)'Read back flux data from ',name(1:lentrim(name))
      return
 101  write(*,*)'Error opening file:',name
      end
