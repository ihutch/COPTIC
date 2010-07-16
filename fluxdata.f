c Initialize the flux data, determining what we are saving and where.
c The objects whose flux is to be tracked are indicated by 
c obj_geom(ofluxtype,i). If this is zero, it is not tracked.
c The number of ibins in each (2) of the surface dimensions is indicated
c by obj_geom(ofn1[/2],i), and data space and addresses allocated.
c The uniform-spacing bin-positions are calculated and set.
      subroutine fluxdatainit(myid)
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
c-------------------------------------------------
c This initialization ought not to be necessary except for partforce.
      do k=1,nf_maxsteps
         do j=1,nf_obj
            do i=1,ns_ndims
               fieldforce(i,j,k)=0.
               pressforce(i,j,k)=0.
               partforce(i,j,k)=0.
            enddo
            charge_ns(j,k)=0.
         enddo
      enddo
c-------------------------------------------------
c Initialize object number
      mf_obj=0
      do i=1,ngeomobj
         if(obj_geom(ofluxtype,i).eq.0)then
c No flux setting for this object.
         elseif(obj_geom(ofluxtype,i).ge.1
     $           .and. obj_geom(ofluxtype,i).le.nf_quant)then
c Might eventually need more interpretation of fluxtype.
            mf_obj=mf_obj+1
            mf_quant(mf_obj)= obj_geom(ofluxtype,i)
c The mapped object number != object number.
            nf_map(i)=mf_obj
            nf_geommap(mf_obj)=i
            nfluxes=obj_geom(ofn1,i)*obj_geom(ofn2,i)
c There are nfluxes positions for each quantity.
c At present there's no way to prescribe different grids for each
c quantity in the input file. But the data structures could 
c accommodate such difference if necessary. 
            do j=1,mf_quant(mf_obj)
c       explicit for posdim=2 for now.
               nf_dimlens(j,mf_obj,1)=int(obj_geom(ofn1,i))
               nf_dimlens(j,mf_obj,2)=int(obj_geom(ofn2,i))
               nf_posno(j,mf_obj)=nfluxes
               if(j.eq.1 .and. myid.eq.0)
     $              write(*,'(a,i3,a,i3,a,i5,a,i3,a,i3,a,2i3)')
     $              ' Fluxinit of object',i
     $           ,'  type',int(obj_geom(ofluxtype,i))
     $           ,'. ',nf_posno(j,mf_obj),' flux positions:'
     $           ,nf_dimlens(j,mf_obj,1),'x',nf_dimlens(j,mf_obj,2)
     $           ,' Quantities',mf_quant(mf_obj)
            enddo
         else
            write(*,*)'==== Unknown flux type',obj_geom(ofluxtype,i)
            stop
         endif
      enddo

c      write(*,*)'nf_posno=',((nf_posno(j,i),j=1,2),i=1,4)
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
      io=0
      do i=1,ngeomobj
         if(obj_geom(ofluxtype,i).gt.0)then
            io=io+1
            ioff=0
            do j=1,mf_quant(io)
               do i2=1,nf_dimlens(j,io,2)
                  p=3.1415926*(-1.+2.*(i2-0.5)/nf_dimlens(j,io,2))
                  do i1=1,nf_dimlens(j,io,1)
                     c=-1.+2.*(i1-0.5)/nf_dimlens(j,io,1)
                     ip=i1+(i2-1)*int(nf_dimlens(j,io,1))
                     ff_data(nf_address(nf_flux,io,0)+ioff+ip-1)=c
                     ff_data(nf_address(nf_flux,io,-1)+ioff+ip-1)=p
                  enddo
               enddo
               ioff=ioff+nf_posno(j,io)
            enddo
c            write(*,*)'Set ff_data',i,ioff,nf_map(i),io
         endif
      enddo

      end
c******************************************************************
      subroutine nfaddressinit()
      include '3dcom.f'

c General iteration given correct settings of nf_posno. Don't change!
c Zero nums to silence incorrect warnings.
      numdata=0
      numobj=0
      nf_address(1,1,1-nf_posdim)=1
      do k=1-nf_posdim,nf_maxsteps
         if(k.gt.1-nf_posdim)
     $        nf_address(1,1,k)=nf_address(1,1,k-1)+numobj
         numobj=0
         do j=1,mf_obj
            if(j.gt.1)nf_address(1,j,k)=nf_address(1,j-1,k)+numdata
            numdata=0
            do i=1,mf_quant(j)
               if(i.gt.1)nf_address(i,j,k)=
     $              nf_address(i-1,j,k)+nf_posno(i-1,j)
               numdata=numdata+nf_posno(i,j)
c               write(*,*)i,j,k,nf_posno(i,j),nf_address(i,j,k),numdata
            enddo
            numobj=numobj+numdata
         enddo
      enddo
c Check if we might overrun the datasize.
      if(nf_address(nf_quant,mf_obj,nf_maxsteps)+numobj
     $     .gt.nf_datasize)then
         write(*,*)'DANGER: data from',nf_quant,mf_obj,nf_maxsteps,
     $        ' would exceed nf_datasize',nf_datasize
         stop
      else
      endif
      end
c******************************************************************
      subroutine tallyexit(i,idiffreg)
c Document the exit of this particle just happened. 
c Assign the exit to a specific object, and bin on object.
c (If it is a mapped object, decided by objsect.)      
c On entry
c        i is particle number, 
c        idiffreg is the difference between its region and active.
c
      include 'partcom.f'
      include '3dcom.f'

      idiff=abs(idiffreg)
      idp=idiff
c Determine (all) the objects crossed and call objsect for each.
      iobj=0
 1    if(idiff.eq.0) return
      iobj=iobj+1
      idiff=idiff/2
      if(idp.ne.idiff*2)then
         call objsect(i,iobj,ierr)
         if(ierr.ne.0)then
            ireg=insideall(ndims,x_part(1,i))
            r=0.
            r1=0.
            do id=1,3
               r=r+x_part(id,i)**2
               r1=r1+(x_part(id,i)-dt*x_part(id+3,i))**2
            enddo
            r=sqrt(r)
            r1=sqrt(r1)
            write(*,*)'Tallyexit error',ierr,i,iobj,idiffreg
            write(*,*)(x_part(k,i),k=1,3),r,ireg
            write(*,*)(x_part(k,i)-dt*x_part(k+3,i),k=1,3),r1
         endif
      endif
      idp=idp/2
      goto 1
      
      end
c******************************************************************
c****************************************************************
      subroutine objsect(j,iobj,ierr)
c Find the intersection of the last step of particle j with  
c object iobj, and update the positioned-fluxes accordingly.
c ierr is returned: 0 good. 1 no intersection. 99 unknown object.
c
c Currently implemented only for object which is
c ndims-dimensional spheroid of radii rc(ndims), center xc.
c With a equally spaced grid in cos\theta and psi.
      include '3dcom.f'
      include 'partcom.f'

      real xc(npdim),rc(npdim),x1(npdim),x2(npdim)

      ierr=0

c Do nothing for untracked object and report no error.
      if(nf_map(iobj).eq.0)return

      itype=int(obj_geom(otype,iobj))
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
c If we can decide direction of sphere crossing
         dir=D-C
         if(dir.ne.0.)then
            if(B.ge.0. .and. A*C.le.0.) then
               fraction=(-B+sqrt(B*B-A*C))/A
            else
               fraction=(-B-sqrt(B*B-A*C))/A
            endif
c That should exhaust the possibilities.
c
c This code decides which of the nf_posno for this object
c to update corresponding to this crossing, and then update it. 
            infobj=nf_map(iobj)
            z12=(1.-fraction)*x1(3)+fraction*x2(3)
c Example: bin by cos(theta)=x12(3) uniform grid in first nf_dimension. 
            ibin=int(nf_dimlens(nf_flux,infobj,1)*(0.999999*z12+1.)*0.5)
            x12=(1.-fraction)*x1(1)+fraction*x2(1)
            y12=(1.-fraction)*x1(2)+fraction*x2(2)
            psi=atan2(y12,x12)
            jbin=int(nf_dimlens(nf_flux,infobj,2)
     $           *(0.999999*psi/3.1415926+1.)*0.5)
            ijbin=ibin+jbin*nf_dimlens(nf_flux,infobj,1)
            iaddress=ijbin+nf_address(nf_flux,infobj,nf_step)
c D is the final radius -1 [!=0]. So its sign determines where we end.
c            sd=-sign(1.,D)
c Better to use the direction D-C, inwards counts positive:
            sd=-sign(1.,dir)
c Particle Flux.
            ff_data(iaddress)=ff_data(iaddress)+sd
c Perhaps ought to consider velocity interpolation.
            if(mf_quant(infobj).ge.2)then
c Momentum               
               iaddress=ijbin+nf_address(nf_gx,infobj,nf_step)
               ff_data(iaddress)=ff_data(iaddress)+ sd*x_part(4,j)
            endif
            if(mf_quant(infobj).ge.3)then
               iaddress=ijbin+nf_address(nf_gy,infobj,nf_step)
               ff_data(iaddress)=ff_data(iaddress)+ sd*x_part(5,j)
            endif
            if(mf_quant(infobj).ge.4)then
               iaddress=ijbin+nf_address(nf_gz,infobj,nf_step)
               ff_data(iaddress)=ff_data(iaddress)+ sd*x_part(6,j)
            endif
            if(mf_quant(infobj).ge.5)then
c Energy
               iaddress=ijbin+nf_address(nf_heat,infobj,nf_step)
               xx=0.
               do k=1,npdim
                  xx=xx+x_part(3+k,j)**2
               enddo
               ff_data(iaddress)=ff_data(iaddress)+ sd*xx
            endif
c Accumulate the particle force= momentum/time over whole object.
c Normalized to rhoinf
            do id=1,ns_ndims
               partforce(id,infobj,nf_step)=partforce(id,infobj,nf_step)
     $              +sd*x_part(ns_ndims+id,j)/dt/rhoinf
            enddo
c            write(*,'(a,3f10.4,i4,i4)')'Partforce='
c     $           ,(partforce(kk,infobj,nf_step),kk=1,3)
c     $           ,infobj,nf_step
c If the bins were different we would have to recalculate ibin. 
         else
c Did not intersect!
            write(*,*)'C=',C,' D=',D
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
c For rhoinf, dt
      include 'partcom.f'

      sum=0
      do i=1,nf_posno(1,1)
         sum=sum+ff_data(nf_address(1,1,nf_step)+i-1)
      enddo
      flux=sum/(4.*3.14159)/rhoinf/dt
      write(*,'(f6.3,''| '',$)')flux
c      write(*,'(a,f7.0,f8.3)')'Total flux',sum,
c     $     sum/(4.*3.14159)/rhoinf/dt
c      write(*,'(10f7.1)')(ff_data(nf_address(1,1,nf_step)+i-1),
c     $     i=1,nf_posno(1,1))

      end
c***********************************************************************
c Averaging the flux data for quantity abs(iquant)
c over all positions for object ifobj,
c for the steps n1 to n2.
c The positions might be described by more than one dimension, but
c that is irrelevant to the averaging.
c Plot the quantity iquant if positive (not if negative).
c Plotting does not attempt to account for the multidimensionality.
c Rhoinf is returned in rinf.
      subroutine fluxave(n1,n2,ifobj,iquant,rinf)
      integer n1,n2
c      logical lplot
      integer iquant
      include '3dcom.f'
      parameter (nfluxmax=200)
      real flux(nfluxmax),angle(nfluxmax)
      real fluxofstep(nf_maxsteps),step(nf_maxsteps)
      character*30 string

c      if(iquant.le.0)then
c         iq=1
c      else
         iq=abs(iquant)
c      endif
c If quantity asked for is not available, do nothing.
      if(iq.gt.mf_quant(ifobj))return

      if(n1.lt.1)n1=1
      if(n2.gt.nf_step)n2=nf_step
      if(n2-n1.lt.0)then
         write(*,*)'fluxave incorrect limits:',n1,n2
         return
      endif

      do i=1,nf_posno(iq,ifobj)
         flux(i)=0.
      enddo
      tot=0
      do i=1,nf_posno(iq,ifobj)
         do is=n1,n2
            flux(i)=flux(i)+ff_data(nf_address(iq,ifobj,is)+i-1)
         enddo
         tot=tot+flux(i)
         flux(i)=flux(i)/(n2-n1+1)
      enddo
      tdur=0.
      rinf=0.
      do is=1,n2
         fluxstep=0
         do i=1,nf_posno(iq,ifobj)
            fluxstep=fluxstep+ff_data(nf_address(iq,ifobj,is)+i-1)
         enddo
         step(is)=is
         fluxofstep(is)=fluxstep
         if(is.ge.n1)then
            tdur=tdur+ff_dt(is)
            rinf=rinf+ff_rho(is)*ff_dt(is)
         endif
c         write(*,*)'step',is,' fluxofstep',fluxofstep(is)
      enddo
      tot=tot/tdur
      rinf=rinf/tdur

c From here on is non-general and is mostly for testing.
      do i=1,nf_posno(iq,ifobj)
c Here's the assumption that k=0 is angle information.
c We could make this more general by binning everything with the same
c angle together. But instead we plot multiple points.
         angle(i)=ff_data(nf_address(iq,ifobj,0)+i-1)
      enddo
      write(*,'(a,i3,a,i3,a,i4,i4,a,f10.3)')' Average flux quant',
     $     iq,' object',ifobj,' over steps',n1,n2,' All Positions:',tot
      write(*,*)'rhoinf',rinf,'  Average collected per step by posn:'
      write(*,'(10f8.2)')(flux(i),i=1,nf_posno(1,ifobj))

      write(*,*)'Flux density*r^2, normalized to rhoinf'
     $     ,tot/(4.*3.14159)/rinf

      if(iquant.gt.0)then
         write(string,'(''Object '',i3,'' Quantity'',i3)')
     $        nf_geommap(ifobj),iquant
         call autoplot(step,fluxofstep,n2)
         call boxtitle(string)
         call axlabels('step','Spatially-summed flux number')
         call pltend()
         call automark(angle,flux,nf_posno(iq,ifobj),1)
         call boxtitle(string)
         call axlabels('First angle variable',
     $        'Time-averaged flux number')
         call pltend()
      endif

      end
c*******************************************************************
      subroutine writefluxfile(name)
c File name:
      character*(*) name
c Common data containing the BC-object geometric information
      include '3dcom.f'
c Particle common data
      include 'partcom.f'
c Plasma common data
      include 'plascom.f'
c 
      character*(100) charout
c Zero the name first. Very Important!
c Construct a filename that contains many parameters
c Using the routines in strings_names.f
      name=' '
      call nameconstruct(name)
c     np=nbcat(name,'.flx')
      call nbcat(name,'.flx')
c      write(*,*)name
      write(charout,51)debyelen,Ti,vd,rs,phip
 51   format('debyelen,Ti,vd,rs,phip:',5f10.4)

      open(22,file=name,status='unknown',err=101)
      close(22,status='delete')
      open(22,file=name,status='new',form='unformatted',err=101)
c This write sequence must be exactly that read below.
      write(22)charout
      write(22)debyelen,Ti,vd,rs,phip
      write(22)nf_step,mf_quant,mf_obj,(nf_geommap(j),j=1,mf_obj)
      write(22)(ff_rho(k),k=1,nf_step)
      write(22)(ff_dt(k),k=1,nf_step)
      write(22)((nf_posno(i,j),(nf_dimlens(i,j,k),k=1,nf_posdim)
     $     ,i=1,mf_quant(j)),j=1,mf_obj)
      write(22)(((nf_address(i,j,k),i=1,mf_quant(j)),j=1,mf_obj),
     $     k=1-nf_posdim,nf_step+1)
c This and three other places:
c      write(*,*)nf_step,'nf_address',(nf_address(1,1,i),i=1,nf_step+1)
      write(22)(ff_data(i),i=1,nf_address(1,1,nf_step+1)-1)

      write(22) fieldforce,pressforce,partforce,charge_ns
      close(22)

c      write(*,*)'Wrote flux data to ',name(1:lentrim(name))

      return
 101  continue
      write(*,*)'Error opening file:',name
      close(22,status='delete')
      end
c*****************************************************************
      subroutine readfluxfile(name,ierr)
      character*(*) name
      include '3dcom.f'
      include 'plascom.f'
      character*(100) charout

      open(23,file=name,status='old',form='unformatted',err=101)
      read(23)charout
      read(23)debyelen,Ti,vd,rs,phip
      read(23)nf_step,mf_quant,mf_obj,(nf_geommap(j),j=1,mf_obj)
      read(23)(ff_rho(k),k=1,nf_step)
      read(23)(ff_dt(k),k=1,nf_step)
      read(23)((nf_posno(i,j),(nf_dimlens(i,j,k),k=1,nf_posdim)
     $     ,i=1,mf_quant(j)),j=1,mf_obj)
      read(23)(((nf_address(i,j,k),i=1,mf_quant(j)),j=1,mf_obj),
     $     k=1-nf_posdim,nf_step+1)
      read(23)(ff_data(i),i=1,nf_address(1,1,nf_step+1)-1)
c If compatibility with flux files earlier than July 10 is required
c use the following version and ignore the last step.
c      read(23)(((nf_address(i,j,k),i=1,mf_quant(j)),j=1,mf_obj),
c     $     k=1-nf_posdim,nf_step)
c      read(23)(ff_data(i),i=1,nf_address(1,1,nf_step)-1)
      read(23) fieldforce,pressforce,partforce,charge_ns
      close(23)

      write(*,*)'Read back flux data from ',name(1:lentrim(name))
      write(*,*)charout(1:lentrim(charout))
      ierr=0
      return
 101  write(*,*)'Error opening file:',name
      ierr=1
      end
