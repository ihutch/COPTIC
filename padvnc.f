c Particle advancing routine
      subroutine padvnc(ndims,iLs,cij,u)
c If ninjcomp (in partcom) is non-zero, then we are operating in a mode
c where the number of reinjections at each timestep is prescribed.
c Otherwise we are using a fixed number npart of particles.
c
c Number of dimensions: ndims
      integer ndims
c Storage size of the mesh arrays.
c      real cij(2*ndims+1,nx,ny,nz)
c      real u(nx,ny,nz)
      real cij(*),u(*)

c Array structure vectors: (1,nx,nx*ny,nx*ny*nz)
      integer iLs(ndims+1)

c sormesh provides ixnp, xn, the mesh spacings. (+ndims_mesh)
      include 'meshcom.f'
c Alternatively they could be passed, but we'd then need parameter.
      include 'myidcom.f'
      include '3dcom.f'
      external linregion
      logical linregion
c Include this only for testing with Coulomb field.
      include 'plascom.f'
c Local storage
      integer ixp(ndims_mesh)
      real field(ndims_mesh)
      real xfrac(ndims_mesh)
      logical linmesh
c Testing storage
      real fieldp(ndims_mesh)
c Needed only to print out averein:
c      common /reinextra/averein,adeficit
c Make this always last to use the checks.
      include 'partcom.f'

      if(ndims.ne.ndims_mesh)
     $        stop 'Padvnc incorrect ndims number of dimensions'
      ic1=2*ndims+1
      ndimsx2=2*ndims
c Initialize. Set reinjection potential. We start with zero reinjections.
c      write(*,*)'Setting averein in padvnc.',phirein
      call avereinset(phirein)
c      write(*,*)'averein',averein
      phirein=0
      nrein=0
      nlost=0
c      ninner=0
      iocthis=0
      n_part=0
c At most do over all particle slots. But generally we break earlier.
      do i=1,n_partmax
         dtprec=dt
         dtpos=dt
 100     continue
c If this particle slot is occupied.
         if(if_part(i).ne.0)then
c Get its region
            iregion=insideall(ndims,x_part(1,i))
c Subcycle start.
 101        continue
c Use dtaccel for acceleration. May be different from dt if there was
c a reinjection (or collision).
            dtaccel=0.5*(dt+dtprec)
c Check the fraction data is not stupid and complain if it is.
            if(x_part(ndimsx2+1,i).eq.0. .and.
     $           x_part(ndimsx2+2,i).eq.0. .and.
     $           x_part(ndimsx2+3,i).eq.0.) then
               write(*,*)'Zero fractions',i,ioc_part,if_part(i)
     $              ,nrein,ninjcomp
            endif
c---------------------------------
            if(.true.)then
c Get the ndims field components at this point. 
c We only use x_part information for location. So we need to pass
c the region information.
               do idf=1,ndims
                  call getfield(
     $                 ndims,cij(ic1),u,iLs
     $                 ,xn(ixnp(idf)+1)
     $                 ,idf
     $                 ,x_part(ndimsx2+1,i)
     $                 ,imaskregion(iregion),field(idf))
                  if(abs(field(idf)).gt.1.e12)then
                     write(*,*)'Field corruption(?)',idf,field
                  endif
               enddo
c               if(.false.)then
c Testing only, of the few-argument field evaluator.
c                  call fieldatpoint(x_part(1,i),u,cij,iLs,fieldp)
c                  if(fieldp(2).ne.field(2))then
c                     write(*,'(i5,a,6f10.6)')i,' Point-field:',
c     $                    fieldp,((fieldp(j)-field(j)),j=1,ndims)
c                  endif
c               endif
c--------------------------------
c            else
c Testing with pure coulomb field from phip potential at r=1.
c               r2=0.
c               do idf=1,ndims
c                  r2=r2+x_part(idf,i)**2
c               enddo
c               r3=sqrt(r2)**3
c               do idf=1,ndims
c                  field(idf)=x_part(idf,i)*phip/r3
c               enddo
            endif
c--------------------------------
c Accelerate          
            do j=ndims+1,2*ndims
               x_part(j,i)=x_part(j,i)+field(j-3)*dtaccel
            enddo
c Move
            do j=1,ndims
               x_part(j,i)=x_part(j,i)+x_part(j+3,i)*dtpos
            enddo          

            call partlocate(i,ixp,xfrac,inewregion,linmesh)
c Old approach 
c            inewregion=insideall(ndims,x_part(1,i))
c If we crossed a boundary, do tallying.
            if(inewregion.ne.iregion)
     $           call tallyexit(i,inewregion-iregion)
            if(.not.linmesh .or.
     $         .not.linregion(ibool_part,ndims,x_part(1,i)))then
c We left the mesh or region. 
c Reinject if we haven't exhausted complement.
               if(ninjcomp.eq.0 .or. nrein.lt.ninjcomp)then
                  if_part(i)=1
                  call reinject(x_part(1,i),ilaunch)
c Find where we are, since we don't yet know.
                  call partlocate(i,ixp,xfrac,iregion,linmesh)
                  if(.not.linmesh)stop 'Reinject out of region'
                  dtpos=dtpos*ran1(myid)
                  dtprec=0.
                  nlost=nlost+1
                  nrein=nrein+ilaunch
                  phi=getpotential(u,cij,iLs,x_part(2*ndims+1,i)
     $                 ,imaskregion(iregion),2)
                  if(.false.)then
c Testing of potential at point. 
                     phicomp=potentialatpoint(x_part(1,i),u,cij,iLs)
                     if(phicomp.ne.phi)write(*,*)
     $                'getpotential vs potential at point difference:'
     $                    ,phi,phicomp
                  endif
                  phirein=phirein+ilaunch*phi
                  call diaginject(x_part(1,i))
c Complete reinjection by advancing by random remaining.
                  goto 101
               else
                  if_part(i)=0
               endif
            else
c The standard exit point for a particle that is active
               iocthis=i
               n_part=n_part+1
            endif
         elseif(ninjcomp.ne.0.and.nrein.lt.ninjcomp)then
c An unfilled slot. Fill it if we need to.
               if_part(i)=1
               call reinject(x_part(1,i),ilaunch)
c Find where we are, since we don't yet know.
               call partlocate(i,ixp,xfrac,iregion,linmesh)
               if(.not.linmesh)stop 'Reinject out of region'
               dtpos=dtpos*ran1(myid)
               dtprec=0.
               nlost=nlost+1
               nrein=nrein+ilaunch
               phi=getpotential(u,cij,iLs,x_part(2*ndims+1,i)
     $                 ,imaskregion(iregion),2)
               phirein=phirein+ilaunch*phi
               call diaginject(x_part(1,i))
c Silence warning of block jump by jumping outside instead of 101.
               goto 100
         elseif(i.ge.ioc_part)then
c We do not need to reinject new particles, and
c this slot is higher than all previously handled. There are no
c more active particles above it. So break
c            write(*,*)'break',i,ioc_part
            goto 102
         endif
c Special diagnostic orbit tracking:
         if(i.le.norbits.and. if_part(i).ne.0)then
            iorbitlen(i)=iorbitlen(i)+1
            xorbit(iorbitlen(i),i)=x_part(1,i)
            yorbit(iorbitlen(i),i)=x_part(2,i)
            zorbit(iorbitlen(i),i)=x_part(3,i)
         endif
      enddo
 102  continue
      if(ninjcomp.ne.0 .and. nrein.lt.ninjcomp)then
         write(*,*)'WARNING: Exhausted n_partmax=',n_partmax,
     $        '  before ninjcomp=',ninjcomp,' . Increase n_partmax?'
      endif
      ioc_part=iocthis
c      write(*,*)'iocthis=',iocthis,i

c Finished this particle step. Calculate the average reinjection 
c potential
      if(nrein.gt.0)then
         phirein=phirein/nrein
c         if(myid.eq.0)
c         write(*,'(a,f12.6,$)')' Phirein=',phirein
c         write(*,*)' nlost=',nlost,' nrein=',nrein,' ninner=',ninner
         if(phirein.gt.0.)then
c            if(myid.eq.0)write(*,*)'PROBLEM: phirein>0:',phirein
            phirein=0.
         endif
      else
         if(ninjcomp.gt.100)write(*,*)'No reinjections'
      endif

c      write(*,*)'Padvnc',n_part,nrein,ilaunch,ninjcomp,n_partmax
      end
c***********************************************************************
      subroutine partlocate(i,ixp,xfrac,iregion,linmesh)

c Locate the particle numbered i (from common partcom) 
c in the mesh (from common meshcom).
c Return the integer cell-base coordinates in ixp(ndims)
c Return the fractions of cell width at which located in xfrac(ndims)
c Return the region identifier in iregion.
c Return whether the particle is in the mesh in linmesh.
c Store the mesh position into common partcom (x_part).

c meshcom provides ixnp, xn, the mesh spacings. (+ndims_mesh)
      include 'meshcom.f'
      parameter (ndimsx2=ndims_mesh*2)
      integer i,iregion
      integer ixp(ndims_mesh)
      real xfrac(ndims_mesh)
      logical linmesh

      include 'partcom.f'

      linmesh=.true.
      iregion=insideall(ndims_mesh,x_part(1,i))
      do id=1,ndims_mesh
c Offset to start of dimension-id-position array.
         ioff=ixnp(id)
c xn is the position array for each dimension arranged linearly.
c Find the index of xprime in the array xn:
         ix=interp(xn(ioff+1),ixnp(id+1)-ioff,x_part(id,i),xm)
         xfrac(id)=xm-ix
         x_part(ndimsx2+id,i)=xm
         ixp(id)=ix
         if(ix.eq.0)then
            linmesh=.false.
c            write(*,'(a,i7,i3,'' '',6g12.4)')
c     $        ' Outside domain',i,id,(x_part(ii,i),ii=1,6)
         endif
      enddo
      
      end
c***************************************************************
      real function getnullpot()
      getnullpot=0.
      end
