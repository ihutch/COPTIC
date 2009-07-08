c Particle advancing routine
      subroutine padvnc(ndims,cij,u,iLs)
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
c Local storage
      integer ixp(ndims_mesh)
      real field(ndims_mesh)
      real xfrac(ndims_mesh)
c Make this always last to use the checks.
      include 'partcom.f'

      if(ndims.ne.ndims_mesh)
     $        stop 'Padvnc incorrect ndims number of dimensions'
      ic1=2*ndims+1
      ndimsx2=2*ndims
c
      nrein=0
      iocthis=0
c We ought not to need to calculate the iregion, since it should be
c known and if a particle is outside it, it would have been reinjected:
c But for now:
c      iregion=insideall(ndims,x_part(1,1))
c which will give an erroneous answer if 1 is an unfilled slot. So
      iregion=iregion_part
      n_part=0
c At most do over all particle slots. But generally we break earlier.
      do i=1,n_partmax
         dtprec=dt
         dtpos=dt
 100     continue
c If this particle slot is occupied.
         if(if_part(i).ne.0)then
c Find out where we are (if we don't already know).
c Should not be necessary if chargetomesh has been called.
c         call partlocate(i,iLs,iu,ixp,xfrac,iregion)
c         write(*,*)(x_part(ndimsx2+kk,i)-xfrac(kk),kk=1,3)

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
c Get the ndims field components at this point. 
c We only use x_part information for location.
            do idf=1,ndims
               call getfield(
     $              ndims,cij(ic1),u,iLs
     $              ,xn(ixnp(idf)+1)
     $              ,idf
     $              ,x_part(ndimsx2+1,i)
     $              ,iregion,field(idf))
            enddo
c Accelerate          
            do j=4,6
               x_part(j,i)=x_part(j,i)+field(j-3)*dtaccel
            enddo
c Move
            do j=1,3
               x_part(j,i)=x_part(j,i)+x_part(j+3,i)*dtpos
            enddo          

            inewregion=insideall(ndims,x_part(1,i))
            if(inewregion.ne.iregion) then
c We left the region. 
               call tallyexit(i,inewregion-iregion)
c Reinject if we haven't exhausted complement.
               if(ninjcomp.eq.0 .or. nrein.lt.ninjcomp)then
                  call reinject(i,x_part(1,i),nrein)
c Find where we are, since we don't yet know?
c Might not be needed if we insert needed information in reinject,
c which might be less costly. (Should be something other than iregion?)
                  call partlocate(i,iLs,iu,ixp,xfrac,irg)
                  dtpos=dtpos*ran1(myid)
                  dtprec=0.
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
               call reinject(i,x_part(1,i),nrein)
c Find where we are, since we don't yet know?
c Might not be needed if we insert needed information in reinject,
               call partlocate(i,iLs,iu,ixp,xfrac,irg)
               dtpos=dtpos*ran1(myid)
               dtprec=0.
c Complete reinjection by advancing by random remaining.
c               goto 101
c Silence warning of jump to different block by jumping outside instead
c gives the same result as 101.
               goto 100
         elseif(i.ge.ioc_part)then
c We do not need to reinject new particles, and
c this slot is higher than all previously handled. There are no
c more active particles above it. So break
            goto 102
         endif
c 
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
c      write(*,*)'iocthis=',iocthis

      end
c***********************************************************************
c Diagnostic Particle advancing routine used only for testing.
      subroutine padvncdiag(ndims,cij,u,iLs)
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
c Local storage
      integer ixp(ndims_mesh)
      real field(ndims_mesh)
      real xfrac(ndims_mesh)
c Make this always last to use the checks.
      include 'partcom.f'

      if(ndims.ne.ndims_mesh)
     $        stop 'Padvnc incorrect ndims number of dimensions'
      ic1=2*ndims+1
      ndimsx2=2*ndims
c
      nrein=0
      iocthis=0
c We ought not to need to calculate the iregion, since it should be
c known and if a particle is outside it, it would have been reinjected:
c But for now:
      iregion=insideall(ndims,x_part(1,1))
      n_part=0
c At most do over all particle slots. But generally we break earlier.
      do i=1,n_partmax
         dtprec=dt
         dtpos=dt
 100     continue
c If this particle slot is occupied.
         if(if_part(i).ne.0)then
c Find out where we are (if we don't already know).
c Should not be necessary if chargetomesh has been called.
c         call partlocate(i,iLs,iu,ixp,xfrac,iregion)
c         write(*,*)(x_part(ndimsx2+kk,i)-xfrac(kk),kk=1,3)

c Subcycle start.
 101        continue
c Use dtaccel for acceleration. May be different from dt if there was
c a reinjection (or collision).
            dtaccel=0.5*(dt+dtprec)
c Get the ndims field components at this point. 
c We only use x_part information for location.
            do idf=1,ndims
               call getfield(
     $              ndims,cij(ic1),u,iLs
     $              ,xn(ixnp(idf)+1)
     $              ,idf
     $              ,x_part(ndimsx2+1,i)
     $              ,iregion,field(idf))
            enddo
c Accelerate          
            do j=4,6
               x_part(j,i)=x_part(j,i)+field(j-3)*dtaccel
            enddo
c Move
            do j=1,3
               x_part(j,i)=x_part(j,i)+x_part(j+3,i)*dtpos
            enddo          

            inewregion=insideall(ndims,x_part(1,i))

            if(i.lt.10)then
               write(*,*)i,field,iregion,inewregion
            endif

            if(inewregion.ne.iregion) then
c We left the region. 
               call tallyexit(i,inewregion-iregion)
c Reinject if we haven't exhausted complement.
               if(ninjcomp.eq.0 .or. nrein.lt.ninjcomp)then
                  call reinject(i,x_part(1,i),nrein)
c Find where we are, since we don't yet know?
c Might not be needed if we insert needed information in reinject,
c which might be less costly. (Should be something other than iregion?)
                  call partlocate(i,iLs,iu,ixp,xfrac,irg)
                  dtpos=dtpos*ran1(myid)
                  dtprec=0.
                  if(i.lt.100)then
                     write(*,*)'Reinjected',i,(x_part(kk,i),kk=1,9)
c     $                    ,dtpos,dtprec,iu,ixp,xfrac,irg
                  endif

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
               call reinject(i,x_part(1,i),nrein)
c Find where we are, since we don't yet know?
c Might not be needed if we insert needed information in reinject,
               call partlocate(i,iLs,iu,ixp,xfrac,irg)
               dtpos=dtpos*ran1(myid)
               dtprec=0.
c Complete reinjection by advancing by random remaining.
c               goto 101
c Silence warning of jump to different block by jumping outside instead
c gives the same result as 101.
               goto 100
         elseif(i.ge.ioc_part)then
c We do not need to reinject new particles, and
c this slot is higher than all previously handled. There are no
c more active particles above it. So break
            goto 102
         endif
c 
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
c      write(*,*)'iocthis=',iocthis

      end
c***********************************************************************
      subroutine partlocate(i,iLs,iu,ixp,xfrac,iregion)
c Locate the particle numbered i (from common partcom) 
c in the mesh (from common meshcom).
c Return the offset of the base of its cell in iu.
c Return the integer cell-base coordinates in ixp(ndims)
c Return the fractions of cell width at which located in xfrac(ndims)
c Return the region identifier in iregion.
c Store the mesh position into common partcom (x_part).

c meshcom provides ixnp, xn, the mesh spacings. (+ndims_mesh)
      include 'meshcom.f'
      parameter (ndimsx2=ndims_mesh*2)
      integer i,iu,iregion
      integer iLs(ndims_mesh+1)
      integer ixp(ndims_mesh)
      real xfrac(ndims_mesh)

      include 'partcom.f'

      iregion=insideall(ndims_mesh,x_part(1,i))
      iu=0
      do id=1,ndims_mesh
c Offset to start of dimension-id-position array.
         ioff=ixnp(id)
c xn is the position array for each dimension arranged linearly.
c Find the index of xprime in the array xn:
         ix=interp(xn(ioff+1),ixnp(id+1)-ioff,x_part(id,i),xm)
         xfrac(id)=xm-ix
         x_part(ndimsx2+id,i)=xm
         ixp(id)=ix
c should be ix-1
         iu=iu+(ix-1)*iLs(id)
      enddo
      
      end
