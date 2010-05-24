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

c Meshcom provides ixnp, xn, the mesh spacings. (+ndims_mesh)
      include 'meshcom.f'
      include 'myidcom.f'
      include '3dcom.f'
      external linregion
      logical linregion
      include 'plascom.f'
c Collision settings.
      include 'colncom.f'

c Local storage
      parameter (fieldtoolarge=1.e12)
      integer ixp(ndims_mesh)
      real field(ndims_mesh)
      real xfrac(ndims_mesh)
      logical linmesh
      logical lcollided
c Testing storage
c      real fieldp(ndims_mesh)
c Needed only to print out averein:
c      common /reinextra/averein,adeficit
c Make this always last to use the checks.
      include 'partcom.f'

      tisq=sqrt(Ti)
      lcollided=.false.

      if(ndims.ne.ndims_mesh)
     $        stop 'Padvnc incorrect ndims number of dimensions'
c-----------------------------------------------------------------
      ic1=2*ndims+1
      ndimsx2=2*ndims
c Initialize. Set reinjection potential. We start with zero reinjections.
c      write(*,*)'Setting averein in padvnc.',phirein
      call avereinset(phirein)
      phirein=0
      nrein=0
      nlost=0
      iocthis=0
      n_part=0
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c At most do over all particle slots. But generally we break earlier.
      do i=1,n_partmax
         dtremain=0.
         dtprec=dt
         dtpos=dt
c=================================
c Decide nature of this slot
         if(if_part(i).ne.0)then
c Standard occupied slot. Proceed.
         elseif(ninjcomp.ne.0.and.nrein.lt.ninjcomp)then
c An unfilled slot needing to be filled. Jump to reinjection.
            goto 200
         elseif(i.ge.ioc_part)then
c We do not need to reinject new particles, and
c this slot is higher than all previously handled. There are no
c more active particles above it. So break from do loop.
            goto 102
         else
c Unoccupied but there might be slots occupied above so just cycle.
            goto 300
         endif
c================ Occupied Slot Treatment =================
c Occupied slot. Get its region
         iregion=insideall(ndims,x_part(1,i))
 100     continue
c Subcycle restart.
c---------- Subcycling ----------------
c One might have a conditional test to set subcycling for this particle.
         if(subcycle.ne.0)then
            dtc=subcycle*dt
            if(dtc.lt.dtpos)then
c Take sub-step.
               dtremain=dtpos+dtremain-dtc
               dtpos=dtc
            endif
         endif
c---------- Collisions ----------------
         if(colntime.ne.0)then
c Time to the first collision
            dtc=-alog(ran1(myid))*colntime
c            dtc=ran1(myid)*colntime
            if(dtc.lt.dtpos)then
c               write(*,*)'Collision',dtpos,dtc,colntime
c We collided during this step. Do the partial step.  
               lcollided=.true.
               dtremain=dtpos+dtremain-dtc
               dtpos=dtc
            endif
         endif
c--------------------------------------
c Check the fraction data is not stupid and complain if it is.
c Ought not to be necessary, but this is a safety check. 
         if(x_part(ndimsx2+1,i).eq.0. .and.
     $        x_part(ndimsx2+2,i).eq.0. .and.
     $        x_part(ndimsx2+3,i).eq.0.) then
            write(*,*)'Zero fractions',i,ioc_part,if_part(i)
     $           ,nrein,ninjcomp
         endif
c---------------------------------
c Get the ndims field components at this point. 
c We only use x_part information for location. So we need to pass
c the region information.
         do idf=1,ndims
            call getfield(ndims,cij(ic1),u,iLs
     $           ,xn(ixnp(idf)+1),idf,x_part(ndimsx2+1,i)
     $           ,imaskregion(iregion),field(idf))
            if(.not.abs(field(idf)).lt.fieldtoolarge)then
               write(*,*)'Field corruption(?)',i,idf,field
     $              ,(x_part(kk,i),kk=1,3*ndims)
            endif
         enddo
c Example of testing code: the few-argument field evaluator.
c                  call fieldatpoint(x_part(1,i),u,cij,iLs,fieldp)
c                  if(fieldp(2).ne.field(2))then
c                     write(*,'(i5,a,6f10.6)')i,' Point-field:',
c     $                    fieldp,((fieldp(j)-field(j)),j=1,ndims)

c---------------- Particle Moving ----------------
c Use dtaccel for acceleration. May be different from dtpos if there was
c a reinjection or collision.
         dtaccel=0.5*(dtpos+dtprec)
c Accelerate          
         do j=ndims+1,2*ndims
            x_part(j,i)=x_part(j,i)+field(j-3)*dtaccel
         enddo
c Move
         do j=1,ndims
            x_part(j,i)=x_part(j,i)+x_part(j+3,i)*dtpos
         enddo          

         if(lcollided)then
c Treat collided particle at (partial) step end
            call postcollide(i,tisq)
            lcollided=.false.
         endif

         call partlocate(i,ixp,xfrac,inewregion,linmesh)
c---------------------------------
c If we crossed a boundary, do tallying.
         if(inewregion.ne.iregion)
     $        call tallyexit(i,inewregion-iregion)
c------------ Possible Reinjection ----------
         if(.not.linmesh .or.
     $        .not.linregion(ibool_part,ndims,x_part(1,i)))then
c We left the mesh or region.
            if(ninjcomp.eq.0 .or. nrein.lt.ninjcomp)goto 200
c Reinject because we haven't exhausted complement. Else empty slot
            if_part(i)=0
         else
c The standard exit point for a particle that is active
            iocthis=i
            n_part=n_part+1
         endif
c--------------------------------------------
c If we haven't completed this full step, subcycle to do so.
         if(dtremain.gt.0.)then
            dtprec=dtpos
            dtpos=dtremain
            dtremain=0.
            iregion=inewregion
            goto 100
         endif
         goto 300
c================= End of Occupied Slot Treatement ================
c Reinjection:
 200     continue
         if_part(i)=1
         call reinject(x_part(1,i),ilaunch)
         call partlocate(i,ixp,xfrac,iregion,linmesh)
         if(.not.linmesh)then
            write(*,*)'Reinject out of region',i,xfrac
            stop
         endif
c         dtpos=dt*ran1(myid)
c         dtpos=dtpos*ran1(myid)
         dtpos=(dtpos+dtremain)*ran1(myid)
         dtprec=0.
         dtremain=0.
         nlost=nlost+1
         nrein=nrein+ilaunch
         phi=getpotential(u,cij,iLs,x_part(2*ndims+1,i)
     $        ,imaskregion(iregion),2)
         phirein=phirein+ilaunch*phi
         call diaginject(x_part(1,i))
c Restart the rest of the advance
         goto 100
c----------------------------------------------
c Non-injection completion.
 300     continue
c Special diagnostic orbit tracking:
         if(i.le.norbits.and. if_part(i).ne.0)then
            iorbitlen(i)=iorbitlen(i)+1
            xorbit(iorbitlen(i),i)=x_part(1,i)
            yorbit(iorbitlen(i),i)=x_part(2,i)
            zorbit(iorbitlen(i),i)=x_part(3,i)
         endif
      enddo
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c End of cycle through particles.
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
         isz=ixnp(id+1)-ioff
         ix=interp(xn(ioff+1),isz,x_part(id,i),xm)
         xfrac(id)=xm-ix
         x_part(ndimsx2+id,i)=xm
         ixp(id)=ix
c         if(ix.eq.0)then
c This more complete test is necessary and costs perhaps 3% extra time.
c It would be cheaper to change interp to be exclusive of limits.
c         if(ix.eq.0.or.xm.le.1..or.xm.ge.float(isz))then
c Invert the test so that any NAN also rejects particle, recovering
c from problems, one might hope.
         if(.not.(ix.ne.0.and.xm.gt.1..and.xm.lt.float(isz)))then
            linmesh=.false.
         endif
c specific particle test
c         if(i.eq.2298489) write(*,*)i,isz,ix,xm,linmesh
      enddo
      
      end
c********************************************************************
      subroutine postcollide(i,tisq)
      include 'partcom.f'
      include 'colncom.f'
c Get new velocity; reflects neutral maxwellian shifted by vneutral.
      do k=1,npdim
         x_part(npdim+k,i)=tisq*gasdev(idum)
      enddo
c Vneutral is in z-direction.
      x_part(npdim+npdim,i)=x_part(npdim+npdim,i)+vneutral
      end
c***************************************************************
      real function getnullpot()
      getnullpot=0.
      end
