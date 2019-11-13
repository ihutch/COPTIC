      subroutine padvnc(iLs,cij,u,ndiags,psum,diagsum,ispecies,ndiagmax)
! Particle advancing routine.
! If ninjcomp (in partcom) is non-zero, then we are operating in a mode
! where the number of reinjections at each timestep is prescribed.
! Otherwise we are using a fixed number n_part of particles.

! Storage size of the mesh arrays.
!      real cij(2*ndims+1,nx,ny,nz)
!      real u(nx,ny,nz)
      integer ndiags,ispecies
      real cij(*),u(*),psum(*)
      real diagsum(*)
!      real diagsum(64,64,128,8,2)

      include 'ndimsdecl.f'
! Array structure vector: (1,nx,nx*ny,nx*ny*nz)
      integer iLs(ndims+1)

! Meshcom provides ixnp, xn, the mesh spacings. (+ndims)
      include 'meshcom.f'
      include 'myidcom.f'
      include '3dcom.f'
      include 'ptchcom.f'
      external linregion
      logical linregion
      include 'plascom.f'
      include 'creincom.f'
! Collision settings.
      include 'colncom.f'
      include 'facebcom.f'

! Local parameters
      include 'dbgcom.f'
! Lowest particle to print debugging data for.
      integer npr
      parameter (npr=0)
      real fieldtoosmall
      parameter (fieldtoosmall=1.e-3)
      real thetamax
      parameter (thetamax=1.)
      real fieldtoolarge
      parameter (fieldtoolarge=1.e12)
! Local storage
      real fractions(10)
      real rannum
      integer ixp(ndims),ninjadd
      real Efield(ndims),adfield(ndims)
      real xfrac(ndims)
      real xprior(2*ndims)
      real xn1(ndims),xn2(ndims)
      logical linmesh
      logical lcollided,ltlyerr,lbtoolarge
      save adfield,lbtoolarge
      data lbtoolarge/.false./

! Make this always last to use the checks.
      include 'partcom.f'

!-------------------------------------------------------------
      reinspecies=ispecies
      ndeposits=0
      tisq=sqrt(Ts(ispecies))
      lcollided=.false.
      ncollided=0
      ic1=2*ndims+1
      ndimsx2=2*ndims
      theta=abs(dt*Bt*eoverms(ispecies))

!-----------------------------------------------------------------
! Initialize. Set reinjection potential. We start with zero reinjections.
      if(ispecies.eq.1)call cavereinset(phirein)
      phirein=0
      nrein=0
      nlost=0
      if(.not.lnotallp)ninjcompa(ispecies)=0
      iocthis=0
      nparta(ispecies)=0
      nsubc=0
      nwmax=20
      echarge=sign(1.,eoverms(ispecies))
      if(ispecies.eq.2)echarge=echarge*(1.-boltzamp)
      do idf=1,ndims
         adfield(idf)=0.
      enddo
! Set additional injection to zero by default.
      ninjadd=0
      if(ninjcompa(ispecies).ne.0)then
         call ranlux(rannum,1)
         if(rannum.lt.pinjcompa(ispecies))ninjadd=1
      endif
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! At most do over all particle slots for this species. 
! But generally we break earlier.
      numratio=numratioa(ispecies)
      do i=iicparta(ispecies),iicparta(ispecies+1)-1
         do icycle=1,numratio
         dtdone=0.
         dtremain=0.
         dtpos=dt/float(numratio)
!=================================
! Decide nature of this slot
         if(x_part(iflag,i).ne.0)then
! Standard occupied slot. Proceed.
         elseif(ninjcompa(ispecies).ne.0
     $        .and.nrein.lt.ninjcompa(ispecies)+ninjadd)then
! An unfilled slot needing to be filled. Jump to reinjection.
            goto 200
         elseif(i.ge.iocparta(ispecies))then
! We do not need to reinject new particles, and
! this slot is higher than all previously handled. There are no
! more active particles above it. So break from do loop.
            goto 102
         else
! Unoccupied but there might be slots occupied above so just cycle.
            goto 400
         endif
!================ Occupied Slot Treatment =================
! Occupied slot. Get its region
         iregion=insideall(ndims,x_part(1,i))
 100     continue
!--------------------------------------
! Check the fraction data is not stupid and complain if it is.
! Ought not to be necessary, but this is a safety check. 
! Remove after the reimplementation.
!         if(x_part(ndimsx2+1,i).eq.0. .and.
!     $        x_part(ndimsx2+2,i).eq.0. .and.
!     $        x_part(ndimsx2+3,i).eq.0.) then
!            write(*,*)'Zero fractions',i,iocparta(ispecies)
!     $           ,x_part(iflag,i),nrein,ninjcompa(ispecies)
!         endif
!---------------------------------
! adfield must always be zero arriving here regardless of irptch. 
! But we don't set it always because that would be expensive. 
! We only reset it below if it was changed
         irptch=IAND(iregion,iptch_mask)
         if(irptch.ne.0)then
! We are in a point-charge region. Get analytic part of force.
            call getadfield(irptch,adfield,x_part(1,i),2,ierr)
!            if(i.eq.1)write(*,'(i7,'' in ptch region'',i8,6f8.3)')
!     $           i,irptch,(x_part(k,i),k=1,3)
!     $           ,(adfield(k),k=1,3)
         endif
! Get the ndims field components at this point, from mesh potential.
! We only use x_part information for location. So we need to pass
! the region information.
         f2=0.
         r2=0.
         v2=0.
         do idf=1,ndims
! Test for need to get field:
            if(.not.(LPF(idf).and.ixnp(idf+1)-ixnp(idf).eq.3))then
!            if(.true.)then   
               call getfield(cij(ic1),u,iLs
     $              ,xn(ixnp(idf)+1),idf,x_part(ndimsx2+1,i)
     $              ,IAND(iregion,ifield_mask),Efield(idf))
               if(.not.abs(Efield(idf)).lt.fieldtoolarge)then
                  write(*,*)'Field corruption(?) by getfield.'
                  write(*,'(''Particle'',i8,'' dimension='',i2,
     $'' iregion='',i3,'' masked='',i3,'' Field='',3f10.5)')i,idf
     $                 ,iregion ,IAND(iregion,ifield_mask),Efield(idf)
                  write(*,'(''xp:'',9f8.4)')(x_part(kk,i),kk=1,3
     $                 *ndims)
                  write(*,'(''In region'',l2,''  iLs=''
     $,4i8)')linregion(ibool_part,ndims,x_part(1,i)),iLs
                  stop
               endif
            else
               Efield(idf)=0.
            endif
            Efield(idf)=Efield(idf)+adfield(idf)
            f2=f2+Efield(idf)**2
            r2=r2+x_part(idf,i)**2
            v2=v2+x_part(idf+3,i)**2
         enddo
! If needed, reset adfield.
         if(irptch.ne.0)then
            do idf=1,ndims
               adfield(idf)=0.
            enddo
         endif
         f1=sqrt(f2)
! Example of testing code: the few-argument field evaluator.
!                  call fieldatpoint(x_part(1,i),u,cij,iLs,fieldp)
!                  if(fieldp(2).ne.Efield(2))then
!                     write(*,'(i5,a,6f10.6)')i,' Point-field:',
!     $                    fieldp,((fieldp(j)-Efield(j)),j=1,ndims)
!--------------------------------------
         dtc=0.
!---------- Subcycling ----------------
         if(subcycle.ne.0)then
            if(f1*dtpos.gt.subcycle .and. 8*dtpos.gt.dt)then
! Reduce the position-step size.
               dtpos=dtpos/2.
! Subcycling backs up to the start of the current drift-step, and then 
! takes a position corresponding to a drift-step smaller by a factor of 2.
! That means the backed-up position is dtprec/4 earlier. There might be
! an error if dtpost and dtaccel are different. Is -dtaccel correct?
               dtback=-x_part(idtp,i)/4.
               call moveparticle(x_part(1,i),ndims,Bt ,Efield,Bfield
     $              ,vperp,dtback,-dtaccel,driftfields(1,ispecies)
     $              ,eoverms(ispecies))
               dtdone=dtdone+dtback
! The remaining time in step is increased by back up and by step loss.
               dtremain=dtremain-dtback+dtpos
! And the new dtprec is half as large:
               x_part(idtp,i)=x_part(idtp,i)/2.
               nsubc=nsubc+1
               if(i.lt.norbits)then
                  write(*,'(/,a,i3,4f10.5,$)')'Subcycle Set',i,dtpos
     $                 ,x_part(idtp,i),dtremain,dtdone+dtpos+dtremain
               endif
               goto 100
            endif
         endif
!---------- End of subcycling ---------
!---------- Collision Decision ----------------
         if(ispecies.eq.1.and.colntime.ne.0.)then
            if(colpow.ne.0.)then
! The collision time scaled by velocity
!     $              *(Ti/(v2+Ti))**(colpow/2.)
               call ranlux(dtc,1)
               dtc=-alog(dtc)*colntime
     $              *(Ti/(v2+Ti))**(colpow/2.)
            else
! Time to the first collision
               call ranlux(dtc,1)
               dtc=-alog(dtc)*colntime
            endif
            if(dtc.lt.dtpos)then
!               write(*,*)'Collision',dtpos,dtc,colntime
! We collided during this step. Do the partial step.  
               lcollided=.true.
! Set dtremain to how much time will be left after this step.
               dtremain=dtpos+dtremain-dtc
               dtpos=dtc
               ncollided=ncollided+1
            endif
         endif 
!---------- End Collision Decision ------------
!         if(i.lt.norbits)write(*,*)x_part(idtp,i),dtremain,dtdone
!---------------- Particle Advance ----------------
! Use dtaccel for acceleration. May be different from dtpos if there was
! a subcycle, reinjection or collision last step.
         dtaccel=0.5*(dtpos+x_part(idtp,i))
!----- Excessive Acceleration test. Drop particle without any tally. ---
         if(f1*dtaccel.gt.dropaccel)then
            ndropped=ndropped+1
! Why don't we tally this particle? It effectively has been absorbed.
! It carried in some momentum. That disappears. It is accounted for
! as being conveyed to any tallying object inside which it disappears.
!            write(*,*)'Excessive acceleration. Dropping particle',i
! Reinject if we haven't exhausted complement:
            if(ninjcompa(ispecies).eq.0
     $           .or. nrein.lt.ninjcompa(ispecies)+ninjadd)goto 200
! Else empty slot and move to next particle.
            x_part(iflag,i)=0
            goto 400            
         endif
! Accelerate used to be here -----------------
! Here we only include the electric field force. B-field force and vxB
! field is accommodated within moveparticle routine.
         if(Eneutral.ne.0.)then
! Calculate collision force --------
            if(Bt.ne.0.)then
! Only the Enparallel is needed for non-zero Bfield.
! This doesn't seem to be correct unless drift is in z-direction.
               Enpar=Eneutral*Bfield(ndims)
               do j=ndims+1,2*ndims
                  Eadd=Enpar*Bfield(j)
                  x_part(j,i)=x_part(j,i)+Eadd*dtaccel
                  do iob=1,mf_obj
                     igeom=nf_geommap(iob)
                     if(btest(iregion,igeom-1))then
                        colnforce(j,iob,nf_step) =colnforce(j,iob
     $                       ,nf_step)+Eadd*dtaccel
                     endif
                  enddo
               enddo
            else
! Including the extra Eneutral field, which is in the ndims direction.
! But correcting the collisional force appropriately.
               x_part(2*ndims,i)=x_part(2*ndims,i)+Eneutral*dtaccel
               do iob=1,mf_obj
                  igeom=nf_geommap(iob)
                  if(btest(iregion,igeom-1))then
                     colnforce(ndims,iob,nf_step)=colnforce(ndims
     $                    ,iob,nf_step)+Eneutral*dtaccel
                  endif
               enddo
            endif
         endif
! Move ----------------
! Save prior position and velocity.
         do id=1,2*ndims
            xprior(id)=x_part(id,i)
         enddo
         if(abs(theta).lt.thetamax)then
            call moveparticle(x_part(1,i),ndims,Bt ,Efield,Bfield,vperp
     $           ,dtpos,dtaccel,driftfields(1,ispecies)
     $           ,eoverms(ispecies))
         else
            if(.not.lbtoolarge)then
               if(myid.eq.0)write(*,*)'Large field drift motion',Bt,dt
               lbtoolarge=.true.
            endif
            do j=1,ndims            
! E-kick. Why do we do this? And why not in driftparticle?
               x_part(j+ndims,i)=x_part(j+ndims,i)+eoverms(ispecies)
     $              *Efield(j)*dtaccel
            enddo
            call driftparticle(x_part(1,i),ndims,Bt
     $           ,Efield,Bfield,vperp,dtpos,i-iicparta(ispecies))
         endif
         dtdone=dtdone+dtpos
! ---------------- End Particle Advance -------------
! ---------------- Special Conditions -------------
         if(lcollided)then
! Treat collided particle at (partial) step end
            call postcollide(x_part(1,i),tisq,iregion)
            lcollided=.false.
         endif
! Find where we've ended.
         call partlocate(x_part(1,i),ixp,xfrac,inewregion,linmesh)
!---------------------------------
! If we crossed an object boundary, do tallying.
         ltlyerr=.false.
! Discovery of all the objects this step crosses.
         icross=icrossall(xprior,x_part(1,i))
         leftregion=0
         if(icross.ne.0)then
! Since we crossed objects we might have leftregion:
            leftregion=leaveregion(ibool_part,ndims,xprior,x_part(1,i)
     $           ,icross,fractions)
! Here we tell tallyexit not to go past fractions(1) because that is the
! end of this step.
            call tallyexit(xprior,x_part(1,i),icross,ltlyerr,ispecies
     $           ,fractions(1))
         endif
! Diagnostic for leakage. Remove when convinced it is fixed:
         if(.not.linregion(ibool_part,ndims,x_part(1,i))
     $        .and..not.leftregion.ne.0)then
            write(*,'(a,i8,/,6f8.3,2i3,/,3f8.3,i3)')'Particle leakage',i
     $           ,(x_part(kk,i),xprior(kk),kk=1,3),icross,inewregion
     $           ,(x_part(kk,i)*dtpos,kk=4,6)
            call world3contra(ndims,xprior,xn1,1)
            call world3contra(ndims,x_part(1,i),xn2,1)
! This gave crashes in cluster running.
!            idbug=1
!            call srvsect(xn1,xn2,1,icross,f,ids)
!            call srvsectplot(1,xn1,xn2,f)
            idbug=0
            write(*,*)'icross=',icross,xn1,xn2
! Get rid of this leaked particle:
            leftregion=1
         endif
!------------ Possible Reinjection ----------
         if(.not.linmesh.or.leftregion.ne.0)then
!  We left the mesh or region. Test if this was a trapped passthrough.
           if(linmesh.and.inewregion.eq.iregion)npassthrough
     $           =npassthrough+1
            if(ninjcompa(ispecies).eq.0
     $           .or. nrein.lt.ninjcompa(ispecies)+ninjadd)goto 200
! Reinject because we haven't exhausted complement. 
! Else empty slot and move to next particle.
            x_part(iflag,i)=0
            goto 400
         endif
!--------------------------------------------
         x_part(idtp,i)=dtpos
         if(dtremain.gt.0.)then
! We haven't completed this full step; do so.
            if(dtc.ne.dtpos)then
! We did not have a collision, we are purely subcycling, just double the
! timestep for next subcycle. So as not to be over-optimistic.
               dtpos=min(2.*dtpos,dtremain)
!               if(i.lt.norbits)write(*,'(a,f10.5,$)')' Double',dtpos
            else
! If collision, default next step size to a full remaining step.
               dtpos=dtremain
            endif
! dtremain is the time remaining after the next dtpos step.
            dtremain=dtremain-dtpos
            iregion=inewregion
            goto 100
         endif
!-------------------------------------------
! The standard exit point for a particle that is active.
         iocthis=i
         nparta(ispecies)=nparta(ispecies)+1
         goto 300
!================= End of Occupied Slot Treatement ================
!----------- Reinjection treatment -----------
 200     continue
         x_part(iflag,i)=1
!         write(*,*)'Calling reinject',i,ispecies
         call reinject(x_part(1,i),ilaunch,ispecies)
         call partlocate(x_part(1,i),ixp,xfrac,iregion,linmesh)
         if(.not.linmesh)then
            write(*,*)'Reinject out of mesh',i,xfrac
            stop
         endif
         if(.not.linregion(ibool_part,ndims,x_part(1,i)))then
! This situation is benign and not an error if we have a region that
! happens not to cover the entire mesh edge. So don't stop, retry.
!            write(*,*)'Reinject out of region',i,iregion,xfrac
            if(ninjcompa(ispecies).eq.0
     $           .or. nrein.lt.ninjcompa(ispecies)+ninjadd)then
               nrein=nrein+ilaunch
               goto 200
            else
! We've exhausted the injection complement, so don't keep trying.
               x_part(iflag,i)=0
               goto 400
            endif
         endif
         call ranlux(ra,1)
! Indicate that this reinjected particle still needs to be moved a
! distance corresponding to a random fraction of the timestep. 
! That will happen on restart of the rest of the advance.
         dtpos=(dtpos+dtremain)*ra
         x_part(idtp,i)=0.
         dtremain=0.
         dtdone=dt-dtpos
         nlost=nlost+1
         nrein=nrein+ilaunch
         phi=getpotential(u,cij,iLs,x_part(2*ndims+1,i)
     $        ,IAND(iregion,ifield_mask),2)
         phirein=phirein*(1+ilaunch-1)+phi
         call diaginject(x_part(1,i))
! Restart the rest of the advance
         goto 100
!----------- End Reinjection treatment --------
!========================================================
!--------- Completion of particle i treatment ------
 300     continue
! Special diagnostic orbit tracking:
         ii=i-iicparta(ispecies)+1
         if(ii.le.norbits .and.x_part(iflag,i).ne.0)then
            if(dtdone-dt.gt.1.e-4)then
               write(*,*)' Step error? i,dtdone,dt=',i,dtdone,dt
            endif
            iorbitlens(ii,ispecies)=iorbitlens(ii,ispecies)+1
            xorbits(iorbitlens(ii,ispecies),ii,ispecies)=x_part(1,i)
            yorbits(iorbitlens(ii,ispecies),ii,ispecies)=x_part(2,i)
            zorbits(iorbitlens(ii,ispecies),ii,ispecies)=x_part(3,i)
         endif
! This is the charge deposition. Diagnostics are put into the 
! appropriate species place of diagsum. Charge magnitude is 1.
         if(x_part(iflag,i).ne.0)then
            call achargetomesh(i,psum,iLs
     $           ,diagsum(1+iLs(ndims+1)*(ndiagmax+1)*(ispecies-1))
     $           ,ndiags,echarge,xprior)
!            ndeposits=ndeposits+1
         endif
! End of icycle step subloop.
         enddo
!         if(ndeposits.ne.n_part)write(*,*)'ndeposits,nparta',ndeposits
!     $        ,n_part
! End of this particle's advance.
 400     continue
      enddo
      write(*,*)'nrein=',nrein
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 102  continue
! -------------- End of cycle through particles ----------
      if(ispecies.eq.1.and.rhoinf*colntime.ne.0.)then
! Normalize the collision force appropriately.
         call bulknorm(1./(rhoinf*dt))
      endif
      if(ninjcompa(ispecies).ne.0.and.nrein.lt.ninjcompa(ispecies)
     $     +ninjadd)then
         write(*,*)'WARNING: Exhausted n_partmax=',n_partmax,
     $        '  before ninjcomp=',ninjcompa(ispecies)
     $        ,' . Increase n_partmax?'
      endif
      iocparta(ispecies)=iocthis
!      write(*,'(a,5i8)')'  i,iocthis,npart,nrein=',i,iocthis
!     $     ,nparta(ispecies),nrein

! Finished this particle step. Calculate the average reinjection 
! potential
      if(nrein.gt.0)then
         phirein=phirein/nrein
!         if(myid.eq.0)
!         write(*,'(a,f12.6,$)')' Phirein=',phirein
!         write(*,*)' nlost=',nlost,' nrein=',nrein,' ninner=',ninner
         cap=2.*Ti
         if(phirein.gt.cap)then
            if(myid.eq.0)write(*,*)'PROBLEM: phirein>0:',cap,phirein
            phirein=cap
         endif
      else
         if(ninjcompa(ispecies).gt.100)write(*,*)'No reinjections'
      endif

      if(ispecies.eq.1.and.colntime.gt.0)then
         fcollided=float(ncollided)/nparta(ispecies)
      endif

!      if(nsubc.ne.0) write(*,'(a,i6,$,'' '')')' Subcycled:',nsubc
!      write(*,*)'Padvnc',n_part,nrein,ilaunch,ninjcomp,n_partmax
!     $     ,ispecies
      end
!***********************************************************************
!***********************************************************************
      subroutine getadfield(irptch,adfield,xp,isw,ierr)
! Cycle through the nonzero bits of irptch and add up the extra
! potential (isw=1), field (isw=2) or charge (isw=3) contributions.
! adfield should be zero on entry, because it ain't set,
! just incremented.
! ierr returns zero for no error, 1 for too close to particle.

      integer irptch,isw
      include 'ndimsdecl.f'
      real adfield(ndims)
      real xp(ndims)
      include '3dcom.f'
      real xd(ndims)
      im=irptch
      do i=1,31
         if(im.eq.0)return
         if(mod(im,2).ne.0)then
! This bit set. Add analytic field. Unequal radii not allowed.
            p2=0.
            xr=obj_geom(oradius,i)
            do id=1,ndims
               xc=xp(id)-obj_geom(ocenter+id-1,i)
               p2=p2+xc**2
               xd(id)=xc
            enddo
            if(.not.p2.ge.1.e-12)then
! Avoid overflows. Ought to happen only v rarely:
               write(*,*)'Correcting ptch field overflow rp^2=',p2
               p2=1.e-12
               ierr=1
            endif
            rp=sqrt(p2)
            p1=rp/xr
            xr2=xr**2
            p2=p2/xr2
! (There are inconsistencies if radii are not equal. Not allowed.)
!            write(*,*)irptch,isw,xp,xr,p1,p2
            if(isw.eq.2)then
               tfield=(obj_geom(omag,i)/xr)*(1./p2-4.*p1+3.*p2)
!               write(*,*)tfield,obj_geom(omag,i),xr
               do id=1,ndims
                  adfield(id)=adfield(id)+
     $                 tfield*(xd(id)/rp)
               enddo
            elseif(isw.eq.1)then
               tpotl=obj_geom(omag,i)*(1/p1+2.*p2-p1*p2-2.)
               adfield(1)=adfield(1)+tpotl
            elseif(isw.eq.3)then
               tchg=(obj_geom(omag,i)/xr2)*12.*(1.-p1)
               adfield(1)=adfield(1)+tchg
            else
               write(*,*)'getadfield switch error'
               stop
            endif
         endif
         im=im/2
      enddo
      
      end
!***********************************************************************
      subroutine setadfield(ifull,iuds,irptch,lsliceplot)
! Set the values of the potential and charge at the grid
! points that compensate for the analytic field of getadfield.
! Also any temperature-gradient and density-gradient factors.
      implicit none
      logical lsliceplot
! Defines iptch_mask uci, rhoc and dimensions.
! ndims must be same as ndimsp.
! To do the slice plot we need this it include grid decl.
      include 'ndimsdecl.f'
      integer ifull(ndims),iuds(ndims),irptch
      include 'meshcom.f'
      include 'ptchcom.f'
      include 'plascom.f'
! Here a segfault was caused if na_m2 was used.
      real zp(na_m,na_m)
      integer ipoint,ifix
      real dum
      external ucrhoset
      iptch_mask=irptch
      gtt_copy=gtt
      gnt_copy=gnt
!      write(*,*)'Point charges included. Mask:',iptch_mask,gtt
      ipoint=0
      ifix=1
      call mditerarg(ucrhoset,ndims,ifull,iuds,ipoint)
!     $     ,uci,rhoci,iptch_mask,Teci,boltzwt)
      if(lsliceplot)then
         call sliceGweb(ifull,iuds,uci,na_m,zp,
     $        ixnp,xn,ifix,'u!dc!d ptch',dum,dum)
         call sliceGweb(ifull,iuds,rhoci,na_m,zp,
     $        ixnp,xn,ifix,'!Ar!@!dc!d ptch',dum,dum)
        call sliceGweb(ifull,iuds,Teci,na_m,zp,
     $        ixnp,xn,ifix,'T!de!d Profile',dum,dum)
      endif
      end
!**********************************************************************
      subroutine ucrhoset(inc,ipoint,indi,mdims,iLs,iuds)
! This routine passes no mditerarg arguments only uses commons.
! Set uci, rhoci, and Teci arrays to compensate for point charges,
! electron temperature gradients, or density gradients.
! These are then subsequently used in faddu to decide the electron
! density. [They are not used in getadfield.]
      integer inc,ipoint,mdims
      integer iuds(mdims),indi(mdims)
! Commons: For position.
      include 'ndimsdecl.f'
      include 'meshcom.f'
! For debyelen and Tempr gradient.
      include 'plascom.f'
! For deciding the region via iboolpart for setting boltzwt.
      include '3dcom.f'
      include 'partcom.f'
      include 'ptchcom.f'
      logical linregion
      external linregion
! Local storage:
      integer isw,iregion,irptch
      real xp(ndims),adfield(ndims)
! Silence warnings
      irptch=iLs
      irptch=iuds(1)
!      write(*,*)'ucrhoset',ipoint,indi,iuds,ndims
!      if(ipoint.gt.100)stop
! Get grid point position, and irptch.
      do id=1,ndims
         xp(id)=xn(ixnp(id)+1+indi(id))
      enddo
      iregion=insideall(ndims,xp)
      irptch=IAND(iregion,iptch_mask)
! Set boltzwt depending on whether we are in the region or not.
      if(linregion(ibool_part,ndims,xp))then
         boltzwti(ipoint+1)=1.
      else
         boltzwti(ipoint+1)=0.
      endif
! Set the Teci factors.
      Teci(ipoint+1)=1.
      if(gtt.ne.0)then
         do id=1,ndims
            Teci(ipoint+1)=Teci(ipoint+1)+(xp(id)-gp0(id))*gt(id)
         enddo
         if(Teci(ipoint+1).le.0.)then
            write(*,*)'Te gradient too large.',
     $           ' Non-positive Te encountered at',xp
            stop
         endif
      endif
      rhoci(ipoint+1)=0.
      uci(ipoint+1)=0.
! Set the density-gradient compensating potential
      if(gnt.ne.0)then
         do id=1,ndims
            uci(ipoint+1)=uci(ipoint+1)+(xp(id)-gp0(id))*gn(id)
     $           *Teci(ipoint+1)
         enddo
      endif
! Set the uci and rhoci compensation for point charges:
      if(irptch.ne.0)then
! Get uc
         isw=1
         adfield(1)=0.
         call getadfield(irptch,adfield,xp,isw,ierr)
         uci(ipoint+1)=uci(ipoint+1)+adfield(1)
! Get charge
         isw=3
         adfield(1)=0.
         call getadfield(irptch,adfield,xp,isw,ierr)
         rhoci(ipoint+1)=rhoci(ipoint+1)+debyelen**2*adfield(1)
      endif
! Always just increment by 1
      inc=1
!      write(*,*)'ucrhoset return',irptch
      end

!***********************************************************************
      subroutine partlocate(xi,ixp,xfrac,iregion,linmesh)
! Locate the particle xi in the mesh (from common meshcom).
! Return the integer cell-base coordinates in ixp(ndims)
! Return the fractions of cell width at which located in xfrac(ndims)
! Return the region identifier in iregion.
! Return whether the particle is in the mesh in linmesh.
! Store the mesh position into upper ndims of xi.
! If particle is relocated by periodicity, advance nrein.

      integer k,iregion
      logical linmesh
! meshcom provides ixnp, xn, the mesh spacings. (+ndims)
      include 'ndimsdecl.f'
      include 'meshcom.f'
      include 'creincom.f'
      real xi(3*ndims)
      integer ixp(ndims)
      real xfrac(ndims)
      parameter (ndimsx2=ndims*2)
      include 'partcom.f'
      logical ldebug
      parameter (ldebug=.false.)

      linmesh=.true.
      do id=1,ndims
! Offset to start of dimension-id-position array.
         ioff=ixnp(id)
! xn is the position array for each dimension arranged linearly.
         isz=ixnp(id+1)-ioff         
! --------------------------------------------------------------
         if(ipartperiod(id).eq.4)then 
! Pure periodic particles. -----------------------------------
            if(xi(id).le.xmeshstart(id))then ! relocate near xmeshend
! Move the particle by one grid length, to the periodic position. Use
! tiny bit less so that if it starts exactly on boundary, it does not
! end on it. The length is between end mid-cell positions. That has
! been accounted for in xmesh variables.
               xi(id)=xi(id)+(xmeshend(id)-xmeshstart(id))*.999999
            elseif(xi(id).ge.xmeshend(id))then ! relocate near xmeshstart
               xi(id)=xi(id)-(xmeshend(id)-xmeshstart(id))*.999999
            endif
            ix=interp(xn(ioff+1),isz,xi(id),xm)
            if(.not.(xi(id).gt.xmeshstart(id).and.
     $           xi(id).lt.xmeshend(id)).or.ix.eq.0)then
! It's conceivable that a particle might move more than one period,
! in which case correction won't work. Don't repeat. Instead, just
! give up and call it lost but announce the problem. 
               write(*,'(a,2i2,3f10.4)')
     $              ' Particle crosses >meshlength',id,ix
     $              ,(xi(k),k=1,ndims)
               linmesh=.false.
               goto 2
            else
! If every dimension is periodic, increment nrein. (Otherwise not)
               if(.not.lnotallp)nrein=nrein+1
            endif
         elseif(ipartperiod(id).eq.5)then 
! Periodic vreset particles -----------------------------
            index=0
            if(xi(id).lt.xmeshstart(id))then ! relocate near xmeshend
               index=2*id-1     ! index odd for velocity reset
               xbdy=xmeshend(id)
            elseif(xi(id).ge.xmeshend(id))then ! relocate near xmeshstart
               index=2*id       ! index even for velocity reset
               xbdy=xmeshstart(id)
            elseif(xi(id).eq.xmeshstart(id))then
               xi(id)=xi(id)+1.e-6*(xmeshend(id)-xmeshstart(id))
            endif
            if(index.ne.0)then  !choose new velocity and relocate
               call ranlux(ra,1)
               ra=ra*ncrein
               ir=int(ra)
               fr=ra-ir
               vold=xi(ndims+id)
               xi(ndims+id)=hreins(ir,index,reinspecies)
     $              *(1-fr)+hreins(ir+1,index,reinspecies)*fr
! Relocate to lie in the periodic domain, with rounding correction.
               xt=xi(id)-xmeshstart(id)
               xmod=(xmeshend(id)-xmeshstart(id))*.999998
               xi(id)=xmeshstart(id)+modulo(xt,xmod)+1e-6*xmod
! Correct distance inside mesh for vold vs vnew difference.
               vfac=xi(ndims+id)/vold
               if(vfac.gt.0.001.and.vfac.lt.10.)then
                  xi(id)=xbdy+(xi(id)-xbdy)*vfac
               endif
            endif
            ix=interp(xn(ioff+1),isz,xi(id),xm)
! The following trap ought not to be necessary if the rounding correction
! works, but it might still be needed during debugging.
            if(.not.(xi(id).gt.xmeshstart(id).and.
     $           xi(id).lt.xmeshend(id)).or.ix.eq.0)then
               write(*,*)'Particle outside meshlen'
               write(*,*)'id, ix, xi(1  ... 3)    xbdy,xt,xmod,vold,xi'
               write(*,'(2i2,10f9.4)')id,ix
     $              ,(xi(k),k=1,ndims),xbdy,xt,xmod,vold,xi(ndims+id)
               linmesh=.false.
               goto 2
            else
               if(.not.lnotallp)nrein=nrein+1
            endif
         else
! Non-periodic Default --------------------
! Find the index of xprime in the array xn:
            ipi=int((xi(id)-xmeshstart(id))/(xmeshend(id)
     $           -xmeshstart(id))*(ipilen-1.00001)+1)
            if(ipi.lt.1)then
! Outside mesh
               ix=0
               xm=0
! was ipilen, but that seemed to be an error, gave index overrun below.
            elseif(ipi.gt.ipilen-1)then
               ix=0
               xm=isz+1
            else 
! Find the position index
               ix=iposindex(ipi,id)
               if(iposindex(ipi+1,id).eq.ix)then
                  xm=(xi(id)-xn(ioff+ix))/(xn(ioff+ix+1)-xn(ioff+ix))+ix
               else
! Default call returns ix and xm
! Interp costs here were 18% of particle intensive runs.
                  ix=interp(xn(ioff+1),isz,xi(id),xm)
               endif
            endif
! Set linmesh. New sensible version using xmesh
            if(.not.(ix.ne.0.and.xi(id).gt.xmeshstart(id).and.
     $           xi(id).lt.xmeshend(id)))then
               linmesh=.false.
               goto 2
            endif
         endif
! --------------------------------------------------------------
         xfrac(id)=xm-ix
         xi(ndimsx2+id)=xm
         ixp(id)=ix
! specific particle test
!         if(i.eq.2298489) write(*,*)i,isz,ix,xm,linmesh
      enddo
 2    continue
      iregion=insideall(ndims,xi(1))
      end
!********************************************************************
      subroutine postcollide(xi,tisq,iregion)
! Get new velocity; reflects neutral maxwellian shifted by vneutral.
! Update the collisional force by momentum change.
      real tisq
      integer iregion
      include 'ndimsdecl.f'
      real xi(2*ndims)
!      include 'partcom.f'
      include 'colncom.f'
      include '3dcom.f'
! Local variables
      integer icount
      real dv(ndims)
      data icount/0/

      do k=1,ndims
         dv(k)=-xi(ndims+k)
         xi(ndims+k)=tisq*gasdev(0)
         dv(k)=dv(k)+xi(ndims+k)
      enddo
! Vneutral is in z-direction.
      xi(ndims+ndims)=xi(ndims+ndims)+vneutral
      dv(ndims)=dv(ndims)+vneutral
! Now contribute the momentum change to the collisional force.
      do j=1,mf_obj
         if(btest(iregion,nf_geommap(j)-1))then
! Inside object nf_geommap(j) which is a flux measuring object.
!            write(*,*)icount,i,j,'  dv=',dv,dtaccel
            icount=icount+1
            do id=1,ndims
! This force accounting should be subsequently bulknormed by rhoinf
! and by dt.
                  colnforce(id,j,nf_step)=colnforce(id,j,nf_step)+dv(id)
            enddo
         endif
      enddo
      end
!*******************************************************************
      subroutine rotate3(xin,s,c,u)
! Rotate the input 3-vector xin, by angle theta whose sine and cosine
! are inputs about the direction given by direction cosines u
! (normalized axis vector) and return in xin.  This is a clockwise
! (right handed) rotation about positive u, I hope.

      real xin(3),u(3)
      real x1,x2
!      s=sin(theta)
!      c=cos(theta)
      d=1.-c
!      if(d.lt.1.e-5)d=s**2*0.5
! Just written out is about twice as fast:
      x1=(c+d*u(1)*u(1))*xin(1)
     $     +(d*u(1)*u(2)-s*u(3))*xin(2)
     $     +(d*u(1)*u(3)+s*u(2))*xin(3)
      x2=(c+d*u(2)*u(2))*xin(2)
     $     +(d*u(2)*u(3)-s*u(1))*xin(3)
     $     +(d*u(2)*u(1)+s*u(3))*xin(1)
      xin(3)=(c+d*u(3)*u(3))*xin(3)
     $     +(d*u(3)*u(1)-s*u(2))*xin(1)
     $     +(d*u(3)*u(2)+s*u(1))*xin(2)
! This shuffle is necessary because xin and xout are the same storage.
      xin(1)=x1
      xin(2)=x2
! And it is not permitted to make xin and xout different dummy
! arguments, if some calls are done with them being the same
! variables. That violates the fortran standard which requires that no
! aliased locations be written to in dummy arguments.
! Therefore this call assumes they are always the same and avoids 
! any aliases.
      end
!********************************************************************
      subroutine translate3(xin,vin,dt,u,xout)
! Translate in the u direction (cosines) by the component of vin 
! in the u-direction times dt from the input position xin to xout.

      real xin(3),vin(3),u(3),xout(3),dt
      vudt=(u(1)*vin(1)+u(2)*vin(2)+u(3)*vin(3))*dt
      xout(1)=xin(1)+u(1)*vudt
      xout(2)=xin(2)+u(2)*vudt
      xout(3)=xin(3)+u(3)*vudt
      end
!*********************************************************************
      subroutine gyro3(Bt,u,xin,vin,xg,xc)
! Bt in this routine is multiplied by eoverm.
! Find gyro radius and center of particle at xin, velocity vin.

! Given a field Bt
      real Bt
! in direction u, 
      real u(3)
! obtain the gyro radius xg 
      real xg(3)
! of the particle at xin, 
      real xin(3)
! with velocity vin. 
      real vin(3)
! Subtract the gyro radius from xin to give the gyrocenter in xc.
      real xc(3)
! Find the perpendicular velocity, rotate it by -90 degrees, and divide 
! by the field. 
      vu=(u(1)*vin(1)+u(2)*vin(2)+u(3)*vin(3))
      xg(1)=(vin(1)-u(1)*vu)/Bt
      xg(2)=(vin(2)-u(2)*vu)/Bt
      xg(3)=(vin(3)-u(3)*vu)/Bt
!      theta=3.1415917*0.5
      s=1.
      c=0.
      call rotate3(xg,s,c,u)
! Now xg is the gyroradius.
      xc(1)=xin(1)-xg(1)
      xc(2)=xin(2)-xg(2)
      xc(3)=xin(3)-xg(3)
      end
!*****************************************************************
      subroutine bulkforce(xpart,iregion,colntime,vneutral,Eneutral)
! Add on the contribution of this particle to the collision force.
! colntime is the collision time, vneutral the drift velocity, dt time
! step.  The total force is the contained momentum difference over the
! collision time. Needs subsequently to be properly normalized.
      include 'ndimsdecl.f'
      include 'partcom.f'
      include '3dcom.f'
      real xpart(3*ndims),colntime,vneutral
      integer iregion
      if(colntime.gt.0)then
         do i=1,mf_obj
            igeom=nf_geommap(i)
!            if(iregion.ne.0)write(*,*)'iregion',iregion
            if(btest(iregion,igeom-1))then
! Inside object i->igeom which is a flux measuring object.
! The scalar drift velocity is in the z-direction.
               do id=1,ndims
                  vid=0.
                  if(id.eq.ndims)vid=vneutral+Eneutral*colntime
                  colnforce(id,i,nf_step)=colnforce(id,i,nf_step)+
     $                 (vid-xpart(ndims+id))
               enddo
!               write(*,*)'colnforce',i,colnforce(3,i,nf_step),vd
!     $              ,xpart(ndims+3),vid,iregion
            endif
         enddo
      endif
      end
!********************************************************************
      subroutine bulknorm(fac)
! Normalize the bulkforce, multiplying by factor fac
      include 'ndimsdecl.f'
      include '3dcom.f'
      do i=1,mf_obj
         do id=1,ndims
            colnforce(id,i,nf_step)=colnforce(id,i,nf_step)
     $           *fac
         enddo
      enddo
      end
! *******************************************************************
      subroutine moveparticle(xr,ndims,Bt,Efield,Bfield,vperp
     $     ,dtpos,dtaccel,driftfield,eom)
! Accelerate and move the passed particle in the fields B,E for timestep
! dtaccel and dtpos. 
      implicit none
      integer ndims
      real Bt,dtpos,dtaccel,eom
      real xr(2*ndims),Bfield(ndims),Efield(ndims)
      real vperp(ndims),driftfield(ndims)
! Local      
      integer j
      real theta,stheta,ctheta,thacc,sthacc,cthacc
      real xg(3),xc(3)
      real thetatoosmall
      parameter (thetatoosmall=1.e-3)
      real thetamax
      parameter (thetamax=1.)

      theta=eom*Bt*dtpos
      thacc=eom*Bt*dtaccel
      if(abs(theta).eq.0.)then
! B-field-less Kick - Move -----------------------------
         do j=1,ndims
            xr(j+ndims)=xr(j+ndims)+eom*Efield(j)*dtaccel !Kick
            xr(j)=xr(j)+xr(j+ndims)*dtpos                 !Drift
         enddo          
         return
      endif

! Only if magnetic field nonzero ---------------------------------
! Rotation is counterclockwise for ions. We only want to call the 
! trig functions once, otherwise they dominate the cost.
      if(abs(theta).lt.thetatoosmall)then
! Weak B-field. ------------------------------------
! Boris Mover in frame of reference moving with external drift
         do j=1,ndims
! First E half-kick
            xr(j+ndims)=xr(j+ndims)+eom*Efield(j)*dtaccel*0.5
! Change reference to drifting frame
            xr(j+ndims)=xr(j+ndims)-vperp(j)
         enddo
! Rotate the velocity to add the magnetic field acceleration.
! This amounts to a presumption that the magnetic field acceleration
! acts at the mid-point of the translation (drift) rather than at
! the kick between translations.
         sthacc=sin(-thacc)
         cthacc=cos(thacc)
         call rotate3(xr(4),sthacc,cthacc,Bfield)
! Second half-move.
         do j=1,ndims
! Transform back to undrifted frame
            xr(j+ndims)=xr(j+ndims)+vperp(j)
! Second E half-kick
            xr(j+ndims)=xr(j+ndims)+eom*Efield(j)*dtaccel*0.5
! Move particle
            xr(j)=xr(j)+xr(j+ndims)*dtpos
         enddo          
      else
! Strong but finite B-field case -------------------------
! Cyclotronic mover
         do j=1,ndims            
! E-kick
            xr(j+ndims)=xr(j+ndims)+eom*Efield(j)*dtaccel
! Subtract off the perpendicular drift velocity.
            xr(j+ndims)=xr(j+ndims)-vperp(j)
         enddo
! Find the gyro radius and gyrocenter.
         call gyro3(eom*Bt,Bfield,xr(1),xr(4),xg,xc)
! Rotate the velocity and gyro radius.
! Rotation is counterclockwise for ions.
! For unequal timesteps, rotate gyro center by theta, not thacc.
         stheta=sin(-theta)
         ctheta=cos(theta)
! Anything other than theta gives velocity strangeness at edge.
         call rotate3(xr(4),stheta,ctheta,Bfield)
         call rotate3(xg,stheta,ctheta,Bfield)
! Move xc along the B-direction.
         call translate3(xc,xr(4),dtpos,Bfield,xc)
! Add the new gyro center and gyro radius
! And add back the drift velocity.
         do j=1,ndims
! Move the gyro-center perpendicular and add gyro-radius:
            xr(j)=xc(j)+vperp(j)*dtpos+xg(j)
! Add back vperp.
            xr(j+ndims)=xr(j+ndims)+vperp(j)
         enddo
      endif

! End of Move -----------------------------------------------------
      end
!********************************************************************
      subroutine driftparticle(xr,ndims,Bt,Efield,Bfield,vperp
     $     ,dtpos,k)
! Drift particle instead of move for cases of large Bt.
      implicit none
      integer ndims,k
      real Bt,dtpos
      real xr(2*ndims),Bfield(ndims),Efield(ndims)
      real vperp(ndims)
      integer j
! Local
      real EB(3)
      real vp
      integer k1,k2

      do j=1,ndims
         k1=mod(j,ndims)+1
         k2=mod(j+1,ndims)+1
! Explicit ExB velocity plus vperp
         EB(j)=vperp(j)+(Efield(k1)*Bfield(k2)-Efield(k2)*Bfield(k1))/Bt
      enddo
! Drift motion parallel plus EB.
! Set non-drift perp particle velocity to zero.
! Bfield is normalized: i.e. a direction cosine, so parallel velocity
      vp=0.
      do j=1,ndims
! v_parallel
         vp=vp+Bfield(j)*xr(j+ndims)
      enddo
      do j=1,ndims
! Parallel particle and perpendicular drift move:
         xr(j)=xr(j)
     $        +(vp*Bfield(j)+EB(j))*dtpos
! Reset perp velocity to just the drift. Effectively we remove the 
! perpendicular energy forcing the magnetic moment to zero. This
! is necessary because the perpendicular acceleration has been applied
! improperly for the whole timestep.
         xr(j+ndims)=vp*Bfield(j)+EB(j)
      enddo

      end
!*********************************************************************
      subroutine initdriftfield
! Initialize the drift electric field = -vxB.  The reason this can be
! different for different species is that we might want to use it to
! represent, e.g. gravitational drift. At any rate, it is the electric
! field needed to give rise to the specified (ExB) drift for that species.
      include 'ndimsdecl.f'
      include 'plascom.f'
      include 'partcom.f'
! Local
      do j=1,nspecies
         do i=1,ndims
            driftfields(i,j)=-Bt*
     $           (vperps(mod(i,ndims)+1,j)*Bfield(mod(i+1,ndims)+1)
     $           -vperps(mod(i+1,ndims)+1,j)*Bfield(mod(i,ndims)+1))
         enddo
      enddo
      end
!***********************************************************************
      subroutine locateinit()
! Initialize particles loaded from restart files so that they have the
! correct mesh-fraction coordinate in their upper values.
! If they are outside the mesh or particle region, drop them.
! Common data:
      include 'ndimsdecl.f'
      include 'partcom.f'
      include 'meshcom.f'
! Local dummy variables for partlocate.
      real xfrac(ndims)
      integer ixp(ndims)
      logical linmesh,linregion
      external linregion

      do i=1,ioc_part
         call partlocate(x_part(1,i),ixp,xfrac,iregion,linmesh)
!         write(*,*)linmesh,linregion(ibool_part,ndims,x_part(1,i))
         if(.not.linmesh .or. .not.
     $        linregion(ibool_part,ndims,x_part(1,i)))then
            x_part(iflag,i)=0
         endif
      enddo
      nrein=0
!      write(*,*)'Redetermined particle mesh locations (locateinit)'
      end
