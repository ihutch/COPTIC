      subroutine padvnc(iLs,cij,u,ndiags,psum,diagsum,ispecies,ndiagmax)
c Particle advancing routine.
c If ninjcomp (in partcom) is non-zero, then we are operating in a mode
c where the number of reinjections at each timestep is prescribed.
c Otherwise we are using a fixed number n_part of particles.

c Storage size of the mesh arrays.
c      real cij(2*ndims+1,nx,ny,nz)
c      real u(nx,ny,nz)
      integer ndiags,ispecies
      real cij(*),u(*),psum(*)
      real diagsum(*)
c      real diagsum(64,64,128,8,2)

      include 'ndimsdecl.f'
c Array structure vector: (1,nx,nx*ny,nx*ny*nz)
      integer iLs(ndims+1)

c Meshcom provides ixnp, xn, the mesh spacings. (+ndims)
      include 'meshcom.f'
      include 'myidcom.f'
      include '3dcom.f'
      include 'ptchcom.f'
      external linregion
      logical linregion
      include 'plascom.f'
c Collision settings.
      include 'colncom.f'
      include 'facebcom.f'

c Local parameters
c Lowest particle to print debugging data for.
      integer npr
      parameter (npr=0)
      real fieldtoosmall
      parameter (fieldtoosmall=1.e-3)
      real thetamax
      parameter (thetamax=1.)
      real fieldtoolarge
      parameter (fieldtoolarge=1.e12)
c Local storage
      real fractions(10)
      real rannum
      integer ixp(ndims),ninjadd
      real Efield(ndims),adfield(ndims)
      real xfrac(ndims)
      real xprior(2*ndims)
      logical linmesh
      logical lcollided,ltlyerr
      save adfield

c Make this always last to use the checks.
      include 'partcom.f'

c-------------------------------------------------------------
      tisq=sqrt(Ts(ispecies))
      lcollided=.false.
      ncollided=0
      ic1=2*ndims+1
      ndimsx2=2*ndims
      theta=abs(dt*Bt*eoverms(ispecies))

c-----------------------------------------------------------------
c Initialize. Set reinjection potential. We start with zero reinjections.
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
c Set additional injection to zero by default.
      ninjadd=0
      if(ninjcompa(ispecies).ne.0)then
         call ranlux(rannum,1)
         if(rannum.lt.pinjcompa(ispecies))ninjadd=1
      endif
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c At most do over all particle slots for this species. 
c But generally we break earlier.
      numratio=numratioa(ispecies)
      do i=iicparta(ispecies),iicparta(ispecies+1)-1
         do icycle=1,numratio
         dtdone=0.
         dtremain=0.
         dtpos=dt/float(numratio)
c=================================
c Decide nature of this slot
         if(x_part(iflag,i).ne.0)then
c Standard occupied slot. Proceed.
         elseif(ninjcompa(ispecies).ne.0
     $        .and.nrein.lt.ninjcompa(ispecies)+ninjadd)then
c An unfilled slot needing to be filled. Jump to reinjection.
            goto 200
         elseif(i.ge.iocparta(ispecies))then
c We do not need to reinject new particles, and
c this slot is higher than all previously handled. There are no
c more active particles above it. So break from do loop.
            goto 102
         else
c Unoccupied but there might be slots occupied above so just cycle.
            goto 400
         endif
c================ Occupied Slot Treatment =================
c Occupied slot. Get its region
         iregion=insideall(ndims,x_part(1,i))
 100     continue
c--------------------------------------
c Check the fraction data is not stupid and complain if it is.
c Ought not to be necessary, but this is a safety check. 
c Remove after the reimplementation.
         if(x_part(ndimsx2+1,i).eq.0. .and.
     $        x_part(ndimsx2+2,i).eq.0. .and.
     $        x_part(ndimsx2+3,i).eq.0.) then
            write(*,*)'Zero fractions',i,iocparta(ispecies)
     $           ,x_part(iflag,i),nrein,ninjcompa(ispecies)
         endif
c---------------------------------
c adfield must always be zero arriving here regardless of irptch. 
c But we don't set it always because that would be expensive. 
c We only reset it below if it was changed
         irptch=IAND(iregion,iptch_mask)
         if(irptch.ne.0)then
c We are in a point-charge region. Get analytic part of force.
            call getadfield(irptch,adfield,x_part(1,i),2,ierr)
c            if(i.eq.1)write(*,'(i7,'' in ptch region'',i8,6f8.3)')
c     $           i,irptch,(x_part(k,i),k=1,3)
c     $           ,(adfield(k),k=1,3)
         endif
c Get the ndims field components at this point, from mesh potential.
c We only use x_part information for location. So we need to pass
c the region information.
         f2=0.
         r2=0.
         v2=0.
         do idf=1,ndims
c Test for need to get field:
            if(.not.(LPF(idf).and.ixnp(idf+1)-ixnp(idf).eq.3))then
c            if(.true.)then   
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
c If needed, reset adfield.
         if(irptch.ne.0)then
            do idf=1,ndims
               adfield(idf)=0.
            enddo
         endif
         f1=sqrt(f2)
c Example of testing code: the few-argument field evaluator.
c                  call fieldatpoint(x_part(1,i),u,cij,iLs,fieldp)
c                  if(fieldp(2).ne.Efield(2))then
c                     write(*,'(i5,a,6f10.6)')i,' Point-field:',
c     $                    fieldp,((fieldp(j)-Efield(j)),j=1,ndims)
c--------------------------------------
         dtc=0.
c---------- Subcycling ----------------
         if(subcycle.ne.0)then
            if(f1*dtpos.gt.subcycle .and. 8*dtpos.gt.dt)then
c Reduce the position-step size.
               dtpos=dtpos/2.
c Subcycling backs up to the start of the current drift-step, and then 
c takes a position corresponding to a drift-step smaller by a factor of 2.
c That means the backed-up position is dtprec/4 earlier. There might be
c an error if dtpost and dtaccel are different. Is -dtaccel correct?
               dtback=-x_part(idtp,i)/4.
               call moveparticle(x_part(1,i),ndims,Bt ,Efield,Bfield
     $              ,vperp,dtback,-dtaccel,driftfields(1,ispecies)
     $              ,eoverms(ispecies))
               dtdone=dtdone+dtback
c The remaining time in step is increased by back up and by step loss.
               dtremain=dtremain-dtback+dtpos
c And the new dtprec is half as large:
               x_part(idtp,i)=x_part(idtp,i)/2.
               nsubc=nsubc+1
               if(i.lt.norbits)then
                  write(*,'(/,a,i3,4f10.5,$)')'Subcycle Set',i,dtpos
     $                 ,x_part(idtp,i),dtremain,dtdone+dtpos+dtremain
               endif
               goto 100
            endif
         endif
c---------- End of subcycling ---------
c---------- Collision Decision ----------------
         if(ispecies.eq.1.and.colntime.ne.0.)then
            if(colpow.ne.0.)then
c The collision time scaled by velocity
c     $              *(Ti/(v2+Ti))**(colpow/2.)
               call ranlux(dtc,1)
               dtc=-alog(dtc)*colntime
     $              *(Ti/(v2+Ti))**(colpow/2.)
            else
c Time to the first collision
               call ranlux(dtc,1)
               dtc=-alog(dtc)*colntime
            endif
            if(dtc.lt.dtpos)then
c               write(*,*)'Collision',dtpos,dtc,colntime
c We collided during this step. Do the partial step.  
               lcollided=.true.
c Set dtremain to how much time will be left after this step.
               dtremain=dtpos+dtremain-dtc
               dtpos=dtc
               ncollided=ncollided+1
            endif
         endif 
c---------- End Collision Decision ------------
c         if(i.lt.norbits)write(*,*)x_part(idtp,i),dtremain,dtdone
c---------------- Particle Advance ----------------
c Use dtaccel for acceleration. May be different from dtpos if there was
c a subcycle, reinjection or collision last step.
         dtaccel=0.5*(dtpos+x_part(idtp,i))
c----- Excessive Acceleration test. Drop particle without any tally. ---
         if(f1*dtaccel.gt.dropaccel)then
            ndropped=ndropped+1
c Why don't we tally this particle? It effectively has been absorbed.
c It carried in some momentum. That disappears. It is accounted for
c as being conveyed to any tallying object inside which it disappears.
c            write(*,*)'Excessive acceleration. Dropping particle',i
c Reinject if we haven't exhausted complement:
            if(ninjcompa(ispecies).eq.0
     $           .or. nrein.lt.ninjcompa(ispecies)+ninjadd)goto 200
c Else empty slot and move to next particle.
            x_part(iflag,i)=0
            goto 400            
         endif
c Accelerate used to be here -----------------
c Here we only include the electric field force. B-field force and vxB
c field is accommodated within moveparticle routine.
         if(Eneutral.ne.0.)then
c Calculate collision force --------
            if(Bt.ne.0.)then
c Only the Enparallel is needed for non-zero Bfield.
c This doesn't seem to be correct unless drift is in z-direction.
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
c Including the extra Eneutral field, which is in the ndims direction.
c But correcting the collisional force appropriately.
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
c Move ----------------
c Save prior position and velocity.
         do id=1,2*ndims
            xprior(id)=x_part(id,i)
         enddo
         if(abs(theta).lt.thetamax)then
            call moveparticle(x_part(1,i),ndims,Bt ,Efield,Bfield,vperp
     $           ,dtpos,dtaccel,driftfields(1,ispecies)
     $           ,eoverms(ispecies))
         else
            do j=1,ndims            
c E-kick
               x_part(j+ndims,i)=x_part(j+ndims,i)+eoverms(ispecies)
     $              *Efield(j)*dtaccel
            enddo
            call driftparticle(x_part(1,i),ndims,Bt
     $           ,Efield,Bfield,vperp,dtpos,i-iicparta(ispecies))
         endif
         dtdone=dtdone+dtpos
c ---------------- End Particle Advance -------------
c ---------------- Special Conditions -------------
         if(lcollided)then
c Treat collided particle at (partial) step end
            call postcollide(x_part(1,i),tisq,iregion)
            lcollided=.false.
         endif
         call partlocate(x_part(1,i),ixp,xfrac,inewregion,linmesh)
c---------------------------------
c If we crossed an object boundary, do tallying.
         ltlyerr=.false.
c Discovery of all the objects this step crosses.
         icross=icrossall(xprior,x_part(1,i))
         leftregion=0
         if(icross.ne.0)then
c Since we crossed objects we might have leftregion:
            leftregion=leaveregion(ibool_part,ndims,xprior,x_part(1,i)
     $           ,icross,fractions)
c Here we tell tallyexit not to go past fractions(1) because that is the
c end of this step.
            call tallyexit(xprior,x_part(1,i),icross,ltlyerr,ispecies
     $           ,fractions(1))
         endif
c Diagnostic for leakage. Remove when convinced it is fixed:
         if(.not.linregion(ibool_part,ndims,x_part(1,i))
     $        .and..not.leftregion.ne.0)then
            write(*,'(a,i8,6f8.3,i4)')'Particle leakage',i,(x_part(kk,i)
     $           ,kk=1,6),icross
c Get rid of this leaked particle:
            leftregion=1
         endif
c------------ Possible Reinjection ----------
         if(.not.linmesh.or.leftregion.ne.0)then
c  We left the mesh or region. Test if this was a trapped passthrough.
           if(linmesh.and.inewregion.eq.iregion)npassthrough
     $           =npassthrough+1
            if(ninjcompa(ispecies).eq.0
     $           .or. nrein.lt.ninjcompa(ispecies)+ninjadd)goto 200
c Reinject because we haven't exhausted complement. 
c Else empty slot and move to next particle.
            x_part(iflag,i)=0
            goto 400
         endif
c--------------------------------------------
         x_part(idtp,i)=dtpos
         if(dtremain.gt.0.)then
c We haven't completed this full step; do so.
            if(dtc.ne.dtpos)then
c We did not have a collision, we are purely subcycling, just double the
c timestep for next subcycle. So as not to be over-optimistic.
               dtpos=min(2.*dtpos,dtremain)
c               if(i.lt.norbits)write(*,'(a,f10.5,$)')' Double',dtpos
            else
c If collision, default next step size to a full remaining step.
               dtpos=dtremain
            endif
c dtremain is the time remaining after the next dtpos step.
            dtremain=dtremain-dtpos
            iregion=inewregion
            goto 100
         endif
c-------------------------------------------
c The standard exit point for a particle that is active.
         iocthis=i
         nparta(ispecies)=nparta(ispecies)+1
         goto 300
c================= End of Occupied Slot Treatement ================
c----------- Reinjection treatment -----------
 200     continue
         x_part(iflag,i)=1
c         write(*,*)'Calling reinject',i,ispecies
         call reinject(x_part(1,i),ilaunch,ispecies)
         call partlocate(x_part(1,i),ixp,xfrac,iregion,linmesh)
c         if(ispecies.eq.2.and.abs(x_part(5,i)).gt..001
c     $        .and.i.lt.131000)then
c Testing
c            write(*,*)'reinject vynonzero',i
c            write(*,*)(x_part(j,i),j=1,6)
c         endif
         if(.not.linmesh)then
            write(*,*)'Reinject out of mesh',i,xfrac
            stop
         endif
         if(.not.linregion(ibool_part,ndims,x_part(1,i)))then
c This situation is benign and not an error if we have a region that
c happens not to cover the entire mesh edge. So don't stop, retry.
c            write(*,*)'Reinject out of region',i,iregion,xfrac
            if(ninjcompa(ispecies).eq.0
     $           .or. nrein.lt.ninjcompa(ispecies)+ninjadd)then
               nrein=nrein+ilaunch
               goto 200
            else
c We've exhausted the injection complement, so don't keep trying.
               x_part(iflag,i)=0
               goto 400
            endif
         endif
         call ranlux(ra,1)
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
c Restart the rest of the advance
         goto 100
c----------- End Reinjection treatment --------
c========================================================
c--------- Completion of particle i treatment ------
 300     continue
c Special diagnostic orbit tracking:
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
c This is the charge deposition. Diagnostics are put into the 
c appropriate species place of diagsum. Charge magnitude is 1.
         if(x_part(iflag,i).ne.0)then
            call achargetomesh(i,psum,iLs
     $           ,diagsum(1+iLs(ndims+1)*(ndiagmax+1)*(ispecies-1))
     $           ,ndiags,echarge,xprior)
c            if(mod(i,500).eq.0)write(*,*)'achargetomesh',ispecies,i,iLs
c     $           ,ndiags,ndiagmax
c     $           ,diagsum(1+iLs(ndims+1)*(ndiagmax+1)*(ispecies-1)+1)
         endif
c End of icycle step subloop.
         enddo
c End of this particle's advance.
c         if(ispecies.eq.2.and.abs(x_part(ndims+2,i)).gt.0.001
c     $        .and.i.lt.131000)then
c Test if we have a non-zero vy.
c            write(*,*)'vynonzero',i
c            write(*,*)(x_part(k,i),k=1,6)
c         endif
 400     continue
      enddo
      write(*,*)'nrein=',nrein
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 102  continue
c -------------- End of cycle through particles ----------
      if(ispecies.eq.1.and.rhoinf*colntime.ne.0.)then
c Normalize the collision force appropriately.
         call bulknorm(1./(rhoinf*dt))
      endif
      if(ninjcompa(ispecies).ne.0.and.nrein.lt.ninjcompa(ispecies)
     $     +ninjadd)then
         write(*,*)'WARNING: Exhausted n_partmax=',n_partmax,
     $        '  before ninjcomp=',ninjcompa(ispecies)
     $        ,' . Increase n_partmax?'
      endif
      iocparta(ispecies)=iocthis
c      write(*,'(a,5i8)')'  i,iocthis,npart,nrein=',i,iocthis
c     $     ,nparta(ispecies),nrein

c Finished this particle step. Calculate the average reinjection 
c potential
      if(nrein.gt.0)then
         phirein=phirein/nrein
c         if(myid.eq.0)
c         write(*,'(a,f12.6,$)')' Phirein=',phirein
c         write(*,*)' nlost=',nlost,' nrein=',nrein,' ninner=',ninner
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

c      if(nsubc.ne.0) write(*,'(a,i6,$,'' '')')' Subcycled:',nsubc
c      write(*,*)'Padvnc',n_part,nrein,ilaunch,ninjcomp,n_partmax
c     $     ,ispecies
      end
c***********************************************************************
c***********************************************************************
      subroutine getadfield(irptch,adfield,xp,isw,ierr)
c Cycle through the nonzero bits of irptch and add up the extra
c potential (isw=1), field (isw=2) or charge (isw=3) contributions.
c adfield should be zero on entry, because it ain't set,
c just incremented.
c ierr returns zero for no error, 1 for too close to particle.

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
c This bit set. Add analytic field. Unequal radii not allowed.
            p2=0.
            xr=obj_geom(oradius,i)
            do id=1,ndims
               xc=xp(id)-obj_geom(ocenter+id-1,i)
               p2=p2+xc**2
               xd(id)=xc
            enddo
            if(.not.p2.ge.1.e-12)then
c Avoid overflows. Ought to happen only v rarely:
               write(*,*)'Correcting ptch field overflow rp^2=',p2
               p2=1.e-12
               ierr=1
            endif
            rp=sqrt(p2)
            p1=rp/xr
            xr2=xr**2
            p2=p2/xr2
c (There are inconsistencies if radii are not equal. Not allowed.)
c            write(*,*)irptch,isw,xp,xr,p1,p2
            if(isw.eq.2)then
               tfield=(obj_geom(omag,i)/xr)*(1./p2-4.*p1+3.*p2)
c               write(*,*)tfield,obj_geom(omag,i),xr
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
c***********************************************************************
      subroutine setadfield(ifull,iuds,irptch,lsliceplot)
c Set the values of the potential and charge at the grid
c points that compensate for the analytic field of getadfield.
c Also any temperature-gradient and density-gradient factors.
      implicit none
      logical lsliceplot
c Defines iptch_mask uci, rhoc and dimensions.
c ndims must be same as ndimsp.
c To do the slice plot we need this it include grid decl.
      include 'ndimsdecl.f'
      integer ifull(ndims),iuds(ndims),irptch
      include 'meshcom.f'
      include 'ptchcom.f'
      include 'plascom.f'
c Here a segfault was caused if na_m2 was used.
      real zp(na_m,na_m)
      integer ipoint,ifix
      real dum
      external ucrhoset
      iptch_mask=irptch
      gtt_copy=gtt
      gnt_copy=gnt
c      write(*,*)'Point charges included. Mask:',iptch_mask,gtt
      ipoint=0
      ifix=1
      call mditerarg(ucrhoset,ndims,ifull,iuds,ipoint)
c     $     ,uci,rhoci,iptch_mask,Teci,boltzwt)
      if(lsliceplot)then
         call sliceGweb(ifull,iuds,uci,na_m,zp,
     $        ixnp,xn,ifix,'u!dc!d ptch',dum,dum)
         call sliceGweb(ifull,iuds,rhoci,na_m,zp,
     $        ixnp,xn,ifix,'!Ar!@!dc!d ptch',dum,dum)
        call sliceGweb(ifull,iuds,Teci,na_m,zp,
     $        ixnp,xn,ifix,'T!de!d Profile',dum,dum)
      endif
      end
c**********************************************************************
      subroutine ucrhoset(inc,ipoint,indi,mdims,iLs,iuds)
c This routine passes no mditerarg arguments only uses commons.
c Set uci, rhoci, and Teci arrays to compensate for point charges,
c electron temperature gradients, or density gradients.
c These are then subsequently used in faddu to decide the electron
c density. [They are not used in getadfield.]
      integer inc,ipoint,mdims
      integer iuds(mdims),indi(mdims)
c Commons: For position.
      include 'ndimsdecl.f'
      include 'meshcom.f'
c For debyelen and Tempr gradient.
      include 'plascom.f'
c For deciding the region via iboolpart for setting boltzwt.
      include '3dcom.f'
      include 'partcom.f'
      include 'ptchcom.f'
      logical linregion
      external linregion
c Local storage:
      integer isw,iregion,irptch
      real xp(ndims),adfield(ndims)
c Silence warnings
      irptch=iLs
      irptch=iuds(1)
c      write(*,*)'ucrhoset',ipoint,indi,iuds,ndims
c      if(ipoint.gt.100)stop
c Get grid point position, and irptch.
      do id=1,ndims
         xp(id)=xn(ixnp(id)+1+indi(id))
      enddo
      iregion=insideall(ndims,xp)
      irptch=IAND(iregion,iptch_mask)
c Set boltzwt depending on whether we are in the region or not.
      if(linregion(ibool_part,ndims,xp))then
         boltzwti(ipoint+1)=1.
      else
         boltzwti(ipoint+1)=0.
      endif
c Set the Teci factors.
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
c Set the density-gradient compensating potential
      if(gnt.ne.0)then
         do id=1,ndims
            uci(ipoint+1)=uci(ipoint+1)+(xp(id)-gp0(id))*gn(id)
     $           *Teci(ipoint+1)
         enddo
      endif
c Set the uci and rhoci compensation for point charges:
      if(irptch.ne.0)then
c Get uc
         isw=1
         adfield(1)=0.
         call getadfield(irptch,adfield,xp,isw,ierr)
         uci(ipoint+1)=uci(ipoint+1)+adfield(1)
c Get charge
         isw=3
         adfield(1)=0.
         call getadfield(irptch,adfield,xp,isw,ierr)
         rhoci(ipoint+1)=rhoci(ipoint+1)+debyelen**2*adfield(1)
      endif
c Always just increment by 1
      inc=1
c      write(*,*)'ucrhoset return',irptch
      end

c***********************************************************************
      subroutine partlocate(xi,ixp,xfrac,iregion,linmesh)
c Locate the particle xi in the mesh (from common meshcom).
c Return the integer cell-base coordinates in ixp(ndims)
c Return the fractions of cell width at which located in xfrac(ndims)
c Return the region identifier in iregion.
c Return whether the particle is in the mesh in linmesh.
c Store the mesh position into upper ndims of xi.
c If particle is relocated by periodicity, advance nrein.

      integer k,iregion
      logical linmesh
c meshcom provides ixnp, xn, the mesh spacings. (+ndims)
      include 'ndimsdecl.f'
      include 'meshcom.f'
      real xi(3*ndims)
      integer ixp(ndims)
      real xfrac(ndims)
      parameter (ndimsx2=ndims*2)
      include 'partcom.f'

      linmesh=.true.
      do id=1,ndims
c Offset to start of dimension-id-position array.
         ioff=ixnp(id)
c xn is the position array for each dimension arranged linearly.
c Find the index of xprime in the array xn:
         isz=ixnp(id+1)-ioff
         
c An inverse lookup table for non-uniform. --------------------
         ipi=int((xi(id)-xmeshstart(id))/(xmeshend(id)-xmeshstart(id))
     $        *(ipilen-1.00001)+1)
         if(ipi.lt.1)then
c Outside mesh
            ix=0
            xm=0
c was ipilen, but that seemed to be an error, gave index overrun below.
         elseif(ipi.gt.ipilen-1)then
            ix=0
            xm=isz+1
         else
            ix=iposindex(ipi,id)
            if(iposindex(ipi+1,id).eq.ix)then
               xm=(xi(id)-xn(ioff+ix))/(xn(ioff+ix+1)-xn(ioff+ix))+ix
            else
               ix=interp(xn(ioff+1),isz,xi(id),xm)
            endif
         endif
c --------------------------------------------------------------
c Interp costs here were 18% of particle intensive runs.
c         ix=interp(xn(ioff+1),isz,xi(id),xm)
         if(ipartperiod(id).eq.4)then
c In periodic directions, we do not allow particles to be closer to the
c mesh boundary than half a cell, so as to use periodicity consistent
c with the potential periodicity. chargetomesh does additional sums 
c to communicate the extra particle weight periodically.
            fisz=float(isz)-0.5
            fist=float(1)+0.5
            if(.not.(ix.ne.0.and.xm.gt.fist.and.xm.lt.fisz))then
c Move the particle by one grid length, to the periodic position. Use
c tiny bit less so that if it starts exactly on boundary, it does not
c end on it. The length is between end mid-cell positions.
               xgridlen=(xn(ixnp(id+1))+xn(ixnp(id+1)-1)
     $              -(xn(ixnp(id)+1)+xn(ixnp(id)+2)))*0.499999
               if(xm.le.fist)then
                  xi(id)=xi(id)+xgridlen
               elseif(xm.ge.fisz)then
                  xi(id)=xi(id)-xgridlen
               endif
               ix=interp(xn(ioff+1),isz,xi(id),xm)
               if(.not.(xm.gt.fist.and.xm.lt.fisz))then
c It's conceivable that a particle might move more than one period,
c in which case correction won't work. Don't repeat. Instead, just
c give up and call it lost but announce the problem. 
                  write(*,'(a,2i2,3f10.4)')
     $                 ' Particle crosses >meshlength',id,ix
     $                 ,(xi(k),k=1,ndims)
                  linmesh=.false.
                  goto 2
               else
c If every dimension is periodic, increment nrein. (Otherwise not)
                  if(.not.lnotallp)nrein=nrein+1
               endif
            endif
         else
            fist=1.
            if(ipartperiod(id)/64-(ipartperiod(id)/128)*2.eq.1)
     $           fist=fist+0.5
            fisz=float(isz)
            if(ipartperiod(id)/128-(ipartperiod(id)/256)*2.eq.1)
     $           fisz=fisz-0.5
c This rather complete test is necessary and costs perhaps 3% extra time.
c Because we must not allow exactly on boundaries.
            if(.not.(ix.ne.0.and.xm.gt.fist.and.xm.lt.fisz))then
               linmesh=.false.
            endif
         endif
         xfrac(id)=xm-ix
         xi(ndimsx2+id)=xm
         ixp(id)=ix
c specific particle test
c         if(i.eq.2298489) write(*,*)i,isz,ix,xm,linmesh
      enddo
 2    continue
      iregion=insideall(ndims,xi(1))
      return
      end
c********************************************************************
      subroutine postcollide(xi,tisq,iregion)
c Get new velocity; reflects neutral maxwellian shifted by vneutral.
c Update the collisional force by momentum change.
      real tisq
      integer iregion
      include 'ndimsdecl.f'
      real xi(2*ndims)
c      include 'partcom.f'
      include 'colncom.f'
      include '3dcom.f'
c Local variables
      integer icount
      real dv(ndims)
      data icount/0/

      do k=1,ndims
         dv(k)=-xi(ndims+k)
         xi(ndims+k)=tisq*gasdev(0)
         dv(k)=dv(k)+xi(ndims+k)
      enddo
c Vneutral is in z-direction.
      xi(ndims+ndims)=xi(ndims+ndims)+vneutral
      dv(ndims)=dv(ndims)+vneutral
c Now contribute the momentum change to the collisional force.
      do j=1,mf_obj
         if(btest(iregion,nf_geommap(j)-1))then
c Inside object nf_geommap(j) which is a flux measuring object.
c            write(*,*)icount,i,j,'  dv=',dv,dtaccel
            icount=icount+1
            do id=1,ndims
c This force accounting should be subsequently bulknormed by rhoinf
c and by dt.
                  colnforce(id,j,nf_step)=colnforce(id,j,nf_step)+dv(id)
            enddo
         endif
      enddo
      end
c*******************************************************************
      subroutine rotate3(xin,s,c,u)
c Rotate the input 3-vector xin, by angle theta whose sine and cosine
c are inputs about the direction given by direction cosines u
c (normalized axis vector) and return in xin.  This is a clockwise
c (right handed) rotation about positive u, I hope.

      real xin(3),u(3)
      real x1,x2
c      s=sin(theta)
c      c=cos(theta)
      d=1.-c
c      if(d.lt.1.e-5)d=s**2*0.5
c Just written out is about twice as fast:
      x1=(c+d*u(1)*u(1))*xin(1)
     $     +(d*u(1)*u(2)-s*u(3))*xin(2)
     $     +(d*u(1)*u(3)+s*u(2))*xin(3)
      x2=(c+d*u(2)*u(2))*xin(2)
     $     +(d*u(2)*u(3)-s*u(1))*xin(3)
     $     +(d*u(2)*u(1)+s*u(3))*xin(1)
      xin(3)=(c+d*u(3)*u(3))*xin(3)
     $     +(d*u(3)*u(1)-s*u(2))*xin(1)
     $     +(d*u(3)*u(2)+s*u(1))*xin(2)
c This shuffle is necessary because xin and xout are the same storage.
      xin(1)=x1
      xin(2)=x2
c And it is not permitted to make xin and xout different dummy
c arguments, if some calls are done with them being the same
c variables. That violates the fortran standard which requires that no
c aliased locations be written to in dummy arguments.
c Therefore this call assumes they are always the same and avoids 
c any aliases.
      end
c********************************************************************
      subroutine translate3(xin,vin,dt,u,xout)
c Translate in the u direction (cosines) by the component of vin 
c in the u-direction times dt from the input position xin to xout.

      real xin(3),vin(3),u(3),xout(3),dt
      vudt=(u(1)*vin(1)+u(2)*vin(2)+u(3)*vin(3))*dt
      xout(1)=xin(1)+u(1)*vudt
      xout(2)=xin(2)+u(2)*vudt
      xout(3)=xin(3)+u(3)*vudt
      end
c*********************************************************************
      subroutine gyro3(Bt,u,xin,vin,xg,xc)
c Bt in this routine is multiplied by eoverm.
c Find gyro radius and center of particle at xin, velocity vin.

c Given a field Bt
      real Bt
c in direction u, 
      real u(3)
c obtain the gyro radius xg 
      real xg(3)
c of the particle at xin, 
      real xin(3)
c with velocity vin. 
      real vin(3)
c Subtract the gyro radius from xin to give the gyrocenter in xc.
      real xc(3)
c Find the perpendicular velocity, rotate it by -90 degrees, and divide 
c by the field. 
      vu=(u(1)*vin(1)+u(2)*vin(2)+u(3)*vin(3))
      xg(1)=(vin(1)-u(1)*vu)/Bt
      xg(2)=(vin(2)-u(2)*vu)/Bt
      xg(3)=(vin(3)-u(3)*vu)/Bt
c      theta=3.1415917*0.5
      s=1.
      c=0.
      call rotate3(xg,s,c,u)
c Now xg is the gyroradius.
      xc(1)=xin(1)-xg(1)
      xc(2)=xin(2)-xg(2)
      xc(3)=xin(3)-xg(3)
      end
c*****************************************************************
      subroutine bulkforce(xpart,iregion,colntime,vneutral,Eneutral)
c Add on the contribution of this particle to the collision force.
c colntime is the collision time, vneutral the drift velocity, dt time
c step.  The total force is the contained momentum difference over the
c collision time. Needs subsequently to be properly normalized.
      include 'ndimsdecl.f'
      include 'partcom.f'
      include '3dcom.f'
      real xpart(3*ndims),colntime,vneutral
      integer iregion
      if(colntime.gt.0)then
         do i=1,mf_obj
            igeom=nf_geommap(i)
c            if(iregion.ne.0)write(*,*)'iregion',iregion
            if(btest(iregion,igeom-1))then
c Inside object i->igeom which is a flux measuring object.
c The scalar drift velocity is in the z-direction.
               do id=1,ndims
                  vid=0.
                  if(id.eq.ndims)vid=vneutral+Eneutral*colntime
                  colnforce(id,i,nf_step)=colnforce(id,i,nf_step)+
     $                 (vid-xpart(ndims+id))
               enddo
c               write(*,*)'colnforce',i,colnforce(3,i,nf_step),vd
c     $              ,xpart(ndims+3),vid,iregion
            endif
         enddo
      endif
      end
c********************************************************************
      subroutine bulknorm(fac)
c Normalize the bulkforce, multiplying by factor fac
      include 'ndimsdecl.f'
      include '3dcom.f'
      do i=1,mf_obj
         do id=1,ndims
            colnforce(id,i,nf_step)=colnforce(id,i,nf_step)
     $           *fac
         enddo
      enddo
      end
c *******************************************************************
      subroutine moveparticle(xr,ndims,Bt,Efield,Bfield,vperp
     $     ,dtpos,dtaccel,driftfield,eom)
c Accelerate and move the passed particle in the fields B,E for timestep
c dtaccel and dtpos. 
      implicit none
      integer ndims
      real Bt,dtpos,dtaccel,eom
      real xr(2*ndims),Bfield(ndims),Efield(ndims)
      real vperp(ndims),driftfield(ndims)
c Local      
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
c B-field-less Kick - Move -----------------------------
         do j=1,ndims
            xr(j+ndims)=xr(j+ndims)+eom*Efield(j)*dtaccel
            xr(j)=xr(j)+xr(j+ndims)*dtpos
         enddo          
         return
      endif

c Only if magnetic field nonzero ---------------------------------
c Rotation is counterclockwise for ions. We only want to call the 
c trig functions once, otherwise they dominate the cost.
      if(abs(theta).lt.thetatoosmall)then
c Weak B-field. ------------------------------------
c Boris Mover in frame of reference moving with external drift
         do j=1,ndims
c First E half-kick
            xr(j+ndims)=xr(j+ndims)+eom*Efield(j)*dtaccel*0.5
c Change reference to drifting frame
            xr(j+ndims)=xr(j+ndims)-vperp(j)
         enddo
c Rotate the velocity to add the magnetic field acceleration.
c This amounts to a presumption that the magnetic field acceleration
c acts at the mid-point of the translation (drift) rather than at
c the kick between translations.
         sthacc=sin(-thacc)
         cthacc=cos(thacc)
         call rotate3(xr(4),sthacc,cthacc,Bfield)
c Second half-move.
         do j=1,ndims
c Transform back to undrifted frame
            xr(j+ndims)=xr(j+ndims)+vperp(j)
c Second E half-kick
            xr(j+ndims)=xr(j+ndims)+eom*Efield(j)*dtaccel*0.5
c Move particle
            xr(j)=xr(j)+xr(j+ndims)*dtpos
         enddo          
      else
c Strong but finite B-field case -------------------------
c Cyclotronic mover
         do j=1,ndims            
c E-kick
            xr(j+ndims)=xr(j+ndims)+eom*Efield(j)*dtaccel
c Subtract off the perpendicular drift velocity.
            xr(j+ndims)=xr(j+ndims)-vperp(j)
         enddo
c Find the gyro radius and gyrocenter.
         call gyro3(eom*Bt,Bfield,xr(1),xr(4),xg,xc)
c Rotate the velocity and gyro radius.
c Rotation is counterclockwise for ions.
c For unequal timesteps, rotate gyro center by theta, not thacc.
         stheta=sin(-theta)
         ctheta=cos(theta)
c Anything other than theta gives velocity strangeness at edge.
         call rotate3(xr(4),stheta,ctheta,Bfield)
         call rotate3(xg,stheta,ctheta,Bfield)
c Move xc along the B-direction.
         call translate3(xc,xr(4),dtpos,Bfield,xc)
c Add the new gyro center and gyro radius
c And add back the drift velocity.
         do j=1,ndims
c Move the gyro-center perpendicular and add gyro-radius:
            xr(j)=xc(j)+vperp(j)*dtpos+xg(j)
c Add back vperp.
            xr(j+ndims)=xr(j+ndims)+vperp(j)
         enddo
      endif

c End of Move -----------------------------------------------------
      end
c********************************************************************
      subroutine driftparticle(xr,ndims,Bt,Efield,Bfield,vperp
     $     ,dtpos,k)
c Drift particle instead of move for cases of large Bt.
      implicit none
      integer ndims,k
      real Bt,dtpos
      real xr(2*ndims),Bfield(ndims),Efield(ndims)
      real vperp(ndims)
      integer j
c Local
      real EB(3)
      real vp
      integer k1,k2

      do j=1,ndims
         k1=mod(j,ndims)+1
         k2=mod(j+1,ndims)+1
c Explicit ExB velocity plus vperp
         EB(j)=vperp(j)+(Efield(k1)*Bfield(k2)-Efield(k2)*Bfield(k1))/Bt
      enddo
c Drift motion parallel plus EB.
c Set non-drift perp particle velocity to zero.
c Bfield is normalized: i.e. a direction cosine, so parallel velocity
      vp=0.
      do j=1,ndims
c v_parallel
         vp=vp+Bfield(j)*xr(j+ndims)
      enddo
      do j=1,ndims
c Parallel particle and perpendicular drift move:
         xr(j)=xr(j)
     $        +(vp*Bfield(j)+EB(j))*dtpos
c Reset perp velocity to just the drift. Effectively we remove the 
c perpendicular energy forcing the magnetic moment to zero. This
c is necessary because the perpendicular acceleration has been applied
c improperly for the whole timestep.
         xr(j+ndims)=vp*Bfield(j)+EB(j)
      enddo

      end
c*********************************************************************
      subroutine initdriftfield
c Initialize the drift electric field = -vxB.  The reason this can be
c different for different species is that we might want to use it to
c represent, e.g. gravitational drift. At any rate, it is the electric
c field needed to give rise to the specified (ExB) drift for that species.
      include 'ndimsdecl.f'
      include 'plascom.f'
      include 'partcom.f'
c Local
      do j=1,nspecies
         do i=1,ndims
            driftfields(i,j)=-Bt*
     $           (vperps(mod(i,ndims)+1,j)*Bfield(mod(i+1,ndims)+1)
     $           -vperps(mod(i+1,ndims)+1,j)*Bfield(mod(i,ndims)+1))
         enddo
      enddo
      end
c***********************************************************************
      subroutine locateinit()
c Initialize particles loaded from restart files so that they have the
c correct mesh-fraction coordinate in their upper values.
c If they are outside the mesh or particle region, drop them.
c Common data:
      include 'ndimsdecl.f'
      include 'partcom.f'
      include 'meshcom.f'
c Local dummy variables for partlocate.
      real xfrac(ndims)
      integer ixp(ndims)
      logical linmesh,linregion
      external linregion

      do i=1,ioc_part
         call partlocate(x_part(1,i),ixp,xfrac,iregion,linmesh)
c         write(*,*)linmesh,linregion(ibool_part,ndims,x_part(1,i))
         if(.not.linmesh .or. .not.
     $        linregion(ibool_part,ndims,x_part(1,i)))then
            x_part(iflag,i)=0
         endif
      enddo
      nrein=0
c      write(*,*)'Redetermined particle mesh locations (locateinit)'
      end
