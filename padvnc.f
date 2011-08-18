      subroutine padvnc(ndims,iLs,cij,u)
c Particle advancing routine.
c If ninjcomp (in partcom) is non-zero, then we are operating in a mode
c where the number of reinjections at each timestep is prescribed.
c Otherwise we are using a fixed number npart of particles.

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

c Local parameters
c Lowest particle to print debugging data for.
      integer npr
      parameter (npr=0)
      real fieldtoosmall
      parameter (fieldtoosmall=1.e-3)
c Local storage
      parameter (fieldtoolarge=1.e12)
      integer ixp(ndims_mesh)
      real field(ndims_mesh),adfield(ndims_mesh)
      real xfrac(ndims_mesh)
      real xg(ndims_mesh),xc(ndims_mesh)
      logical linmesh
      logical lcollided,ltlyerr
      save adfield
c Testing storage
      real ptot,atot
      data ptot/0./atot/0./

c      real fieldp(ndims_mesh)
c Needed only to print out averein:
c      common /reinextra/averein,adeficit
c Make this always last to use the checks.
      include 'partcom.f'

      tisq=sqrt(Ti)
      lcollided=.false.

      if(ndims.ne.ndims_mesh.or. ndims.ne.3)
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
      nsubc=0
      nwmax=20
      do idf=1,ndims
         adfield(idf)=0.
      enddo
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c At most do over all particle slots. But generally we break earlier.
      do i=1,n_partmax
         dtremain=0.
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
c call getptchfield.
         irptch=IAND(iregion,iptch_mask)
         if(irptch.ne.0)then
c adfield must always be zero arriving here. But we don't set it
c always because that would be expensive.
c We are in a point-charge region. Get analytic part of force.
            call getadfield(ndims,irptch,adfield,x_part(1,i),2)
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
            call getfield(ndims,cij(ic1),u,iLs
     $           ,xn(ixnp(idf)+1),idf,x_part(ndimsx2+1,i)
     $           ,IAND(iregion,ifield_mask),field(idf))
            if(.not.abs(field(idf)).lt.fieldtoolarge)then
               write(*,*)'Field corruption(?)',i,idf,iregion,field
     $              ,(x_part(kk,i),kk=1,3*ndims)
     $              ,linregion(ibool_part,ndims,x_part(1,i)),iLs
               stop
c               call partlocate(i,ixp,xfrac,inewregion,linmesh)
c               write(*,*)'inewregion,linmesh',inewregion,linmesh
            endif
            if(i.eq.1)then
c               write(*,*)
c     $             idf,irptch,' field,adfield',field(idf),adfield(idf)
c     $              ,adfield(idf)/field(idf),field(idf)
c     $              /x_part(idf,i)
            endif
            field(idf)=field(idf)+adfield(idf)
            f2=f2+field(idf)**2
            r2=r2+x_part(idf,i)**2
            v2=v2+x_part(idf+3,i)**2
         enddo
c If needed, reset adfield.
         if(irptch.ne.0)then
            do idf=1,ndims
               adfield(idf)=0
            enddo
         endif
         f1=sqrt(f2)
c Example of testing code: the few-argument field evaluator.
c                  call fieldatpoint(x_part(1,i),u,cij,iLs,fieldp)
c                  if(fieldp(2).ne.field(2))then
c                     write(*,'(i5,a,6f10.6)')i,' Point-field:',
c     $                    fieldp,((fieldp(j)-field(j)),j=1,ndims)
c--------------------------------------
c Subcycle restart.
c---------- Subcycling ----------------
         if(subcycle.ne.0)then            
c One might have a conditional test to set subcycling for this particle
c such as this automatic one based on avoiding excessive acceleration.
c            dtc=dt/max(1.,anint(f1*dt/subcycle))
c            dtc=.1000
c For testing we just set directly.
            dtc=subcycle
            if(dtc.lt.dtpos)then
c Take sub-step.
               nsubc=nsubc+1
c               write(*,*)dtpos,dtremain,dtc
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
c---------------- Particle Moving ----------------
c Use dtaccel for acceleration. May be different from dtpos if there was
c a subcycle, reinjection or collision last step.
         dtaccel=0.5*(dtpos+dtprec(i))
         ptot=ptot+dtpos
         atot=atot+dtaccel
c Accelerate ----------
         do j=ndims+1,2*ndims
            x_part(j,i)=x_part(j,i)+field(j-ndims)*dtaccel
         enddo
c Move ----------------
         if(Bt.eq.0.)then
            do j=1,ndims
               x_part(j,i)=x_part(j,i)+x_part(j+ndims,i)*dtpos
            enddo          
         elseif(Bt.lt.Btinf)then
c Magnetic field non-zero but finite
            theta=Bt*dtpos
c Rotation is counterclockwise for ions. We only want to call the 
c trig functions once, otherwise they dominate the cost.
            stheta=sin(-theta)
            ctheta=cos(theta)
c            if(i.eq.1)write(*,*)'theta=',theta,fieldtoosmall
            if(theta.lt.fieldtoosmall)then
c Weak B-field. Advance using summed accelerations. First half-move
               do j=1,ndims
                  x_part(j,i)=x_part(j,i)+x_part(j+ndims,i)*dtpos*0.5
               enddo          
c Rotate the velocity to add the magnetic field acceleration.
c This amounts to a presumption that the magnetic field acceleration
c acts at the mid-point of the translation (drift) rather than at
c the kick between translations.
               call rotate3(x_part(4,i),stheta,ctheta,Bfield)
c Second half-move.
               do j=1,ndims
                  x_part(j,i)=x_part(j,i)+x_part(j+ndims,i)*dtpos*0.5
               enddo          
            else
c Strong but finite B-field case:
c Subtract off the perpendicular drift velocity.
               do j=ndims+1,2*ndims
                  x_part(j,i)=x_part(j,i)-vperp(j-ndims)
               enddo
c Find the gyro radius and gyrocenter.
               call gyro3(Bt,Bfield,x_part(1,i),x_part(4,i),xg,xc)
c Rotate the velocity and gyro radius.
               call rotate3(x_part(4,i),stheta,ctheta,Bfield)
               call rotate3(xg,stheta,ctheta,Bfield)
c Move xc along the B-direction.
               call translate3(xc,x_part(4,i),dtpos,Bfield,xc)
c Add the new gyro center and gyro radius
c And add back the drift velocity.
               do j=1,ndims
c Move the gyro-center perpendicular and add gyro-radius:
                  x_part(j,i)=xc(j)+vperp(j)*dtpos+xg(j)
                  x_part(j+ndims,i)=x_part(j+ndims,i)+vperp(j)
               enddo
            endif
 801        format(a,6f10.6)
            if(i.le.npr)write(*,801)'x_part2',(x_part(k,i),k=1,6)
         else
c Infinite magnetic field; i.e. one-dimensional motion plus a
c perpendicular steady drift velocity. Set non-drift perp particle
c velocity to zero.
c Bfield is normalized: i.e. a direction cosine.
            vp=0.
            do j=1,ndims
               vp=vp+Bfield(j)*x_part(j+ndims,i)
            enddo
            do j=1,ndims
c Parallel particle and perpendicular drift move:
               x_part(j,i)=x_part(j,i)+(vp*Bfield(j)+vperp(j))*dtpos
c Zero perp velocity.
               x_part(j+ndims,i)=vp*Bfield(j)
            enddo
         endif
c End of Move ---------

         if(lcollided)then
c Treat collided particle at (partial) step end
            call postcollide(i,tisq)
            lcollided=.false.
         endif

         call partlocate(i,ixp,xfrac,inewregion,linmesh)
c---------------------------------
c If we crossed a boundary, do tallying.
         ltlyerr=.false.
         if(inewregion.ne.iregion)
     $        call tallyexit(i,inewregion-iregion,ltlyerr)
c------------ Possible Reinjection ----------
         if(ltlyerr .or. .not.linmesh .or.
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
         dtprec(i)=dtpos
c If we haven't completed this full step, do so.
         if(dtremain.gt.0.)then
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
            write(*,*)'Reinject out of mesh',i,xfrac
            stop
         endif
         if(.not.linregion(ibool_part,ndims,x_part(1,i)))then
c This situation is benign and not an error if we have a region that
c happens not to cover the entire mesh edge. So don't stop, retry.
c            write(*,*)'Reinject out of region',i,iregion,xfrac
c            stop
            goto 200
         endif
         dtpos=(dtpos+dtremain)*ran1(myid)
         dtprec(i)=0.
         dtremain=0.
         nlost=nlost+1
         nrein=nrein+ilaunch
         phi=getpotential(u,cij,iLs,x_part(2*ndims+1,i)
     $        ,IAND(iregion,ifield_mask),2)
c This version multiply-weights a relaunched case:
c         phirein=phirein+ilaunch*phi
c But has to compensate for final division by nrein.
c It's probably better to count the relaunches as average:
         phirein=phirein*(1+ilaunch-1)+phi
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
         cap=2.*Ti
         if(phirein.gt.cap)then
c            if(myid.eq.0)write(*,*)'PROBLEM: phirein>0:',phirein
            phirein=cap
         endif
      else
         if(ninjcomp.gt.100)write(*,*)'No reinjections'
      endif

c      if(nsubc.ne.0) write(*,*)'Subcycled:',nsubc
c      write(*,*)'Padvnc',n_part,nrein,ilaunch,ninjcomp,n_partmax
      end

c***********************************************************************
      subroutine getadfield(ndims,irptch,adfield,xp,isw)
c Cycle through the nonzero bits of irptch and add up the extra
c potential (isw=1), field (isw=2) or charge (isw=3) contributions.
c adfield should be zero on entry, because it ain't set,
c just incremented.

      integer ndims,irptch,isw
      real adfield(ndims)
      real xp(ndims)
      include '3dcom.f'
      real xd(ns_ndims)
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
            if(p2.lt.1.e-12)then
c Avoid overflows:
               p2=1.e-12
               write(*,*)'ptch field overflow corrected'
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
      subroutine setadfield(ndimsp,ifull,iuds,irptch,lsliceplot)
c Set the values of the potential and charge at the grid
c points that compensate for the analytic field of getadfield.

      integer ndimsp
      integer ifull(ndimsp),iuds(ndimsp),irptch
      logical lsliceplot
c Defines iptch_copy uci, rhoc and dimensions.
      include 'griddecl.f'
      include 'ptchcom.f'
c ndims must be same as ndimsp.
c To do the slice plot we need:
      include 'meshcom.f'
      real zp(na_m,na_m)
      integer ipoint
      external ucrhoset
      iptch_copy=irptch
c      write(*,*)'Point charges included. Mask:',iptch_copy
      ipoint=0
      ifix=1
      call mditerarg(ucrhoset,ndimsp,ifull,iuds,ipoint
     $     ,uci,rhoci,iptch_copy,idum)
      if(lsliceplot)then
         call sliceGweb(ifull,iuds,uci,na_m,zp,
     $        ixnp,xn,ifix,'u!dc!d ptch')
         call sliceGweb(ifull,iuds,rhoci,na_m,zp,
     $        ixnp,xn,ifix,'!Ar!@!dc!d ptch')
      endif
      end
c**********************************************************************
      subroutine ucrhoset(inc,ipoint,indi,ndims,iuds,
     $     uci,rhoci,iptch_copy)
c Set uci and rhoci when point charges are present.
      integer inc,ipoint,ndims,indi(ndims)
      integer iuds(ndims)
      real uci(*),rhoci(*)
c Commons: For position.
      include 'meshcom.f'
c For debyelen
      include 'plascom.f'
c Local storage:
      integer isw,iregion,irptch
      real xp(ndims_mesh),adfield(ndims_mesh)
c Silence warnings
      irptch=iuds(1)
c      write(*,*)'ucrhoset',ipoint,indi,iuds,ndims
c      if(ipoint.gt.100)stop
c Get grid point position, and irptch.
      do id=1,ndims
         xp(id)=xn(ixnp(id)+1+indi(id))
      enddo
      iregion=insideall(ndims,xp)
      irptch=IAND(iregion,iptch_copy)
      if(irptch.ne.0)then
c Get uc
         isw=1
         adfield(1)=0.
         call getadfield(ndims,irptch,adfield,xp,isw)
         uci(ipoint+1)=adfield(1)
c Get charge
         isw=3
         adfield(1)=0.
         call getadfield(ndims,irptch,adfield,xp,isw)
         rhoci(ipoint+1)=debyelen**2*adfield(1)
      else
         uci(ipoint+1)=0.
         rhoci(ipoint+1)=0.
      endif
c Always just increment by 1
      inc=1
c      write(*,*)'ucrhoset return',irptch
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

      integer i,iregion
      logical linmesh
c meshcom provides ixnp, xn, the mesh spacings. (+ndims_mesh)
      include 'meshcom.f'
      integer ixp(ndims_mesh)
      real xfrac(ndims_mesh)
      parameter (ndimsx2=ndims_mesh*2)

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
c Get new velocity; reflects neutral maxwellian shifted by vneutral.

      integer i
      real tisq
      include 'partcom.f'
      include 'colncom.f'
      do k=1,npdim
         x_part(npdim+k,i)=tisq*gasdev(idum)
      enddo
c Vneutral is in z-direction.
      x_part(npdim+npdim,i)=x_part(npdim+npdim,i)+vneutral
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
