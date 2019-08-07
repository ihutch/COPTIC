!*********************************************************************
! Wrapper:
      subroutine reinject(xr,ilaunch,ispecies)
      real xr(*)
      include 'colncom.f'
      include 'ndimsdecl.f'
      include 'partcom.f'

! Test whether the distribution is separable. It MUST be if pinit.f
! has been replaced by something that does not initialize colreinject.
      if(notseparable(ispecies).eq.0)then
! Cartreinject is the version to be used when the distribution is 
! separable in the directions of the cartesian axes.
         call cartreinject(xr,ilaunch,caverein,ispecies)
      else
! colreinject started as the collisional reinjection scheme. 
! It is based upon sampling from a 3-D distribution function.
! It can therefore be used for situations where the distribution
! is not separable in coordinate directions.
         call colreinject(xr,ipartperiod,caverein,ispecies)
! Only one launch allowed here (and actually in the other too).
         ilaunch=1
      endif
      end
**********************************************************************
! Initialize whether we are initialized!
      block data creinset
      include 'ndimsdecl.f'
      include 'creincom.f'
      data lreininit/.false./
      end
**********************************************************************
! Cartesian reinjection routine.
      subroutine cartreinject(xr,ilaunch,caverein,ispecies)
      implicit none
      integer mdims
      parameter (mdims=3)
      real xr(3*mdims)
      real caverein
      integer ilaunch
      integer ispecies
      real ra,fc,sg,fr,x2,x1
      integer index,k,iother,ir
      real ranlenposition
      external ranlenposition

      include 'ndimsdecl.f'
! Random choice data for these routines:
      include 'creincom.f'
! Plasma common data
      include 'plascom.f'
! Mesh data
      include 'meshcom.f'
! Particle data to define the species of interest

      if(.not.lreininit)call cinjinit()
!      write(*,*)'Returned from cinjinit'
!----------------------------------------
! Pick the face from which to reinject.
      call ranlux(ra,1)
      call invtfunc(gintreins(0,ispecies),7,ra,fc)
! Returns fc between 1+ and 7-.
      index=int(fc)
! Make idrein 1-3, and sg +-1. 
! Plus/minus refers to the direction of velocity.
! So + means plane xstart, and - xend relative to a normally ordered mesh.
      idrein=(index+1)/2
! Just less than +-1.
      sg=-(2*mod(index,2)-1)*0.999999
! Position in this dimension (either end), sg ensures (just) inside the
! mesh.
      xr(idrein)=0.5*(xmeshstart(idrein)*(1.+sg)+
     $     xmeshend(idrein)*(1.-sg))
!----------------------------------------
! Pick the velocities and positions parallel to this face:
! Doing position and velocity simultaneously requires velocity distrib
! to be independent of position.
      do k=1,2
         iother=mod(idrein+k-1,3)+1
! Position: Ensure we never quite reach the mesh edge:
! Or Accounting for density gradients:
         xr(iother)=ranlenposition(iother)
! Velocity
         call ranlux(ra,1)
         ra=ra*ncrein
         ir=int(ra)
         fr=ra-ir
         if(fr.gt.1. .or. fr.lt.0.)then
! This should never happen.
            write(*,*)'Creinject fraction wrong',ra,ir,fr
         endif
! velocity, linear interpolation:
         xr(mdims+abs(iother))= preins(ir,iother,ispecies)*(1-fr)
     $        +preins(ir+1,iother,ispecies)*fr
      enddo
!----------------------------------------
! Pick the velocity perpendicular to this face:
      call ranlux(ra,1)
      ra=ra*ncrein
      ir=int(ra)
      fr=ra-ir
!      write(*,*)ra,ir,index,fr,ncrein
      xr(mdims+idrein)=hreins(ir,index,ispecies)*(1-fr)+hreins(ir+1
     $     ,index,ispecies)*fr
! In this version of reinject we never try more than one launch
      ilaunch=1
!----------------------------------------
! Correct the total energy for caverein. 
      x2=0.
      do k=1,mdims
         x2=x2+xr(mdims+k)**2
      enddo
      x1=(x2+2.*max(0.,-caverein))/x2
      x1=sqrt(x1)
      do k=1,mdims
         xr(mdims+k)=xr(mdims+k)*x1
      enddo

      end
!*********************************************************************
      subroutine cinjinit()
! Cartesian reinjection initialization for drift velocity, vds, in the
! z-direction and maxwellians of width given by Ts

      include 'ndimsdecl.f'
      include 'plascom.f'
      include 'meshcom.f'
      include 'creincom.f'
      include 'myidcom.f'
      include 'partcom.f'
      external ffdrein
      external fvdrein
      external ff1crein
      external fv1crein
      real fv1crein,ff1crein,fvdrein,ffdrein
      parameter (bdys=6.)
      real thetamax
      parameter (thetamax=1.)
      common /species/ispec

      do ispec=1,nspecies
      theta=Bt*eoverms(ispec)*dt
      xc=0.
! For three coordinate directions
      do id=1,3
! Width of velocity distribution
         xw=sqrt(Ts(ispec)*abs(eoverms(ispec)))
!    For two ends (and hence velocity polarities)
         do i2=1,2
! idrein determines the sign of velocity. i2 odd => idrein negative.
            idrein=id*(2*i2-3)
            index=2*(id-1)+i2
            if(abs(theta).lt.thetamax)then
! Set the inverse cumulative probability fn hrein, and flux grein
               call cumprob(ffdrein,xw,xc,
     $           ncrein,hreins(0,index,ispec),greins(index,ispec),myid)
            else
! Infinite-Bt case 1-d projection.
               xc=vpars(ispec)*Bfield(id)+vperps(id,ispec)
               xw=xw*Bfield(id)
!               write(*,*)'Calling cumprob',id,xc,xw
               call cumprob(ff1crein,xw,xc,
     $           ncrein,hreins(0,index,ispec),greins(index,ispec),myid)
            endif
! Zero grein on absorbing faces.
            ip=ipartperiod(id)
! Ordering in grein etc is (+,-) for each dim. Upper face/Lower face. 
! That's the opposite of what I adopted for ipartperiod because
! it corresponds to negative/positive velocity. Pity!
            ip=ip/2**(2-i2)
            ip=ip-2*(ip/2)
            if(ip.eq.1)greins(index,ispec)=0.
            if(gnt.ne.0)then
! Correcting for uniform density scale-length.
               if(idrein.gt.0)then
                  dengfac=exp((min(xmeshstart(id),xmeshend(id))
     $                 -gp0(id))*gn(id))
               else
                  dengfac=exp((max(xmeshstart(id),xmeshend(id))
     $                 -gp0(id))*gn(id))
               endif
               do k=0,1
                  ii=mod(abs(id)+k,ndims)+1
                  if(gn(ii).ne.0)then
                     dengfac=dengfac*(
     $                   exp((max(xmeshstart(ii),xmeshend(ii))-gp0(ii))
     $                   *gn(ii))-
     $                   exp((min(xmeshstart(ii),xmeshend(ii))-gp0(ii))
     $                   *gn(ii)) )/
     $                   (gn(ii)*abs(xmeshend(ii)-xmeshstart(ii)))
                  else
                  endif
               enddo
               greins(index,ispec)=greins(index,ispec)*dengfac
            endif
! Kludge fix of ends to avoid negative velocity injections.
            if(idrein.gt.0)then
               if(hreins(0,index,ispec).lt.0.)hreins(0,index,ispec)=0.
            else
               if(hreins(ncrein,index,ispec).gt.0.)hreins(ncrein,index
     $              ,ispec)=0.
            endif
!            write(*,*)index,(hreins(kk,index,ispec),kk=0,5)
!     $           ,(hreins(kk,index,ispec),kk=ncrein-4,ncrein)
         enddo
         idrein=id
         if(abs(theta).lt.thetamax)then
            call cumprob(fvdrein,xw,xc,
     $           ncrein,preins(0,id,ispec),gdummy,myid)
         else
            call cumprob(fv1crein,xw,xc,
     $           ncrein,preins(0,id,ispec),gdummy,myid)
         endif
      enddo
!
!      grein(1)=1.*grein(1)
!      write(*,*)'grein arbirtrary adjustment',greins
      gtot=0.
! Alternative general-dimension fcarea calculation:
      do i=1,ndims
         fcarea(i)=1.
! Set area tiny for period directions
         if(lnotallp.and.(ipartperiod(i).eq.4.or.ipartperiod(i).eq.5)
     $        )fcarea(i)=1.e-6
         do j=1,ndims-1
            id=mod(i+j-1,ndims)+1
            fcarea(i)=fcarea(i)*abs(xmeshend(id)-xmeshstart(id))
         enddo
!         write(*,*)'fcarea(',i,')=',fcarea(i)
      enddo
      do id=1,6
         gtot=gtot+greins(id,ispec)*fcarea((id+1)/2)
      enddo
      gintreins(0,ispec)=-0.0000005
      do id=1,6
         gintreins(id,ispec)=gintreins(id-1,ispec) +
     $        1.000001*greins(id,ispec)*fcarea((id+1)/2)/gtot
      enddo
      if(.not.gintreins(6,ispec).gt.1.)write(*,*)'gintrein problem!'
      enddo
      if(myid.eq.0)then
      write(*,'(a,3i2,a)')'Injection initialization ipartperiod'
     $        ,ipartperiod,'  greins,gintrein:'
      write(*,*)'greins is the total flux across each of 6 faces'
      write(*,'(6f10.5)')greins
      write(*,'(7f10.5)')gintrein
      endif
      lreininit=.true.
      
      end
!**********************************************************************
! Given a function f(x), whose integral from x=-\infty to +infty exists,
! obtain the definite integral g(x) = \int_-\infty^x f(x) dx.
! Return the array h_j (j=0,K) equally spaced on g, such that 
!    g(h_j)=g(\infty) (j/K)
! Thus h_j is the inverse of the cumulative probability distribution,
! g/g(\infty) based on f(x). Also return g(\infty).
! On entry an estimate of the width of f is provided in xw
!                  and of the center of f in xc
! The range is adjusted automatically to give a decently accurate integral.
! It is possible for the adjustment process to fail if xc is in error by, 
! or xw is, a large factor times the actual width.
! The most efficient usage is for xw to be slightly larger than the 
! points at which f is 10^-2 times its peak.
! f must be declared external in the calling routine.
! If ginfty is returned as zero, then the routine failed with cumulative
! probability too small to be significant. 
! If myid .gt. 0 then don't print warning messages. 
! If myid .lt. 0 print extra warnings.
      subroutine cumprob(f,xw,xc,K,h,ginfty,myid)
      external f
      real xw,xc
      integer K
      real h(0:K)

      parameter (tiny=1.e-14)
! internal storage
      integer m
      parameter (m=2000)
      real fv(m),g(m),x(m)

      if(K.gt.m .and. myid.le.0)
     $     write(*,*)'cumprob warning: too fine a grid ',K

!      write(*,*)'cumprob entry',xw,xc,K

! If we return prematurely, it is with h's=0.
      do i=0,K
         h(i)=0.
      enddo

      xr=abs(xw)
      if(xr.eq.0)xr=1.
      x0=xc-2.*xr
      x1=xc+2.*xr
      icount=0

! Cycle over integrations.
 1    continue
      icount=icount+1
      xd=(x1-x0)/m
      g(1)=0.
      x(1)=x0
      fv(1)=f(x0)
! Form the integral \int_x0^{x0+m*xd}
      do i=2,m
         x(i)=x0+i*xd
         fv(i)=f(x(i))
         g(i)=g(i-1)+xd*(fv(i)+fv(i-1))*.5
      enddo
!      write(*,*)'Integration limits',x0,x1,' Value=',g(m)

! Decide whether we have covered enough range.
! Criteria are: dg0/gm < 1/(10.K.m) ; dg(m/4)/gm > 1/(10.K.m)
! and similarly for top end.  dg=fv*xd
! If dg0 is violated, double distance from middle.
! If dg(m/4) is violated, halve it. 

      ii=1
      ia=m
      dgc=g(m)/(10.*K*m*xd)
      isearch=0

! Find acceptable end points, for this integral
! allowing only 20 adjustment (per cycle)
      do j=1,20
         isearch=isearch+1
         ll=0
         lr=0

!         write(*,'(a,2i5,3g12.4)')'ii,ia,x0,x1,gm',ii,ia,x0,x1,g(m)
         ii2=(3*ii+ia)/4
         ia2=(ii+3*ia)/4
         if(fv(ii).ge.dgc)then
! We cover not enough -x
            ii=ii-(ii2-ii)
         elseif(fv(ii2).lt.dgc)then
! We cover too much -x         
            ii=ii2
         else
            ll=1
         endif

         if(fv(ia).ge.dgc)then
! We cover too little +x
            ia=ia+(ia-ia2)
         elseif(fv(ia2).lt.dgc)then
! We cover too much +x
            ia=ia2
         else
            lr=1
         endif
! If both ends were acceptable this iteration
         if(ll.eq.1 .and. lr.eq.1)then
! If we haven't moved the ends break; we are done.
            if(isearch.eq.1) goto 2
! Else reintegrate with these end points (which should be acceptable). 
            x0=x0+(ii-1)*xd
            x1=x1+(ia-m)*xd
!            write(*,*)'New ends',x0,x1
            goto 1
         endif

         if(ii.lt.1 .or. ia.gt.m)then
! Range is too short. Need to reintegrate from scratch
! Double or triple the range.
            dx=x1-x0
            if(ii.lt.1) x0=x0-dx
            if(ia.gt.m) x1=x1+dx
            if(icount.gt.20)then
!               write(*,'(a,/,a,g12.4,a,g12.4,a)')
!     $              'Integrate too high icount.',' Starting range'
!     $              ,xc,'+-',xw,' may be too wide.'
               ginfty=0.
               return
!               stop
            endif
            goto 1
         endif

         if(abs(ia-ii).lt.4)then
! Range is too big.
            x1=x0+(ii+2)*xd
            x0=x0+(ii-2)*xd
            if(icount.gt.20)then
!               write(*,*)'Integrate count too high.'
               ginfty=0.
               return
            endif
            goto 1
         endif
      enddo
      if(myid.le.0)write(*,*) "Cumprob exhausted end adjustment steps"
      ginfty=0.
      return
 2    continue

      if(ii.eq.1 .and. ia.eq.m)then
! We have the right range. Finish.
         if(.not.g(m).gt.tiny)then
            if(myid.le.0)write(*,*)'Cumprob total too small',g(m)
     $           ,' set to zero',' xw=',xw,' xc=',xc
            ginfty=0.
            return
         endif
         do i=0,K
! solve g(h_i)=gm*i/K with linear interpolation.
            gt=(.999998*(i/float(K))+.000001)*g(m)
            ihi=interp(g,m,gt,hi)
            if(ihi.eq.0)then
               write(*,*)'Cumprob interp error',i,K,m,gt,g(m)
               write(*,*)'Integration range used:',x0,x1
               write(*,*)(g(kk),kk=1,5),(g(kk),kk=m-4,m)
               write(*,*)'Probably the value of tiny',tiny,' is too big'
               stop
            endif
            hf=hi-ihi
            h(i)=x(ihi)*(1.-hf)+x(ihi+1)*hf
!            if(i.le.1.or.i.ge.(K-1))write(*,*)'i,ihi,h(i)',i,ihi,h(i)
         enddo
      else
         goto 1
      endif
!      write(*,*)h
!      call yautoplot(g,m)
!      call pltend()
!      call yautoplot(h(0),K+1)
!      call pltend()
      ginfty=g(m)
      end
!******************************************************************
! Routines for reinjection calculations with ions drifting relative
! to neutrals, driven by external force (e.g. E-field).
! Direct replacements for the fvcrein, ffcrein functions.
!*******************************************************************
      real function fvdrein(v)
! Return the probability distribution for reinjection,
! in coordinate direction idrein (signed) in creincom,
! from a drift-distribution shifted by vds (in plascom),
! whose direction cosines are vdrift(3) (in plascom),
! colliding with neutrals of velocity vneutral (in colncom).
! If v is normalized by sqrt(ZT_e/m_i), then Tsi is the ratio T_i/ZT_e.
! Because the drift distribution is not separable except in the directions
! perpendicular and parallel to the drift, only degenerate Maxwellian
! (vds-vneutral=0) non-z cases are allowed so far.
! Diamagnetic drift is implemented as a Maxwellian shift.
      implicit none
      real v,vn,u,ud,fvcx
      external fvcx
      include 'ndimsdecl.f'
      include 'plascom.f'
      include 'creincom.f'
      include 'colncom.f'
      real vdia
      integer id2,id3
      integer ispec
      common /species/ispec

      ud=vds(ispec)-vneutral
      if(ispec.ne.1)ud=0.
      vdia=0.
! Are there background diamagnetic drifts? Density gradient. 
      if(gnt.ne.0.)then
! Gradient is always perpendicular to B.
         id2=mod(abs(idrein),3)+1
         id3=mod(abs(idrein)+1,3)+1
! v_dia = -(T_i \nabla n \times B / qnB^2)
         vdia=- (gn(id2)*Bfield(id3)-gn(id3)*Bfield(id2))*Ts(ispec)/Bt
! The approximation is to represent the perpendicular distribution by
! a Maxwellian shifted by the diamagnetic drift.
      endif
      vn=sqrt(2.*Ts(ispec)*abs(eoverms(ispec)))
      if(vdrifts(3,ispec).eq.1.)then
! Z-drift cases. (equiv old)
         if(abs(idrein).eq.3)then
            u=(v-vds(ispec)+ud-vdia)/vn
            ud=ud/vn
            fvdrein=fvcx(u,ud)
            fvdrein=fvdrein/vn
         else
            fvdrein=exp(-((v-vdia)/vn)**2)/(vn*sqrt(3.1415926))
         endif
      else
! Non-z
         if(ud.ne.0.)then
            write(*,*)'Non-z collisional drift not implemented.'
     $           ,' Aborting.'
            stop
         else
! Maxwellian non-z drift.
            fvdrein= exp(-((v-vds(ispec)*vdrifts(abs(idrein),ispec)
     $           -vdia)/vn)**2) /(vn*sqrt(3.1415926))
         endif
      endif
      end
!*******************************************************************
      real function ffdrein(v)
! Return the one-way differential flux distribution as a fn of velocity
! from a drift distribution (in dimension 3).
! In the positive or negative direction, determined by idrein's sign.
! Diamagnetic drift is implemented as a Maxwellian shift.
      implicit none
      real v,vn,u,ud,fvcx
      external fvcx
      include 'ndimsdecl.f'
      include 'plascom.f'
      include 'creincom.f'
      include 'colncom.f'
      real vdia
      integer id2,id3
      integer ispec
      common /species/ispec

      ud=vds(ispec)-vneutral
      if(ispec.ne.1)ud=0.
      vdia=0.
! Are there background diamagnetic drifts? Density gradient. 
      if(gnt.ne.0.)then
! Gradient is always perpendicular to B.
         id2=mod(abs(idrein),3)+1
         id3=mod(abs(idrein)+1,3)+1
! v_dia = -(T_i \nabla n \times B / qnB^2)
         vdia=- (gn(id2)*Bfield(id3)-gn(id3)*Bfield(id2))*Ts(ispec)/Bt
! The approximation is to represent the perpendicular distribution by
! a Maxwellian shifted by the diamagnetic drift.
      endif
!      if(vdia.ne.0..and.abs(v-1.).lt.0.0005)then
!         write(*,*)'v,ffdrein,idrein,vdia',v,ffdrein,idrein,vdia,id2,id3
!      endif
      vn=sqrt(2.*Ts(ispec)*abs(eoverms(ispec)))
      if(vdrifts(3,ispec).eq.1.)then
! Z-drift cases. (equiv old)
         if(abs(idrein).eq.3)then
! In z-direction use appropriate drift distribution.
            u=(v-vds(ispec)+ud-vdia)/vn
            ud=ud/vn
            ffdrein=fvcx(u,ud)
            ffdrein=abs(v)*ffdrein/vn
         else
            ffdrein=abs(v)*exp(-((v-vdia)/vn)**2)/(vn*sqrt(3.1415926))
         endif
      else
! Non-z
         if(ud.ne.0.)then
            write(*,*)'Non-z collisional drift not implemented.'
     $           ,' Aborting.'
            stop
         else
! Maxwellian non-z drift.
            ffdrein= exp(-((v-vds(ispec)*vdrifts(abs(idrein),ispec)
     $           -vdia)/vn)**2) *abs(v)/(vn*sqrt(3.1415926))
         endif
      endif
! Don't count particles going the wrong way:
      if(int(sign(1.,v)).ne.sign(1,idrein))ffdrein=0.

      end
!****************************************************************
! FVCX function for 1-d drifting CX distribution.
      function fvcx(u,ud)
      real u,ud,v,vw,fvcx
! Return the normalized distribution function f(u)=v_n f(v) for constant
! cx collision frequency at a value of normalized velocity u=v/v_n, when
! the normalized drift velocity is ud= (a/\nu_c) /v_n, with v_n =
! sqrt(2T_n/m). a is acceleration, nu_c collision freq.  This is the
! solution of the steady Boltzmann equation.
      if(ud.lt.0.) then
         v=-u
         vw=-ud
      else
         v=u
         vw=ud
      endif
      if(vw.eq.0.)then
         carg=20.
         earg=100
      else
         carg=0.5/vw-v
         earg=(0.5/vw)**2-v/vw
      endif
      if(carg.gt.10)then
! asymptotic form for large exp argument (small vw):
!  exp(-v^2)/[sqrt(\pi)(1-2 v_d v)]:
         fvcx=exp(-v**2)/1.77245385/(1.-2.*vw*v)
      elseif(carg.gt.-5.)then
         fvcx=exp(-v**2)*experfcc(carg)*0.5/vw
      else
!         fvcx=exp(earg)*erfcc(carg)*0.5/vw
         fvcx=exp(earg)/vw
      endif
!      write(*,*)'fvcx:vw,v,earg,fvcx',vw,v,earg,fvcx
      if(.not.fvcx.ge.0) then
         write(*,*)'fvcx error. u=',u,' ud=',ud,' f=',fvcx,carg
         fvcx=0.
      endif
      end
!**********************************************************************
      real function ff1crein(v)
! This is the flux function for 1-D motion in the direction of Bfield.
!
! Return the flux from a maxwellian for dimension idrein (in creincom)
! In the positive or negative direction, determined by idrein's sign.
! Maxwellian is shifted by vdj=Bfield(j)*vpar+vperp(j).
! If v is normalized by sqrt(ZT_e/m_i), then Tsi is the ratio T_i/ZT_e.
! But here, the thermal spread along Bfield must be projected into
! the coordinate direction idrein. 

      include 'ndimsdecl.f'
      include 'plascom.f'
      include 'creincom.f'
      real vt2min,argmax
      parameter (vt2min=1.e-6,argmax=12.)
      integer ispec
      common /species/ispec
      
      if(int(sign(1.,v)).eq.sign(1,idrein))then
         Tsi=Ts(ispec)*abs(eoverms(ispec))
         j=abs(idrein)
         vs=Bfield(j)*vpars(ispec)+vperps(j,ispec)
         vt2=2.*Tsi*Bfield(j)**2
         if(vt2.gt.vt2min)then
            arg=-(v-vs)**2/vt2
            scale=1./Bfield(j)
         else
            arg=-(v-vs)**2/vt2min
            scale=sqrt(2.*Tsi/vt2min)
         endif
         if(abs(arg).gt.argmax)then
            ff1crein=0.
         else
!            if(Bfield(j).eq.0.)write(*,*)'v,vs,arg',v,vs,arg
! Corrected the scaling to be consistent:
            ff1crein=scale*abs(v)*exp(arg)/sqrt(2.*Tsi*3.1415926)
         endif
      else
         ff1crein=0.
      endif
      end
!**********************************************************************
      real function fv1crein(v)
! This is the 1-d probability distribution projected in coordinate
! direction idrein (in crein).
!
! Return the probability distribution value.

      include 'ndimsdecl.f'
      include 'plascom.f'
      include 'creincom.f'
      real vt2min,argmax
      parameter (vt2min=1.e-5,argmax=12.)
      integer ispec
      common /species/ispec

      j=abs(idrein)
      vs=Bfield(j)*vpars(ispec)+vperps(j,ispec)
      Tsi=Ts(ispec)*abs(eoverms(ispec))
      vt2=2.*Tsi*Bfield(j)**2
      if(vt2.gt.vt2min)then
         arg=-(v-vs)**2/vt2
         scale=1./Bfield(j)
      else
         arg=-(v-vs)**2/vt2min
         scale=sqrt(2.*Tsi/vt2min)
      endif
      if(abs(arg).gt.argmax)then
         fv1crein=0.
      else
         fv1crein=scale*exp(arg)/sqrt(2.*Tsi*3.1415926)
      endif

      end
!*********************************************************************
      subroutine rhoinfcalc(dtin)
! Obtain the rhoinf to be used in calculating the electron shielding,
! based upon the number and average potential of the reinjections.
      include 'ndimsdecl.f'
      include 'plascom.f'
! No time-averaging for now.
! Use particle information for initializing.
      include 'partcom.f'
      include 'meshcom.f'
      include 'creincom.f'
      include 'myidcom.f'
! Local variables
      real volume,flux
      real cfactor

      cfactor=1.
      volume=1.
      flux=0.
! In reality this calculation is needed only once.
      do i=1,ndims
         fcarea(i)=1.
         if(lnotallp.and.(ipartperiod(i).eq.4.or.ipartperiod(i).eq.5)
     $        )fcarea(i)=1.e-6
         do j=1,ndims-1
            id=mod(i+j-1,ndims)+1
            fcarea(i)=fcarea(i)*(xmeshend(id)-xmeshstart(id))
         enddo
!         write(*,*)'fcarea(',i,')=',fcarea(i)
! Multiply the face area by the factor relative to unshifted maxwellian
! that specifies the flux for this face, and add to total. 
         flux=flux+(grein(2*i-1)+grein(2*i))*fcarea(i)
         volume=volume*(xmeshend(i)-xmeshstart(i))
      enddo

      if(ninjcomp.ne.0.or..not.lnotallp)then
! Fixed injection rate implies fixed rhoinf. Set it only once.
! .not.lnotallp is the all periodic particles case. 
         if(rhoinf.eq.0)then
!            if(.true.)then
            rhoinf=numprocs*(ninjcompa(1)+pinjcompa(1))/(dtin*flux)
            chi=0.
!            write(*,*)'Direct periodic',numprocs,n_part
!     $           ,numprocs*n_part/volume
!            write(*,*)'*************flux,rhoinf=',flux,rhoinf
         else
!            write(*,*)'rhoinf=',rhoinf,ninjcomp
         endif
      else
         if(nrein.ge.10)then
! Calculate rhoinf from nrein if there are enough.
! Correct approximately for edge potential depression (OML).
! Use the first species drift velocity as for calculating rhoinf.
            chi=max(crelax*(-phirein/Ti)+(1.-crelax)*chi,0.)
            cfactor=smaxflux(vds(1)/sqrt(2.*Ti),chi)
     $           /smaxflux(vds(1)/sqrt(2.*Ti),0.)
            rhoinf=(nrein/(dtin*cfactor*flux))
!         write(*,*)nrein,dtin,phirein,chi,cfactor,flux,rhoinf
         else
            if(rhoinf.lt.1.e-4)then
! Approximate initialization
               rhoinf=numprocs*n_part/volume
               write(*,*)'Rhoinf in rhoinfcalc approximated as',rhoinf
     $              ,numprocs,n_part
            endif
! Else just leave it alone.
         endif
      endif
      if(myid.eq.0.and..false.)then
         write(*,*)
         write(*,'(a,2i8,10f9.4)') 'Ending rhoinfcalc',nrein,n_part
     $        ,rhoinf ,phirein,chi,cfactor,dtin
      endif
!,flux

      end
!*********************************************************************
      subroutine ninjcalc(dtin)
! Given ripernode, decide the number of reinjections per step ninjcomp
! for average edge potential. This is only called at initialization.
! Particle information
      include 'ndimsdecl.f'
      include 'plascom.f'
      include 'partcom.f'
      include 'meshcom.f'
      include 'creincom.f'
      real volume,flux
!
      idebug=0
      do ispecies=1,nspecies
         pinjcompa(ispecies)=0.
! An ion species that has n_part set needs no ninjcalc.
         if(nparta(ispecies).eq.0)then
            if(.not.lreininit)call cinjinit()
! Calculate ninjcomp from ripernode
            volume=1.
            flux=0.
            ripn=ripernode
            do i=1,ndims
               fcarea(i)=1.
! We don't correct area here, because we now count every relocation as
! a reinjection.
               if(lnotallp.and.(ipartperiod(i).eq.4.
     $              or.ipartperiod(i).eq.5))fcarea(i)=1.e-6
               do j=1,ndims-1
                  id=mod(i+j-1,ndims)+1
                  fcarea(i)=fcarea(i)*(xmeshend(id)-xmeshstart(id))
               enddo
! greins contains the reinjection flux per face for each species when set.
               flux=flux+(greins(2*i-1,ispecies)
     $              +greins(2*i,ispecies))*fcarea(i)
               if(gn(i).eq.0.)then
                  volume=volume*(xmeshend(i)-xmeshstart(i))
               else
! Corrected effective volume for density gradient cases.
                  volume=volume*(
     $                   exp((max(xmeshstart(i),xmeshend(i))-gp0(i))
     $                   *gn(i))-
     $                   exp((min(xmeshstart(i),xmeshend(i))-gp0(i))
     $                   *gn(i)) )/
     $                   (gn(i))
               endif
            enddo
            if(ispecies.gt.1)ripn=nparta(1)/volume
            fpinj=ripn*dtin*flux/numratioa(ispecies)
            ninjcompa(ispecies)=int(fpinj)
! Partial reinjection is indicated by this fraction.
            pinjcompa(ispecies)=fpinj-int(fpinj)
            nparta(ispecies)=int(ripn*volume/numratioa(ispecies))
      if(idebug.gt.0)then      
      write(*,*)'ispecies,ripn,dtin,cfactor,flux,nparta,ninjcomp,pinj',
     $        ',volume'
      write(*,'(i2,2f8.3,2f10.3,2i8,f8.4,f8.1)')ispecies,ripn,dtin
     $     ,cfactor,flux,nparta(ispecies),ninjcompa(ispecies)
     $     ,pinjcompa(ispecies),volume
      endif
            if(n_part.gt.n_partmax)then
               write(*,*)'ERROR. Too many particles required.'
               write(*,101)ripn,nparta(ispecies),n_partmax
 101           format('ripernode=',f8.2,'  needs n_part=',i9
     $              ,'  which exceeds n_partmax=',i9)
               stop
            endif
         endif
      enddo
      end
!********************************************************************
      real function smaxflux(uc,chi)
!     Return the total flux to a unit radius sphere from a unit density
!     maxwellian distribution shifted by velocity
      real uc
!     normalized to sqrt(2T/m), in a spherically symmetric potential
!     having a value on the sphere normalized to Ti of minus
      real chi
      real eps,pi
      data eps/1.e-3/pi/3.1415927/
      erf=1.-erfcc(uc)
      sqpi=sqrt(pi)
      if(abs(uc).lt.eps) then
         erfbyu=(2./sqpi)*(1.-uc**2 /3.)
      else
         erfbyu=erf/uc
      endif
      smaxflux=pi*sqrt(2.)*(uc*erf +(0.5+chi)*erfbyu + exp(-uc**2)/sqpi)
      end
!********************************************************************
      subroutine geominit(myid)
! Null
! Silence warnings
      i=myid
      end
!********************************************************************
      subroutine cavereinset(phi)
! Null
      include 'reincom.f'
      include 'ndimsdecl.f'
      include 'partcom.f'
       caverein=crelax*phi+(1.-crelax)*caverein
!      write(*,*)'CAvereinset',caverein
      end
!********************************************************************
!********************************************************************
! Reinjection based upon sample distribution for collisional
! distributions. 
      subroutine colreinject(xr,ipartperiod,cdummy,ispecies)
      implicit none
! Collisional distribution data.
      include 'ndimsdecl.f'
      include 'cdistcom.f'
      real xr(3*ndims)
      real cdummy
      integer ipartperiod(ndims),ispecies
! Reinjection data needed for idrein only, needed for creintest only.
      include 'creincom.f'
      include 'meshcom.f'
! fcarea is in partcom but also has ipartperiod which clashes with arg.
!      include 'partcom.f'
! Local data
      real ra,face,fr,rx
      integer i,id,k,iother,ip

! Choose the normal-dimension for reinjection, from cumulative dist.
 2    continue
      call ranlux(ra,1)
      do i=1,ndims
         id=i
         if(ra.lt.cdistcums(i+1,ispecies))goto 1
      enddo
 1    continue
      idrein=id
! Grab a random v-sample weighted by the normal-dimension mod-v:
      call ranlux(ra,1)
      ra=ra*fxvcols(ncdists(ispecies)+1,id,ispecies)
      call invtfunc(fxvcols(1,id,ispecies),ncdists(ispecies)+1,ra,rx)
! Now the integer part of rx is the index of the chosen particle.
      ip=int(rx)
      do i=1,ndims
         xr(ndims+i)=vcols(i,ip,ispecies)
      enddo
!      if(ispecies.eq.2.and.abs(xr(ndims+2)).gt..001)then
!         write(*,*)'colreinject vynonzero',rx,ip,ra
!         write(*,*)(xr(i),i=1,6)
!         stop
!      endif

! Determine face
      if(ipartperiod(id).ge.3)then
         write(*,*)'Got erroneous periodic dimension',id
         face=0.
         goto 2
      elseif(xr(ndims+id).ge.0.)then
         face=-.999999
         if(ipartperiod(id).eq.1.)goto 1
! Injecting at a periodic/absorbing face. Get a different particle.
      else
         face=.999999
         if(ipartperiod(id).eq.2.)goto 1
      endif
! Got Satisfactory particle. 

! Determine position.
! Position in this dimension (either end), (just) inside the mesh.
      xr(id)=0.5*(xmeshstart(id)*(1.-face)+ xmeshend(id)*(1.+face))
! Pick the velocities and positions parallel to this face:
      do k=1,2
         iother=mod(id+k-1,3)+1
! Position: Ensure we never quite reach the mesh edge:
         call ranlux(ra,1)
         fr=ra*0.999999+.0000005
         xr(iother)=xmeshstart(iother)*(1-fr)+xmeshend(iother)*fr
      enddo
!      if(xr(5).gt..001)then
!         write(*,*)'colreinject ',ra,rx,ip,id,ipartperiod(id)
!         write(*,*)(xr(k),k=1,6)
!      endif
      end

!********************************************************************
      subroutine colreinit(myid,ispecies)
! Initialize and normalize the cdistflux factors from the
! Collisional distribution data, presumed already calculated by pinit.
! Based upon the ipartperiod settings.
      include 'ndimsdecl.f'
      include 'cdistcom.f'
      include 'partcom.f'
      if(ncdists(ispecies).eq.0)return
      ctot=0.
      do i=1,ndims
!         if(myid.eq.0)write(*,*)ipartperiod(i),cdistflux(i)
         if(ipartperiod(i).eq.3.or.
     $        ipartperiod(i).eq.4.or.ipartperiod(i).eq.5)
     $        cdistfluxs(i,ispecies)=0.
! Normalize the cdistflux and cdistcum to face area.
         cdistfluxs(i,ispecies)=cdistfluxs(i,ispecies)*fcarea(i)
! Quiet warnings by explict recognition that cdistfluxs is real*8.
         ctot=ctot+real(cdistfluxs(i,ispecies))
      enddo
      if(ctot.ne.0.)then
         cdistcums(1,ispecies)=0.
         do i=1,ndims
            cdistfluxs(i,ispecies)=cdistfluxs(i,ispecies)/ctot
            cdistcums(i+1,ispecies)=cdistcums(i,ispecies)
     $           +real(cdistfluxs(i,ispecies))
         enddo
! Avoid rounding problems.
         cdistcums(ndims+1,ispecies)=1.
      else
         if(myid.eq.0)write(*,*
     $        )'PROBLEM. colreinject: No reinjection faces'
      endif

! Now evaluate the cumulative distribution in 3 normal-directions.
      do id=1,ndims
         fxvcols(1,id,ispecies)=0.
         do i=1,ncdists(ispecies)
            fxvcols(i+1,id,ispecies)=fxvcols(i,id,ispecies)
     $           +abs(vcols(id,i,ispecies))
         enddo
      enddo
! There's a resolution issue in that a million steps can hardly
! be resolved by single precision. However, it is unlikely that
! substantial statistical distortion will occur. If it did we could
! use double precision.
!      if(myid.eq.0)write(*,'(a,7f8.4)')' colreinit completed',cdistcum
!     $     ,cdistflux
!      if(myid.eq.0)write(*,'(a,3f16.4)')' fxvcol:',(fxvcol(ncdist+1,j),j
!     $     =1,3)
      end
!********************************************************************
      real function ranlenposition(id)
! Return a random fractional position in the coordinate direction id,
! accounting for the density scale length if present
! or for nonuniform background if present.
      include 'ndimsdecl.f'
      include 'meshcom.f'
      include 'creincom.f'
      include 'plascom.f'
      include 'partcom.f'
      real expsa(ndims),expsi(ndims)
      logical lfirst
      data lfirst/.true./
      save lfirst,expsa,expsi

      if(lfirst)then
         do i=1,ndims
            g=gn(i)
            s0=gp0(i)
            si=min(xmeshstart(i),xmeshend(i))-s0
            sa=max(xmeshstart(i),xmeshend(i))-s0
            expsa(i)=exp(g*sa)
            expsi(i)=exp(g*si)
!            write(*,*)'ranlenposition',i,expsa(i),expsi(i),g
         enddo
         lfirst=.false.
      endif
      g=gn(id)
 1    call ranlux(P,1)
      if(abs(g).ne.0)then
!         write(*,*)'Nonuniform plasma ranlenposition'
         sp=gp0(id)+alog(P*expsa(id)+(1.-P)*expsi(id))/g
      else
         sp=(1.-P)*xmeshstart(id)+P*xmeshend(id)
         if(bgmax(id).gt.0.)then
! Nonuniform initialization using a rejection scheme.
            call ranlux(Q,1)
            if(Q.gt.(1.+bgofx(sp,id))/(1.+bgmax(id)))then
!               write(*,*)'Nonuniform pinit',Q,bgmax(id)
               goto 1
            endif
         endif
      endif
      ranlenposition=sp*0.999999+.0000005
      if(.false.)then
         write(*,*)' ranlenpos',id,P,sp,g,expsi,expsa
         write(*,*)'Ranlenpos error',id,P,sp,g,ranlenposition
         write(*,*)expsi,expsa
      endif
      end
!********************************************************************
! Now Obsolete below here.
!********************************************************************
!********************************************************************
!********************************************************************
!********************************************************************
