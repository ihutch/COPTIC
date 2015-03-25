c*********************************************************************
c Wrapper:
      subroutine reinject(xr,ilaunch,ispecies)
      real xr(*)
      include 'colncom.f'
      include 'ndimsdecl.f'
      include 'partcom.f'

c This test ought to be just whether the distribution is separable.
c And probably it ought to be decided just once.
c      if(colntime.eq.0. .or. colpow.eq.0.)then
      if(notseparable(ispecies).eq.0)then
c Cartreinject is the version to be used when the distribution is 
c separable in the directions of the cartesian axes.
         call cartreinject(xr,ilaunch,caverein,ispecies)
      else
c colreinject started as the collisional reinjection scheme. 
c It is based upon sampling from a 3-D distribution function.
c It can therefore be used for situations where the distribution
c is not separable in coordinate directions.
         call colreinject(xr,ipartperiod,caverein,ispecies)
c Only one launch allowed here (and actually in the other too).
         ilaunch=1
      endif
      end
**********************************************************************
c Initialize whether we are initialized!
      block data creinset
      include 'ndimsdecl.f'
      include 'creincom.f'
      data lreininit/.false./
      end
**********************************************************************
c Cartesian reinjection routine.
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
c Random choice data for these routines:
      include 'creincom.f'
c Plasma common data
      include 'plascom.f'
c Mesh data
      include 'meshcom.f'
c Particle data to define the species of interest

      if(.not.lreininit)call cinjinit()
c      write(*,*)'Returned from cinjinit'
c----------------------------------------
c Pick the face from which to reinject.
      call ranlux(ra,1)
      call invtfunc(gintreins(0,ispecies),7,ra,fc)
c Returns fc between 1+ and 7-.
      index=int(fc)
c Make idrein 1-3, and sg +-1. 
c Plus/minus refers to the direction of velocity.
c So + means plane xstart, and - xend relative to a normally ordered mesh.
      idrein=(index+1)/2
c Just less than +-1.
      sg=-(2*mod(index,2)-1)*0.999999
c Position in this dimension (either end), sg ensures (just) inside the
c mesh.
      xr(idrein)=0.5*(xmeshstart(idrein)*(1.+sg)+
     $     xmeshend(idrein)*(1.-sg))
c----------------------------------------
c Pick the velocities and positions parallel to this face:
c Doing position and velocity simultaneously requires velocity distrib
c to be independent of position.
      do k=1,2
         iother=mod(idrein+k-1,3)+1
c Position: Ensure we never quite reach the mesh edge:
c Or Accounting for density gradients:
         xr(iother)=ranlenposition(iother)
c Velocity
         call ranlux(ra,1)
         ra=ra*ncrein
         ir=int(ra)
         fr=ra-ir
         if(fr.gt.1. .or. fr.lt.0.)then
c This should never happen.
            write(*,*)'Creinject fraction wrong',ra,ir,fr
         endif
c velocity, linear interpolation:
         xr(mdims+abs(iother))= preins(ir,iother,ispecies)*(1-fr)
     $        +preins(ir+1,iother,ispecies)*fr
      enddo
c----------------------------------------
c Pick the velocity perpendicular to this face:
      call ranlux(ra,1)
      ra=ra*ncrein
      ir=int(ra)
      fr=ra-ir
c      write(*,*)ra,ir,index,fr,ncrein
      xr(mdims+idrein)=hreins(ir,index,ispecies)*(1-fr)+hreins(ir+1
     $     ,index,ispecies)*fr
c In this version of reinject we never try more than one launch
      ilaunch=1
c----------------------------------------
c Correct the total energy for caverein. 
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
c*********************************************************************
      subroutine cinjinit()
c Cartesian reinjection initialization for drift velocity, vd, in the
c z-direction and maxwellians of width given by Ts

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
c For three coordinate directions
      do id=1,3
c Width of velocity distribution
         xw=sqrt(Ts(ispec))*abs(eoverms(ispec))
c    For two ends (and hence velocity polarities)
         do i2=1,2
c idrein determines the sign of velocity. i2 odd => idrein negative.
            idrein=id*(2*i2-3)
            index=2*(id-1)+i2
            if(abs(theta).lt.thetamax)then
c Set the inverse cumulative probability fn hrein, and flux grein
               call cumprob(ffdrein,xw,xc,
     $           ncrein,hreins(0,index,ispec),greins(index,ispec),myid)
            else
c Infinite-Bt case 1-d projection.
               xc=vpars(ispec)*Bfield(id)+vperps(id,ispec)
               xw=xw*Bfield(id)
c               write(*,*)'Calling cumprob',id,xc,xw
               call cumprob(ff1crein,xw,xc,
     $           ncrein,hreins(0,index,ispec),greins(index,ispec),myid)
            endif
c Zero grein on absorbing faces.
            ip=ipartperiod(id)
c Ordering in grein etc is (+,-) for each dim. Upper face/Lower face. 
c That's the opposite of what I adopted for ipartperiod because
c it corresponds to negative/positive velocity. Pity!
            ip=ip/2**(2-i2)
            ip=ip-2*(ip/2)
            if(ip.eq.1)greins(index,ispec)=0.
            if(gnt.ne.0)then
c Correcting for uniform density scale-length.
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
c Kludge fix of ends to avoid negative velocity injections.
            if(idrein.gt.0)then
               if(hreins(0,index,ispec).lt.0.)hreins(0,index,ispec)=0.
            else
               if(hreins(ncrein,index,ispec).gt.0.)hreins(ncrein,index
     $              ,ispec)=0.
            endif
c            write(*,*)index,(hreins(kk,index,ispec),kk=0,5)
c     $           ,(hreins(kk,index,ispec),kk=ncrein-4,ncrein)
         enddo
         idrein=id
         if(theta.lt.thetamax)then
            call cumprob(fvdrein,xw,xc,
     $           ncrein,preins(0,id,ispec),gdummy,myid)
         else
            call cumprob(fv1crein,xw,xc,
     $           ncrein,preins(0,id,ispec),gdummy,myid)
         endif
      enddo
c
c      grein(1)=1.*grein(1)
c      write(*,*)'grein arbirtrary adjustment',greins
      gtot=0.
c Alternative general-dimension fcarea calculation:
      do i=1,ndims
         fcarea(i)=1.
c Set area tiny for period directions
         if(lnotallp.and.ipartperiod(i).eq.4)fcarea(i)=1.e-6
         do j=1,ndims-1
            id=mod(i+j-1,ndims)+1
            fcarea(i)=fcarea(i)*abs(xmeshend(id)-xmeshstart(id))
         enddo
c         write(*,*)'fcarea(',i,')=',fcarea(i)
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
c      write(*,'(a,3i2,a)')'ipartperiod',ipartperiod,'  greins,gintrein:'
c      write(*,'(6f10.5)')greins
c      write(*,'(7f10.5)')gintrein

      lreininit=.true.
      
      end
c**********************************************************************
c Given a function f(x), whose integral from x=-\infty to +infty exists,
c obtain the definite integral g(x) = \int_-\infty^x f(x) dx.
c Return the array h_j (j=0,K) equally spaced on g, such that 
c    g(h_j)=g(\infty) (j/K)
c Thus h_j is the inverse of the cumulative probability distribution,
c g/g(\infty) based on f(x). Also return g(\infty).
c On entry an estimate of the width of f is provided in xw
c                  and of the center of f in xc
c The range is adjusted automatically to give a decently accurate integral.
c It is possible for the adjustment process to fail if xc is in error by, 
c or xw is, a large factor times the actual width.
c The most efficient usage is for xw to be slightly larger than the 
c points at which f is 10^-2 times its peak.
c f must be declared external in the calling routine.
c If ginfty is returned as zero, then the routine failed with cumulative
c probability too small to be significant. 
c If myid .gt. 0 then don't print warning messages. 
c If myid .lt. 0 print extra warnings.
      subroutine cumprob(f,xw,xc,K,h,ginfty,myid)
      external f
      real xw,xc
      integer K
      real h(0:K)

      parameter (tiny=1.e-14)
c internal storage
      integer m
      parameter (m=2000)
      real fv(m),g(m),x(m)

      if(K.gt.m .and. myid.le.0)
     $     write(*,*)'cumprob warning: too fine a grid ',K

c      write(*,*)'cumprob entry',xw,xc,K

c If we return prematurely, it is with h's=0.
      do i=0,K
         h(i)=0.
      enddo

      xr=abs(xw)
      if(xr.eq.0)xr=1.
      x0=xc-2.*xr
      x1=xc+2.*xr
      icount=0

c Cycle over integrations.
 1    continue
      icount=icount+1
      xd=(x1-x0)/m
      g(1)=0.
      x(1)=x0
      fv(1)=f(x0)
c Form the integral \int_x0^{x0+m*xd}
      do i=2,m
         x(i)=x0+i*xd
         fv(i)=f(x(i))
         g(i)=g(i-1)+xd*(fv(i)+fv(i-1))*.5
      enddo
c      write(*,*)'Integration limits',x0,x1,' Value=',g(m)

c Decide whether we have covered enough range.
c Criteria are: dg0/gm < 1/(10.K.m) ; dg(m/4)/gm > 1/(10.K.m)
c and similarly for top end.  dg=fv*xd
c If dg0 is violated, double distance from middle.
c If dg(m/4) is violated, halve it. 

      ii=1
      ia=m
      dgc=g(m)/(10.*K*m*xd)
      isearch=0

c Find acceptable end points, for this integral
c allowing only 20 adjustment (per cycle)
      do j=1,20
         isearch=isearch+1
         ll=0
         lr=0

c         write(*,'(a,2i5,3g12.4)')'ii,ia,x0,x1,gm',ii,ia,x0,x1,g(m)
         ii2=(3*ii+ia)/4
         ia2=(ii+3*ia)/4
         if(fv(ii).ge.dgc)then
c We cover not enough -x
            ii=ii-(ii2-ii)
         elseif(fv(ii2).lt.dgc)then
c We cover too much -x         
            ii=ii2
         else
            ll=1
         endif

         if(fv(ia).ge.dgc)then
c We cover too little +x
            ia=ia+(ia-ia2)
         elseif(fv(ia2).lt.dgc)then
c We cover too much +x
            ia=ia2
         else
            lr=1
         endif
c If both ends were acceptable this iteration
         if(ll.eq.1 .and. lr.eq.1)then
c If we haven't moved the ends break; we are done.
            if(isearch.eq.1) goto 2
c Else reintegrate with these end points (which should be acceptable). 
            x0=x0+(ii-1)*xd
            x1=x1+(ia-m)*xd
c            write(*,*)'New ends',x0,x1
            goto 1
         endif

         if(ii.lt.1 .or. ia.gt.m)then
c Range is too short. Need to reintegrate from scratch
c Double or triple the range.
            dx=x1-x0
            if(ii.lt.1) x0=x0-dx
            if(ia.gt.m) x1=x1+dx
            if(icount.gt.20)then
c               write(*,'(a,/,a,g12.4,a,g12.4,a)')
c     $              'Integrate too high icount.',' Starting range'
c     $              ,xc,'+-',xw,' may be too wide.'
               ginfty=0.
               return
c               stop
            endif
            goto 1
         endif

         if(abs(ia-ii).lt.4)then
c Range is too big.
            x1=x0+(ii+2)*xd
            x0=x0+(ii-2)*xd
            if(icount.gt.20)then
c               write(*,*)'Integrate count too high.'
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
c We have the right range. Finish.
         if(.not.g(m).gt.tiny)then
            if(myid.le.0)write(*,*)'Cumprob total too small',g(m)
     $           ,' set to zero',' xw=',xw,' xc=',xc
            ginfty=0.
            return
         endif
         do i=0,K
c solve g(h_i)=gm*i/K with linear interpolation.
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
c            if(i.le.1.or.i.ge.(K-1))write(*,*)'i,ihi,h(i)',i,ihi,h(i)
         enddo
      else
         goto 1
      endif
c      write(*,*)h
c      call yautoplot(g,m)
c      call pltend()
c      call yautoplot(h(0),K+1)
c      call pltend()
      ginfty=g(m)
      end
c******************************************************************
c Routines for reinjection calculations with ions drifting relative
c to neutrals, driven by external force (e.g. E-field).
c Direct replacements for the fvcrein, ffcrein functions.
c*******************************************************************
      real function fvdrein(v)
c Return the probability distribution for reinjection,
c in coordinate direction idrein (signed) in creincom,
c from a drift-distribution shifted by vd (in plascom),
c whose direction cosines are vdrift(3) (in plascom),
c colliding with neutrals of velocity vneutral (in colncom).
c If v is normalized by sqrt(ZT_e/m_i), then Tsi is the ratio T_i/ZT_e.
c Because the drift distribution is not separable except in the directions
c perpendicular and parallel to the drift, only degenerate Maxwellian
c (vd-vneutral=0) non-z cases are allowed so far.
c Diamagnetic drift is implemented as a Maxwellian shift.
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

      ud=vd-vneutral
      if(ispec.ne.1)ud=0.
      vdia=0.
c Are there background diamagnetic drifts? Density gradient. 
      if(gnt.ne.0.)then
c Gradient is always perpendicular to B.
         id2=mod(abs(idrein),3)+1
         id3=mod(abs(idrein)+1,3)+1
c v_dia = -(T_i \nabla n \times B / qnB^2)
         vdia=- (gn(id2)*Bfield(id3)-gn(id3)*Bfield(id2))*Ti/Bt
c The approximation is to represent the perpendicular distribution by
c a Maxwellian shifted by the diamagnetic drift.
      endif
      vn=sqrt(2.*Ts(ispec)*abs(eoverms(ispec)))
      if(vdrift(3).eq.1.)then
c Z-drift cases. (equiv old)
         if(abs(idrein).eq.3)then
            u=(v-vd+ud-vdia)/vn
            ud=ud/vn
            fvdrein=fvcx(u,ud)
            fvdrein=fvdrein/vn
         else
            fvdrein=exp(-((v-vdia)/vn)**2)/(vn*sqrt(3.1415926))
         endif
      else
c Non-z
         if(ud.ne.0.)then
            write(*,*)'Non-z collisional drift not implemented.'
     $           ,' Aborting.'
            stop
         else
c Maxwellian non-z drift.
            fvdrein=exp(-((v-vd*vdrift(abs(idrein))-vdia)/vn)**2)
     $           /(vn*sqrt(3.1415926))
         endif
      endif
      end
c*******************************************************************
      real function ffdrein(v)
c Return the one-way differential flux distribution as a fn of velocity
c from a drift distribution (in dimension 3).
c In the positive or negative direction, determined by idrein's sign.
c Diamagnetic drift is implemented as a Maxwellian shift.
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

      ud=vd-vneutral
      if(ispec.ne.1)ud=0.
      vdia=0.
c Are there background diamagnetic drifts? Density gradient. 
      if(gnt.ne.0.)then
c Gradient is always perpendicular to B.
         id2=mod(abs(idrein),3)+1
         id3=mod(abs(idrein)+1,3)+1
c v_dia = -(T_i \nabla n \times B / qnB^2)
         vdia=- (gn(id2)*Bfield(id3)-gn(id3)*Bfield(id2))*Ti/Bt
c The approximation is to represent the perpendicular distribution by
c a Maxwellian shifted by the diamagnetic drift.
      endif
c      if(vdia.ne.0..and.abs(v-1.).lt.0.0005)then
c         write(*,*)'v,ffdrein,idrein,vdia',v,ffdrein,idrein,vdia,id2,id3
c      endif
      vn=sqrt(2.*Ts(ispec)*abs(eoverms(ispec)))
      if(vdrift(3).eq.1.)then
c Z-drift cases. (equiv old)
         if(abs(idrein).eq.3)then
c In z-direction use appropriate drift distribution.
            u=(v-vd+ud-vdia)/vn
            ud=ud/vn
            ffdrein=fvcx(u,ud)
            ffdrein=abs(v)*ffdrein/vn
         else
            ffdrein=abs(v)*exp(-((v-vdia)/vn)**2)/(vn*sqrt(3.1415926))
         endif
      else
c Non-z
         if(ud.ne.0.)then
            write(*,*)'Non-z collisional drift not implemented.'
     $           ,' Aborting.'
            stop
         else
c Maxwellian non-z drift.
            ffdrein=exp(-((v-vd*vdrift(abs(idrein))-vdia)/vn)**2)
     $           *abs(v)/(vn*sqrt(3.1415926))
         endif
      endif
c Don't count particles going the wrong way:
      if(int(sign(1.,v)).ne.sign(1,idrein))ffdrein=0.

      end
c****************************************************************
c FVCX function for 1-d drifting CX distribution.
      function fvcx(u,ud)
      real u,ud,v,vd,fvcx
c Return the normalized distribution function f(u)=v_n f(v) for constant
c cx collision frequency at a value of normalized velocity u=v/v_n, when
c the normalized drift velocity is ud= (a/\nu_c) /v_n, with v_n =
c sqrt(2T_n/m). a is acceleration, nu_c collision freq.  This is the
c solution of the steady Boltzmann equation.
      if(ud.lt.0.) then
         v=-u
         vd=-ud
      else
         v=u
         vd=ud
      endif
      if(vd.eq.0.)then
         carg=20.
         earg=100
      else
         carg=0.5/vd-v
         earg=(0.5/vd)**2-v/vd
      endif
      if(carg.gt.10)then
c asymptotic form for large exp argument (small vd):
c  exp(-v^2)/[sqrt(\pi)(1-2 v_d v)]:
         fvcx=exp(-v**2)/1.77245385/(1.-2.*vd*v)
      elseif(carg.gt.-5.)then
         fvcx=exp(-v**2)*experfcc(carg)*0.5/vd
      else
c         fvcx=exp(earg)*erfcc(carg)*0.5/vd
         fvcx=exp(earg)/vd
      endif
c      write(*,*)'fvcx:vd,v,earg,fvcx',vd,v,earg,fvcx
      if(.not.fvcx.ge.0) then
         write(*,*)'fvcx error. u=',u,' ud=',ud,' f=',fvcx,carg
         fvcx=0.
      endif
      end
c**********************************************************************
      real function ff1crein(v)
c This is the flux function for 1-D motion in the direction of Bfield.
c
c Return the flux from a maxwellian for dimension idrein (in creincom)
c In the positive or negative direction, determined by idrein's sign.
c Maxwellian is shifted by vdj=Bfield(j)*vpar+vperp(j).
c If v is normalized by sqrt(ZT_e/m_i), then Tsi is the ratio T_i/ZT_e.
c But here, the thermal spread along Bfield must be projected into
c the coordinate direction idrein. 

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
         vs=Bfield(j)*vpar+vperp(j)
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
c            if(Bfield(j).eq.0.)write(*,*)'v,vs,arg',v,vs,arg
c Corrected the scaling to be consistent:
            ff1crein=scale*abs(v)*exp(arg)/sqrt(2.*Tsi*3.1415926)
         endif
      else
         ff1crein=0.
      endif
      end
c**********************************************************************
      real function fv1crein(v)
c This is the 1-d probability distribution projected in coordinate
c direction idrein (in crein).
c
c Return the probability distribution value.

      include 'ndimsdecl.f'
      include 'plascom.f'
      include 'creincom.f'
      real vt2min,argmax
      parameter (vt2min=1.e-5,argmax=12.)
      integer ispec
      common /species/ispec

      j=abs(idrein)
      vs=Bfield(j)*vpar+vperp(j)
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
c*********************************************************************
      subroutine rhoinfcalc(dtin)
c Obtain the rhoinf to be used in calculating the electron shielding,
c based upon the number and average potential of the reinjections.
      include 'ndimsdecl.f'
      include 'plascom.f'
c No time-averaging for now.
c Use particle information for initializing.
      include 'partcom.f'
      include 'meshcom.f'
      include 'creincom.f'
      include 'myidcom.f'
      real volume,flux
      real cfactor

      volume=1.
      flux=0.

      do i=1,ndims
         fcarea(i)=1.
         if(lnotallp.and.ipartperiod(i).eq.4)fcarea(i)=1.e-6
         do j=1,ndims-1
            id=mod(i+j-1,ndims)+1
            fcarea(i)=fcarea(i)*(xmeshend(id)-xmeshstart(id))
         enddo
c         write(*,*)'fcarea(',i,')=',fcarea(i)
c Multiply the face area by the factor relative to unshifted maxwellian
c that specifies the flux for this face, and add to total. 
         flux=flux+(grein(2*i-1)+grein(2*i))*fcarea(i)
         volume=volume*(xmeshend(i)-xmeshstart(i))
      enddo
c      if(nrein.ne.0)then
c Better to use a significant number to avoid bias at low reinjections.
      if(nrein.ge.10)then
c Calculate rhoinf from nrein if there are enough.
c Correct approximately for edge potential depression (OML).
c         chi=min(-phirein/Ti,0.5)
         chi=max(crelax*(-phirein/Ti)+(1.-crelax)*chi,0.)
         cfactor=smaxflux(vd/sqrt(2.*Ti),chi)
     $        /smaxflux(vd/sqrt(2.*Ti),0.)
         rhoinf=(nrein/(dtin*cfactor*flux))
c         write(*,*)nrein,dtin,phirein,chi,cfactor,flux,rhoinf
      else
         if(rhoinf.lt.1.e-4)then
c Approximate initialization
            rhoinf=numprocs*n_part/volume
            write(*,*)'Rhoinf in rhoinfcalc approximated as',rhoinf
     $           ,numprocs,n_part
         endif
c Else just leave it alone.
      endif
c      write(*,*)
c      if(myid.eq.0)write(*,'(a,2i8,10f9.4)')
c     $ 'Ending rhoinfcalc',nrein,n_part,rhoinf
c     $     ,phirein,chi,cfactor,dtin
c,flux
      end
c*********************************************************************
      subroutine ninjcalc(dtin)
c Given ripernode, decide the number of reinjections per step ninjcomp
c for average edge potential. This is only called at initialization.
c Particle information
      include 'ndimsdecl.f'
      include 'plascom.f'
      include 'partcom.f'
      include 'meshcom.f'
      include 'creincom.f'
      real volume,flux
c 
      do ispecies=1,nspecies
         pinjcompa(ispecies)=0.
c An ion species that has n_part set needs no ninjcalc.
         if(nparta(ispecies).eq.0)then
            if(.not.lreininit)call cinjinit()
c Calculate ninjcomp from ripernode
            volume=1.
            flux=0.
            ripn=ripernode
            do i=1,ndims
               fcarea(i)=1.
c We don't correct area here, because we now count every relocation as
c a reinjection.
               if(lnotallp.and.ipartperiod(i).eq.4)fcarea(i)=1.e-6
               do j=1,ndims-1
                  id=mod(i+j-1,ndims)+1
                  fcarea(i)=fcarea(i)*(xmeshend(id)-xmeshstart(id))
               enddo
c greins contains the reinjection flux per face for each species when set.
               flux=flux+(greins(2*i-1,ispecies)
     $              +greins(2*i,ispecies))*fcarea(i)
               volume=volume*(xmeshend(i)-xmeshstart(i))
            enddo
            if(ispecies.gt.1)ripn=nparta(1)/volume
            fpinj=ripn*dtin*flux/numratioa(ispecies)
            ninjcompa(ispecies)=int(fpinj)
c Partial reinjection is indicated by this fraction.
            pinjcompa(ispecies)=fpinj-int(fpinj)

            nparta(ispecies)=int(ripn*volume/numratioa(ispecies))

c      write(*,*)'ispecies,ripn,dtin,cfactor,flux,nparta,ninjcomp'
c      write(*,*) ispecies,ripn,dtin,cfactor,flux
c     $     ,nparta(ispecies),ninjcompa(ispecies)
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
c********************************************************************
      real function smaxflux(uc,chi)
c     Return the total flux to a unit radius sphere from a unit density
c     maxwellian distribution shifted by velocity
      real uc
c     normalized to sqrt(2T/m), in a spherically symmetric potential
c     having a value on the sphere normalized to Ti of minus
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
c********************************************************************
      subroutine geominit(myid)
c Null
c Silence warnings
      i=myid
      end
c********************************************************************
      subroutine cavereinset(phi)
c Null
      include 'reincom.f'
      include 'ndimsdecl.f'
      include 'partcom.f'
       caverein=crelax*phi+(1.-crelax)*caverein
c      write(*,*)'CAvereinset',caverein
      end
c********************************************************************
c********************************************************************
c Reinjection based upon sample distribution for collisional
c distributions. 
      subroutine colreinject(xr,ipartperiod,cdummy,ispecies)
      implicit none
c Collisional distribution data.
      include 'ndimsdecl.f'
      include 'cdistcom.f'
      real xr(3*ndims)
      real cdummy
      integer ipartperiod(ndims),ispecies
c Reinjection data needed for idrein only, needed for creintest only.
      include 'creincom.f'
      include 'meshcom.f'
c fcarea is in partcom but also has ipartperiod which clashes with arg.
c      include 'partcom.f'
c Local data
      real ra,face,fr,rx
      integer i,id,k,iother,ip

c Choose the normal-dimension for reinjection, from cumulative dist.
 2    continue
      call ranlux(ra,1)
      do i=1,ndims
         id=i
         if(ra.lt.cdistcums(i+1,ispecies))goto 1
      enddo
 1    continue
      idrein=id
c Grab a random v-sample weighted by the normal-dimension mod-v:
      call ranlux(ra,1)
      ra=ra*fxvcols(ncdists(ispecies)+1,id,ispecies)
      call invtfunc(fxvcols(1,id,ispecies),ncdists(ispecies)+1,ra,rx)
c Now the integer part of rx is the index of the chosen particle.
      ip=int(rx)
      do i=1,ndims
         xr(ndims+i)=vcols(i,ip,ispecies)
      enddo
c      if(ispecies.eq.2.and.abs(xr(ndims+2)).gt..001)then
c         write(*,*)'colreinject vynonzero',rx,ip,ra
c         write(*,*)(xr(i),i=1,6)
c         stop
c      endif

c Determine face
      if(ipartperiod(id).ge.3)then
         write(*,*)'Got erroneous periodic dimension',id
         face=0.
         goto 2
      elseif(xr(ndims+id).ge.0.)then
         face=-.999999
         if(ipartperiod(id).eq.1.)goto 1
c Injecting at a periodic/absorbing face. Get a different particle.
      else
         face=.999999
         if(ipartperiod(id).eq.2.)goto 1
      endif
c Got Satisfactory particle. 

c Determine position.
c Position in this dimension (either end), (just) inside the mesh.
      xr(id)=0.5*(xmeshstart(id)*(1.-face)+ xmeshend(id)*(1.+face))
c Pick the velocities and positions parallel to this face:
      do k=1,2
         iother=mod(id+k-1,3)+1
c Position: Ensure we never quite reach the mesh edge:
         call ranlux(ra,1)
         fr=ra*0.999999+.0000005
         xr(iother)=xmeshstart(iother)*(1-fr)+xmeshend(iother)*fr
      enddo
c      if(xr(5).gt..001)then
c         write(*,*)'colreinject ',ra,rx,ip,id,ipartperiod(id)
c         write(*,*)(xr(k),k=1,6)
c      endif
      end

c********************************************************************
      subroutine colreinit(myid,ispecies)
c Initialize and normalize the cdistflux factors from the
c Collisional distribution data, presumed already calculated by pinit.
c Based upon the ipartperiod settings.
      include 'ndimsdecl.f'
      include 'cdistcom.f'
      include 'partcom.f'
      if(ncdists(ispecies).eq.0)return
      ctot=0.
      do i=1,ndims
c         if(myid.eq.0)write(*,*)ipartperiod(i),cdistflux(i)
         if(ipartperiod(i).ge.3)cdistfluxs(i,ispecies)=0.
c Normalize the cdistflux and cdistcum to face area.
         cdistfluxs(i,ispecies)=cdistfluxs(i,ispecies)*fcarea(i)
c Quiet warnings by explict recognition that cdistfluxs is real*8.
         ctot=ctot+real(cdistfluxs(i,ispecies))
      enddo
      if(ctot.ne.0.)then
         cdistcums(1,ispecies)=0.
         do i=1,ndims
            cdistfluxs(i,ispecies)=cdistfluxs(i,ispecies)/ctot
            cdistcums(i+1,ispecies)=cdistcums(i,ispecies)
     $           +real(cdistfluxs(i,ispecies))
         enddo
c Avoid rounding problems.
         cdistcums(ndims+1,ispecies)=1.
      else
         if(myid.eq.0)write(*,*
     $        )'PROBLEM. colreinject: No reinjection faces'
      endif

c Now evaluate the cumulative distribution in 3 normal-directions.
      do id=1,ndims
         fxvcols(1,id,ispecies)=0.
         do i=1,ncdists(ispecies)
            fxvcols(i+1,id,ispecies)=fxvcols(i,id,ispecies)
     $           +abs(vcols(id,i,ispecies))
         enddo
      enddo
c There's a resolution issue in that a million steps can hardly
c be resolved by single precision. However, it is unlikely that
c substantial statistical distortion will occur. If it did we could
c use double precision.
c      if(myid.eq.0)write(*,'(a,7f8.4)')' colreinit completed',cdistcum
c     $     ,cdistflux
c      if(myid.eq.0)write(*,'(a,3f16.4)')' fxvcol:',(fxvcol(ncdist+1,j),j
c     $     =1,3)
      end
c********************************************************************
      real function ranlenposition(id)
c Return a random fractional position in the coordinate direction id,
c accounting for the density scale length if present
      include 'ndimsdecl.f'
      include 'meshcom.f'
      include 'creincom.f'
      include 'plascom.f'
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
c            write(*,*)'ranlenposition',i,expsa(i),expsi(i),g
         enddo
         lfirst=.false.
      endif
      g=gn(id)
      call ranlux(P,1)
      if(abs(g).ne.0)then
c         write(*,*)'Nonuniform ranlenposition'
         sp=gp0(id)+alog(P*expsa(id)+(1.-P)*expsi(id))/g
      else
         sp=(1.-P)*xmeshstart(id)+P*xmeshend(id)
      endif
      ranlenposition=sp*0.999999+.0000005
      if(.false.)then
         write(*,*)' ranlenpos',id,P,sp,g,expsi,expsa
         write(*,*)'Ranlenpos error',id,P,sp,g,ranlenposition
         write(*,*)expsi,expsa
      endif
      end
c********************************************************************
c Now Obsolete below here.
c********************************************************************
c********************************************************************
c********************************************************************
c********************************************************************
