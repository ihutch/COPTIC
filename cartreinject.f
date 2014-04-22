c*********************************************************************
c Wrapper:
      subroutine reinject(xr,ilaunch,cdummy)
      include 'colncom.f'
      include 'ndimsdecl.f'
      include 'partcom.f'
      if(colntime.eq.0. .or. colpow.eq.0.)then
         call cartreinject(xr,ilaunch,caverein)
      else
         call colreinject(xr,ipartperiod,caverein)
c Only one launch allowed here (and actually in the other too).
         ilaunch=1
      endif
      end
**********************************************************************
c Initialize whether we are initialized!
      block data creinset
      include 'creincom.f'
      data lreininit/.false./
      end
**********************************************************************
c Cartesian reinjection routine.
      subroutine cartreinject(xr,ilaunch,caverein)
      implicit none
      integer mdims
      parameter (mdims=3)
      real xr(3*mdims)
      real caverein
      integer ilaunch
      real ra,fc,sg,fr,x2,x1
      integer index,k,iother,ir
      real ran1
      external ran1
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
      ra=ran1(0)
      call invtfunc(gintrein,7,ra,fc)
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
c         fr=ran1(0)*0.999999+.0000005
c         if(fr.le.0 .or. fr.ge.1.)write(*,*)'fr=',fr
c         xr(iother)=xmeshstart(iother)*(1-fr)+xmeshend(iother)*fr
c Or Accounting for density gradients:
         xr(iother)=ranlenposition(iother)
c Velocity
         ra=ran1(0)*ncrein
         ir=int(ra)
         fr=ra-ir
         if(fr.gt.1. .or. fr.lt.0.)then
c This should never happen.
            write(*,*)'Creinject fraction wrong',ra,ir,fr
         endif
c velocity, linear interpolation:
         xr(mdims+abs(iother))=
     $        prein(ir,iother)*(1-fr)+prein(ir+1,iother)*fr
c velocity, cubic/parabolic interpolation:
c Works very badly. Don't use.
c         xr(mdims+abs(iother))=yinterp4pt(prein(ir-1,iother),fr)
      enddo
c----------------------------------------
c Pick the velocity perpendicular to this face:
      ra=ran1(0)*ncrein
      ir=int(ra)
      fr=ra-ir
c      write(*,*)ra,ir,index,fr,ncrein
      xr(mdims+idrein)=hrein(ir,index)*(1-fr)+hrein(ir+1,index)*fr
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
c Cartesian reinjection initialization
c drift velocity, vd, in the z-direction
c maxwellians of width given by Ti

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

      xc=0.
      xw=sqrt(Ti)
c For three coordinate directions
      do id=1,3
c    For two ends (and hence velocity polarities)
         do i2=1,2
c idrein determines the sign of velocity. i2 odd => idrein negative.
            idrein=id*(2*i2-3)
            index=2*(id-1)+i2
            if(Bt.lt.Btinf)then
c Set the inverse cumulative probability fn hrein, and flux grein
               call cumprob(ffdrein,xw,xc,
     $           ncrein,hrein(0,index),grein(index),myid)
            else
c Infinite-Bt case 1-d projection.
               xc=vpar*Bfield(id)+vperp(id)
               xw=sqrt(Ti)*Bfield(id)
               call cumprob(ff1crein,xw,xc,
     $           ncrein,hrein(0,index),grein(index),myid)
            endif
c Zero grein on absorbing faces.
            ip=ipartperiod(id)
c Ordering in grein etc is (+,-) for each dim. Upper face/Lower face. 
c That's the opposite of what I adopted for ipartperiod because
c it corresponds to negative/positive velocity. Pity!
            ip=ip/2**(2-i2)
            ip=ip-2*(ip/2)
            if(ip.eq.1)grein(index)=0.
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
               grein(index)=grein(index)*dengfac
            endif
c Kludge fix of ends to avoid negative velocity injections.
            if(idrein.gt.0)then
               if(hrein(0,index).lt.0.)hrein(0,index)=0.
            else
               if(hrein(ncrein,index).gt.0.)hrein(ncrein,index)=0.
            endif
c            write(*,*)index,(hrein(kk,index),kk=0,5)
c     $           ,(hrein(kk,index),kk=ncrein-4,ncrein)
         enddo
         idrein=id
         if(Bt.lt.Btinf)then
            call cumprob(fvdrein,xw,xc,
     $           ncrein,prein(0,id),gdummy,myid)
         else
            call cumprob(fv1crein,xw,xc,
     $           ncrein,prein(0,id),gdummy,myid)
         endif
      enddo
c
c      write(*,*)'grein',grein
      gtot=0.
c Alternative general-dimension fcarea calculation:
      do i=1,ndims
         fcarea(i)=1.
         if(lnotallp.and.ipartperiod(i).eq.4)fcarea(i)=1.e-6
         do j=1,ndims-1
            id=mod(i+j-1,ndims)+1
            fcarea(i)=fcarea(i)*abs(xmeshend(id)-xmeshstart(id))
         enddo
c         write(*,*)'fcarea(',i,')=',fcarea(i)
      enddo
      do id=1,6
         gtot=gtot+grein(id)*fcarea((id+1)/2)
      enddo
      gintrein(0)=-0.0000005
      do id=1,6
         gintrein(id)=gintrein(id-1) +
     $        1.000001*grein(id)*fcarea((id+1)/2)/gtot
      enddo
      if(.not.gintrein(6).gt.1.)write(*,*)'gintrein problem!'
c      write(*,*)'ipartperiod',ipartperiod,' grein',grein
c      write(*,*)'gintrein',gintrein
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
c If v is normalized by sqrt(ZT_e/m_i), then Ti is the ratio T_i/ZT_e.
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
      if(vdrift(1).eq.0)then
c Z-drift cases. (equiv old)
         if(abs(idrein).eq.3)then
            vn=sqrt(2.*Tneutral)
            ud=(vd-vneutral)/vn
            u=(v-vneutral-vdia)/vn
            fvdrein=fvcx(u,ud)
            fvdrein=fvdrein/vn
         else
            vn=sqrt(2.*Ti)
            fvdrein=exp(-((v-vdia)/vn)**2)/(vn*sqrt(3.1415926))
         endif
      else
c Non-z
         ud=(vd-vneutral)
         if(ud.ne.0.)then
            write(*,*)'Non-z collisional drift not implemented.'
     $           ,' Aborting.'
            stop
         else
c Maxwellian non-z drift.
            vn=sqrt(2.*Ti)
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
      if(vdrift(1).eq.0.and.vdrift(2).eq.0.)then
c Z-drift cases. (equiv old)
         if(abs(idrein).eq.3)then
c In z-direction use appropriate drift distribution.
            vn=sqrt(2.*Tneutral)
            ud=(vd-vneutral)/vn
            u=(v-vneutral-vdia)/vn
            ffdrein=fvcx(u,ud)
            ffdrein=abs(v)*ffdrein/vn
         else
            vn=sqrt(2.*Ti)
            ffdrein=abs(v)*exp(-((v-vdia)/vn)**2)/(vn*sqrt(3.1415926))
         endif
      else
c Non-z
         ud=(vd-vneutral)
         if(ud.ne.0.)then
            write(*,*)'Non-z collisional drift not implemented.'
     $           ,' Aborting.'
            stop
         else
c Maxwellian non-z drift.
            vn=sqrt(2.*Ti)
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
c****************************************************************
c This is exp(X^2)*erfc(X)
      FUNCTION expERFCC(X)
      Z=ABS(X)      
      T=1./(1.+0.5*Z)
      expERFCC=T*EXP(-1.26551223+T*(1.00002368+T*(.37409196+
     *    T*(.09678418+T*(-.18628806+T*(.27886807+T*(-1.13520398+
     *    T*(1.48851587+T*(-.82215223+T*.17087277)))))))))
      IF (X.LT.0.) expERFCC=2.*exp(z**2)-expERFCC
      END
c*********************************************************************
c**********************************************************************
      real function ff1crein(v)
c This is the flux function for 1-D motion in the direction of Bfield.
c
c Return the flux from a maxwellian for dimension idrein (in creincom)
c In the positive or negative direction, determined by idrein's sign.
c Maxwellian is shifted by vdj=Bfield(j)*vpar+vperp(j).
c If v is normalized by sqrt(ZT_e/m_i), then Ti is the ratio T_i/ZT_e.
c But here, the thermal spread along Bfield must be projected into
c the coordinate direction idrein. 

      include 'ndimsdecl.f'
      include 'plascom.f'
      include 'creincom.f'
      real vt2min,argmax
      parameter (vt2min=1.e-6,argmax=12.)
      
      if(int(sign(1.,v)).eq.sign(1,idrein))then
         j=abs(idrein)
         vs=Bfield(j)*vpar+vperp(j)
         vt2=2.*Ti*Bfield(j)**2
         if(vt2.gt.vt2min)then
            arg=-(v-vs)**2/vt2
            scale=1./Bfield(j)
         else
            arg=-(v-vs)**2/vt2min
            scale=sqrt(2.*Ti/vt2min)
         endif
         if(abs(arg).gt.argmax)then
            ff1crein=0.
         else
c            if(Bfield(j).eq.0.)write(*,*)'v,vs,arg',v,vs,arg
c Corrected the scaling to be consistent:
            ff1crein=scale*abs(v)*exp(arg)/sqrt(2.*Ti*3.1415926)
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
      
      j=abs(idrein)
      vs=Bfield(j)*vpar+vperp(j)
      vt2=2.*Ti*Bfield(j)**2
      if(vt2.gt.vt2min)then
         arg=-(v-vs)**2/vt2
         scale=1./Bfield(j)
      else
         arg=-(v-vs)**2/vt2min
         scale=sqrt(2.*Ti/vt2min)
      endif
      if(abs(arg).gt.argmax)then
         fv1crein=0.
      else
         fv1crein=scale*exp(arg)/sqrt(2.*Ti*3.1415926)
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
c The sum of the grein for this dimension ought to be equal to 
c fonefac*sqrt(2.*Ti/\pi). And generally it is, unless B=\infty.
c         ffc=fonefac(i)
c         ratio=(grein(2*i-1)+grein(2*i))
c     $    /(ffc+1.e-8)/sqrt(2.*Ti/3.1415926)
c         write(*,*)'i,fonefac(i),grein',i,ffc ,grein(2*i-1),grein(2*i)
c     $        ,ratio
c         flux=flux+fonefac(i)*fcarea(i)*sqrt(2.*Ti/3.1415926)
c This change to get rid of fonefac produces tiny rounding differences.
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
c         write(*,*)nrein,dtin,chi,cfactor,flux,rhoinf
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
c      write(*,'(a,2i8,10f9.4)')
c     $ 'Ending rhoinfcalc',nrein,n_part,rhoinf
c     $     ,phirein,chi,cfactor,dtin
c,flux
      end
c*********************************************************************
      subroutine ninjcalc(dtin)
c Given ripernode, decide the number of reinjections per step ninjcomp
c for average edge potential. This is only called at initialization.
c No time-averaging for now.
c Particle information
      include 'ndimsdecl.f'
      include 'plascom.f'
      include 'partcom.f'
      include 'meshcom.f'
      include 'creincom.f'
      real volume,flux
      real cfactor
c 
      if(.not.lreininit)call cinjinit()
c Calculate ninjcomp from ripernode
      volume=1.
      flux=0.
      do i=1,ndims
         fcarea(i)=1.
c We don't correct area here, because we now count every relocation as
c a reinjection.
         if(lnotallp.and.ipartperiod(i).eq.4)fcarea(i)=1.e-6
         do j=1,ndims-1
            id=mod(i+j-1,ndims)+1
            fcarea(i)=fcarea(i)*(xmeshend(id)-xmeshstart(id))
         enddo
c         flux=flux+fonefac(i)*fcarea(i)*sqrt(2.*Ti/3.1415926)
         flux=flux+(grein(2*i-1)+grein(2*i))*fcarea(i)
         volume=volume*(xmeshend(i)-xmeshstart(i))
      enddo
c Correct approximately for edge potential depression (OML).
      chi=crelax*min(-phirein/Ti,0.5)
      cfactor=smaxflux(vd/sqrt(2.*Ti),chi)/smaxflux(vd/sqrt(2.*Ti),0.)
      ninjcomp=nint(ripernode*dtin*cfactor*flux)
c      write(*,*)'ripernode,dtin,cfactor,flux,ninjcomp',ripernode,dtin
c     $     ,cfactor,flux,ninjcomp
      nrein=ninjcomp*numprocs
      if(ninjcomp.le.0)ninjcomp=1
      n_part=int(ripernode*volume)
      rhoinf=ripernode*numprocs
      if(n_part.gt.n_partmax)then
         write(*,*)'ERROR. Too many particles required.'
         write(*,101)rhoinf,n_part,n_partmax
 101     format('rhoinf=',f8.2,'  needs n_part=',i9
     $        ,'  which exceeds n_partmax=',i9)
         stop
      endif
c      write(*,*)'Ending ninjcalc',rhoinf,nrein,n_part

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
      subroutine colreinject(xr,ipartperiod,cdummy)
      implicit none
c Collisional distribution data.
      include 'ndimsdecl.f'
      include 'cdistcom.f'
      real xr(3*ndims)
      real cdummy
      integer ipartperiod(ndims)
c Reinjection data needed for idrein only, needed for creintest only.
      include 'creincom.f'
      include 'meshcom.f'
c Local data
      real ra,ran1,face,fr,rx
      integer i,id,k,iother
      external ran1

c Choose the normal-dimension for reinjection, from cumulative dist.
 2    ra=ran1(1)
      do i=1,ndims
         id=i
         if(ra.lt.cdistcum(i+1))goto 1
      enddo
 1    continue
      idrein=id
c Grab a random v-sample. Needs to be changed to reflect not just a 
c random grab,
c      call colvget(xr(ndims+1))
c but a grab weighted by the normal-dimension mod-v:
      ra=ran1(1)*fxvcol(ncdist+1,id)
      call invtfunc(fxvcol(1,id),ncdist+1,ra,rx)
c Now the integer part of rx is the index of the chosen particle.
      k=int(rx)
      do i=1,ndims
         xr(ndims+i)=v_col(i,k)
      enddo

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
         fr=ran1(0)*0.999999+.0000005
         xr(iother)=xmeshstart(iother)*(1-fr)+xmeshend(iother)*fr
      enddo
      end

c********************************************************************
      subroutine colreinit(myid)
c Initialize and normalize the cdistflux factors from the
c Collisional distribution data, presumed already calculated by pinit.
c Based upon the ipartperiod settings.
      include 'ndimsdecl.f'
      include 'cdistcom.f'
      include 'partcom.f'
      if(ncdist.eq.0)return
      ctot=0.
      do i=1,ndims
c         if(myid.eq.0)write(*,*)ipartperiod(i),cdistflux(i)
         if(ipartperiod(i).ge.3)cdistflux(i)=0.
         ctot=ctot+cdistflux(i)
      enddo
      if(ctot.ne.0.)then
         cdistcum(1)=0.
         do i=1,ndims
            cdistflux(i)=cdistflux(i)/ctot
            cdistcum(i+1)=cdistcum(i)+cdistflux(i)
         enddo
c Avoid rounding problems.
         cdistcum(ndims+1)=1.
      else
         if(myid.eq.0)write(*,*
     $        )'PROBLEM. colreinject: No reinjection faces'
      endif

c Now evaluate the cumulative distribution in 3 normal-directions.
      do id=1,ndims
         fxvcol(1,id)=0.
         do i=1,ncdist
            fxvcol(i+1,id)=fxvcol(i,id)+abs(v_col(id,i))
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
c            write(*,*)i,expsa(i),expsi(i),g
         enddo
         lfirst=.false.
      endif
      g=gn(id)
      P=ran1(0)
      if(abs(g).ne.0)then
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
