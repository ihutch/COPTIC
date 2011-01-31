c**********************************************************************
c Cartesian reinjection routine.
      subroutine reinject(xr,ilaunch)
      parameter (mdims=3)
      real xr(3*mdims)
      integer ilaunch

c Random choice data for these routines:
      include 'creincom.f'
c Plasma common data
      include 'plascom.f'
      include 'reincom.f'
c Mesh data
      include 'meshcom.f'
      logical lfirst
      save lfirst
      data lfirst/.true./

      if(lfirst) then
c      write(*,*)'Calling cinjinit'      
         call cinjinit()
         lfirst=.false.
c      write(*,*)'Returned from cinjinit'
      endif

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
c Position on start or end face, sg ensures (just) inside the mesh.
      xr(idrein)=0.5*(xmeshstart(idrein)*(1.+sg)+
     $     xmeshend(idrein)*(1.-sg))

c Pick the velocities and positions parallel to this face:
c Doing position and velocity simultaneously requires velocity distrib
c to be independent of position.
      do k=1,2
         iother=mod(idrein+k-1,3)+1
c Position: Ensure we never quite reach the mesh edge:
         fr=ran1(0)*0.999999+.0000005
         xr(iother)=xmeshstart(iother)*(1-fr)+xmeshend(iother)*fr
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

c      write(*,*)'Picking perpendicular velocities'
c Pick the velocity perpendicular to this face:
      ra=ran1(0)*ncrein
      ir=int(ra)
      fr=ra-ir
c      write(*,*)ra,ir,index,fr,ncrein
      xr(mdims+idrein)=hrein(ir,index)*(1-fr)+hrein(ir+1,index)*fr
c In this version of reinject we never try more than one launch
      ilaunch=1
c Correct the total energy for averein. 
      x2=0.
      do k=1,mdims
         x2=x2+xr(mdims+k)**2
      enddo
      x1=(x2+2.*max(0.,-averein))/x2
      x1=sqrt(x1)
      do k=1,mdims
         xr(mdims+k)=xr(mdims+k)*x1
      enddo

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
      subroutine cumprob(f,xw,xc,K,h,ginfty)
      external f
      real xw,xc
      integer K
      real h(0:K)

      parameter (tiny=1.e-14)
c internal storage
      integer m
      parameter (m=2000)
      real fv(m),g(m),x(m)

      if(K.gt.m)write(*,*)'cumprob warning: too fine a grid ',K

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
c Form \int_x0^{x0+m*xd}
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
      write(*,*) "Cumprob exhausted end adjustment steps"
      ginfty=0.
      return
 2    continue

      if(ii.eq.1 .and. ia.eq.m)then
c We have the right range. Finish.
         if(.not.g(m).gt.tiny)then
            write(*,*)'Cumprob total too small',g(m),' set to zero'
     $           ,' xw=',xw,' xc=',xc
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
c*********************************************************************
      subroutine cinjinit()
c Cartesian reinjection initialization
c drift velocity, vd, in the z-direction
c maxwellians of width given by Ti

      include 'plascom.f'
      include 'meshcom.f'
      include 'creincom.f'
      external ffcrein
      external fvcrein
      parameter (bdys=6.)
      real area(3)

c testing only
c      parameter (nt=1000)
c      real yt(nt)

      do id=1,3
         do i2=1,2
c idrein determines the sign of velocity. id odd => idrein negative.
            idrein=id*(2*i2-3)
            index=2*(id-1)+i2
            call cumprob(ffcrein,0.,0.,
     $           ncrein,hrein(0,index),grein(index))
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
         call cumprob(fvcrein,0.,0.,
     $           ncrein,prein(0,id),gdummy)
      enddo
c      write(*,*)'grein',grein
      gtot=0.
      do id=1,3
         i2=mod(id,3)+1
         i3=mod(id+1,3)+1
         area(id)=abs((xmeshend(i2)-xmeshstart(i2))
     $        *(xmeshend(i3)-xmeshstart(i3)))
      enddo
      do id=1,6
         gtot=gtot+grein(id)*area((id+1)/2)
      enddo
      gintrein(0)=-0.0000005
      do id=1,6
         gintrein(id)=gintrein(id-1) +
     $        1.000001*grein(id)*area((id+1)/2)/gtot
      enddo
      if(.not.gintrein(6).gt.1.)write(*,*)'gintrein problem!'
c      write(*,*)'gintrein',gintrein

      end
c**********************************************************************
      real function ffcrein(v)
c Return the flux from a maxwellian for dimension idrein (in creincom)
c In the positive or negative direction, determined by idrein's sign.
c If abs(idrein) == 3, then maxwellian is shifted by vd (in plascom).
c If v is normalized by sqrt(ZT_e/m_i), then Ti is the ratio T_i/ZT_e.

      include 'plascom.f'
      include 'creincom.f'

      if(int(sign(1.,v)).eq.sign(1,idrein))then
c      if(idrein*v.gt.0)then
         if(abs(idrein).eq.3)then
            ffcrein=abs(v)*exp(-(v-vd)**2/(2.*Ti))
         else
            ffcrein=abs(v)*exp(-v**2/(2.*Ti))
         endif
      else
         ffcrein=0.
      endif
      end
c**********************************************************************
      real function fvcrein(v)
c Return the probability distribution 
c from a maxwellian for dimension idrein (in creincom)
c If abs(idrein) == 3, then maxwellian is shifted by vd (in plascom).
c If v is normalized by sqrt(ZT_e/m_i), then Ti is the ratio T_i/ZT_e.

      include 'plascom.f'
      include 'creincom.f'

      if(abs(idrein).eq.3)then
         fvcrein=exp(-(v-vd)**2/(2.*Ti))
      else
         fvcrein=exp(-v**2/(2.*Ti))
      endif
      end
c*********************************************************************
      subroutine rhoinfcalc(dtin)
c Obtain the rhoinf to be used in calculating the electron shielding,
c based upon the number and average potential of the reinjections.
      include 'plascom.f'
c No time-averaging for now.
c Use particle information for initializing.
      include 'partcom.f'
      include 'meshcom.f'
      real area(ndims_mesh),volume,flux
      real a,td,cfactor
      real chi
      save chi
      data chi/0./

      volume=1.
      flux=0.
      do i=1,ndims_mesh
         area(i)=1.
         do j=1,ndims_mesh-1
            id=mod(i+j-1,ndims_mesh)+1
            area(i)=area(i)*(xmeshend(id)-xmeshstart(id))
         enddo
         a=area(i)*sqrt(2.*Ti/3.1415926)
         if(i.eq.ndims_mesh)then
c Assume vd is in the last dimension
            td=vd/sqrt(2.*Ti)
            a=a*(exp(-td**2)+
     $           0.5*sqrt(3.1415926)*td*(erfcc(-td)-erfcc(td)))
         endif
         flux=flux+a
         volume=volume*(xmeshend(i)-xmeshstart(i))
      enddo
c      if(nrein.ne.0)then
c Better to use a significant number to avoid bias at low reinjections.
      if(nrein.ge.10)then
c Calculate rhoinf from nrein if there are enough.
c Correct approximately for edge potential depression (OML).
c         chi=min(-phirein/Ti,0.5)
         chi=crelax*(-phirein/Ti)+(1.-crelax)*chi
         cfactor=smaxflux(vd/sqrt(2.*Ti),chi)
     $        /smaxflux(vd/sqrt(2.*Ti),0.)
         rhoinf=(nrein/(dtin*cfactor*flux))
      else
         if(rhoinf.lt.1.e-4)then
c Approximate initialization
            rhoinf=numprocs*n_part/(volume)
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
      include 'plascom.f'
c No time-averaging for now.
c Particle information
      include 'partcom.f'
      include 'meshcom.f'
      real area(ndims_mesh),volume,flux
      real a,td,cfactor
c 
c Calculate ninjcomp from ripernode
      volume=1.
      flux=0.
      do i=1,ndims_mesh
         area(i)=1.
         do j=1,ndims_mesh-1
            id=mod(i+j-1,ndims_mesh)+1
            area(i)=area(i)*(xmeshend(id)-xmeshstart(id))
         enddo
         a=area(i)*sqrt(2.*Ti/3.1415926)
         if(i.eq.ndims_mesh)then
c Assume vd is in the last dimension
            td=vd/sqrt(2.*Ti)
            a=a*(exp(-td**2)+
     $           0.5*sqrt(3.1415926)*td*(erfcc(-td)-erfcc(td)))
         endif
         flux=flux+a
         volume=volume*(xmeshend(i)-xmeshstart(i))
      enddo
c Correct approximately for edge potential depression (OML).
      chi=crelax*min(-phirein/Ti,0.5)
      cfactor=smaxflux(vd/sqrt(2.*Ti),chi)/smaxflux(vd/sqrt(2.*Ti),0.)
      ninjcomp=nint(ripernode*dtin*cfactor*flux)
      nrein=ninjcomp*numprocs
      if(ninjcomp.le.0)ninjcomp=1
      n_part=ripernode*volume
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
      subroutine avereinset(phi)
c Null
      include 'reincom.f'
      include 'partcom.f'
c      averein=phi
      averein=crelax*phi+(1.-crelax)*averein
      end
c********************************************************************
c Cubic interpolation. Not used.
      real function yinterp4pt(y,xi)
c Given four points on an equally spaced mesh: y(1-4)
c and a fractional mesh position xi between the middle two,
c so that xi=x-x(2), [and we form et=1-xi],
c return the interpolated y value, cubically interpolated.
      real y(4)
      real xi

      et=1.-xi
      a= -et*xi*y(1)+(2.+2.*xi-3.*xi**2)*y(2)
      b= -xi*et*y(4)+(2.+2.*et-3.*et**2)*y(3)
      yinterp4pt=0.5*(et*a+xi*b)
      end
