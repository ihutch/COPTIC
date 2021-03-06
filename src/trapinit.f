!****************************************************************************
! Optional main to run test program.
!      call untrappeddentest
!      call untrapcumtest
!      end
!***********************************************************************
! Initializing particles with a trapped hole region.  
! 
! The 1-d parallel direction is that of B if it is coordinate-aligned,
! or else dimension 1 (x).
!
! Hole parameters [switch-mnemonic,defaults] are 
! holepsi, the peak hole potential [psi,0] 
! holeum, the drift speed of the hole relative to the distribution [u,0]
! holelen, the parallel length of the hole [l,4*Debyelen]
! holeeta, the requested power of u determining trapped distribution shape
!          2 gives parabolic. <=0 gives flattop [h,2]
! holepow, the power governing trapped transverse temperature variation [p,1]
! holerad, the transverse hole radius [r,0]
! holespeed, the rest-frame hole speed [derived: =holeum+vds]
!
! The hole is created by adjusting the particle density consistent with
! a potential profile that is sech^4(x/4l) in the parallel direction.
! Adjustment is achieved by a rejection scheme based upon the analytic
! density variation for specified peak potential psi, which gives the
! density as a function of phi. When the hole is transversely uniform,
! no transverse density gradient is anyway induced by rejection.
!
! However, to avoid generating a transverse density gradient when the
! hole is transverse-localized, the rejection scheme retries until
! successful, but with the transverse position unchanged. Therefore the
! average density has no gradient in the transverse direction. Ideally
! one wants the external density to have no gradient, so this is close
! but not quite correct.
!
! If Tperp is different from T, then if holepow=0 this is interpreted as
! the transverse temperature of only the trapped region and the
! untrapped Tperp is set back equal to T. If holepow!=0 all Tperps are set.
!
      subroutine trapinit(sprior)
      implicit none
! Common data:
      include 'ndimsdecl.f'
      include 'partcom.f'
      include 'plascom.f'
      include 'myidcom.f'
      include 'colncom.f'
      include 'cdistcom.f'
      real sprior
      external linregion,ranlenposition,gasdev,interp
      real ranlenposition,gasdev
      integer interp
      logical linregion
! Local dummy variables for partlocate.
      real xfrac(ndims)
      integer ixp(ndims)
      logical linmesh
! Local declarations
      integer i,i1,islotmax,ispecies,j,iregion,k,ko,nonzero,ntries,id,iv
      integer iaccept,ireject,jj,ntrapcount
      real theta,tisq,denmax,dentot,vr,v,vs,Ttrans,eopsi,Ttr,vid,r2
! Factor by which we leave additional space for increased number of 
! particles if not fixed:
      real slotsurplus
      real thetamax
      parameter (thetamax=1.)
! Arguments of untrapcum
      real um,eta,psi,phi
      integer nt,nin,np,nm,ntot
      parameter (nin=400,nt=16,ntot=nin+nt)
      real fpa(nin),fma(nin),ua(nin)
      real cump(-ntot:ntot),cumv(-ntot:ntot)
! Internal diagnostics of the initialized distribution
! Collect particles in bins from -vdrange to +vdrange as initialized.
      integer idiag,npbins,ibin
      parameter (npbins=50)
      real vdrange,xcount(npbins,ndims),vdist(npbins),emax
      
      idiag=0
      vdrange=4.
      emax=0.7
      if(idiag.ne.0)then
         do i=1,npbins
            vdist(i)=-vdrange+2.*vdrange*(i-1.)/(npbins-1.)
            do j=1,ndims
               xcount(i,j)=0.
            enddo
         enddo
      endif
      npassthrough=0
      iaccept=0
      ireject=0
      ntrapcount=0
! silence a spurious warning.
      r2=0.
      eopsi=0.
!-----------------------------------------------------------------
      i1=1
! Point to the bottom of the particle stack for start of species 1.
      iicparta(1)=1
! Decide whether the external distribution is separable.
      do ispecies=1,nspecies
! Conveniently here initialize distribution numbers.
         ncdists(ispecies)=0
         if(ninjcompa(ispecies).gt.0)then
            slotsurplus=1.3
         else
            slotsurplus=1.
         endif
         ntries=0
         Eneutral=0.
! Set tperp to a tiny number if species is infinitely magnetized.
         theta=Bt*eoverms(ispecies)*dt
         if(abs(theta).gt.thetamax)Tperps(ispecies)=1.e-24
         notseparable(ispecies)=0
         Ttrans=Tperps(ispecies)
         if(colntime.ne.0..and.ispecies.eq.1)then
            if(myid.eq.0)write(*,*)'Collisional holes not implemented'
     $           ,colntime
            stop
         elseif(Tperps(ispecies).gt.1.e-20
     $           .and.Tperps(ispecies).ne.Ts(ispecies)
     $           .and.holepow.ne.0.)then
! Interpret Tperp unequal to T as being the trapped-only particle 
! transverse temperature. 
            if(myid.eq.0)write(*,'(a,f8.4,a,f8.4,a,i3,/,a,f8.4)')
     $           ' Tperp',Tperps(ispecies),' !=Ts',Ts(ispecies)
     $           ,' interpreted as trapped only. Species',ispecies
     $           ,' Resetting passing Tperp as isotropic.'
     $           ,Ts(ispecies)
            Tperps(ispecies)=Ts(ispecies)
         else
! Count the number of non-zero vdrift components. If it is more than one
! then vdrift is not along a coordinate axis => nonseparable.
            nonzero=0
            do k=1,ndims
               if(vdrifts(k,ispecies).ne.0)nonzero=nonzero+1
            enddo
            if(nonzero.gt.1)then
               if(myid.eq.0)
     $   write(*,*)'Non-separable oblique vdrift holes not implemented'
               stop
            endif
         endif
! ---------- Finished notseparable checks ------------------
! Scale the number of particles for higher slots to the ions(?)
         nparta(ispecies)=nparta(1)/numratioa(ispecies)
!     $        *sqrt(abs(eoverms(1)/eoverms(ispecies)))
         islotmax=nparta(ispecies)+iicparta(ispecies)-1
         if((islotmax+1)*slotsurplus.gt.n_partmax)then
            write(*,*)'Too few particle slots',n_partmax
     $           ,' max, to accommodate',islotmax
            write(*,*)'multiplied by slotsurplus factor'
     $           ,slotsurplus,int(islotmax*slotsurplus)
            stop
         endif
! The following will actually be overridden later:
         tisq=sqrt(Ts(ispecies)*abs(eoverms(ispecies)))
! Determine the trapping dimension
         id=0
         do i=1,ndims
            if(Bfield(i).ne.0.)then
               if(id.eq.0)then
                  id=i
               else
                  write(*,*)'More than one Bfield component nonzero',i
     $                 ,id,Bfield(i)
                  id=0
               endif
            endif
         enddo
! Default x
         if(id.eq.0)id=1
! Calculate holespeed using vds component in projection dimension.
! Holeum is minus f drift speed relative to hole so holeum=-vds+holespeed:
         holespeed=holeum+vds(ispecies)*vdrift(id)
         if(myid.eq.0.and.ispecies.eq.hspecies)
     $        write(*,'(2a,i2,a,f7.3,a,f7.3,/,a)')' Hole particle'
     $        ,' initialization. Trapping id',id,'. Speed',holespeed
     $        ,' psi',holepsi
     $        ,' Please wait...'
! ----------------------------- Actual Particle Setting ----------
         do i=iicparta(ispecies)+i1-1,islotmax
            x_part(iflag,i)=0
            r2=0.
 1          continue
            ntries=ntries+1
! Position choice including density gradients.
            if(x_part(iflag,i).eq.1)then  ! Retry
! For transverse-localized hole it is essential not to move the retried
! particle in the transverse direction. 
               x_part(id,i)=ranlenposition(id)
            else  ! First time get all three positions.
               call ranlen3position(x_part(1,i),id)
               do j=0,ndims-2 ! calculate transverse radius
                  r2=r2+x_part(mod(id+j,ndims)+1,i)**2
               enddo
            endif
            x_part(iflag,i)=1
            if(.not.linregion(ibool_part,ndims,x_part(1,i)))then
!     If we are not in the plasma region, try again if we are doing fixed
!     particle number else just set this slot empty.
               x_part(iflag,i)=0
               if(ninjcompa(ispecies).ne.0)then
! The only place a slot is left empty.
               else
                  goto 1
               endif
            endif
! Shifted Gaussians um is in units sqrt(2T/m), holeum in sqrt(T/m)
            um=holeum/sqrt(2.)
            psi=holepsi
            if(holerad.gt.0)then
! Scale the potential in transverse direction.
! I wonder if this is allowable or if there's a scaling problem when
! the hole is not uniform in transverse direction.
               psi=psi*max((1.-r2/holerad**2),0.)
!                  write(*,*)'psi',psi,'r2',r2
            endif
            do jj=1,ndims
! Start with the id direction, then the others.
               j=mod(id+jj-2,ndims)+1
               if(j.eq.id.and.psi.gt.0..and.ispecies.eq.hspecies)then
! Parallel distribution. Get potential. Decide if to reject.
! Decide velocity component.
                  call getholepotl(psi,holelen,phi,x_part(id,i),id)
! Total electron density based upon the sech^4 potential variation
                  dentot=1.+ phi-(5./4./sqrt(psi))*phi**1.5
! Value of dentot at its max (when phi/psi=(8/15)^2)
                  denmax=1+(64./45.)*psi
! Reject if random*denmax >= dentot at this phi:
                  call ranlux(vr,1)
                  vr=vr*denmax
!                  write(*,*)'dentot,denmax,vr',dentot,denmax,vr
                  if(vr.lt.dentot)then
! Accept
                     iaccept=iaccept+1
                     eta=holeeta
! Calculate the probability distribution for both untrapped and trapped
! velocities in this direction at this phi.
                     call untrapcum(phi,um,psi,nin,np,nm,fpa,fma,ua,nt
     $                    ,cump,cumv,eta)
! Now cump contains the offset cumulative probability function at cumv.
! It begins at -nm-nt. We decide the velocity using the rejection random.
! But I don't see why this is correct. At this point, vr is randomly
! uniformly distributed between 0 and dentot; shift it to cump minimum.
!                     vs=vr+cump(-nm-nt)
! Find the position in the array where cump=vs. Why is cump normalized
! such that this is the right choice? This seems an error prior to Nov
! 17. It was not obvious in results because actually the cump range is
! close to 1 and denmax is close to 1, so the selection is not too
! bad. Also no moving holes were investigated.  We ought to scale vs as
! follows.
                     vs=vr*(cump(np+nt)-cump(-nm-nt))/dentot
     $                    +cump(-nm-nt)
! However, printing out the following shows that the new vr scaling
! is only 1.0002 times the old, so the correction is negligible.
!                     if(mod(i,1000).eq.0)write(*,'(3f8.4)')
!     $                    (cump(np+nt)-cump(-nm-nt))/dentot,dentot
                     iv=interp(cump(-nm-nt),np+nm+2*nt+1,vs,v)
                     if(iv.eq.0)then
                        write(*,*)'vr overflow',cump(-nm-nt),vr,vs
     $                       ,cump(np+nt),denmax,dentot,cump(np+nt)
     $                       -cump(-nm-nt)
                        stop
                     endif 
                     v=v-iv
! Make iv the absolute (not relative) index of the cumulative arrays.
                     iv=iv-1-nm-nt
! Interpolate back the velocity.
                     vid=(cumv(iv)*(1-v)+cumv(iv+1)*(v))
! Units of cumv are sqrt(2T/m). Add holespeed in sqrt(T/m) units.
                     x_part(ndims+id,i)=vid*sqrt(2.)+holespeed
! Now we have selected the particle velocity 
! Calculate the id energy in the hole frame and hence the transverse
! temperature for cases where it is anisotropic.  Here charge is
! presumed negative, psi positive (since several places the sqrt of psi
! is taken).
                     Ttr=Tperps(ispecies)
                     if(phi.gt.0)then
! eopsi=1 at psi and zero at sepx.
                        eopsi=(phi-vid**2)/psi
                        if(eopsi.lt.0.)eopsi=0.
                        if(holepow.ne.0.)then
                           eopsi=eopsi**holepow
                           Ttr=(Ttrans*eopsi
     $                          +Tperps(ispecies)*(1.-eopsi))
!                        write(*,*)'phi=',phi,' eopsi=',eopsi,' Ttr=',Ttr
                        endif
                        if(eopsi.gt.0.)ntrapcount=ntrapcount+1
                     endif
                     tisq=sqrt(Ttr*abs(eoverms(ispecies)))
                  else
! Reject. Try again.
                     ireject=ireject+1
                     goto 1
                  endif
               else
! Transverse directions unaffected by hole potential.
                  x_part(ndims+j,i)=tisq*gasdev(myid)
     $                 + vds(ispecies)*vdrift(j)
               endif
            enddo
! Selective velocity diagnostics.
            if(idiag.ne.0.and.eopsi.gt.emax)then
!            if(idiag.ne.0.and.phi.gt..5*psi)then
               do j=1,ndims
                  ibin=1+nint((npbins-1)*(x_part(ndims+j,i)+vdrange)
     $                 /(2.*vdrange))
                  if(ibin.lt.1)ibin=1
                  if(ibin.gt.npbins)ibin=npbins
                  xcount(ibin,j)=xcount(ibin,j)+1
               enddo
            endif
! The previous timestep length.
            x_part(idtp,i)=0.
! Initialize the mesh fraction data in x_part.
            call partlocate(x_part(1,i),ixp,xfrac,iregion,linmesh)
! This test rejects particles exactly on mesh boundary:
            if(.not.linmesh)goto 1
! 2          continue
         enddo
         if(idiag.ne.0)then
! Plot the diagnostic histograms
            call multiframe(3,1,2)
            call autoplot(vdist,xcount(1,1),npbins)
            call axlabels(' ','v!d1!d')
            call autoplot(vdist,xcount(1,2),npbins)
            call axlabels(' ','v!d2!d')
            call autoplot(vdist,xcount(1,3),npbins)
            call axlabels(' ','v!d3!d')
            call pltend()
         endif
!         write(*,*)'Accept,Reject,Ntrap',iaccept,ireject,ntrapcount
!------------------------------- End of Actual Particle Setting --
! The maximum used slot for this species
         iocparta(ispecies)=i-1
! Start of next slot-set may give a gap for overflow.
         iicparta(ispecies+1)=int(i*slotsurplus)
! Zero the overflow slots' flag
         do i=iocparta(ispecies)+1,iicparta(ispecies+1)-1
            x_part(iflag,i)=0
         enddo
         if(myid.eq.0)then
            write(*,101)ispecies,nprocs,
     $           iocparta(ispecies)-iicparta(ispecies)+1,ntries
 101        format(' Initialized species',i2,i4,'x',i7
     $           ,' ntries=',i7,$)
            if(nspecies.gt.0)write(*,'(a,3i8)')' Slots'
     $           ,iicparta(ispecies),iocparta(ispecies)
!     $           ,ninjcompa(ispecies)
         endif
! Initialize orbit tracking
         do ko=1,norbits
            iorbitlen(ko)=0
         enddo
         
! Don't shift for special particle the subsequent sections.
         i1=1
      enddo
! Set flag of unused slots to 0
      do i=iicparta(ispecies),n_partmax
         x_part(iflag,i)=0
      enddo
! Allow the last species to fill the array:
      iicparta(ispecies)=n_partmax+1
      end
!****************************************************************************
!****************************************************************************
!*********************************************************************
! Untrapped density function
!  untrappedden(\phi,u_m) = {2\over \sqrt{\pi}}\int_{\sqrt{\phi}}^\infty
!     \exp(-u^2+\phi-u_m^2)\cosh(2u_m\sqrt{u^2-\phi}) du.
!   = {1\over\sqrt{\pi}} \int_{\sqrt{\phi}}^\infty \sum_{\pm}
!             \exp(-[\pm\sqrt{u^2-\phi}-u_m]^2)du
! written as f = \int g du.

! This is the density (normalized to the background density) of
! untrapped particles having a Maxwellian background shifted by a speed
! -um, measured in units of sqrt(2T/m), at a potential energy -phi
! measured in units of T/e.

! The integration is carried out on uniform u-grid, initially with
! velocity spacing dui. The spacing is halved until relative difference
! from the prior integral is less than a bit more than df.  [Some sort
! of Gaussian integration would probably be faster but tricky because of
! the infinite range.]

! If nin is >0, on entry, then return the integrated distributions
! for positive and negative velocity in fpa,fma, arrays of (max) size
! nin. Also return the corresponding speed (magnitude) in ua.
! Typically nin=300 ought to be enough to avoid precision compromise.
! On exit np and nm label the upper bounds above which the integrals are 
! essentially complete. There is no data above the larger of np,nm.
! If nin is zero on entry then subsequent arguments may be omitted 
! (but warnings may arise from the compiler).

      real function untrappedden(phi,um,nin,fpa,fma,ua,np,nm)
      real phi,um
      integer nin,np,nm
      real fpa(0:nin-1),fma(0:nin-1),ua(0:nin-1)

      parameter (dui=.04,df=1.e-4,nppar=10000)
!      data iwarn/0/

      nmax=nppar
      sqpi=sqrt(3.1415926)
! silence warnings; not really needed.
      f=0.
      fp=0.
      fm=0.
      v=0.
      f1=0
      g1=1.
      g1p=g1
      g1m=g1
      gp=0.
      gm=0.

      um2=um**2
      du=dui
      if(phi.lt.0)stop 'Negative phi in untrappedden call not allowed'
! Iterate over steps sizes
      do k=1,5
! Integrate
         do i=0,nmax
            u=sqrt(phi)+i*du
            u2=u**2
            if(i.eq.0)then
               g=exp(-um2)
               gp=g
               gm=g
               g1=g
               if(g.eq.0)stop 'untrappedden error um2 overflow'
               f=0.
               fp=0.
               fm=0.
            else
               gp=exp(-(sqrt(u2-phi)+abs(um))**2)
               gm=exp(-(sqrt(u2-phi)-abs(um))**2)
               g=0.5*(gp+gm)
               fp=fp+(u-v)*(gp+g1p)*0.5
               fm=fm+(u-v)*(gm+g1m)*0.5
            endif
            if(.not.g.ge.0 .or. .not.g.lt.1.e30)then
! Error trap
               write(*,*)'u,um,phi,g error',u,um,phi,g
               g=0.
            endif
!            write(*,'(10f8.4)')f,fp,fm,g,gp,gm
            f=fp+fm
            g1p=gp
            g1m=gm
! If new contributions are negligible, break
            if(g.lt.g1*df)goto 1
            v2=u2
            v=u
            if(gp.lt.g1*df.and.gm.lt.g1*df)
     $           write(*,*)'g,gp,gm,g1',g,gp,gm,g1
            if(i.lt.nin)then
! Return the cumulative distributions.
               fpa(i)=fp/sqpi
               fma(i)=fm/sqpi
               ua(i)=u
               if(gp.gt.g1*df)np=i
               if(gm.gt.g1*df)nm=i
            endif
         enddo
         write(*,*)'untrappedden exhausted number of steps np',i
 1       continue
! Converged when the totals don't differ by a bit more than df. Break
         if(abs(f1-f).le.1.5*df)then
!            write(*,*)'Converged',k,i,f-f1,g
            goto 2
         elseif(nin.gt.0.and.i.gt.nin/2)then
! Next time we will over-run, so break and warn.
            write(*,'(a,2i4,a,2f8.4,g11.3)')
     $           'Precision compromised by nin',nin,i
     $           ,' u,f,Df=',u,f,f-f1
            goto 2
         endif
         f1=f
         du=du/2.
      enddo
      write(*,*)'untrappedden exhausted step iterations',k,f1,f,du
      write(*,*)'Steps, du, u, f, um',i,du,u,f,um
      write(*,*)u,g,f
 2    continue

      untrappedden=f/sqpi

!      if(iwarn.eq.0)then
!         call autoplot(ua,fma,nm)
!         call axlabels('u','fma')
!         call polyline(ua,fpa,np)
!         call pltend()
!         iwarn=iwarn+1
!      endif
!      if(nin.ne.0)write(*,*)'i,np,nm,fpa(np),fma(nm),ua(np)'
!     $     ,i,np,nm,fpa(np),fma(nm),ua(np)

      end
!*********************************************************************
      subroutine untrappeddentest
      integer nphi,num
      real phimax
      parameter (nphi=1000,num=5,phimax=1.,ustep=.2,fbinit=1.)
      real phi(0:nphi-1),f(0:nphi-1,0:num-1)
!      real ef(0:nphi-1)
      real fiu(0:nphi-1,0:num-1),fb(0:nphi-1,0:num-1)
      real Vhatp(0:nphi-1,0:num-1)
      real refden(0:nphi-1,0:num-1),Vhatf(0:nphi-1,0:num-1)
      real fb1(0:num-1),uma(0:num-1)
      real sq(0:nphi),us(0:nphi)
      parameter (nu=300)
      real fpa(nu),fma(nu),ua(nu)
      character*20 string
      real Kcap

      phistep=phimax/(nphi-1)
      denom=0.

! Do over values of um index km.
      write(*,*)'km,  i,  phi(i),   f(i,km),    fiu(i,km),'
     $     ,' fb(i,km),   Vhatp,   sum,   denom,'
      do km=0,num-1
! To diagnose asymmetric cases only:    +ustep
         um=(km)*ustep 
         expum2=exp(-um**2)
         Kcap=2.*sqrt(3.1415926)/phistep**1.5/expum2
         write(*,*)'phistep=',phistep,' Kcap=',Kcap
         fb(0,km)=fbinit
         do i=0,nphi-1
            sq(0)=0.
            phi(i)=i*phistep
            us(i)=sqrt(phi(i))
! Get the untrapped electron density at this potential and drift.
            nin=0
            nin=nu
            f(i,km)=untrappedden(phi(i),um,nin,fpa,fma,ua,idum,idum)
! Add to reference flat-top trapped density to get total reference
            refden(i,km)=2.*expum2*sqrt(phi(i)/3.1415926)+f(i,km)
            sum=0.
            if(i.eq.0)then           
               fiu(i,km)=0.
               Vhatp(i,km)=0.
               Vhatf(i,km)=0.
               fb(i,km)=fbinit
            else
! Integrate a step in \phi to find untrapped electron integrated density.
               fiu(i,km)=fiu(i-1,km)
     $              +(phi(i)-phi(i-1))*(f(i,km)+f(i-1,km))*0.5
! And Vhatp
               Vhatp(i,km)=Vhatp(i-1,km)
     $              +(phi(i)-phi(i-1))*(1.-(f(i,km)+f(i-1,km))*0.5)
! And minus Vhatf, which is the reference density including immobile ions.
               Vhatf(i,km)=Vhatf(i-1,km) -(phi(i)-phi(i-1))
     $              *(1.-(refden(i,km)+refden(i-1,km))*0.5)
! Now we have the fiu (Vhatp) value.
! Step the solution of the integral equation. 
!               l=i
! Store sqrt(i)
               sq(i)=sqrt(float(i))
               fack=2.
               do k=0,i
                  if(k.eq.i)fack=1.
                  md=0
                  do m=0,k-1
                     n=i+m-k
                     sum=sum+fack*(sq(m+1)-sq(m-md))*(fb(n,km))
                     md=1
                  enddo
! Here we need to calculate the trapped density from this integral.
! and store it in trapden.
               enddo
               denom=(Kcap*Vhatp(i,km)-sq(i)-sq(i-1))
               fb(i,km)=(sum)/denom
            endif
 301        format(2i4,5f10.3,4e10.2)
            if(i.lt.10)
     $           write(*,301)km,i,phi(i),f(i,km),fiu(i,km),fb(i,km)
     $           ,Vhatp(i,km),sum,denom
     $           
         enddo

         if(km.eq.0)then
            call pfset(3)
            call autoplot(phi,f,nphi)
            call axlabels('!Af!@ [e/T!de!d]','n!dp!d [/n!db!d]')
         else
            call color(km)
!      call dashset(2)
            write(*,*)km,nphi
            call polyline(phi,f(0,km),nphi)
         endif
         string='u!dm!d='
         call fwrite(um,iwidth,2,string(8:))
         call jdrwstr(wx2nx(phi(nphi-1)),wy2ny(f(nphi-1,km))
     $           ,string,-1.)
      enddo
      call pltend()

      call pltinit(0.,phi(nphi-1),.85,1.6)
      call axis()
      call axlabels('!Af!@ [e/T!de!d]','n!df!d [/n!db!d]')
      do km=0,num-1
         if(km.ne.0)call color(km)
         call polyline(phi,refden(0,km),nphi)
         string='u!dm!d='
         call fwrite((km)*ustep,iwidth,2,string(8:))
         call jdrwstr(wx2nx(phi(nphi-1)),wy2ny(refden(nphi-1,km))
     $           ,string,-1.)
      enddo
      call pltend()
      call pltinit(0.,phi(nphi-1),0.,phi(nphi-1))
      call axis()
      call axlabels('!Af!@ [e/T!de!d]','!AJ!@ !p!o^!o!qn!dp!d d!Af!@')
      do km=0,num-1
         if(km.ne.0)call color(km)
         call polyline(phi,fiu(0,km),nphi)
         string='u!dm!d='
         call fwrite((km)*ustep,iwidth,2,string(8:))
         call jdrwstr(wx2nx(phi(nphi-1)),wy2ny(fiu(nphi-1,km))
     $           ,string,-1.)
      enddo
      call pltend()
      call pltinit(0.,phi(nphi-1),-.1,.1)
      call scalewn(1.e-2,phi(nphi-1),1.e-5,1.,.true.,.true.)
      call axis()
      call axis2()
      call axlabels('!Ay!@ [e/T!de!d]'
     $     ,'!p!o^!o!qV!df!d= 1- !AJ!@!p!o^!o!qn!df!d d!Af!@')
      call winset(.true.)
      do km=0,num-1
         if(km.ne.0)then 
            call color(km)
            call dashset(km)
         endif
         call polyline(phi,Vhatf(0,km),nphi)
         string='u!dm!d='
         call fwrite((km)*ustep,iwidth,2,string(8:))
         xg=.1
         yg=.9-.05*km
         if(Vhatf(nphi-1,km).gt.0.)call legendline(xg,yg,0,string)
      enddo
      call dashset(0)
      call pltend()
! f_b vs us plot.
      call accisinit()
      write(*,*)us(1),us(nphi-1),fb(1,0),fb(nphi-1,num-1)
      call scalewn(0.,us(nphi-1),1.,min(1.e6,abs(fb(nphi-1,0))),
     $     .false.,.true.)
      call axis()
      call axis2()
      call axlabels('u','f!db!d(u)/f!db!d(0)')
      do km=0,num-1
         if(km.ne.0)then 
            call color(km)
            call dashset(km)
         endif
         call polyline(us(1),fb(1,km),nphi-1)
         string='u!dm!d='
         call fwrite((km)*ustep,iwidth,2,string(8:))
         xg=.1
         yg=.9-.05*km
         call legendline(xg,yg,0,string)
!         call jdrwstr(wx2nx(us(nphi-1)),wy2ny(fb(nphi-1,km))
!     $           ,string,-1.)
         uma(km)=(km)*ustep
!         fb1(km)=fb(nphi-1,km)
      enddo
      call dashset(0)
      call pltend()
      usmax=1.
      call pltinit(0.,usmax,0.,3.)
      call axis()
      call axlabels('um','fb1')
      do iu=1,(nphi-1)/5,5
         do km=0,num-1
            usm=us(iu)
            fb1(km)=alog(fb(iu,km))/(usm)
!            write(*,'(a,2i3,2f8.3)')'iu, km, um, fb1='
!     $           ,iu,km,uma(km),fb1(km)
!      call autoplot(uma(0),fb1(0),num)
         enddo
         call polyline(uma(0),fb1(0),num)
      enddo
      do km=0,num-1
         uma(km)=km*usmax/num
! The following fit more or less captures the slope of the lines calculated
! in fb1 above.
         fb1(km)=2.54*(1-(uma(km)/0.93)**2)
      enddo
      call color(6)
      call polyline(uma(0),fb1(0),num)
      call pltend()
      end
!*********************************************************************
!*********************************************************************
      subroutine untrapcumtest
      integer nin,np,nm
      parameter (nin=400,nt=16,ntot=nin+nt)
      real fpa(nin),fma(nin),ua(nin)
      real cump(-ntot:ntot),cumv(-ntot:ntot)

      idebug=2
      phi=1.
      psi=phi
      um=0.3
      eta=2.
      call untrapcum(phi,um,psi,nin,np,nm,fpa,fma,ua,nt,cump,cumv,eta)

      if(idebug.gt.1)then
         call autoplot(ua,fma,nm)
         call axlabels('u','fma')
         call polyline(ua,fpa,np)
         call pltend()
      endif
      maxp=nm+nt
      maxm=np+nt
      call autoplot(cumv(-maxm),cump(-maxm),maxm+maxp+1)
      call axis2()
      call axlabels('u','cumulative distribution')
      call polymark(cumv(-nt-1),cump(-nt-1),2*nt+3,1)
      call pltend()
      end
!*********************************************************************
      subroutine untrapcum(phi,um,psi,nin,np,nm,fpa,fma,ua,nt,cump,cumv
     $     ,eta)
! Calculate the unnormalized cumulative probability distribution for a
! shifted Maxwellian at potential phi including a section for trapped
! particles in the middle. The trapped region contains 2*nt points in
! addition to that at index zero, which has velocity zero and is the
! zero of the bipolar distribution. If nt=0, then the trapped distribution
! is flat, giving a straight cumulative distribution section.
! eta is the power of u determining the shape of the trapped distribution.
! 2 gives parabolic. <=0 gives flattop.

      real phi,um,psi
      integer nin,np,nm
      real fpa(nin),fma(nin),ua(nin)
      real cump(-nin-nt:nin+nt),cumv(-nin-nt:nin+nt)

      idebug=0

      sqphi=sqrt(phi)
      den=untrappedden(phi,um,nin,fpa,fma,ua,np,nm)
! Separatrix value of distribution function.
      sepf=exp(-(um)**2)/sqrt(3.1415926)
      untrapden=(fpa(np)+fma(nm))
! Total electron density based upon the sech^4 potential variation
      dentot=1.+ phi-(5./4./sqrt(psi))*phi**1.5
! If the trapped distribution is sepf times 
!      (1+\gamma |u/sqrt{phi}|**\eta)/(1+\gamma)
! then the integrated trapped distribution is 
!      u (1+\gamma |u/sqrt{phi}|**\eta / (\eta+1))/(1+\gamma)
! and the integral 0 to sqrt{phi} is 
!      sqphi (1+\gamma/(\eta+1))/(1+\gamma)
! Thus 2*sepf*theabove = trapden and we must solve theabove=trapden/2*sepf
! to find gamma. Write r=ratio= trapden/(1*sqphi*sepf). Then 
! gamma = (1-r)/(r-1/(1+eta))
! Negative gamma is unacceptable, but occurs if r<1/(1+eta). Therefore
! eta is adjusted to prevent this.

! That enables us to relate \gamma (the depth parameter) to the required
! deficit, if \eta is specified. The special case \eta=2 has the
! property that for distributions that are functions of total energy,
! the distribution shape has the same parabolic form at all potentials
! (although different widths). We do not here suppose that f is a
! function of total energy.
      holeheight=-1
      if(eta.gt.0)then
! Deficit calculation
         trapden=dentot-untrapden
         ratio=min(1.,trapden/(2.*sepf*sqphi))
         if(ratio.lt.1./(1.+eta))then
            if(idebug.ne.0)write(*,'(a,4f8.4)')
     $           'Trapden impossible. phi,trapden,dentot,ratio='
     $           ,phi,trapden,dentot,ratio
            ratio=max(.01,ratio)
            eta=1/ratio-.9999
         endif
         gamma=(1-ratio)/(ratio-1./(1.+eta))
      else
! Flattop calculation
         trapden=sepf*sqphi*2.
         flattop=untrapden+trapden
! It gives total density greater than 1.
!         write(*,*)'Flattop: untrapden',untrapden,' trapden',trapden
!         write(*,*)'dentot',dentot,' flattop',flattop
         gamma=0.
      endif
      holeheight=1/(1.+gamma)

      if(idebug.gt.0)then
         write(*,*)'phi=',phi,' um=',um
         write(*,*)'np,nm=',np,nm,'fp,fm=',fpa(np),fma(nm)
         write(*,*)'Sepx speed',ua(1),'Sqrt phi',sqphi
         write(*,*)'sepf*2*sqphi',sepf*2*sqphi,' ratio',ratio
         write(*,*)'flattop-untrapden',flattop-untrapden
      endif

! Add enough to t2 to prevent vr overflow.
      t2=0.5*trapden+.0001
      dv=sqphi/(nt+1.)
! Trapped part (which is antisymmetric)
      do i=0,nt
         cumv(i)=i*dv
         cumv(-i)=-i*dv
         cump(i)=cumv(i)*(1+gamma*abs(cumv(i)/sqphi)**eta/(eta+1.))
     $        /(1.+gamma) *sepf
         cump(-i)=-cump(i)
      enddo
! Untrapped part. 
      do i=nt+1,nt+nin
         if(i.le.nm)then
            cump(-i)=-fma(i-nt)-t2
            cumv(-i)=-ua(i-nt)
         else
            cump(-i)=-fma(nm-nt)-t2
            cumv(-i)=-ua(nm-nt)
         endif
         if(i.le.np)then
            cump(i)=fpa(i-nt)+t2
            cumv(i)=ua(i-nt)
         else
            cump(i)=fpa(np-nt)+t2
            cumv(i)=ua(np-nt)
         endif
      enddo

      end
!*************************************************************************
!*************************************************************************
      subroutine getholepotl(psi,holelen,phi,x,id)
      integer id
      real psi,holelen,phi,x
      include 'ndimsdecl.f'
      include 'plascom.f'
      include 'meshcom.f'

      center=0.5*(xmeshstart(id)+xmeshend(id))
      arg=(x-center)/(debyelen*holelen)
      phi=psi/cosh(arg)**4
      if(.not.phi.ge.0.)write(*,*)'psi,arg,phi',psi,arg,phi

      end
!********************************************************************
      subroutine ranlen3position(x,id)
! Return a random fractional position in all 3 coordinate directions,
! accounting for the density scale length if present
! or for nonuniform background if present. id is the parallel dimension.
      include 'ndimsdecl.f'
      include 'meshcom.f'
      include 'creincom.f'
      include 'plascom.f'
      include 'partcom.f'
      real x(idtp)
      integer id
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
         enddo
         lfirst=.false.
      endif
      do j=1,ndims
         if(j.eq.id.or.x(iflag).ne.1)then
! Select new transverse position only the first time.
            g=gn(j)
            call ranlux(P,1)
            if(abs(g).ne.0)then
!         write(*,*)'Nonuniform plasma ranlenposition'
               sp=gp0(j)+alog(P*expsa(j)+(1.-P)*expsi(j))/g
            else
 1             sp=(1.-P)*xmeshstart(j)+P*xmeshend(j)
               if(bgmax(j).gt.0.)then
! Nonuniform initialization using a rejection scheme.
                  call ranlux(Q,1)
                  if(Q.gt.(1.+bgofx(sp,j))/(1.+bgmax(j)))then
!               write(*,*)'Nonuniform pinit',Q,bgmax(j)
                     call ranlux(P,1)
                     goto 1
                  endif
               endif
            endif
            x(j)=sp*0.999999+.0000005
         endif
      enddo
      x(iflag)=1
      end
!********************************************************************
