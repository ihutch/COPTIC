!**********************************************************************
!     L(u) + f(u) = q(x,y,...), 
!     where L is a second order elliptical differential operator 
!     represented by a difference stencil of specified coefficients,
!     f is some additional function, and q is the "charge density".
! f is exp(u) here for Boltzmann electrons and densities normalized
! to unity at infinity.
      real function faddu(u,fprime,index)
      real u,fprime
      integer index
      include 'ndimsdecl.f'
! In order to access point-charge information we need:
      include '3dcom.f'
! Needed for ptchcom.f:
      include 'griddecl.f'
      include 'ptchcom.f'
      real ubig,um
      parameter (ubig=40.)
!
      if(iptch_mask.eq.0 .and. gtt_copy.eq.0. .and. gnt_copy.eq.0)then
         fprime=exp(u)
         faddu=fprime
      else
! Need to compensate for point charges and/or Te-gradient.
         um=(u+uci(index))/Teci(index)
         if(abs(um).gt.ubig)um=sign(ubig,um)
         fprime=exp(um)
         faddu=fprime-rhoci(index)
      endif
      end
!************************************************************************

!**********************************************************************
!     L(u) + f(u) = q(x,y,...), 
!     where L is a second order elliptical differential operator 
!     represented by a difference stencil of specified coefficients,
!     f is some additional function, and q is the "charge density".
! f is exp(u) here for Boltzmann electrons and densities normalized
! to unity at infinity.
! This version applies the additional factor weighted by the value
! of boltzwt on the indexed mesh. boltzwt is set to zero outside the
! particle region where there ought to be no electrons.
      real function fadcomp(u,fprime,index)
      implicit none
      real u,fprime
      integer index
      include 'ndimsdecl.f'
! In order to access point-charge information we need:
      include '3dcom.f'
! Needed for ptchcom.f:
      include 'griddecl.f'
      include 'ptchcom.f'
      real ubig,um
      parameter (ubig=40.)
!
      if(gtt_copy.ne.0. .or. gnt_copy.ne.0)then
! Need to compensate for point charges and/or Te-gradient.
         um=(u+uci(index))/Teci(index)
         if(abs(um).gt.ubig)um=sign(ubig,um)*boltzsign
      else
         um=u*boltzsign
      endif
      fprime=boltzsign*boltzamp*boltzwti(index)*exp(um)
      if(iptch_mask.ne.0)then 
! We access rhoci only if we know we need to.
         fadcomp=fprime-rhoci(index)
      else
         fadcomp=fprime
      endif
      end
!************************************************************************
!**********************************************************************
!     L(u) + f(u) = q(x,y,...), 
!     where L is a second order elliptical differential operator 
!     represented by a difference stencil of specified coefficients,
!     f is some additional function, and q is the "charge density".
! Linearized f with some small amplitude helps smooth long wavelength
! modes when electrons are used.
! to unity at infinity.
! Obsolete/Not used.
      real function faddlin(u,fprime,index)
      real u,fprime
      integer index
      include 'ndimsdecl.f'
! In order to access point-charge information we need:
      include '3dcom.f'
! Needed for ptchcom.f:
      include 'griddecl.f'
      include 'ptchcom.f'
      real fbig,um
      parameter (fbig=1.)
!
      fprime=boltzamp
      if(iptch_mask.eq.0 .and. gtt_copy.eq.0. .and. gnt_copy.eq.0)then
         faddlin=u*boltzamp
      else
! Need to compensate for point charges and/or Te-gradient.
         um=(u+uci(index))/Teci(index)
         faddlin=um*boltzamp-rhoci(index)
      endif
      if(abs(faddlin).gt.fbig)faddlin=sign(fbig,faddlin)
      end
!************************************************************************
