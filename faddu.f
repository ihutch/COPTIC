c**********************************************************************
c     L(u) + f(u) = q(x,y,...), 
c     where L is a second order elliptical differential operator 
c     represented by a difference stencil of specified coefficients,
c     f is some additional function, and q is the "charge density".
c f is exp(u) here for Boltzmann electrons and densities normalized
c to unity at infinity.
      real function faddu(u,fprime,index)
      real u,fprime
      integer index
      include 'ndimsdecl.f'
c In order to access point-charge information we need:
      include '3dcom.f'
c Needed for ptchcom.f:
      include 'griddecl.f'
      include 'ptchcom.f'
      real ubig,um
      parameter (ubig=40.)
c
      if(iptch_mask.eq.0 .and. gtt_copy.eq.0. .and. gnt_copy.eq.0)then
         fprime=exp(u)
         faddu=fprime
      else
c Need to compensate for point charges and/or Te-gradient.
         um=(u+uci(index))/Teci(index)
         if(abs(um).gt.ubig)um=sign(ubig,um)
         fprime=exp(um)
         faddu=fprime-rhoci(index)
      endif
      end
c************************************************************************

c**********************************************************************
c     L(u) + f(u) = q(x,y,...), 
c     where L is a second order elliptical differential operator 
c     represented by a difference stencil of specified coefficients,
c     f is some additional function, and q is the "charge density".
c f is exp(u) here for Boltzmann electrons and densities normalized
c to unity at infinity.
c This version applies the additional factor weighted by the value
c of boltzwt on the indexed mesh. boltzwt is set to zero outside the
c particle region where there ought to be no electrons.
      real function fadcomp(u,fprime,index)
      implicit none
      real u,fprime
      integer index
      include 'ndimsdecl.f'
c In order to access point-charge information we need:
      include '3dcom.f'
c Needed for ptchcom.f:
      include 'griddecl.f'
      include 'ptchcom.f'
      real ubig,um
      parameter (ubig=40.)
c
      if(gtt_copy.ne.0. .or. gnt_copy.ne.0)then
c Need to compensate for point charges and/or Te-gradient.
         um=(u+uci(index))/Teci(index)
         if(abs(um).gt.ubig)um=sign(ubig,um)*boltzsign
      else
         um=u*boltzsign
      endif
      fprime=boltzsign*boltzamp*boltzwti(index)*exp(um)
      if(iptch_mask.ne.0)then 
c We access rhoci only if we know we need to.
         fadcomp=fprime-rhoci(index)
      else
         fadcomp=fprime
      endif
      end
c************************************************************************
c**********************************************************************
c     L(u) + f(u) = q(x,y,...), 
c     where L is a second order elliptical differential operator 
c     represented by a difference stencil of specified coefficients,
c     f is some additional function, and q is the "charge density".
c Linearized f with some small amplitude helps smooth long wavelength
c modes when electrons are used.
c to unity at infinity.
c Obsolete/Not used.
      real function faddlin(u,fprime,index)
      real u,fprime
      integer index
      include 'ndimsdecl.f'
c In order to access point-charge information we need:
      include '3dcom.f'
c Needed for ptchcom.f:
      include 'griddecl.f'
      include 'ptchcom.f'
      real fbig,um
      parameter (fbig=1.)
c
      fprime=boltzamp
      if(iptch_mask.eq.0 .and. gtt_copy.eq.0. .and. gnt_copy.eq.0)then
         faddlin=u*boltzamp
      else
c Need to compensate for point charges and/or Te-gradient.
         um=(u+uci(index))/Teci(index)
         faddlin=um*boltzamp-rhoci(index)
      endif
      if(abs(faddlin).gt.fbig)faddlin=sign(fbig,faddlin)
      end
c************************************************************************
