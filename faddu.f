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
c If the overflow trap above is working, then this ought not to be needed.
      if(.not.faddu.lt.1.e20)then
         write(*,*)'ERROR: faddu singularity'
     $        ,u,um,uci(index),fprime,faddu,index
         stop
      endif
      end
c************************************************************************
c**********************************************************************
c This version returns the Boltzmann derivative but zero actual function
c value. That is supposed to make the iteration linearized using the 
c Boltzmann response but allow explicit electrons.
      real function fnodensity(u,fprime,index)
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
c If the overflow trap above is working, then this ought not to be needed.
      if(.not.faddu.lt.1.e20)then
         write(*,*)'ERROR: faddu singularity'
     $        ,u,um,uci(index),fprime,faddu,index
         stop
      endif
      fnodensity=0.
      end
c************************************************************************

