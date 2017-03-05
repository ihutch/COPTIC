! Needs ndimsdecl.f for parameter nspeciesmax.
! Parameters and random integrals for cartesian reinjection
! from a separable distribution function.
! The six values correspond to +-x,+-y,+-z.
! Have to use a pretty large array for the inverse cumulative probs.
      integer ncrein
      parameter (ncrein=1000)
! The inverse cumulative flux probability distribution 
! for each direction (face), such that vchoice(id)=hrein(random,id).
      real hreins(0:ncrein,6,nspeciesmax)
! The total flux for each direction unnormalized:
      real greins(6,nspeciesmax)
! The inverse cumulative velocity probability distribution
! (not weighted by velocity) for the three (not 6) dimensions,
! used for orthogonal velocity choices. vchoice(id)=prein(random,id).
      real preins(0:ncrein,3,nspeciesmax)
! The cumulative distribution of flux among the six faces. 
      real gintreins(0:6,nspeciesmax)
      integer idreins(nspeciesmax)
! Whether reinjection is initialized or not
      logical lreininit
      real hrein(0:ncrein,6),grein(6),prein(0:ncrein,3)
      real gintrein(0:6)
      integer idrein
      equivalence (hrein,hreins),(grein,greins),(prein,preins)
      equivalence (gintrein,gintreins),(idrein,idreins)

      common /crein/hreins,greins,preins,gintreins,idreins,lreininit
