c Parameters and random integrals for cartesian reinjection
c from a separable distribution function.
c The six values correspond to +-x,+-y,+-z.
c Have to use a pretty large array for the inverse cumulative probs.
      integer ncrein
      parameter (ncrein=1000)
c The inverse cumulative flux probability distribution 
c for each direction (face), such that vchoice(id)=hrein(random,id).
      real hrein(0:ncrein,6)
c The total flux for each direction unnormalized:
      real grein(6)
c The inverse cumulative velocity probability distribution
c (not weighted by velocity) for the three (not 6) dimensions,
c used for orthogonal velocity choices. vchoice(id)=prein(random,id).
      real prein(0:ncrein,3)
c The cumulative distribution of flux among the six faces. 
      real gintrein(0:6)
      integer idrein

      common /crein/hrein,grein,prein,gintrein,idrein
