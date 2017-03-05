! Common for storing collisional distribution
! Requires ndimsdecl.f to give ndims
      integer ncdistmax
      parameter (ncdistmax=1000000)
      integer ncdist
      real vcol(ndimsmax,ncdistmax)
      real fxvcol(ncdistmax+1,ndims)
      real*8 cdistflux(ndimsmax)
      real cdistcum(ndimsmax+1)
! Multispecies
      integer ncdists(nspeciesmax)
      real vcols(ndimsmax,ncdistmax,nspeciesmax)
      real fxvcols(ncdistmax+1,ndims,nspeciesmax)
      real*8 cdistfluxs(ndimsmax,nspeciesmax)
      real cdistcums(ndimsmax+1,nspeciesmax)
      equivalence (vcol,vcols)
      equivalence (cdistflux,cdistfluxs)
      equivalence (ncdist,ncdists)
      equivalence (fxvcol,fxvcols)
      equivalence (cdistcum,cdistcums)
      common /cdistcom/ncdists,vcols,cdistfluxs,cdistcums,fxvcols
