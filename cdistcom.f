c Common for storing collisional distribution
c Requires ndimsdecl.f to give ndims
      integer ncdistmax
      parameter (ncdistmax=1000000)
      integer ncdist
      real v_col(ndimsmax,ncdistmax)
      real fxvcol(ncdistmax+1,ndims)
      real cdistflux(ndimsmax),cdistcum(ndimsmax+1)
c Multispecies
      integer ncdists(nspeciesmax)
      real vcols(ndimsmax,ncdistmax,nspeciesmax)
      real cdistfluxs(ndimsmax,nspeciesmax)
      equivalence (v_col,vcols),(cdistflux,cdistfluxs)
      equivalence (ncdist,ncdists)
      common /cdistcom/ncdists,vcols,cdistfluxs,cdistcum,fxvcol
