c Common for storing collisional distribution
c Requires ndimsdecl.f to give ndims
      integer ncdistmax
      parameter (ncdistmax=1000000)
      integer ncdist
      real v_col(ndimsmax,ncdistmax)
      real fxvcol(ncdistmax+1,ndims)
      real cdistflux(ndimsmax),cdistcum(ndimsmax+1)
      common /cdistcom/ncdist,v_col,cdistflux,cdistcum,fxvcol
