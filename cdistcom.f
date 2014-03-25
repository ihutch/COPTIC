c Common for storing collisional distribution
c Requires ndimsdecl.f to give ndims
      integer ncdistmax
      parameter (ncdistmax=1000000)
      integer ncdist
      real v_col(ndims,ncdistmax)
      real fxvcol(ncdistmax+1,ndims)
      real cdistflux(ndims),cdistcum(ndims+1)
      common /cdistcom/ncdist,v_col,cdistflux,cdistcum,fxvcol
