c Common for storing collisional distribution
      integer ncdistmax,nc_ndims
      parameter (ncdistmax=1000000,nc_ndims=3)
      integer ncdist
      real v_col(nc_ndims,ncdistmax)
      real fxvcol(ncdistmax+1,nc_ndims)
      real cdistflux(nc_ndims),cdistcum(nc_ndims+1)
      common /cdistcom/ncdist,v_col,cdistflux,cdistcum,fxvcol
