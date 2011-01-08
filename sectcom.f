c Data structures for collecting intersections with objects

      integer sc_ndims,sc_npts,sc_ipt
      parameter (sc_ndims=3,sc_npts=10000)
      real x_sc(sc_ndims,2,sc_npts)
      integer iob_sc(sc_npts)

      common /sccom/x_sc,iob_sc,sc_ipt
