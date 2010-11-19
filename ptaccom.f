c The following makes this dependent on having meshcom loaded already.
      integer nptdiag,mdims
      parameter (nptdiag=400,mdims=ndims_mesh)
      real fv(nptdiag,mdims),cumfv(0:nptdiag,mdims)
      real px(nptdiag,mdims)
      real vdiag(nptdiag,mdims)
      real xdiag(nptdiag,mdims)
      common /cartdiag/fv,px,vdiag,xdiag,cumfv
