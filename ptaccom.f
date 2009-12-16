      parameter (ndiag=200,mdims=3)
      real xr(3*mdims)
      real fv(ndiag,mdims)
      real px(ndiag,mdims)
      real diagv(ndiag)
      real diagx(ndiag,mdims)
      common /cartdiag/fv,px,diagv,diagx
