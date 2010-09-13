      parameter (ndiag=100,mdims=3)
      real xr(3*mdims)
      real fv(ndiag,mdims)
      real px(ndiag,mdims)
      real diagv(ndiag,mdims),cumdiagv(ndiag,mdims)
      real diagx(ndiag,mdims),cumdiagx(ndiag,mdims)
      common /cartdiag/fv,px,diagv,diagx,cumdiagv,cumdiagx
