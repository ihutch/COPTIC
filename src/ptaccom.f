! Dependent on having ndimsdecl.f loaded already, defining ndimsmax.
      integer nptdiagmax,nptdiag
      parameter (nptdiagmax=400)
      real fv(nptdiagmax,ndimsmax),cumfv(0:nptdiagmax,ndimsmax)
      real px(nptdiagmax,ndimsmax)
      real vdiag(nptdiagmax,ndimsmax)
      real xdiag(nptdiagmax,ndimsmax)
      integer nfvaccum,ivproj
      common /cartdiag/fv,px,vdiag,xdiag,cumfv,nfvaccum,ivproj,nptdiag
! In this section there is an assumption that we are in 3 dimensions.
      integer nsub_i,nsub_j,nsub_k
! These determine the number of sub-divisions of the three directions
! into which the mesh is divided for the purpose of particle diagnostics
! However, they functionally only set defaults and the total nsub_tot.
      parameter (nsub_i=15,nsub_j=15,nsub_k=15)
      integer nsub_tot
      parameter (nsub_tot=nsub_i*nsub_j*nsub_k)
      integer isfull(ndimsmax),isuds(ndimsmax)

! Data for summed nonuniform binning to reduce the size of the output
! If nsbins.eq.nptdiag, then uniform binning will be used for the summed
! bins. Behavior for nsbins>nptdiag is uncertain. Avoid that.
      integer nsbins
      parameter (nsbins=64)
      integer ibinmap(nptdiagmax,ndimsmax)
      real vsbin(nsbins,ndimsmax),csbin(nsbins,ndimsmax)
      real fsv(nsbins,ndimsmax)
      real vhbin(0:nsbins,ndimsmax)
      real fvx(nsbins,ndimsmax,nsub_tot)
      real f2vx(nsbins,nsbins,ndimsmax,nsub_tot)
      real denfvx(nsub_tot)
      real vtkudata(nsub_i+1,nsub_j+1,nsub_k+1,0:nsbins,2*ndimsmax)
! ibinmap is the map from uniform to combined bins
! vsbin is the center velocity of the combined bins
! csbin is the number of fine bins in each combined bin
! fsv is the sum of fv in each combined bin during the initial
!     accumulation and bin calculation.
! vhbin is the histogram boundaries of the combined bins
! fvx is the 1-D distribution in 3 directions at each sub-position.
! f2vx is the 2-D distribution normal to 3 directions ditto.
! All of these must be common to all processes.

      common /subdiag/ibinmap,isfull,isuds,vsbin,csbin,vhbin,fsv,fvx
     $     ,denfvx,f2vx,vtkudata
