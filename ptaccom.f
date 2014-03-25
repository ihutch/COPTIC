c Dependent on having meshcom loaded already, defining ndims.
      integer nptdiag,mdims
      parameter (nptdiag=400,mdims=ndims)
      real fv(nptdiag,mdims),cumfv(0:nptdiag,mdims)
      real px(nptdiag,mdims)
      real vdiag(nptdiag,mdims)
      real xdiag(nptdiag,mdims)
      integer nfvaccum,ivproj
      common /cartdiag/fv,px,vdiag,xdiag,cumfv,nfvaccum,ivproj
c In this section there is an assumption that we are in 3 dimensions.
      integer nsub_i,nsub_j,nsub_k
c These determine the number of sub-divisions of the three directions
c into which the mesh is divided for the purpose of particle diagnostics
c However, they functionally only set defaults and the total nsub_tot.
      parameter (nsub_i=15,nsub_j=15,nsub_k=15)
      integer nsub_tot
      parameter (nsub_tot=nsub_i*nsub_j*nsub_k)
      integer isfull(mdims),isuds(mdims)

c Data for summed nonuniform binning to reduce the size of the output
c If nsbins.eq.nptdiag, then uniform binning will be used for the summed
c bins. Behavior for nsbins>nptdiag is uncertain. Avoid that.
      integer nsbins
      parameter (nsbins=32)
      integer ibinmap(nptdiag,mdims)
      real vsbin(nsbins,mdims),csbin(nsbins,mdims)
      real fsv(nsbins,mdims)
      real vhbin(0:nsbins,mdims)
      real fvx(nsbins,mdims,nsub_tot)
      real denfvx(nsub_tot)
      real vtkudata(nsub_i+1,nsub_j+1,nsub_k+1,0:nsbins,2*mdims)
c ibinmap is the map from uniform to combined bins
c vsbin is the center velocity of the combined bins
c csbin is the number of fine bins in each combined bin
c fsv is the sum of fv in each combined bin during the initial
c     accumulation and bin calculation.
c vhbin is the histogram boundaries of the combined bins
c All of these must be common to all processes.
      common /subdiag/ibinmap,isfull,isuds,vsbin,csbin,vhbin,fsv,fvx
     $     ,denfvx,vtkudata
