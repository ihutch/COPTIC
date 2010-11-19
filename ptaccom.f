c The following makes this dependent on having meshcom loaded already.
      integer nptdiag,mdims
      parameter (nptdiag=400,mdims=ndims_mesh)
      real fv(nptdiag,mdims),cumfv(0:nptdiag,mdims)
      real px(nptdiag,mdims)
      real vdiag(nptdiag,mdims)
      real xdiag(nptdiag,mdims)
      integer nfvaccum
      common /cartdiag/fv,px,vdiag,xdiag,cumfv,nfvaccum

c------------Former fvgriddecl.f------------
c Declarations of the spatial grid on which particle distribution
c diagnostics are to be accumulated.

c In this section there is an assumption that we are in 3 dimensions.
      integer nsub_i,nsub_j,nsub_k
      parameter (nsub_i=10,nsub_j=10,nsub_k=10)
      integer Lsi,Lsj,nsub_tot
      parameter (Lsi=nsub_i,Lsj=Lsi*nsub_j,nsub_tot=Lsj*nsub_k)
      integer isfull(mdims)

      integer nsbins
      parameter (nsbins=32)
      integer ibinmap(nptdiag,mdims)
      real vsbin(nsbins,mdims),csbin(nsbins,mdims)
      real fsv(nsbins,mdims)
      real vhbin(0:nsbins,mdims)
      real fvx(nsbins,mdims,nsub_tot)
      real denfvx(nsub_tot)
c ibinmap is the map from uniform to combined bins
c vsbin is the center velocity of the combined bins
c csbin is the number of fine bins in each combined bin
c fsv is the sum of fv in each combined bin during the initial
c     accumulation and bin calculation.
c vhbin is the histogram boundaries of the combined bins
c All of these must be common to all processes.
      common /subdiag/ibinmap,isfull,vsbin,csbin,vhbin,fsv,fvx,denfvx
