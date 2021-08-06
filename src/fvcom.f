c Particle Distribution function parameters for non-Maxwellian
c Requires ndimsdecl.f to define nspeciesmax
      integer ncmax,nofv
      parameter (ncmax=4,nofv=400)
c Number of Gaussians:
      integer nc(nspeciesmax)
c Gaussian v-shift, vthermal[sqrt(T/m)], density, store.
      real vsc(ncmax,nspeciesmax),vtc(ncmax,nspeciesmax)
      real dcc(ncmax,nspeciesmax),vscs(ncmax,nspeciesmax)
      real vofv(nofv),fofv(nofv),Pfofv(nofv)
      logical lfv
      common/fvcom/lfv,nc,vsc,vtc,dcc,vscs,vofv,fofv,Pfofv

