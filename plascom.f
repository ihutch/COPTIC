c Common data describing the plasma.
c Requires ndimsdecl.f to be preloaded to give parameters
c    ndims, nspeciesmax
c
c Contains the following which refer to the first species only:
c lambda_{De}, Ti/Te, Drift velocity, region spherical boundary radius.
c phiprobe, ratio of ion mass(in proton masses) to Z 
c Btotal,Value taken as \infty,
      real debyelen,Ti,vd,rs,pi,phip,eoverm,Bt,Btinf
      parameter (pi=3.1415927,Btinf=1.e3)
c Bfield direction,vparallel,vperp, vdrift=direction cosines of v.
      real Bfield(ndims),vpar,vperp(ndims),vdrift(ndims)
c Reference point of Te, Te gradient components, gtt total gradient.
      real gp0(ndims),gt(ndims),gtt
c Reference point of ne gradient is same as Te gradient.
c ne gradient components, total gradient:
      real gn(ndims),gnt
c
c Those used to be directly stored but now some are just equivalences to
c the first species of the following array forms.
      real Ts(nspeciesmax),vds(nspeciesmax),eoverms(nspeciesmax)
      real vpars(nspeciesmax), Tperps(nspeciesmax)
      real vperps(ndims,nspeciesmax),vdrifts(ndims,nspeciesmax)
      equivalence (Ti,Ts),(vd,vds),(eoverm,eoverms)
      equivalence (vpar,vpars),(vperp,vperps),(vdrift,vdrifts)

      common/plascom/debyelen,Ts,Tperps,vds,rs,phip,eoverms,Bt,Bfield
     $     ,vpars,vperps,vdrifts,gp0,gt,gtt,gn,gnt
