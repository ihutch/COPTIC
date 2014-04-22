c Common data describing the plasma.
c Requires ndimsdecl.f to be preloaded to give parameters
c    ndims, nspeciesmax
c
c Contains the following which refer to the first species only:
c lambda_{De}, Ti/Te, Drift velocity, region spherical boundary radius.
c phiprobe, ratio of ion mass(in proton masses) to Z 
c Btotal,Value taken as \infty,
      real debyelen,Ti,vd,rs,pi,phip,rmtoz,Bt,Btinf
      parameter (pi=3.1415927,Btinf=1.e3)
c Bfield direction,vparallel,vperp, vdrift=direction cosines of v.
      real Bfield(ndims),vpar,vperp(ndims),vdrift(ndims)
c Reference point of Te, Te gradient components, gtt total gradient.
      real gp0(ndims),gt(ndims),gtt
c Reference point of ne gradient is same as Te gradient.
c ne gradient components, total gradient:
      real gn(ndims),gnt
c
c Those used to be directly stored:
c    common/plascom/debyelen,Ti,vd,rs,phip,rmtoz,Bt,Bfield,vpar,vperp
c     $     ,vdrift,gp0,gt,gtt,gn,gnt
c but now some are just equivalences to the first species of the 
c following array forms.
      real Ts(nspeciesmax),vds(nspeciesmax),rmtozs(nspeciesmax)
      real vpars(nspeciesmax)
      real vperps(ndims,nspeciesmax),vdrifts(ndims,nspeciesmax)
      equivalence (Ti,Ts),(vd,vds),(rmtoz,rmtozs)
      equivalence (vpar,vpars),(vperp,vperps),(vdrift,vdrifts)

      common/plascom/debyelen,Ts,vds,rs,phip,rmtozs,Bt,Bfield,vpars
     $     ,vperps,vdrifts,gp0,gt,gtt,gn,gnt
