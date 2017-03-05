! Common data describing the plasma.
! Requires ndimsdecl.f to be preloaded to give parameters
!    ndims, nspeciesmax
!
! Contains the following which refer to the first species only:
! lambda_{De}, Ti/Te, Drift velocity, region spherical boundary radius.
! phiprobe, ratio of ion mass(in proton masses) to Z 
! Btotal,Value taken as \infty,
      real debyelen,Ti,vd,rs,pi,phip,eoverm,Bt,Btinf
      parameter (pi=3.1415927,Btinf=1.e3)
! Bfield direction,vparallel,vperp, vdrift=direction cosines of v.
      real Bfield(ndims),vpar,vperp(ndims),vdrift(ndims)
! Reference point of Te, Te gradient components, gtt total gradient.
      real gp0(ndims),gt(ndims),gtt
! Reference point of ne gradient is same as Te gradient.
! ne gradient components, total gradient:
      real gn(ndims),gnt
!
! Those used to be directly stored but now some are just equivalences to
! the first species of the following array forms.
      real Ts(nspeciesmax),vds(nspeciesmax),eoverms(nspeciesmax)
      real vpars(nspeciesmax), Tperps(nspeciesmax)
      real vperps(ndims,nspeciesmax),vdrifts(ndims,nspeciesmax)
      real driftfields(ndims,nspeciesmax)
      equivalence (Ti,Ts),(vd,vds),(eoverm,eoverms)
      equivalence (vpar,vpars),(vperp,vperps),(vdrift,vdrifts)

      common/plascom/debyelen,Ts,Tperps,vds,rs,phip,eoverms,Bt,Bfield
     $     ,vpars,vperps,vdrifts,driftfields,gp0,gt,gtt,gn,gnt
