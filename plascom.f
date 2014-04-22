c Common data describing the plasma.
c Requires ndimsdecl.f to be preloaded to give parameter ndims
c lambda_{De}, Ti/Te, Drift velocity, region spherical boundary radius.
c phiprobe, ratio of ion mass(in proton masses) to Z 
c Btotal,Value taken as \infty,
      real debyelen,Ti,vd,rs,pi,phip,rmtoz,Bt,Btinf
      parameter (pi=3.1415927,Btinf=1.e3)
c Bfield directoin,vparallel,vperp,vdrift direction cosines of v.
      real Bfield(ndims),vpar,vperp(ndims),vdrift(ndims)
c Reference point of Te, Te gradient components, gtt total gradient.
      real gp0(ndims),gt(ndims),gtt
c Reference point of ne gradient is same as Te gradient.
c ne gradient components, total gradient:
      real gn(ndims),gnt

      common/plascom/debyelen,Ti,vd,rs,phip,rmtoz,Bt,Bfield,vpar,vperp
     $     ,vdrift,gp0,gt,gtt,gn,gnt
