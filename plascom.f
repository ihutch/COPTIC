c Common data describing the plasma.
c lambda_{De}, Ti/Te, Drift velocity, region spherical boundary radius.
c phiprobe, ratio of ion mass(in proton masses) to Z 
c Btotal,Value taken as \infty,
      real debyelen,Ti,vd,rs,pi,phip,rmtoz,Bt,Btinf
      integer nplasdims
      parameter (pi=3.1415927,nplasdims=3,Btinf=1.e3)
c Bfield directoin,vparallel,vperp,vdrift direction cosines of v.
      real Bfield(nplasdims),vpar,vperp(nplasdims),vdrift(nplasdims)
c Reference point of Te, Te gradient components, gtt total gradient.
      real gp0(nplasdims),gt(nplasdims),gtt
c Reference point of ne gradient is same as Te gradient.
c ne gradient components, total gradient:
      real gn(nplasdims),gnt

      common/plascom/debyelen,Ti,vd,rs,phip,rmtoz,Bt,Bfield,vpar,vperp
     $     ,vdrift,gp0,gt,gtt,gn,gnt
