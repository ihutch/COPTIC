c Common data describing the plasma.
c lambda_{De}, Ti/Te, Drift velocity, region spherical boundary radius.
c phiprobe, ratio of ion mass(in proton masses) to Z 
      real debyelen,Ti,vd,rs,pi,phip,rmtoz,Bt
      integer nplasdims
      parameter (pi=3.1415927,nplasdims=3)
      real Bfield(nplasdims)

      common/plascom/debyelen,Ti,vd,rs,phip,rmtoz,Bt,Bfield
