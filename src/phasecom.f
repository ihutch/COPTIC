c Phase space accumulation variables
      integer npsx,npsv,npsbuf
      integer ipsftri
      parameter (npsx=200,npsv=50,npsbuf=10000)
      real psfxv(npsx,npsv,nspeciesmax),psx(npsx),psn(npsx,nspeciesmax)
     $     ,psv(npsv,nspeciesmax)
      real psvmax,psvmin,psxmax,psxmin,psfmax
      common /phasespace/psfxv,psvmax,psvmin,psxmax,psxmin,psx,psv
     $     ,psn,psfmax,ipsftri
