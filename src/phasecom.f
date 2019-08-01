c Phase space accumulation variables
      integer npsx,npsv,npsbuf
      integer ipsftri
      parameter (npsx=200,npsv=50,npsbuf=10000)
      real psfxv(npsx,npsv),psx(npsx),psn(npsx),psv(npsv)
      real psvmax,psvmin,psxmax,psxmin,psfmax
      common /phasespace/psfxv,psvmax,psvmin,psxmax,psxmin,psx,psv
     $     ,psn,psfmax,ipsftri
