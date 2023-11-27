c Phase space accumulation variables
      integer npsx,npsv,npsbuf
      integer ipsftri
      parameter (npsx=200,npsv=50,npsbuf=10000)
      real psfxv(npsx,npsv,nspeciesmax),psx(npsx),psn(npsx,nspeciesmax)
     $     ,psv(npsv,nspeciesmax),finfofv(npsv,nspeciesmax)
      real psvmax(nspeciesmax),psvmin(nspeciesmax),psxmax,psxmin,psfmax
      common /phasespace/psfxv,psvmax,psvmin,psxmax,psxmin,psx,psv
     $     ,psn,psfmax,ipsftri,finfofv
