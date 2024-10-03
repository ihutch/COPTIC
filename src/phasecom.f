c Phase space accumulation variables
      integer npsx,npsv,npsbuf
      integer ipsftri,ilogspec,ipsversion
      parameter (npsx=200,npsv=50,npsbuf=10000)
      logical lsideplot
      real psfxv(npsx,npsv,nspeciesmax),psx(npsx),psn(npsx,nspeciesmax)
     $     ,psv(npsv,nspeciesmax),finfofv(npsv,nspeciesmax)
     $     ,psvave(npsx,nspeciesmax)
      real psvmax(nspeciesmax),psvmin(nspeciesmax),psxmax,psxmin
      real finfmax(nspeciesmax),psfmax(nspeciesmax)
      common /phasespace/psfxv,psvmax,psvmin,psxmax,psxmin,psx,psv,psn
     $     ,psvave,psfmax,ipsftri,finfofv,finfmax,ipsversion,lsideplot
     $     ,ilogspec
