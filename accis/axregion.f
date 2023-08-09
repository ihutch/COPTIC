c*********************************************************************
      subroutine axregion(nxmin,nxmax,nymin,nymax)
c Set the normalized position of the plotting axis window.
      real nxmin,nxmax,nymin,nymax
      include 'plotcom.h'
      if(nxmin.eq.0.and.nxmax.eq.0)then
         naxmin=0.31
         naxmax=0.91
         naymin=0.1
         naymax=0.7
      else
         naxmin=nxmin
         naxmax=nxmax
         naymin=nymin
         naymax=nymax
      endif
      naxpt=nxmin
      naypt=nymin
      end

