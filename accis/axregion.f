c*********************************************************************
      subroutine axregion(nxmin,nxmax,nymin,nymax)
c Set the normalized position of the plotting axis window.
      real nxmin,nxmax,nymin,nymax
      include 'plotcom.h'
      naxmin=nxmin
      naxmax=nxmax
      naymin=nymin
      naymax=nymax
      naxpt=nxmin
      naypt=nymin
      return
      end

