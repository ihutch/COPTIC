c***************************************************************************
	subroutine winset(wsw)
	logical wsw
      include 'plotcom.h'
      if(wsw)then
      call truncf(naxmin,naxmax,naymin,naymax)
      else
      call truncf(0.,0.,0.,0.)
      endif
      return
      end
