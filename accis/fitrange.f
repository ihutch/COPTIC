cd 
c************************************************************************
c Obsolete version.
      subroutine fitrangeold(xmin,xmax,ntics,nxfac,xfac,xtic,xt1st
     $     ,xtlast)
c
c  To fit a suitable axis range for reasonable scales.
c  Inputs:
c       xmin, xmax: range to be fitted.
c       ntics: (maximum) number of divisions (tics) to fit to it.
c  Outputs:
c       nxfac: power of ten by which the range is scaled.
c       xfac: 10**nxfac. World-value=xfac*axis-label.
c       xtic: Tic-spacing in world units.
c       xt1st: The integer multiple of xtic closest to xmin 
c             lying outside the range (xmin,xmax).
c       xtlast: The integer multiple of xtic closest to xmax
c             lying outside the range (xmin,xmax).
      real span,nsfac,sfac

      span=(xmax-xmin)
      if(xmax.eq.0. .and. xmin.eq.0)then
         write(*,*)'Fitrange warning. xmin=xmax=0'
         nxfac=0
         span=1.
      else
         nxfac=nint(log10(max(abs(xmin),abs(xmax)))-0.4999999)
      endif
      xfac=10.**nxfac
      if(ntics.le.0)then
         write(*,'('' ntics<=0'')')
         return
      endif
      xtic=span/ntics
      nsfac=nint(log10(0.099999*abs(xtic))+0.500001)
      sfac=10.**nsfac
      xtic=abs(xtic)/sfac
      if(xtic.lt.1.)then
         write(*,'('' Fitrange error 1. xtic='',f16.7)')xtic
      elseif(xtic.le.2.)then
         xtic=2.
c A prior version used just .le.3. here which favors xtic=5, but leads
c to a ratcheting up with successive calls, which is unsatisfactory.
c Therefore suppose that if we are exactly 4 it is because we did an
c earlier fitrange.
      elseif(xtic.le.3 .or. (xtic-4.).lt.0.0001)then
         xtic=4.
      elseif(xtic.le.5.)then
         xtic=5.
      elseif(xtic.le.10.0001)then
         xtic=10.
      else
         write(*,'('' Fitrange error NAN. Range:'',2g11.4)')xmin,xmax
         xtic=1.
         nxfac=0
         sfac=1.
         xfac=1.
      endif
      xtic=xtic*sfac
      xtic=sign(xtic,span)
c Use span to cope with error cases:
      xtlast=xtic*anint((xmin+span)/xtic+0.49999)
c Change xt1st the last thing, since it might be xmin, itself.
      xt1st=xtic*anint(xmin/xtic-0.49999)
c      write(*,*)'xtic,sfac,xt1st,xtlast',xtic,sfac,xt1st,xtlast,ntics
      return
      end
c********************************************************************

c************************************************************************
c New version examines all nice increments for the shortest range.
c
      subroutine fitrange(xmin,xmax,ntics,nxfac,xfac,xtic,xt1st,xtlast)
c
c  To fit a suitable axis range for reasonable scales.
c  Inputs:
c       xmin, xmax: range to be fitted.
c       ntics: (maximum) number of divisions (tics) to fit to it.
c  Outputs:
c       nxfac: power of ten by which the range is scaled.
c       xfac: 10**nxfac. World-value=xfac*axis-label.
c       xtic: Tic-spacing in world units.
c       xt1st: The integer multiple of xtic closest to xmin 
c             lying outside the range (xmin,xmax).
c       xtlast: The integer multiple of xtic closest to xmax
c             lying outside the range (xmin,xmax).
      real span,sfac
      parameter (npos=7)
      integer incpos(npos),nsfac
      data incpos/1,2,4,5,10,20,40/

      span=(xmax-xmin)
      if(xmax.eq.0. .and. xmin.eq.0)then
         write(*,*)'Fitrange error. xmin=xmax=0'
         nxfac=0
         span=1.
      else
         nxfac=nint(log10(max(abs(xmin),abs(xmax)))-0.4999999)
      endif
      xfac=10.**nxfac
      if(ntics.le.0)then
         write(*,'('' ntics<=0'',i8)')ntics
         return
      endif
      if(span.eq.0)then
         write(*,'('' Fitrange error span. Range:'',2g11.4)')xmin,xmax
         span=max(1.e-6,max(abs(xmin),abs(xmax)))
      endif
      xtic=span/ntics
      nsfac=nint(log10(0.099999*abs(xtic))+0.500001)
      sfac=10.**nsfac
      xtic=abs(xtic)/sfac
      if(xtic.lt.1.)then
         write(*,'('' Fitrange error 1. xtic='',f16.7)')xtic
      elseif(.not.xtic.lt.10.0001)then
         write(*,'('' Fitrange error NAN. Range:'',2g11.4)')xmin,xmax
         xtic=1.
         nxfac=0
         sfac=1.
         xfac=1.
      else
c Choose the increment
         iret=0
 201     fspan=1.e31
         ichoice=0
         do i=1,npos
            xt=incpos(i)*sfac
            xt=sign(xt,span)
            n2=nint((xmin+span)/xt-0.49999)
            n1=nint(xmin/xt+0.49999)
            atr=abs((n2-n1)*xt)
            if(iret.eq.1)write(*,'(2i3,4f7.3,3i4)')i,incpos(i),atr,xt
     $           ,sfac,fspan,n1,n2,ntics
c            if(atr.lt.fspan .and. abs(n2-n1).le.ntics)then
            if(atr.lt.fspan .and. abs(n2-n1).le.ntics
     $           .and. abs(n2-n1).ge.ntics/2)then
               ichoice=i
               xtic=incpos(i)
               fspan=atr*.9999
            endif
         enddo
         if(ichoice.eq.0 .and. iret.eq.0)then
            write(*,*)'Fitrange choice error',ichoice,nsfac
     $        ,sfac,xmin,xmax,span,ntics,xtic
            write(*,*)'i,incpos(i),atr,xt,sfac,fspan,n1,n2,ntics'
            iret=1
            goto 201
         endif
      endif
c Don't do this final setting till the end.
      xtic=sign(xtic,span)*sfac
c Use span to cope with error cases:
      xtlast=xtic*anint((xmin+span)/xtic+0.49999)
c Change xt1st the last thing, since it might be xmin, itself.
      xt1st=xtic*anint(xmin/xtic-0.49999)
c      write(*,*)'xtic,sfac,xt1st,xtlast',xtic,sfac,xt1st,xtlast,ntics
      return
      end
