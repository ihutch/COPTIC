c*********************************************************************
c	Set character width and height in NORMAL units. 0 resets.
      subroutine charsize(wd, ht)
      real wd,ht
      include 'plotcom.h'
      chrswdth=wd
      if(abs(wd) .le. 0.)chrswdth=0.015
      chrshght=ht
      if(abs(ht) .le. 0.)chrshght=0.015
      return
      end
c*********************************************************************
c    Set character angle in degrees.
	subroutine charangl(theta)
	real theta,tr
      include 'plotcom.h'
	tr=3.14159*theta/180.
	chrscos=cos(tr)
	chrssin=sin(tr)
	end
c*********************************************************************
c     Get character angle in degrees.
      subroutine getcangl(theta)
      real theta
      include 'plotcom.h'
      theta=atan2(chrssin,chrscos)
      theta=theta*180./3.14159
      end
c*************************************************************************
c********************************************************************
      integer function istlen(str,ilmax)
c Return the active length of string str, length ilmax, ignoring 
c blanks on the end.
      integer ilmax
      character*(*) str
      do 1 istlen=ilmax,1,-1
	 if(str(istlen:istlen) .ne. ' ') return
    1 continue
      end

c**************************************************************************
c Convert spaces at the end of string to 0 so it is a C-string termination.
      subroutine termchar(string)
      character*(*) string
      do 1 i=len(string),1,-1
         if(ichar(string(i:i)) .ne. ichar(' ')) goto 2
         string(i:i)=char(0)
 1    continue
 2    end
c********************************************************************
      integer function istpos(str,ilmin,ilmax,ch)
c Return the first position of character ch in string str,
c searching from position ilmin to ilmax.
      integer ilmin,ilmax
      character*(*) str
      character*1 ch
      if(ilmax.gt.ilmin) then
	 do 1 istpos=ilmin,ilmax
	    if(str(istpos:istpos) .eq. ch) return
    1	 continue
      else
	 do 2 istpos=ilmin,ilmax,-1
	    if(str(istpos:istpos) .eq. ch) return
    2	 continue
      endif
      istpos=0
      end
c******************************//**************************************
c Subroutine for decently formatted floating write.
c Returns total  width of   string.
      subroutine fwrite(x,width,point,string)
      real x,xx
      integer width,point,iax,neg,ibx,i
      character*(*) string
c      character*8 sformat
      character*1 st(20)
c Use g format if too long.
      xx=x*10**point
      if(abs(xx).gt.1.e9)then
	 write(string,'(g12.6)')x
	 width=12
	 return
      endif
c Else standard minimum f.point
      iax=abs(xx)
      if(abs(xx)-iax.gt.0.5)iax=iax+1
      width=0
      do 2 i=1,point
	 width=width+1
	 ibx=iax/10
	 st(width)=char(48+iax-10*ibx)
	 iax=ibx
    2 continue
      width=width+1
      st(width)='.'
    1 width=width+1
      ibx=iax/10
      st(width)=char(48+iax-10*ibx)
      iax=ibx
      if(iax .ne. 0) goto 1
      neg=0
      if(x.lt.0) then
	 string(1:1)='-'
	 neg=1
	 width=width+1
      endif
      do 3 i=1+neg,width
	 string(i:i)=st(width-i+1)
    3 continue
c Put a null at the end so we can treat as a C-string.
      string(width+1:width+1)=char(0)
      return
      end

c********************************************************************
c Subroutine for decently formatted integer write.
c Returns total  width of string. 
c Avoid write statements. Version 4. Nearly 4 timesfaster than version 3.
      subroutine iwrite(i,width,string)
      integer i,ai,aj
      integer width,nc,neg
      character*(*) string
      character*(1) st(20)

      ai=abs(i)
      nc=0
    1 nc=nc+1
       aj=ai/10
       st(nc)=char(48+ai-10*aj)
       ai=aj
      if(ai .ne. 0) goto 1
      width=nc
      neg=0
      if(i .lt. 0) then
	 neg=1
	 string(1:1)='-'
      endif
      do 2 nc=1,width
	 string(nc+neg:nc+neg)=st(width+1-nc)
    2 continue
      if(neg .eq.1) width=width+1
c Put a null at the end so we can treat as a C-string.
      string(width+1:width+1)=char(0)
      return
      end

