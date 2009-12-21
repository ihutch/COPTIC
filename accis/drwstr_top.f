c*******************************************************************
      subroutine drwstrdo( px, py, str1,drw,width)
c  Draw string code. drw.ne.0 means actually draw, else just return length.
c  Draw a string str1 from the point px,py, (right justified).  */
c  Return the final normalized position in px,py.
c  Sep 92 revision to allow use of only printable characters (VAX set)
c  Also removal of the ^A and ^B fonts switches. All controls are now
c  via special switch: \ or ! (or ^\) with subsequent character used for
c  control as follows:
c	@	to normal font 1.
c	A	to font 1 (math)
c	B	to font 2 (italic)
c	D	Toggle subscript mode
c	U	Toggle superscript mode
      real px,py,width
      character*(*) str1
      integer drw
      include 'plotcom.h'
      integer j,n1,cwidth,isw,lenstr
      real dd
      integer STDWDTH,offset
      parameter (STDWDTH=21)
      save

      lenstr=len(str1)
      width=0
c Reset to standard font at start of each string.
      offset=0
      if(drw.ne.0) call vecn(px,py,0)
      j=1
    1 continue
	if(j.gt.lenstr)goto 2
c Now making everything a switch not momentary.
      n1=ichar(str1(j:j))
c Fix for f2c characters over 128.
      if(n1 .lt. 0) n1=n1+256
c Jan 97 ! switch upgrade.
      if(n1.eq.28.or.n1.eq.92.or.n1.eq.33.or.n1.eq.1.or.n1.eq.2)then
c Special switch: ctrl-\, or \. Get the next and handle it.
         if(n1.eq.1.or.n1.eq.2)then
c Grandfathering of old style. ctrl-A, ctrl-B
            n1=n1+96
	 else	
	    j=j+1
            n1=ichar(str1(j:j))
c If we ended jump out Aug 98.
            if(n1.eq.0) goto 2
	 endif
         call spstr(n1,offset,isw)
         j=j+isw
         if(isw.eq.1) then
c Switch operated; adjust position and start on next character.
            px=crsrx
	    py=crsry
	    goto 1
         endif
c No valid switch character found after \. Interpret as plain.
         n1=ichar(str1(j:j))
c Fixed for f2c/gcc upper characters.
	 if(n1 .lt. 0) n1=n1+256
      endif
      if(n1.eq.0) goto 2
      n1=n1+offset
c
      call drwchar(n1,px,py,drw,cwidth)
      dd=chrswdth*cwidth/STDWDTH
      width=width+dd
      px=px+chrscos*dd
      py=py+chrssin*dd
      if(drw.ne.0)call vecn(px,py,0)
c  Terminate after end of string. 
c   crsrx crsry contain end crsr. */
      j=j+1
      if(j.le.300) goto 1
    2 continue
      return
      end
c*********************************************************************
c         Special string switch handler   */
      subroutine spstr( n1, offset,isw)
      integer n1,sw,offset,isw
      include 'plotcom.h'
      real  dx,dy,height,width,sgn
      integer su
      save
      data su/0/
      isw=1
      sw=n1
c  sw=tolower(sw)
      if(sw.lt.91.and.sw.gt.63)sw=sw+32
      if(sw.eq.96)then
c Standard font: (\@)
         offset=0
      elseif(sw.eq.97) then
c Font shift: \A.
         offset=128
      elseif(sw.eq.98) then
c Font shift: \B
         offset=256
      elseif(sw.eq.92.or.sw.eq.33)then
c Not a switch literal: \\, or \!, or !!.
	 isw=0
      elseif(sw.eq.100.or.sw.eq.117)then
c /* Toggle super/sub-script mode  */
	 if(su.eq.0)then
	    if(sw.eq.ichar('u'))then
	       sgn=1.
	    elseif(sw.eq.ichar('d'))then
	       sgn=-1.
	    endif
	    dx=-chrshght*.6*chrssin*sgn
	    dy= chrshght*.6*chrscos*sgn
	    crsrx=crsrx+dx
	    crsry=crsry+dy
	    height=chrshght
	    width=chrswdth
	    chrshght=0.7*height
	    chrswdth=0.7*width
	    su=1
	 else
	    crsrx=crsrx-dx
	    crsry=crsry-dy
	    chrshght=height
	    chrswdth=width
	    su=0
	 endif
      else
c Not a valid switch
	 isw=-1
      endif
      return
      end

c************************************************************************
      subroutine drwchar(n1,px,py,drw,width)
c Primitive: draw (if drw.ne.0) character at px,py.  Return width.
      integer n1,width,drw
      real px,py
      integer base,xci,xca,ud,i,n2,xc,yc
      real xcn,ycn,xn,yn,ch,cw
      INCLUDE 'plotcom.h'
      integer STDWDTH
      parameter (STDWDTH=21)
      i=0
      base=chrsaddr(n1+1)
      i=i+1
      xci=ichar(chrsfont(base+i))-97
      i=i+1
      xca=ichar(chrsfont(base+i))-97
      width=xca-xci
      ch=chrshght/STDWDTH
      cw=chrswdth/STDWDTH
      if(drw.ne.0)then
	 ud=0
	 do 3 n2=1,300
	    i=i+1
	    xc=ichar(chrsfont(base+i))-97
	    i=i+1
	    yc=ichar(chrsfont(base+i))-97
c F2C fix. 25 May 96.
	    if(yc .lt. -97) yc=yc+256
	    if(xc.ne.-64)then
	       ycn=-yc*ch
	       xcn=(xc-xci)*cw +chrsslnt*ycn
	       xn=xcn*chrscos - ycn*chrssin
	       yn=xcn*chrssin + ycn*chrscos
	       call vecn((px+xn),(py+yn),ud)
	       ud=1
	    else
	       ud=0
	    endif
	    if((xc.eq.-64).and.(yc.eq.-64)) goto 4
    3	 continue
    4	 continue
      endif
      end
c*********************************************************************
      subroutine jdrwstr(pex, pey, str1, js)
c   Draw a string (norm-units) justified relative to px,py per parameter 
c   js :  0 => centered, -1. => right-just., +1. => left-just.
c   Do not change the external positions pex and pey.
      real pex,pey,js
      character*(*) str1
      INCLUDE 'plotcom.h'
      real dd,px,py
      real wstr,width

      px=pex
      py=pey
      width=wstr(str1)
c/* Offset and draw */
      dd=0.5*(js-1.)*width
      px=px+chrscos*dd
      py=py+chrssin*dd
      call drwstr(px,py,str1)
      return
      end

c********************************************************************/
      subroutine drcstr(str1)
c draw a string from current position. Leave at end of string.
      character*(*) str1
      include 'plotcom.h'
      real px,py
      px=crsrx
      py=crsry
      call drwstr(px,py,str1)
      return
      end
c********************************************************************/
      subroutine drwstr(px,py,str1)
c draw a string from norm position (px,py). leave at end of string.
      character*(*) str1
      real px,py,width
      call drwstrdo(px,py,str1,1,width)
      return
      end
c*********************************************************************
      function wstr(str1)
c return the normalized length of the string str1 using current widths.
      character*(*) str1
      real x,y,width
      real wstr
      call drwstrdo(x,y,str1,0,width)
      wstr=width
      end
c*******************************************************************
