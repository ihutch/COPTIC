c*******************************************************************
      subroutine drwstrdo( px, py, str1,drw,width)
c  Draw string code. drw.ne.0 means actually draw, else just return length.
c  Draw a string str1 from the point px,py, (right justified).  */
c  Return the final normalized position in px,py.
c  Sep 92 revision to allow use of only printable characters (VAX set)
c  Also removal of the ^A and ^B fonts switches. All controls are now
c  via special switch: ! with subsequent character used for
c  control as follows:
c       @       to normal font 1.
c       A       to font 1 (math)
c       B       to font 2 (italic)
c       R       to font 3 (roman)
c       E       to font 4 (English gothic)
c       G       to font 5 (German gothic)
c       d/D     Toggle subscript mode
c       u       Toggle superscript mode
c       p       Save this position
c       q       Return to saved position (and save current position)
c       o       Toggle Overbar mode.
c       n/N     Toggle Underbar mode
      real px,py,width
      character*(*) str1
      integer drw
      include 'plotcom.h'
      integer j,n1,cwidth,isw,lenstr
      real dd
      integer STDWDTH,offset,offprev
      parameter (STDWDTH=21)
c The following statement is to ensure that the block data program 
c fontdatais linked in when a program is compiled using character fonts.
      external fontdata
      save

      lenstr=len(str1)
      width=0
c Reset to standard font at start of each string.
      offset=0
      if(pfPS.eq.1 .and. pfsw.ge.2)call PSsetfont(offset/128)
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
            if(pfPS.eq.1 .and. pfsw.ge.2)call PSsetfont(offset/128)
            goto 1
         endif
c No valid switch character found after \. Interpret as plain.
         n1=ichar(str1(j:j))
c Fixed for f2c/gcc upper characters.
         if(n1 .lt. 0) n1=n1+256
      endif
      if(n1.eq.0) goto 2
      if(pfPS.eq.1 .and. pfsw.ge.2 .and. drw.ne.0)then
         if(offset.eq.128)then
            n2=iPSsymsub(n1)
         else
            n2=n1
         endif
         if(n2.eq.40.or.n2.eq.41.or.n2.eq.92)then
            call PSchardrw(char(92)//char(n2))
         else
            if(char(n2).eq.'-')then
c Make the PS minus sign the version from symbol font. It's longer.
               call PSsetfont(1)
               call PSchardrw(char(n2))
               call PSsetfont(offset/128)
            else
               call PSchardrw(char(n2))
            endif
         endif
      endif
      n1=n1+offset
      call drwchar(n1,px,py,drw,cwidth)
      dd=chrswdth*cwidth/STDWDTH
      width=width+dd
      px=px+chrscos*dd
      py=py+chrssin*dd
      if(drw.ne.0)then 
         if(pfPS.eq.0) then
            call vecn(px,py,0)
         else
c Use spaces to accommodate width variation between sysfont and Hershey
c Also if we are approaching a special character, force a space to ensure
c that the subscript or superscript etc is well aligned wrt its parent.
            if(n1.eq.ichar(' ').or.str1(j+2:j+2).eq.'!')then
c This call repositions after each PS character using Hershey width.
               call vecn(px,py,0)
            else
c This call does not. It uses the internal PS width.
               call vecnnops(px,py,0)
            endif
         endif
      endif
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
      real scrsx,scrsy,sx,sy
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
      elseif(sw.eq.114) then
c Font shift: \R
         offset=384
      elseif(sw.eq.101) then
c Font shift: \E
         offset=512
      elseif(sw.eq.103) then
c Font shift: \G
         offset=640
      elseif(sw.eq.92.or.sw.eq.33)then
c Not a switch literal: \\, or \!, or !!.
         isw=0
      elseif(sw.eq.100.or.sw.eq.117.or.sw.eq.111.or.sw.eq.110)then
c /* Toggle sUper/sub-script mode or Overbar/uNderbar mode */
         if(su.eq.0)then
            if(sw.eq.ichar('u').or.sw.eq.ichar('o'))then
               sgn=1.
            elseif(sw.eq.ichar('d').or.sw.eq.ichar('n'))then
               sgn=-1.
            endif
            if(sw.eq.100.or.sw.eq.117)then
               dxc=0.6
               hs=0.7
            elseif(sw.eq.111.or.sw.eq.110)then
               dxc=1.
               hs=1.
            endif
            dx=-chrshght*dxc*chrssin*sgn
            dy= chrshght*dxc*chrscos*sgn
            crsrx=crsrx+dx
            crsry=crsry+dy
c Needed to ensure that the ps fonts move down.
            call vecn(crsrx,crsry,0)
            height=chrshght
            width=chrswdth
            chrshght=hs*height
            chrswdth=hs*width
            if(pfPS.eq.1 .and. pfsw.ge.2) call PSsetfont(offset/128)
            su=1
         else
            crsrx=crsrx-dx
            crsry=crsry-dy
            call vecn(crsrx,crsry,0)
            chrshght=height
            chrswdth=width
            if(pfPS.eq.1 .and. pfsw.ge.2)call PSsetfont(offset/128)
            su=0
         endif
      elseif(sw.eq.112)then
c save position
         scrsx=crsrx
         scrsy=crsry
      elseif(sw.eq.113)then
c return to saved position
         sx=crsrx
         sy=crsry
         call vecn(scrsx,scrsy,0)
         scrsx=sx
         scrsy=sy
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
               if(pfPS.eq.0) then
                  call vecn((px+xn),(py+yn),ud)
               else
                  call vecnnops((px+xn),(py+yn),ud)
               endif
               ud=1
            else
               ud=0
            endif
            if((xc.eq.-64).and.(yc.eq.-64)) goto 4
    3    continue
    4    continue
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
      subroutine jdrcstr(str1,js)
c draw a justified string relative to the current position.
      character*(*) str1
      real js
      include 'plotcom.h'
      real px,py,dd,wstr,width
      external wstr
      px=crsrx
      py=crsry
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
c Changed May 2001 to leave px,py unchanged.
      character*(*) str1
      real px,py,width
      real pix,piy
      pix=px
      piy=py
      call drwstrdo(pix,piy,str1,1,width)
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
