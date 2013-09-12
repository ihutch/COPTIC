c Routines in support of using internal fonts rather than vector drawing
c of characters. Postscript here.
c**********************************************************************
      subroutine PSfontdefine(fname,fref)
c Define a font whose PS name is fname, scaled appropriately as fref.
      character*(*) fname, fref
      integer ipscale
      parameter (iunit=12)
      data htstnd/.015/wdstnd/.015/
      save htstnd,wdstnd
      ipscale=16
      shift=-.3
      if(fname(2:6).eq.'Times')then
         ipscale=18
         shift=-.26
      endif
c The scaling of PS fonts is different from hershey size.
      call abufstring(' ',iunit)
      call abufstring(fref,iunit)
      call abufstring(' { ',iunit)
      call abufstring(fname,iunit)
      call abufstring(' findfont ',iunit)
c Translate the font so that its reference point is same as Hershey.
      call abufstring(' 0. ',iunit)
      call fbufwrt(shift,2,iunit)
      call abufstring(' matrix translate makefont ',iunit)
c scale the font to approximately correct size
      crh=htstnd*1000*ipscale
c      crw=1.1*wdstnd*1000*ipscale
      crw=1.*wdstnd*1000*ipscale
      call fbufwrt(crw,4,iunit)
      call abufstring(' ',iunit)
      call fbufwrt(crh,4,iunit)
      call abufstring(' matrix scale makefont ',iunit)
      call abufstring('setfont } def ',iunit)
c      call abufstring(fref,iunit)
c      call abufstring(' exch definefont pop ',iunit)
      call crbufwrt(iunit)
      end
c**********************************************************************
      subroutine PSfontinit()
      include 'plotcom.h'
c Define the three Accis fonts to be the standard ones.
      call PSfontdefine('/Helvetica','/Accisfont0')
      call PSfontdefine('/Symbol','/Accisfont1')
      call PSfontdefine('/Times-Italic','/Accisfont2')
      psini=1
      end
c**********************************************************************
      subroutine PSscalefont(wd,ht)
c Scale the current font to correspond to normalized width wd and height ht.
c If there's no change, do nothing. So it is safe to call this often in the
c form call scalefont(chrswdth,chrshght).
c Zero width or height resets to default.
      real wd,ht
      integer iunit
      real htstnd,wdstnd
      integer ipoint
      parameter (ipoint=2)
      parameter (iunit=12)
c      data htstnd/.015/wdstnd/.015/
      data htstnd/.015/wdstnd/.015/
      save htstnd,wdstnd
      if(ht.eq.0.)ht=.015
      if(wd.eq.0.)wd=.015
      if(ht.ne.htstnd .or. wd.ne.wdstnd)then
         htrt=ht/htstnd
         wdrt=wd/wdstnd
         call abufwrt(' currentfont ',13,iunit)
c scale the font to approximately correct size
         call fbufwrt(wdrt,ipoint,iunit)
         call abufwrt(' ',1,iunit)
         call fbufwrt(htrt,ipoint,iunit)
         call abufwrt(' matrix scale makefont setfont ',31,iunit)
      endif
      end
c***********************************************************************
      subroutine PSsetfont(ifont)
c Set to font 0,1,2 accoring to ifont. Then scale and rotate if nec.
      integer ifont,iunit
      include 'plotcom.h'
      character*(12) cl
      parameter (iunit=12)      
      cl='Accisfont'
      call iwrite(ifont,iwidth,cl(10:))
      call abufstring(cl(1:10),iunit)
      call abufstring(' ',iunit)
c      call abufstring(' findfont setfont ',iunit)
      call PSscalefont(chrswdth,chrshght)
c Rotate font to the current character angle if necessary
      if(chrssin.ne.0..or.chrscos.ne.1.)then
         call getcangl(angle)
         call PSfontangle(angle)
      endif
      end
c***********************************************************************
      subroutine PSstartfont(fname)
c Find a named PS font and scale it to current charsize. PS
      character*(*) fname
      include 'plotcom.h'
      integer ipscale
      integer iunit
      parameter (iunit=12)      
      ipscale=16
      if(fname(2:6).eq.'Times')ipscale=18
c The scaling of PS fonts is different from hershey size.
      call abufwrt(' ',1,iunit)
      call abufstring(fname,iunit)
      call abufstring(' findfont ',iunit)
c Translate the font so that its reference point is same as Hershey.
      call abufstring(' 0. -.3 matrix translate makefont ',iunit)
c scale the font to approximately correct size
      crh=chrshght*1000*ipscale
c      crw=1.1*chrswdth*1000*ipscale
      crw=1.*chrswdth*1000*ipscale
      call fbufwrt(crw,4,iunit)
      call abufstring(' ',iunit)
      call fbufwrt(crh,4,iunit)
      call abufstring(' matrix scale makefont setfont ',iunit)
c Rotate font to the current character angle if necessary
      if(chrssin.ne.0..or.chrscos.ne.1.)then
         call getcangl(angle)
         call PSfontangle(angle)
      endif
      end
c***********************************************************************
      subroutine vecnnops(nx,ny,ud)
c Defeat writing to the PS file. For text when using hard fonts.
      real nx,ny
      integer ud
      include 'plotcom.h'
      pfin=pfsw
      pfsw=0
      call vecn(nx,ny,ud)
      pfsw=pfin
      end
c***********************************************************************
      subroutine abufstring(string,iunit)
      character*(*) string
      integer iunit
c add the full length of string to output buffer
      ilen=len(string)
      if(ilen.gt.0) then
         call abufwrt(string,ilen,iunit)
      else
         write(*,*) 'Incorrect ilen', ilen,'  in abufstring'
      endif
      end
c***********************************************************************
      subroutine crbufwrt(iunit)
c Do an end of line on the output string buffer.
      integer sblen,iunit
      character*80 sbuf
      common /wbuf/sblen,sbuf
      write(iunit,*)sbuf(1:sblen-1)
      sblen=1
      end
c***********************************************************************
      subroutine PSchardrw(char)
      character*(*) char
c Ghastly kludge because concatenation of *(*) does not work.
      character*(256) clocal
      parameter (iunit=12)
      il=len(char)
      clocal(1:1)='('
      clocal(2:2+il-1)= char(1:len(char))
      clocal(2+il:9+il)=') show '
c Have to do this in one write so that it does not get split.
      call abufwrt(clocal,il+9,iunit)
      end
c**********************************************************************
      subroutine PSfontangle(angle)
c Rotate the standard font *to* angle in degrees, 
      real angle
      integer iunit
      parameter (ipoint=3)
      parameter (iunit=12)
      if(angle.ne.0.)then
         call abufstring(' currentfont ',iunit)
         call fbufwrt(angle,ipoint,iunit)
         call abufstring(' matrix rotate makefont setfont ',iunit)
      endif
      end
c**********************************************************************
      function iPSsymsub(ic)
c Substitute the type1 symbol font position of an accis symbol character
      integer iPSsymsub,ic
      integer isym(0:127)
      data isym/32*0,
     $     32,124,171,35,36,184,181,162,206,214,42,177,201,45,215,214,
     $     79,196,197,183,168,170,210,211,224,169,92,165,163,186,179,204
     $     ,182,192,66,88,68,36,70,71,35,193,242,209,76,77,199,176,
     $     80,81,194,83,200,161,86,87,180,89,198,91,185,93,172,174,
     $     94,97,98,120,100,101,102,103,104,105,166,107,108,109,110,101,
     $     112,113,114,115,116,117,218,119,99,121,122,175,124,173,187,32
     $     /

      iPSsymsub=isym(ic)
      end
c*********************************************************************
      subroutine pfPSset(ips)
c Set the PSfonts on or off.      
      include 'plotcom.h'
      integer ipfpsprev
      data ipfpsprev/0/
      save
c Need to decide if we are ready to write to file.
      if(psini.eq.2) call PSfontinit()
      itemp=pfPS
      if(ips.eq.3)then
         pfPS=ipfpsprev
      else
         pfPS=ips
      endif
      ipfpsprev=itemp
      end
c*********************************************************************
      function ipfps()
      include 'plotcom.h'
      ipfps=pfPS
      end
