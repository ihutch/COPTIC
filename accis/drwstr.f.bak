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
      integer STDWDTH,offset,offprev
      parameter (STDWDTH=21)
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
            call PSchardrw(char(n2))
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
            call vecnnops(px,py,0)
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
c Needed to ensure that the ps fonts move down.
            call vecn(crsrx,crsry,0)
	    height=chrshght
	    width=chrswdth
	    chrshght=0.7*height
	    chrswdth=0.7*width
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
c Be careful editing this file. Emacs breaks it.
      block data fonts
      character*19000 ca
      integer*2 a(384)
      include 'plotcom.h'
      character*2 crlf
            parameter (crlf='')
      equivalence (ca,chrsfont)
      equivalence (a,chrsaddr)
      data chrscos,chrssin,chrsslnt,chrswdth
     $   ,chrshght/ 1.,0.,0.,.015,.015 /
      data (a(j),j=  1,  6)/ 8197, 8197, 8197, 8197, 8197, 8197/
      data (a(j),j=  7, 12)/ 8197, 8197, 8197, 8197, 8197, 8197/
      data (a(j),j= 13, 18)/ 8197, 8197, 8197, 8197, 8197, 8197/
      data (a(j),j= 19, 24)/ 8197, 8197, 8197, 8197, 8197, 8197/
      data (a(j),j= 25, 30)/ 8197, 8197, 8197, 8197, 8197, 8197/
      data (a(j),j= 31, 36)/ 8197, 8197, 2443, 2451, 2505, 2555/
      data (a(j),j= 37, 42)/ 2597, 2705, 2775, 2923, 2951, 3005/
      data (a(j),j= 43, 48)/ 3059, 3143, 3181, 3235, 3257, 3295/
      data (a(j),j= 49, 54)/ 3317, 3407, 3437, 3511, 3613, 3655/
      data (a(j),j= 55, 60)/ 3767, 3897, 3927, 4069, 4199, 4269/
      data (a(j),j= 61, 66)/ 4355, 4369, 4407, 4421, 4543, 4661/
      data (a(j),j= 67, 72)/ 4707, 4801, 4883, 4953, 5013, 5061/
      data (a(j),j= 73, 78)/ 5155, 5205, 5227, 5273, 5323, 5357/
      data (a(j),j= 79, 84)/ 5415, 5461, 5547, 5607, 5709, 5783/
      data (a(j),j= 85, 90)/ 5875, 5911, 5965, 5999, 6057, 6095/
      data (a(j),j= 91, 96)/ 6135, 6181, 6211, 6223, 6253, 6275/
      data (a(j),j= 97,102)/ 6287, 6341, 6419, 6497, 6567, 6645/
      data (a(j),j=103,108)/ 6723, 6777, 6879, 6935, 6989, 7043/
      data (a(j),j=109,114)/ 7093, 7115, 7205, 7261, 7339, 7417/
      data (a(j),j=115,120)/ 7495, 7543, 7649, 7687, 7743, 7777/
      data (a(j),j=121,126)/ 7835, 7873, 7913, 7959, 8045, 8057/
      data (a(j),j=127,132)/ 8143, 8197,13243,13243,13243,13243/
      data (a(j),j=133,138)/13243,13243,13243,13243,13243,13243/
      data (a(j),j=139,144)/13243,13243,13243,13243,13243,13243/
      data (a(j),j=145,150)/13243,13243,13243,13243,13243,13243/
      data (a(j),j=151,156)/13243,13243,13243,13243,13243,13243/
      data (a(j),j=157,162)/13243,13243,13243,13243, 8205, 8213/
      data (a(j),j=163,168)/ 8235, 8257, 8299, 8381, 8417, 8465/
      data (a(j),j=169,174)/ 8483, 8521, 8547, 8571, 8595, 8627/
      data (a(j),j=175,180)/ 8651, 8681, 8707, 8749, 8767, 8783/
      data (a(j),j=181,186)/ 8859, 8919, 8959, 8989, 9039, 9057/
      data (a(j),j=187,192)/ 9085, 9127, 9185, 9211, 9235, 9261/
      data (a(j),j=193,198)/ 9293, 9387, 9485, 9603, 9681, 9717/
      data (a(j),j=199,204)/ 9747, 9849, 9883, 9925, 9997,10067/
      data (a(j),j=205,210)/10103,10139,10239,10271,10305,10353/
      data (a(j),j=211,216)/10471,10607,10653,10685,10757,10823/
      data (a(j),j=217,222)/10915,10933,11021,11127,11161,11185/
      data (a(j),j=223,228)/11219,11247,11275,11293,11379,11499/
      data (a(j),j=229,234)/11595,11689,11759,11849,11911,11981/
      data (a(j),j=235,240)/12017,12093,12157,12209,12271,12325/
      data (a(j),j=241,246)/12381,12431,12523,12591,12667,12705/
      data (a(j),j=247,252)/12773,12791,12885,12937,13011,13079/
      data (a(j),j=253,258)/13107,13125,13153,13243,18873,18873/
      data (a(j),j=259,264)/18873,18873,18873,18873,18873,18873/
      data (a(j),j=265,270)/18873,18873,18873,18873,18873,18873/
      data (a(j),j=271,276)/18873,18873,18873,18873,18873,18873/
      data (a(j),j=277,282)/18873,18873,18873,18873,18873,18873/
      data (a(j),j=283,288)/18873,18873,18873,18873,18873,18873/
      data (a(j),j=289,294)/13251,13259,13299,13317,13351,13439/
      data (a(j),j=295,300)/13509,13625,13647,13693,13739,13763/
      data (a(j),j=301,306)/13801,13823,13845,13863,13885,13975/
      data (a(j),j=307,312)/14011,14101,14207,14233,14317,14415/
      data (a(j),j=313,318)/14481,14613,14711,14739,14773,14793/
      data (a(j),j=319,324)/14831,14851,14925,15043,15085,15173/
      data (a(j),j=325,330)/15247,15313,15363,15409,15499,15559/
      data (a(j),j=331,336)/15589,15637,15697,15731,15797,15845/
      data (a(j),j=337,342)/15935,15995,16123,16213,16289,16327/
      data (a(j),j=343,348)/16383,16419,16473,16521,16567,16605/
      data (a(j),j=349,354)/16635,16647,16677,16699,16711,16733/
      data (a(j),j=355,360)/16817,16895,16951,17041,17099,17175/
      data (a(j),j=361,366)/17267,17335,17393,17463,17537,17579/
      data (a(j),j=367,372)/17689,17769,17839,17929,18001,18059/
      data (a(j),j=373,378)/18121,18163,18243,18301,18389,18479/
      data (a(j),j=379,384)/18567,18635,18721,18733,18819,18873/
      data ca(00001:00042)/'16439 data bytes after 24 index lines, fro'/
      data ca(00043:00084)/'m files sserif.dat and greek.dat and ital.'/
      data ca(00085:00087)/'dat'/
      data ca(00088:00089)/crlf/
      data ca(00090:00131)/'  5756  5756  5756  5756  5756  5756  5756'/
      data ca(00132:00173)/'  5756  5756  5756  5756  5756  5756  5756'/
      data ca(00174:00185)/'  5756  5756'/
      data ca(00186:00187)/crlf/
      data ca(00188:00229)/'  5756  5756  5756  5756  5756  5756  5756'/
      data ca(00230:00271)/'  5756  5756  5756  5756  5756  5756  5756'/
      data ca(00272:00283)/'  5756  5756'/
      data ca(00284:00285)/crlf/
      data ca(00286:00327)/'     2    10    64   114   156   264   334'/
      data ca(00328:00369)/'   482   510   564   618   702   740   794'/
      data ca(00370:00381)/'   816   854'/
      data ca(00382:00383)/crlf/
      data ca(00384:00425)/'   876   966   996  1070  1172  1214  1326'/
      data ca(00426:00467)/'  1456  1486  1628  1758  1828  1914  1928'/
      data ca(00468:00479)/'  1966  1980'/
      data ca(00480:00481)/crlf/
      data ca(00482:00523)/'  2102  2220  2266  2360  2442  2512  2572'/
      data ca(00524:00565)/'  2620  2714  2764  2786  2832  2882  2916'/
      data ca(00566:00577)/'  2974  3020'/
      data ca(00578:00579)/crlf/
      data ca(00580:00621)/'  3106  3166  3268  3342  3434  3470  3524'/
      data ca(00622:00663)/'  3558  3616  3654  3694  3740  3770  3782'/
      data ca(00664:00675)/'  3812  3834'/
      data ca(00676:00677)/crlf/
      data ca(00678:00719)/'  3846  3900  3978  4056  4126  4204  4282'/
      data ca(00720:00761)/'  4336  4438  4494  4548  4602  4652  4674'/
      data ca(00762:00773)/'  4764  4820'/
      data ca(00774:00775)/crlf/
      data ca(00776:00817)/'  4898  4976  5054  5102  5208  5246  5302'/
      data ca(00818:00859)/'  5336  5394  5432  5472  5518  5604  5616'/
      data ca(00860:00871)/'  5702  5756'/
      data ca(00872:00873)/crlf/
      data ca(00874:00915)/' 10802 10802 10802 10802 10802 10802 10802'/
      data ca(00916:00957)/' 10802 10802 10802 10802 10802 10802 10802'/
      data ca(00958:00969)/' 10802 10802'/
      data ca(00970:00971)/crlf/
      data ca(00972:01013)/' 10802 10802 10802 10802 10802 10802 10802'/
      data ca(01014:01055)/' 10802 10802 10802 10802 10802 10802 10802'/
      data ca(01056:01067)/' 10802 10802'/
      data ca(01068:01069)/crlf/
      data ca(01070:01111)/'  5764  5772  5794  5816  5858  5940  5976'/
      data ca(01112:01153)/'  6024  6042  6080  6106  6130  6154  6186'/
      data ca(01154:01165)/'  6210  6240'/
      data ca(01166:01167)/crlf/
      data ca(01168:01209)/'  6266  6308  6326  6342  6418  6478  6518'/
      data ca(01210:01251)/'  6548  6598  6616  6644  6686  6744  6770'/
      data ca(01252:01263)/'  6794  6820'/
      data ca(01264:01265)/crlf/
      data ca(01266:01307)/'  6852  6946  7044  7162  7240  7276  7306'/
      data ca(01308:01349)/'  7408  7442  7484  7556  7626  7662  7698'/
      data ca(01350:01361)/'  7798  7830'/
      data ca(01362:01363)/crlf/
      data ca(01364:01405)/'  7864  7912  8030  8166  8212  8244  8316'/
      data ca(01406:01447)/'  8382  8474  8492  8580  8686  8720  8744'/
      data ca(01448:01459)/'  8778  8806'/
      data ca(01460:01461)/crlf/
      data ca(01462:01503)/'  8834  8852  8938  9058  9154  9248  9318'/
      data ca(01504:01545)/'  9408  9470  9540  9576  9652  9716  9768'/
      data ca(01546:01557)/'  9830  9884'/
      data ca(01558:01559)/crlf/
      data ca(01560:01601)/'  9940  9990 10082 10150 10226 10264 10332'/
      data ca(01602:01643)/' 10350 10444 10496 10570 10638 10666 10684'/
      data ca(01644:01655)/' 10712 10802'/
      data ca(01656:01657)/crlf/
      data ca(01658:01699)/' 16432 16432 16432 16432 16432 16432 16432'/
      data ca(01700:01741)/' 16432 16432 16432 16432 16432 16432 16432'/
      data ca(01742:01753)/' 16432 16432'/
      data ca(01754:01755)/crlf/
      data ca(01756:01797)/' 16432 16432 16432 16432 16432 16432 16432'/
      data ca(01798:01839)/' 16432 16432 16432 16432 16432 16432 16432'/
      data ca(01840:01851)/' 16432 16432'/
      data ca(01852:01853)/crlf/
      data ca(01854:01895)/' 10810 10818 10858 10876 10910 10998 11068'/
      data ca(01896:01937)/' 11184 11206 11252 11298 11322 11360 11382'/
      data ca(01938:01949)/' 11404 11422'/
      data ca(01950:01951)/crlf/
      data ca(01952:01993)/' 11444 11534 11570 11660 11766 11792 11876'/
      data ca(01994:02035)/' 11974 12040 12172 12270 12298 12332 12352'/
      data ca(02036:02047)/' 12390 12410'/
      data ca(02048:02049)/crlf/
      data ca(02050:02091)/' 12484 12602 12644 12732 12806 12872 12922'/
      data ca(02092:02133)/' 12968 13058 13118 13148 13196 13256 13290'/
      data ca(02134:02145)/' 13356 13404'/
      data ca(02146:02147)/crlf/
      data ca(02148:02189)/' 13494 13554 13682 13772 13848 13886 13942'/
      data ca(02190:02231)/' 13978 14032 14080 14126 14164 14194 14206'/
      data ca(02232:02243)/' 14236 14258'/
      data ca(02244:02245)/crlf/
      data ca(02246:02287)/' 14270 14292 14376 14454 14510 14600 14658'/
      data ca(02288:02329)/' 14734 14826 14894 14952 15022 15096 15138'/
      data ca(02330:02341)/' 15248 15328'/
      data ca(02342:02343)/crlf/
      data ca(02344:02385)/' 15398 15488 15560 15618 15680 15722 15802'/
      data ca(02386:02427)/' 15860 15948 16038 16126 16194 16280 16292'/
      data ca(02428:02439)/' 16378 16432'/
      data ca(02440:02441)/crlf/
      data ca(02442:02447)/'kvYi!!'/
      data ca(02448:02449)/crlf/
      data ca(02450:02491)/'zv\gaUacbc!aaUbUbc!aag`h`iajbjcichbgag!aah'/
      data ca(02492:02501)/'aibibhah!!'/
      data ca(02502:02503)/crlf/
      data ca(02504:02545)/'àvXj]U\V\\!a]V\\!a]U^V\\!afUeVe\!afVe\!afU'/
      data ca(02546:02551)/'gVe\!!'/
      data ca(02552:02553)/crlf/
      data ca(02554:02593)/'•sYi^V^m!adUdl!a\]f[!a\^f\!a\ffd!a\gfe!!'/
      data ca(02594:02595)/crlf/
      data ca(02596:02637)/'vXkaQanbn!aaQbQbn!afXhXfVcU`U]V[X[Z\\]]ea'/
      data ca(02638:02679)/'fbgdgffhci`i^h]g!afXeWcV`V]W\X\Z]\e`gbhdhf'/
      data ca(02680:02701)/'ghficj`j]i[g]g!aggdi!!'/
      data ca(02702:02703)/crlf/
      data ca(02704:02745)/'ørUmjUXj!a]U_W_Y^[\\Z\XZXXYV[U]U_VbWeWhVjU'/
      data ca(02746:02771)/'!afcddcfchejgjiijgjehcfc!!'/
      data ca(02772:02773)/crlf/
      data ca(02774:02815)/'~vUnk\i\g]f_decgbh`i\iZhYfYdZb[a`^b\cZcXbV'/
      data ca(02816:02857)/'`U_U]V\X\Z]]_`dfgiijkj!ak\k]i]g^!ah]g_eedg'/
      data ca(02858:02899)/'bi`j\jZiYhXfXdYb[``]a\bZbXaV!abW`V_V]W!a^V'/
      data ca(02900:02919)/']X]Z^]``efghiikikj!!'/
      data ca(02920:02921)/crlf/
      data ca(02922:02947)/'áv]fbUaVa\!abVa\!abUcVa\!!'/
      data ca(02948:02949)/crlf/
      data ca(02950:02991)/'ÅvZhdQbS`V^Z]_]c^h`lbodqeq!adQeQcSaV_Z^_^c'/
      data ca(02992:03001)/'_halcoeq!!'/
      data ca(03002:03003)/crlf/
      data ca(03004:03045)/'ÇvZh]Q_SaVcZd_dcchal_o]q^q!a]Q^Q`SbVdZe_ec'/
      data ca(03046:03055)/'dhbl`o^q!!'/
      data ca(03056:03057)/crlf/
      data ca(03058:03099)/'ÉvYiaU`Vb`aa!aaUaa!aaUbV``aa!a\X]Xe^f^!a\X'/
      data ca(03100:03139)/'f^!a\X\Yf]f^!afXeX]^\^!afX\^!afXfY\]\^!!'/
      data ca(03140:03141)/crlf/
      data ca(03142:03177)/'ÖvUnaXaibi!aaXbXbi!aY`j`ja!aY`Yaja!!'/
      data ca(03178:03179)/crlf/
      data ca(03180:03221)/'wv\gcibjaj`i`hagbgchckbm`n!aahaibibhah!abj'/
      data ca(03222:03231)/'ck!acibm!!'/
      data ca(03232:03233)/crlf/
      data ca(03234:03253)/'ÑvUnY`j`ja!aY`Yaja!!'/
      data ca(03254:03255)/crlf/
      data ca(03256:03291)/'vv\gag`h`iajbjcichbgag!aahaibibhah!!'/
      data ca(03292:03293)/crlf/
      data ca(03294:03313)/'ÄvVmjQXqYq!ajQkQYq!!'/
      data ca(03314:03315)/crlf/
      data ca(03316:03357)/'lvWk`U]V[YZ^Za[f]i`jbjeigfhah^gYeVbU`U!a^V'/
      data ca(03358:03399)/'\Y[^[a\f^i!a]h`ibieh!adiffgag^fYdV!aeWbV`V'/
      data ca(03400:03403)/']W!!'/
      data ca(03404:03405)/crlf/
      data ca(03406:03433)/'mvWk]Y_XbUbj!a]Y]Z_YaWajbj!!'/
      data ca(03434:03435)/crlf/
      data ca(03436:03477)/'nvWk[Z[Y\W]V_UcUeVfWgYg[f]d`[j!a[Z\Z\Y]W_V'/
      data ca(03478:03507)/'cVeWfYf[e]c`Zj!a[ihihj!aZjhj!!'/
      data ca(03508:03509)/crlf/
      data ca(03510:03551)/'ovWk\UgU`^!a\U\VfV!afU_^!a`]b]e^g`hchdggei'/
      data ca(03552:03593)/'bj_j\i[hZf[f!a_^b^e_gb!ac^f`gcgdfgci!ageeh'/
      data ca(03594:03609)/'bi_i\h[f!a^i[g!!'/
      data ca(03610:03611)/crlf/
      data ca(03612:03651)/'pvWkdXdjej!aeUej!aeUZeie!adX[e!a[didie!!'/
      data ca(03652:03653)/crlf/
      data ca(03654:03695)/'qvWk\U[^!a]V\]!a\UfUfV!a]VfV!a\]_\b\e]g_hb'/
      data ca(03696:03737)/'hdggeibj_j\i[hZf[f!a[^\^^]b]e^ga!ac]f_gbgd'/
      data ca(03738:03763)/'fgci!ageehbi_i\h[f!a^i[g!!'/
      data ca(03764:03765)/crlf/
      data ca(03766:03807)/'rvWkeVfXgXfVcUaU^V\Y[^[c\g^iajbjeigghdhcg`'/
      data ca(03808:03849)/'e^b]a]^^\`!afWcVaV^W!a_V]Y\^\c]g`i!a\e^hai'/
      data ca(03850:03891)/'biehge!acifggdgcf`c^!agbe_b^a^^_\b!a`^]`\c'/
      data ca(03892:03893)/'!!'/
      data ca(03894:03895)/crlf/
      data ca(03896:03923)/'svWkZUhU^j!aZUZVgV!agU]j^j!!'/
      data ca(03924:03925)/crlf/
      data ca(03926:03967)/'tvWk_U\V[X[Z\\]]_^c_e`fagcgffhci_i\h[f[c\a'/
      data ca(03968:04009)/']`__c^e]f\gZgXfVcU_U!a]V\X\Z]\_]c^e_gahchf'/
      data ca(04010:04051)/'ghficj_j\i[hZfZc[a]__^c]e\fZfXeV!afWcV_V\W'/
      data ca(04052:04065)/'!a[g^i!adigg!!'/
      data ca(04066:04067)/crlf/
      data ca(04068:04109)/'uvWkf_daab`b]a[_Z\Z[[X]V`UaUdVfXg\gaffdiaj'/
      data ca(04110:04151)/'_j\i[g\g]i!af\e_ba!af]d`aa`a]`[]!a_a\_[\[['/
      data ca(04152:04193)/'\X_V!a[Z]W`VaVdWfZ!abVeXf\faefci!adhai_i\h'/
      data ca(04194:04195)/'!!'/
      data ca(04196:04197)/crlf/
      data ca(04198:04239)/'xv\ga\`]`^a_b_c^c]b\a\!aa]a^b^b]a]!aag`h`i'/
      data ca(04240:04265)/'ajbjcichbgag!aahaibibhah!!'/
      data ca(04266:04267)/crlf/
      data ca(04268:04309)/'yv\ga\`]`^a_b_c^c]b\a\!aa]a^b^b]a]!acibjaj'/
      data ca(04310:04351)/'`i`hagbgchckbm`n!aahaibibhah!abjck!acibm!!'/
      data ca(04352:04353)/crlf/
      data ca(04354:04365)/'°rUmiXYaij!!'/
      data ca(04366:04367)/crlf/
      data ca(04368:04403)/'ÜvUnY\j\j]!aY\Y]j]!aYdjdje!aYdYeje!!'/
      data ca(04404:04405)/crlf/
      data ca(04406:04417)/'¢rUmYXiaYj!!'/
      data ca(04418:04419)/crlf/
      data ca(04420:04461)/'{vXk[Z[Y\W]V`UcUfVgWhYh[g]f^d_a`!a[Z\Z\Y]W'/
      data ca(04462:04503)/'`VcVfWgYg[f]d^a_!a\X_V!adVgX!ag\c_!aa_acbc'/
      data ca(04504:04539)/'b_!aag`h`iajbjcichbgag!aahaibibhah!!'/
      data ca(04540:04541)/crlf/
      data ca(04542:04583)/'¡rTof]e[cZ`Z^[]\\_\b]d_ebeddeb!a`Z^\]_]b^d'/
      data ca(04584:04625)/'_e!afZebedgeiekcl`l^k[jYhWfVcU`U]V[WYYX[W^'/
      data ca(04626:04657)/'WaXdYf[h]i`jcjfihhig!agZfbfdge!!'/
      data ca(04658:04659)/crlf/
      data ca(04660:04701)/'•tWkaUYj!aaXZjYj!aaXhjij!aaUij!a\dfd!a[ege'/
      data ca(04702:04703)/'!!'/
      data ca(04704:04705)/crlf/
      data ca(04706:04747)/'¶tWk[U[j!a\V\i!a[UcUfVgWhYh\g^f_c`!a\VcVfW'/
      data ca(04748:04789)/'gYg\f^c_!a\_c_f`gahchfghficj[j!a\`c`fagcgf'/
      data ca(04790:04797)/'fhci\i!!'/
      data ca(04798:04799)/crlf/
      data ca(04800:04841)/'ßtWliZhXfVdU`U^V\X[ZZ]Zb[e\g^i`jdjfihgie!a'/
      data ca(04842:04879)/'iZhZgXfWdV`V^W\Z[][b\e^h`idifhggheie!!'/
      data ca(04880:04881)/crlf/
      data ca(04882:04923)/'®tWl[U[j!a\V\i!a[UbUeVgXhZi]ibheggeibj[j!a'/
      data ca(04924:04949)/'\VbVeWfXgZh]hbgefgehbi\i!!'/
      data ca(04950:04951)/crlf/
      data ca(04952:04993)/'©tXk\U\j!a]V]i!a\UhU!a]VhVhU!a]_c_c`!a]`c`'/
      data ca(04994:05009)/'!a]ihihj!a\jhj!!'/
      data ca(05010:05011)/crlf/
      data ca(05012:05053)/'™tXj\U\j!a]V]j\j!a\UhU!a]VhVhU!a]_c_c`!a]`'/
      data ca(05054:05057)/'c`!!'/
      data ca(05058:05059)/crlf/
      data ca(05060:05101)/'´tWliZhXfVdU`U^V\X[ZZ]Zb[e\g^i`jdjfihgieia'/
      data ca(05102:05143)/'da!aiZhZgXfWdV`V^W]X\Z[][b\e]g^h`idifhgghe'/
      data ca(05144:05151)/'hbdbda!!'/
      data ca(05152:05153)/crlf/
      data ca(05154:05195)/'¨tVlZUZj!aZU[U[jZj!ahUgUgjhj!ahUhj!a[_g_!a'/
      data ca(05196:05201)/'[`g`!!'/
      data ca(05202:05203)/crlf/
      data ca(05204:05223)/'≠t]faUajbj!aaUbUbj!!'/
      data ca(05224:05225)/crlf/
      data ca(05226:05267)/'ÆtYjeUeedhbi`i^h]e\e!aeUfUfeehdibj`j^i]h\e'/
      data ca(05268:05269)/'!!'/
      data ca(05270:05271)/crlf/
      data ca(05272:05313)/'ØtWl[U[j\j!a[U\U\j!aiUhU\a!aiU\b!a_^hjij!a'/
      data ca(05314:05319)/'`^ij!!'/
      data ca(05320:05321)/crlf/
      data ca(05322:05353)/'∞tXi\U\j!a\U]U]i!a]ihihj!a\jhj!!'/
      data ca(05354:05355)/crlf/
      data ca(05356:05397)/'±tUmYUYj!aZZZjYj!aZZaj!aYUag!aiUag!ahZaj!a'/
      data ca(05398:05411)/'hZhjij!aiUij!!'/
      data ca(05412:05413)/crlf/
      data ca(05414:05455)/'≤tVlZUZj!a[X[jZj!a[Xhj!aZUgg!agUgg!agUhUhj'/
      data ca(05456:05457)/'!!'/
      data ca(05458:05459)/crlf/
      data ca(05460:05501)/'≥tVl_U]V[XZZY]YbZe[g]i_jcjeiggheibi]hZgXeV'/
      data ca(05502:05543)/'cU_U!a`V]W[ZZ]Zb[e]h`ibiehgehbh]gZeWbV`V!!'/
      data ca(05544:05545)/crlf/
      data ca(05546:05587)/'¥tWk[U[j!a\V\j[j!a[UdUfVgWhYh\g^f_d`\`!a\V'/
      data ca(05588:05603)/'dVfWgYg\f^d_\_!!'/
      data ca(05604:05605)/crlf/
      data ca(05606:05647)/'µtVl_U]V[XZZY]YbZe[g]i_jcjeiggheibi]hZgXeV'/
      data ca(05648:05689)/'cU_U!a`V]W[ZZ]Zb[e]h`ibiehgehbh]gZeWbV`V!a'/
      data ca(05690:05705)/'bgglhl!abgcghl!!'/
      data ca(05706:05707)/crlf/
      data ca(05708:05749)/'∂tWk[U[j!a\V\j[j!a[UcUfVgWhYh\g^f_c`\`!a\V'/
      data ca(05750:05779)/'cVfWgYg\f^c_\_!aa`gjhj!ab`hj!!'/
      data ca(05780:05781)/crlf/
      data ca(05782:05823)/'∑tWkhXfVcU_U\VZXZZ[\\]^^c`eafbgdggfhci_i]h'/
      data ca(05824:05865)/'\gZg!ahXfXeWcV_V\W[X[Z\\^]c_e`gbhdhgficj_j'/
      data ca(05866:05871)/'\iZg!!'/
      data ca(05872:05873)/crlf/
      data ca(05874:05907)/'∏tYjaVaj!abVbjaj!a[UhUhV!a[U[VhV!!'/
      data ca(05908:05909)/crlf/
      data ca(05910:05951)/'πtVlZUZd[g]i`jbjeigghdhU!aZU[U[d\g]h`ibieh'/
      data ca(05952:05961)/'fggdgUhU!!'/
      data ca(05962:05963)/crlf/
      data ca(05964:05995)/'∫tWkYUaj!aYUZUag!aiUhUag!aiUaj!!'/
      data ca(05996:05997)/crlf/
      data ca(05998:06039)/'ªtTnVU\j!aVUWU\g!aaU\g!aaX\j!aaXfj!aaUfg!a'/
      data ca(06040:06053)/'lUkUfg!alUfj!!'/
      data ca(06054:06055)/crlf/
      data ca(06056:06091)/'ºtWkZUgjhj!aZU[Uhj!ahUgUZj!ahU[jZj!!'/
      data ca(06092:06093)/crlf/
      data ca(06094:06131)/'ΩtXkZUa_ajbj!aZU[Ub_!aiUhUa_!aiUb_bj!!'/
      data ca(06132:06133)/crlf/
      data ca(06134:06175)/'ætWkgUZj!ahU[j!aZUhU!aZUZVgV!a[ihihj!aZjhj'/
      data ca(06176:06177)/'!!'/
      data ca(06178:06179)/crlf/
      data ca(06180:06207)/'èrZh^Q^q!a_Q_q!a^QeQ!a^qeq!!'/
      data ca(06208:06209)/crlf/
      data ca(06210:06219)/'ÑgZhZUhm!!'/
      data ca(06220:06221)/crlf/
      data ca(06222:06249)/'êrZhcQcq!adQdq!a]QdQ!a]qdq!!'/
      data ca(06250:06251)/crlf/
      data ca(06252:06271)/'ßrVlYca^ic!aYca_ic!!'/
      data ca(06272:06273)/crlf/
      data ca(06274:06283)/'bcWkXiji!!'/
      data ca(06284:06285)/crlf/
      data ca(06286:06327)/'|v\gcUaV`X`[a\b\c[cZbYaY`Z!aaZa[b[bZaZ!aaV'/
      data ca(06328:06337)/'`Z!a`XaY!!'/
      data ca(06338:06339)/crlf/
      data ca(06340:06381)/'âuWkf\fjgj!af\g\gj!af_d]b\_\]][_ZbZd[g]i_j'/
      data ca(06382:06415)/'bjdifg!af_b]_]]^\_[b[d\g]h_ibifg!!'/
      data ca(06416:06417)/crlf/
      data ca(06418:06459)/'äuWk[U[j\j!a[U\U\j!a\_^]`\c\e]g_hbhdggeicj'/
      data ca(06460:06493)/'`j^i\g!a\_`]c]e^f_gbgdfgehci`i\g!!'/
      data ca(06494:06495)/crlf/
c 'c' editted
      data ca(06496:06537)/'ãuXjg_e]c\`\^]\_[b[d\g^i`jcjeigh!ag_g`e^c]'/
      data ca(06538:06563)/'`]^^]_\b\d]g^h`iciehfhgh!!'/
      data ca(06564:06565)/crlf/
      data ca(06566:06607)/'åuWkfUfjgj!afUgUgj!af_d]b\_\]][_ZbZd[g]i_j'/
      data ca(06608:06641)/'bjdifg!af_b]_]]^\_[b[d\g]h_ibifg!!'/
      data ca(06642:06643)/crlf/
c 'e' editted.
      data ca(06644:06685)/'çuXj\cgcg`f^e]c\`\^]\_[b[d\g^i`jcjeigh!a\b'/
      data ca(06686:06719)/'fbf`e^c]`]^^]_\b\d]g^h`iciehfhgh!!'/
      data ca(06720:06721)/crlf/
      data ca(06722:06763)/'éu[ifUdUbVaYajbj!afUfVdVbW!acVbYbj!a^\e\e]'/
      data ca(06764:06773)/'!a^\^]e]!!'/
      data ca(06774:06775)/crlf/
      data ca(06776:06817)/'èuWkg\f\fkendobp`p^o]n[n!ag\gkfndpbq_q]p[n'/
      data ca(06818:06859)/'!af_d]b\_\]][_ZbZd[g]i_jbjdifg!af_b]_]]^\_'/
      data ca(06860:06875)/'[b[d\g]h_ibifg!!'/
      data ca(06876:06877)/crlf/
      data ca(06878:06919)/'êuWk[U[j\j!a[U\U\j!a\`_]a\d\f]g`gj!a\`_^a]'/
      data ca(06920:06931)/'c]e^f`fjgj!!'/
      data ca(06932:06933)/crlf/
      data ca(06934:06975)/'ëu]faU`V`WaXbXcWcVbUaU!aaVaWbWbVaV!aa\ajbj'/
      data ca(06976:06985)/'!aa\b\bj!!'/
      data ca(06986:06987)/crlf/
      data ca(06988:07029)/'íu]faU`V`WaXbXcWcVbUaU!aaVaWbWbVaV!aa\aqbq'/
      data ca(07030:07039)/'!aa\b\bq!!'/
      data ca(07040:07041)/crlf/
      data ca(07042:07083)/'ìuWj[U[j\j!a[U\U\j!ag\f\\f!ag\\g!a_cejgj!a'/
      data ca(07084:07089)/'`bgj!!'/
      data ca(07090:07091)/crlf/
      data ca(07092:07111)/'îu]faUajbj!aaUbUbj!!'/
      data ca(07112:07113)/crlf/
      data ca(07114:07155)/'ïuRqV\VjWj!aV\W\Wj!aW`Z]\\_\a]b`bj!aW`Z^\]'/
      data ca(07156:07197)/'^]`^a`ajbj!ab`e]g\j\l]m`mj!ab`e^g]i]k^l`lj'/
      data ca(07198:07201)/'mj!!'/
      data ca(07202:07203)/crlf/
      data ca(07204:07245)/'ñuWk[\[j\j!a[\\\\j!a\`_]a\d\f]g`gj!a\`_^a]'/
      data ca(07246:07257)/'c]e^f`fjgj!!'/
      data ca(07258:07259)/crlf/
      data ca(07260:07301)/'óuXk`\^]\_[b[d\g^i`jcjeigghdhbg_e]c\`\!a`]'/
      data ca(07302:07335)/'^^]_\b\d]g^h`iciehfggdgbf_e^c]`]!!'/
      data ca(07336:07337)/crlf/
      data ca(07338:07379)/'òuWk[\[q\q!a[\\\\q!a\_^]`\c\e]g_hbhdggeicj'/
      data ca(07380:07413)/'`j^i\g!a\_`]c]e^f_gbgdfgehci`i\g!!'/
      data ca(07414:07415)/crlf/
      data ca(07416:07457)/'ôuWkf\fqgq!af\g\gq!af_d]b\_\]][_ZbZd[g]i_j'/
      data ca(07458:07491)/'bjdifg!af_b]_]]^\_[b[d\g]h_ibifg!!'/
      data ca(07492:07493)/crlf/
      data ca(07494:07535)/'öuZh^\^j_j!a^\_\_j!a_b`_b]d\g\!a_b``b^d]g]'/
      data ca(07536:07539)/'g\!!'/
      data ca(07540:07541)/crlf/
      data ca(07542:07583)/'õuYjg_f]c\`\]]\_]a_bddfe!aedfffgei!afhci`i'/
      data ca(07584:07625)/']h!a^i]g\g!ag_f_e]!af^c]`]]^!a^]]_^a!a]`_a'/
      data ca(07626:07645)/'dcfdgfggficj`j]i\g!!'/
      data ca(07646:07647)/crlf/
      data ca(07648:07683)/'úu\gaUajbj!aaUbUbj!a^\e\e]!a^\^]e]!!'/
      data ca(07684:07685)/crlf/
      data ca(07686:07727)/'ùuWk[\[f\i^jajciff!a[\\\\f]h_iaichff!af\fj'/
      data ca(07728:07739)/'gj!af\g\gj!!'/
      data ca(07740:07741)/crlf/
      data ca(07742:07773)/'ûuYi[\aj!a[\\\ah!ag\f\ah!ag\aj!!'/
      data ca(07774:07775)/crlf/
      data ca(07776:07817)/'üuUmX\]j!aX\Y\]g!aa\]g!aa_]j!aa_ej!aa\eg!a'/
      data ca(07818:07831)/'j\i\eg!aj\ej!!'/
      data ca(07832:07833)/crlf/
      data ca(07834:07869)/'†uXj[\fjgj!a[\\\gj!ag\f\[j!ag\\j[j!!'/
      data ca(07870:07871)/crlf/
      data ca(07872:07909)/'°uYi[\aj!a[\\\ah!ag\f\ah]q!ag\aj^q]q!!'/
      data ca(07910:07911)/crlf/
      data ca(07912:07953)/'¢uXje][j!ag\]i!a[\g\!a[\[]e]!a]igigj!a[jgj'/
      data ca(07954:07955)/'!!'/
      data ca(07956:07957)/crlf/
      data ca(07958:07999)/'ërZhcQaR`S_U_W`YaZb\b^``!aaR`T`VaXbYc[c]b_'/
      data ca(08000:08041)/'^abccecgbiaj`l`nap!a`bbdbfah`i_k_m`oapcq!!'/
      data ca(08042:08043)/crlf/
      data ca(08044:08053)/'ïr]eaQaq!!'/
      data ca(08054:08055)/crlf/
      data ca(08056:08097)/'írZh_QaRbScUcWbYaZ`\`^b`!aaRbTbVaX`Y_[_]`_'/
      data ca(08098:08139)/'da`c_e_g`iajblbnap!abb`d`fahbickcmboap_q!!'/
      data ca(08140:08141)/crlf/
      data ca(08142:08183)/'¶rUmXdXbY_[^]^__cbecgcibj`!aXbY`[_]__`cced'/
      data ca(08184:08193)/'gdicj`j^!!'/
      data ca(08194:08195)/crlf/
      data ca(08196:08201)/'ôfaa!!'/
      data ca(08202:08203)/crlf/
      data ca(08204:08209)/'õfYi!!'/
      data ca(08210:08211)/crlf/
      data ca(08212:08231)/'≠t]faUajbj!aaUbUbj!!'/
      data ca(08232:08233)/crlf/
      data ca(08234:08253)/'‹bRp[jVdld!agXl^V^!!'/
      data ca(08254:08255)/crlf/
      data ca(08256:08295)/'tsYi]U]g!ae[em!a]]e[!a]^e\!a]fed!a]gee!!'/
      data ca(08296:08297)/crlf/
      data ca(08298:08339)/'pcUmjViWjXkWkVjUiUhUfVeWcZa`_d^f\iZjXjWiWg'/
      data ca(08340:08377)/'XfZf\g_iajdjfihg!aei_ea`c[eWfV!ad_^c!!'/
      data ca(08378:08379)/crlf/
      data ca(08380:08413)/'ùrTnaX`YaZbYaX!aXaja!aah`iajbiah!!'/
      data ca(08414:08415)/crlf/
      data ca(08416:08457)/'•rUnjfhffedca_`^^]\]Z^Y`YbZd\e^e`dacd_f]h\'/
      data ca(08458:08461)/'j\!!'/
      data ca(08462:08463)/crlf/
      data ca(08464:08479)/'àr]eaU`\!abU`\!!'/
      data ca(08480:08481)/crlf/
      data ca(08482:08517)/'¥rUmiYbY^Z\[Z]Y`YbZe\g^hbiii!aYaea!!'/
      data ca(08518:08519)/crlf/
      data ca(08520:08543)/'ØrTjW\[\ah!aZ\aj!ajQaj!!'/
      data ca(08544:08545)/crlf/
      data ca(08546:08567)/'Øg\fa[ag!a\^fd!af^\d!!'/
      data ca(08568:08569)/crlf/
      data ca(08570:08591)/'ôrUmaYaj!aYaia!aYjij!!'/
      data ca(08592:08593)/crlf/
      data ca(08594:08623)/'≤rUmYY`YdZf[h]i`ibhefgdh`iYi!!'/
      data ca(08624:08625)/crlf/
      data ca(08626:08647)/'örUmaYaj!aYYiY!aYaia!!'/
      data ca(08648:08649)/crlf/
      data ca(08650:08677)/'ms]e```bbbb```!a``bb!ab``b!!'/
      data ca(08678:08679)/crlf/
      data ca(08680:08703)/'ªrPqS\X\ah!aW]aj!aqIaj!!'/
      data ca(08704:08705)/crlf/
      data ca(08706:08745)/'®gZh`Z][[]Z`Zb[e]g`hbheggehbh`g]e[bZ`Z!!'/
      data ca(08746:08747)/crlf/
      data ca(08748:08763)/'©g[g[[[gggg[[[!!'/
      data ca(08764:08765)/crlf/
      data ca(08766:08779)/'™gZhaYZeheaY!!'/
      data ca(08780:08781)/crlf/
      data ca(08782:08823)/'≤g]e`]^^]`]b^d`ebeddebe`d^b]`]!a^`^b!a___c'/
      data ca(08824:08855)/'!a`^`d!aa^ad!ab^bd!ac_cc!ad`db!!'/
      data ca(08856:08857)/crlf/
      data ca(08858:08899)/'≥g]e]]]eeee]]]!a^^^d!a_^_d!a`^`d!aa^ad!ab^'/
      data ca(08900:08915)/'bd!ac^cd!ad^dd!!'/
      data ca(08916:08917)/crlf/
      data ca(08918:08955)/'¥g\fa[\dfda[!aa^^c!aa^dc!aaa`c!aaabc!!'/
      data ca(08956:08957)/crlf/
      data ca(08958:08985)/'¨gYiaX_^Y^^b\hadfhdbi^c^aX!!'/
      data ca(08986:08987)/crlf/
      data ca(08988:09029)/'∏g[ga[]fg_[_efa[!aaaa[!aaa[_!aaa]f!aaaef!a'/
      data ca(09030:09035)/'aag_!!'/
      data ca(09036:09037)/crlf/
      data ca(09038:09053)/'´g[gaW[aakgaaW!!'/
      data ca(09054:09055)/crlf/
      data ca(09056:09081)/'πgahaZah!aaZh]a`!ab\e]b^!!'/
      data ca(09082:09083)/crlf/
      data ca(09084:09123)/'ƒfTnaX`YaZbYaX!aXhWiXjYiXh!ajhiijjkijh!!'/
      data ca(09124:09125)/crlf/
      data ca(09126:09167)/'ærUnkbjdhefeddcc`__^]][]Y^X`XbYd[e]e_d`cc_'/
      data ca(09168:09181)/'d^f]h]j^k`kb!!'/
      data ca(09182:09183)/crlf/
      data ca(09184:09207)/'£rUmiUY\ic!aYeie!aYjij!!'/
      data ca(09208:09209)/crlf/
      data ca(09210:09231)/'†rTnX\j\!aXaja!aXfjf!!'/
      data ca(09232:09233)/crlf/
      data ca(09234:09257)/'§rUmYUi\Yc!aYeie!aYjij!!'/
      data ca(09258:09259)/crlf/
      data ca(09260:09289)/'∞rUmiYbY^Z\[Z]Y`YbZe\g^hbiii!!'/
      data ca(09290:09291)/crlf/
      data ca(09292:09333)/'πrXkgaf^e]c\a\^]\`[c[f\h]i_jajdifggdh_hZgW'/
      data ca(09334:09375)/'fVdUaU_V^W^X_X_W!aa\_]]`\c\g]i!aajciegfdg_'/
      data ca(09376:09383)/'gZfWdU!!'/
      data ca(09384:09385)/crlf/
      data ca(09386:09427)/'}qWkZX[Zgfhhhj!a[[gg!aZXZZ[\ghhj!a__[cZeZg'/
      data ca(09428:09469)/'[iZj!aZe\i!a[c[e\g\iZj!abbg]!aeXe[f]h]h[fZ'/
      data ca(09470:09481)/'eX!aeXf[h]!!'/
      data ca(09482:09483)/crlf/
      data ca(09484:09525)/'∆rYiaU`WaYbWaU!aaUac!aa_`abeag`ebaa_!aacaq'/
      data ca(09526:09567)/'!aam`oaqboam!a[\]]_\][[\!a[\g\!ac\e]g\e[c\'/
      data ca(09568:09599)/'!a[j]k_j]i[j!a[jgj!acjekgjeicj!!'/
      data ca(09600:09601)/crlf/
      data ca(09602:09643)/'ÿpVlZTYY!aiThY!a^]]b!ae]db!aZfYk!aifhk!aZV'/
      data ca(09644:09677)/'hV!aZWhW!a^_d_!a^`d`!aZhhh!aZihi!!'/
      data ca(09678:09679)/crlf/
      data ca(09680:09713)/'ŒpWkaUYj!aaUij!aaXhj!aZihi!aYjij!!'/
      data ca(09714:09715)/crlf/
      data ca(09716:09743)/'«rXkgUgj!aZUgU!a__g_!aZjgj!!'/
      data ca(09744:09745)/crlf/
      data ca(09746:09787)/'ﬂpWlaUaj!abUbj!a_Z\[[\Z^Za[c\d_edegdhciai^'/
      data ca(09788:09829)/'h\g[dZ_Z!a_Z][\\[^[a\c]d_e!adefdgchah^g\f['/
      data ca(09830:09845)/'dZ!a^UeU!a^jej!!'/
      data ca(09846:09847)/crlf/
      data ca(09848:09879)/'ÕpXj]U]j!a^U^j!aZUiUi[hU!aZjaj!!'/
      data ca(09880:09881)/crlf/
      data ca(09882:09921)/'ssYi^V^m!adUdl!a\]f[!a\^f\!a\ffd!a\gfe!!'/
      data ca(09922:09923)/crlf/
      data ca(09924:09965)/'ﬂtXhfWdZb_`d_f]i[j!ah[f]c^`^^]][]Y^W`VdUhU'/
      data ca(09966:09993)/'fWeYc_ae`g^i[jYjXiXgYfZgYh!!'/
      data ca(09994:09995)/crlf/
      data ca(09996:10037)/'ºrUmjRiSjTkSkRjQhQfRdTcVbYa]_i^m]o!aeSdUcY'/
      data ca(10038:10063)/'ae`i_l^n\pZqXqWpWoXnYoXp!!'/
      data ca(10064:10065)/crlf/
      data ca(10066:10099)/'∫rWkYUaj!aZUah!aiUaj!aYUiU!aZVhV!!'/
      data ca(10100:10101)/crlf/
      data ca(10102:10135)/'’pWkaUZj!aaUhj!aaXgj!aXj^j!adjjj!!'/
      data ca(10136:10137)/crlf/
      data ca(10138:10179)/'ec[h^_a_b_e_!aa_aS!a_UaSbRb_!aha\a!a^f]f]d'/
      data ca(10180:10221)/'^d`cccdeffehcibj^k]m]o]p!a^n_nbpepeofnfm!a'/
      data ca(10222:10235)/'eobo!adhefdd!!'/
      data ca(10236:10237)/crlf/
      data ca(10238:10267)/'≥rUmYiYbZ^[\]Z`YbYeZg\h^ibii!!'/
      data ca(10268:10269)/crlf/
      data ca(10270:10301)/'ärZh`U^V]X]Z^\`]b]d\eZeXdVbU`U!!'/
      data ca(10302:10303)/crlf/
      data ca(10304:10345)/'⁄pUmZUZj!a[U[j!agUgj!ahUhj!aWUkU!aWj^j!adj'/
      data ca(10346:10349)/'kj!!'/
      data ca(10350:10351)/crlf/
      data ca(10352:10393)/'“pVl`U]V[XZZY^YaZe[g]i`jbjeiggheiai^hZgXeV'/
      data ca(10394:10435)/'bU`U!a`U^V\X[ZZ^Za[e\g^i`j!abjdifggehah^gZ'/
      data ca(10436:10467)/'fXdVbU!a^\^c!ad\dc!a^_d_!a^`d`!!'/
      data ca(10468:10469)/crlf/
      data ca(10470:10511)/'huUmdVcWbY`^^d]f[iYj!acWbZ`b_e^g\iYjWjViVg'/
      data ca(10512:10553)/'WfXgWh!a^[]]\^Z^Y]Y[ZY\W^VaUfUiVjXjZi\h]e^'/
      data ca(10554:10595)/'a^!afUhViXiZh\g]e^!aa^d_eafhgj!aa^c_daehgj'/
      data ca(10596:10603)/'hjjilg!!'/
      data ca(10604:10605)/crlf/
      data ca(10606:10647)/'‹pWlZUa_Yj!aYU`_!aYUhUi[gU!aZigi!aYjhjidgj'/
      data ca(10648:10649)/'!!'/
      data ca(10650:10651)/crlf/
      data ca(10652:10681)/'±rUmYYY`Zd[f]h`ibiehgfhdi`iY!!'/
      data ca(10682:10683)/crlf/
      data ca(10684:10725)/'ﬁpXkZZZX[V\U^U_V`Xa\aj!aZX\V^V`X!aiZiXhVgU'/
      data ca(10726:10753)/'eUdVcXb\bj!aiXgVeVcX!a^jej!!'/
      data ca(10754:10755)/crlf/
      data ca(10756:10797)/'≈rYiaU`WaYbWaU!aaUaq!aa``caqbca`!a[\]]_\]['/
      data ca(10798:10819)/'[\!a[\g\!ac\e]g\e[c\!!'/
      data ca(10820:10821)/crlf/
      data ca(10822:10863)/'bqVlYgZj^j\fZbY_Y[ZX\V_UcUfVhXi[i_hbffdjhj'/
      data ca(10864:10905)/'ig!a\f[cZ_Z[[X]V_U!acUeVgXh[h_gcff!aZi]i!a'/
      data ca(10906:10911)/'eihi!!'/
      data ca(10912:10913)/crlf/
      data ca(10914:10929)/'õrVlZZhh!ahZZh!!'/
      data ca(10930:10931)/crlf/
      data ca(10932:10973)/'aqVmaUaj!abUbj!aX\Y[[\\`]b^c`d!aY[Z\[`\b]c'/
      data ca(10974:11015)/'`dcdfcgbh`i\j[!acdecfbg`h\j[k\!a^UeU!a^jej'/
      data ca(11016:11017)/'!!'/
      data ca(11018:11019)/crlf/
      data ca(11020:11061)/'ΩrUmjRiSjTkSkRjQhQfRdTcVbYa]_i^m]o!aeSdUcY'/
      data ca(11062:11103)/'ae`i_l^n\pZqXqWpWoXnYoXp!a`Z][[]Z`Zb[e]g`h'/
      data ca(11104:11123)/'bheggehbh`g]e[bZ`Z!!'/
      data ca(11124:11125)/crlf/
      data ca(11126:11157)/'≈sXj\:\a\à!a]:]a]à!a\:g:!a\àgà!!'/
      data ca(11158:11159)/crlf/
      data ca(11160:11181)/'ürTnhXZj!aX^j^!aXdjd!!'/
      data ca(11182:11183)/crlf/
      data ca(11184:11215)/'∆sXje:eaeà!af:fafà!a[:f:!a[àfà!!'/
      data ca(11216:11217)/crlf/
      data ca(11218:11243)/'∑rTn[_Xa[c!a^\Ya^f!aYaja!!'/
      data ca(11244:11245)/crlf/
      data ca(11246:11271)/'µrTng_jagc!ad\iadf!aXaia!!'/
      data ca(11272:11273)/crlf/
      data ca(11274:11289)/'¬fUma\aj!aXjjj!!'/
      data ca(11290:11291)/crlf/
      data ca(11292:11333)/'ØqVm`\]][_ZaYdYgZi]j_jaidffch_i\!a`\^]\_[a'/
      data ca(11334:11375)/'ZdZg[i]j!a`\b\d]e_gghiij!ab\c]d_fggiijjj!!'/
      data ca(11376:11377)/crlf/
      data ca(11378:11419)/'∞qVkcU`V^X\\[_ZcYiXq!acUaV_X]\\_[cZiYq!acU'/
      data ca(11420:11461)/'eUgVhWhZg\f]c^_^!aeUgWgZf\e]c^!a_^c_eafcff'/
      data ca(11462:11495)/'ehdiaj_j]i\h[e!a_^b_daecefdhciaj!!'/
      data ca(11496:11497)/crlf/
      data ca(11498:11539)/'ºqXicUaV`W`XaYdZgZ!adZ`[^\]^]`_bbcec!adZa['/
      data ca(11540:11581)/'_\^^^``bbc!abc^d\e[g[i]kbmcncpaq_q!abc_d]e'/
      data ca(11582:11591)/'\g\i^kbm!!'/
      data ca(11592:11593)/crlf/
      data ca(11594:11635)/'≤qXke]c\a\^]\`[c[f\h]i_jajdiffgcg`f^bYaWaU'/
      data ca(11636:11677)/'bTdTfUhW!aa\_]]`\c\g]i!aajcieffcf_e]cZbXbV'/
      data ca(11678:11685)/'cUeUhW!!'/
      data ca(11686:11687)/crlf/
      data ca(11688:11729)/'≥qXjg_e]c\_\]]]__abb!a_\^]^_`abb!abb]c[e[g'/
      data ca(11730:11755)/'\i_jbjdifg!abb^c\e\g]i_j!!'/
      data ca(11756:11757)/crlf/
      data ca(11758:11799)/'jrVldU^q!aeU]q!a`\\]Z_YbYeZg\i_jbjfihgidia'/
      data ca(11800:11841)/'h_f]c\`\!a`\]][_ZbZe[g]i_j!abjeigghdhag_e]'/
      data ca(11842:11845)/'c\!!'/
      data ca(11846:11847)/crlf/
      data ca(11848:11889)/'±qWkX_Z]\\^\`]a^babeai^q!aY^[]_]a^!ai\h_ga'/
      data ca(11890:11907)/'bh_m]q!ah\g_fabh!!'/
      data ca(11908:11909)/crlf/
      data ca(11910:11951)/'µqVlW`X^Z\]\^]^_]c[j!a\\]]]_\cZj!a]c__a]c\'/
      data ca(11952:11977)/'e\g]h^hagfdq!ae\g^gaffcq!!'/
      data ca(11978:11979)/crlf/
      data ca(11980:12013)/'∑q[ga\_c^g^i_jbjdhef!ab\`c_g_i`j!!'/
      data ca(12014:12015)/crlf/
      data ca(12016:12057)/'ÃqZiiVhWiXjWjViUgUeVdWcYb\_j^n]p!agUeWdYc]'/
      data ca(12058:12089)/'af`j_m^o]p[qYqXpXoYnZoYp!a^\h\!!'/
      data ca(12090:12091)/crlf/
      data ca(12092:12133)/'∏qWk]\Yj!a^\Zj!ag\h]i]h\f\d]`a^b\b!a^b`cbi'/
      data ca(12134:12153)/'cj!a^b_caibjdjfihf!!'/
      data ca(12154:12155)/crlf/
      data ca(12156:12197)/'πqWkZU\U^V_W`Yfggihj!a\U^W_Yegfihjij!aa\Yj'/
      data ca(12198:12205)/'!aa\Zj!!'/
      data ca(12206:12207)/crlf/
      data ca(12208:12249)/'∫qUl\\Vq!a]\Wq!a\_[e[h]j_jaicged!ag\dgdiej'/
      data ca(12250:12267)/'hjjhkf!ah\egeifj!!'/
      data ca(12268:12269)/crlf/
      data ca(12270:12311)/'ªqWk]\[j!a^\]b\g[j!ah\g`ed!ai\h_gaedcf`h^i'/
      data ca(12312:12321)/'[j!aZ\^\!!'/
      data ca(12322:12323)/crlf/
      data ca(12324:12365)/'hrXif]d\a\^]\_[b[e\h]i`jcjei!aa\_]]_\b\e]h'/
      data ca(12366:12377)/'^i`j!a\cdc!!'/
      data ca(12378:12379)/crlf/
      data ca(12380:12421)/'æqVl_][j!a_]\j!ae]ej!ae]fj!aX_Z]]\j\!aX_Z^'/
      data ca(12422:12427)/']]j]!!'/
      data ca(12428:12429)/crlf/
      data ca(12430:12471)/'irXkcU`V^Y][\^[c[g\i^j`jcieffdgah\hXgVeUcU'/
      data ca(12472:12513)/'!acUaV_Y^[]^\c\g]i^j!a`jbidfedfag\gXfVeU!a'/
      data ca(12514:12519)/']_f_!!'/
      data ca(12520:12521)/crlf/
      data ca(12522:12563)/'øqWj[e\h]i_jajdiffgcg`f^e]c\a\^]\`[cWq!aaj'/
      data ca(12564:12587)/'cieffcf_e]!aa\_]]`\cXq!!'/
      data ca(12588:12589)/crlf/
      data ca(12590:12631)/'¿qWlj\`\]][`ZcZf[h\i^j`jcieffcf`e^d]b\!a`\'/
      data ca(12632:12663)/'^]\`[c[g\i!a`jbidfece_d]!ad]j]!!'/
      data ca(12664:12665)/crlf/
      data ca(12666:12701)/'¡qWkb]_j!ab]`j!aY_[]^\i\!aY_[^^]i]!!'/
      data ca(12702:12703)/crlf/
      data ca(12704:12745)/'¬qWkX`Y^[\^\_]__]e]h_j!a]\^]^_\e\h]i_j`jci'/
      data ca(12746:12769)/'eggdhah^g\f]g^ha!agdh^!!'/
      data ca(12770:12771)/crlf/
      data ca(12772:12787)/'“bWkY_aei_adY_!!'/
      data ca(12788:12789)/crlf/
      data ca(12790:12831)/'∆qUlY`[^^]]\[]Y`XcXfYiZj\j^i`fac!aXfYhZi\i'/
      data ca(12832:12873)/'^h`f!a`c`faibjdjfihfici`h]g\f]h^i`!a`fahbi'/
      data ca(12874:12881)/'difhhf!!'/
      data ca(12882:12883)/crlf/
      data ca(12884:12925)/'ƒqXjZ\\\^]__dnepfq!a\\]]^_cndpfqhq!ai\h^fa'/
      data ca(12926:12933)/'\lZoYq!!'/
      data ca(12934:12935)/crlf/
      data ca(12936:12977)/'≈qUldU^q!aeU]q!aV`W^Y\\\]]]_\d\g^iaichfehb'/
      data ca(12978:13007)/'!a[\\]\_[d[g\i^jajcieggdhbj\!!'/
      data ca(13008:13009)/crlf/
      data ca(13010:13051)/'¥qXjcUaV`W`XaYdZiZiYfZb\_^\a[d[f\h_jblcncp'/
      data ca(13052:13075)/'bq`q_p!ad[`^]a\d\f]h_j!!'/
      data ca(13076:13077)/crlf/
      data ca(13078:13103)/'∏rYi_gajcg!a\daifd!aaXai!!'/
      data ca(13104:13105)/crlf/
      data ca(13106:13121)/'¡fZh^Q^j!adQdj!!'/
      data ca(13122:13123)/crlf/
      data ca(13124:13149)/'∂rYi_[aXc[!a\^aYf^!aaYaj!!'/
      data ca(13150:13151)/crlf/
      data ca(13152:13193)/'ÿbUmXbY`Y_[^]^__cbecgcibj`!aicgdedcc_`]_[_'/
      data ca(13194:13235)/'Y`!aXhYfYe[d]d_echeigiihjf!aiigjejci_f]e[e'/
      data ca(13236:13239)/'Yf!!'/
      data ca(13240:13241)/crlf/
      data ca(13242:13247)/'ôfaa!!'/
      data ca(13248:13249)/crlf/
      data ca(13250:13255)/'ùvYi!!'/
      data ca(13256:13257)/crlf/
      data ca(13258:13295)/'¨v\gdUcVab!adVab!adUeVab!a_h^i_j`i_h!!'/
      data ca(13296:13297)/crlf/
      data ca(13298:13313)/'πv]fdUb\!aeUb\!!'/
      data ca(13314:13315)/crlf/
      data ca(13316:13347)/'ªvZibU`V_X_Z`\b]d]f\gZgXfVdUbU!!'/
      data ca(13348:13349)/crlf/
      data ca(13350:13391)/'±vWlcQ[n!ahQ`n!aiYhZi[jZjYiWhVeUaU^V\X\Z]\'/
      data ca(13392:13433)/'^]eagc!a\Z^\e`fagcgffheibj^j[iZhYfYeZd[eZf'/
      data ca(13434:13435)/'!!'/
      data ca(13436:13437)/crlf/
      data ca(13438:13479)/'ørUmjUXj!a]U_W_Y^[\\Z\XZXXYV[U]U_VbWeWhVjU'/
      data ca(13480:13505)/'!afcddcfchejgjiijgjehcfc!!'/
      data ca(13506:13507)/crlf/
      data ca(13508:13549)/'∞vTnk]j^k_l^l]k\j\h]f_ag_i]jZjWiVgVeWcXbZa'/
      data ca(13550:13591)/'__a^c\dZdXcVaU_V^X^[_a`dbgdifjhjihig!aZjXi'/
      data ca(13592:13621)/'WgWeXcYb__!a^[_``cbfdhfihiih!!'/
      data ca(13622:13623)/crlf/
      data ca(13624:13643)/'Øv\gdWcVdUeVeWdYb[!!'/
      data ca(13644:13645)/crlf/
      data ca(13646:13687)/'≥vZiiQeTbW`Z^^]c]g^l_o`q!aeTbX`\__^d^i_n`q'/
      data ca(13688:13689)/'!!'/
      data ca(13690:13691)/crlf/
      data ca(13692:13733)/'¥vYhbQcSdVe[e_ddbh`k]nYq!abQcTdYd^ccbf`j]n'/
      data ca(13734:13735)/'!!'/
      data ca(13736:13737)/crlf/
      data ca(13738:13759)/'µvYjcUca!a^Xh^!ahX^^!!'/
      data ca(13760:13761)/crlf/
      data ca(13762:13797)/'ÖvUnaXaibi!aaXbXbi!aY`j`ja!aY`Yaja!!'/
      data ca(13798:13799)/crlf/
      data ca(13800:13819)/'©v\g_j^i_h`i`j_l]n!!'/
      data ca(13820:13821)/crlf/
      data ca(13822:13841)/'ÑvUnY`j`ja!aY`Yaja!!'/
      data ca(13842:13843)/crlf/
      data ca(13844:13859)/'®v\g_h^i_j`i_h!!'/
      data ca(13860:13861)/crlf/
      data ca(13862:13881)/'ÄvVmjQXqYq!ajQkQYq!!'/
      data ca(13882:13883)/crlf/
      data ca(13884:13925)/'ûvWlcU`V^X\[[^ZbZe[h\i^j`jcieggdhai]iZhWgV'/
      data ca(13926:13967)/'eUcU!acUaV_X][\^[b[e\h^j!a`jbidgfdgah]hZgW'/
      data ca(13968:13971)/'eU!!'/
      data ca(13972:13973)/crlf/
      data ca(13974:14007)/'üvWlcY^j!aeU_j!aeUbX_Z][!adX`Z][!!'/
      data ca(14008:14009)/crlf/
      data ca(14010:14051)/'†vWl^Y_Z^[]Z]Y^W_VbUeUhViXiZh\f^c`_b\dZfXj'/
      data ca(14052:14093)/'!aeUgVhXhZg\e^_b!aYhZg\gaidifhgf!a\gajdjfi'/
      data ca(14094:14097)/'gf!!'/
      data ca(14098:14099)/crlf/
      data ca(14100:14141)/'°vWl^Y_Z^[]Z]Y^W_VbUeUhViXiZh\e^b_!aeUgVhX'/
      data ca(14142:14183)/'hZg\e^!a`_b_e`fagcgffheibj^j[iZhYfYeZd[eZf'/
      data ca(14184:14203)/'!ab_d`eafcffehdibj!!'/
      data ca(14204:14205)/crlf/
      data ca(14206:14229)/'¢vWlgVaj!ahUbj!ahUYdid!!'/
      data ca(14230:14231)/crlf/
      data ca(14232:14273)/'£vWl`U[_!a`UjU!a`VeVjU!a[_\^_]b]e^f_gagdfg'/
      data ca(14274:14313)/'diaj^j[iZhYfYeZd[eZf!ab]d^e_fafdegciaj!!'/
      data ca(14314:14315)/crlf/
      data ca(14316:14357)/'§vWlhXgYhZiYiXhVfUcU`V^X\[[^ZbZf[h\i^jajdi'/
      data ca(14358:14399)/'fggegbf`e_c^`^^_\a[c!acUaV_X][\^[b[g\i!aaj'/
      data ca(14400:14411)/'ciegfefae_!!'/
      data ca(14412:14413)/crlf/
      data ca(14414:14455)/'•vWl]U[[!ajUiXg[ba`d_f^j!ag[aa_d^f]j!a\X_U'/
      data ca(14456:14477)/'aUfX!a]W_VaVfXhXiWjU!!'/
      data ca(14478:14479)/crlf/
      data ca(14480:14521)/'¶vWlbU_V^W]Y]\^^`_c_g^h]i[iXhVeUbU!abU`V_W'/
      data ca(14522:14563)/'^Y^\_^`_!ac_f^g]h[hXgVeU!a`_\`ZbYdYgZi]jaj'/
      data ca(14564:14605)/'eifhgfgcfae`c_!a`_]`[bZdZg[i]j!aajdiehfffb'/
      data ca(14606:14609)/'e`!!'/
      data ca(14610:14611)/crlf/
      data ca(14612:14653)/'ßvWlh\g^e`ca`a^`]_\]\Z]X_VbUeUgVhWiYi]hagd'/
      data ca(14654:14695)/'egci`j]j[iZgZf[e\f[g!a^`]^]Z^X`VbU!agVhXh]'/
      data ca(14696:14707)/'gafddgbi`j!!'/
      data ca(14708:14709)/crlf/
      data ca(14710:14735)/'™v\gb\a]b^c]b\!a_h^i_j`i!!'/
      data ca(14736:14737)/crlf/
      data ca(14738:14769)/'´v\gb\a]b^c]b\!a_j^i_h`i`j_l]n!!'/
      data ca(14770:14771)/crlf/
      data ca(14772:14789)/'°rUmiXYaijiiZaiW!!'/
      data ca(14790:14791)/crlf/
      data ca(14792:14827)/'ÜvUnY\j\j]!aY\Y]j]!aYdjdje!aYdYeje!!'/
      data ca(14828:14829)/crlf/
      data ca(14830:14847)/'¢rUmYXiaYjYkhaYW!!'/
      data ca(14848:14849)/crlf/
      data ca(14850:14891)/'≠vWl^Y_Z^[]Z]Y^W_VbUfUiVjXjZi\h]b_```baccc'/
      data ca(14892:14921)/'!afUhViXiZh\g]e^!a_h^i_j`i_h!!'/
      data ca(14922:14923)/crlf/
      data ca(14924:14965)/'¡rTof]e[cZ`Z^[]\\_\b]d_ebeddeb!a`Z^\]_]b^d'/
      data ca(14966:15007)/'_e!afZebedgeiekcl`l^k[jYhWfVcU`U]V[WYYX[W^'/
      data ca(15008:15039)/'WaXdYf[h]i`jcjfihhig!agZfbfdge!!'/
      data ca(15040:15041)/crlf/
      data ca(15042:15081)/'cqWkdUWj!adUej!acWdj!a[ddd!aUj[j!aajgj!!'/
      data ca(15082:15083)/crlf/
      data ca(15084:15125)/'dqUm^UXj!a_UYj!a[UfUiVjXjZi]h^e_!afUhViXiZ'/
      data ca(15126:15167)/'h]g^e_!a\_e_g`hbhdggeiajUj!ae_f`gbgdfgdiaj'/
      data ca(15168:15169)/'!!'/
      data ca(15170:15171)/crlf/
      data ca(15172:15213)/'eqWliWjWkUj[jYiWhVfUcU`V^X\[[^ZbZe[h\i_jbj'/
      data ca(15214:15243)/'difgge!acUaV_X][\^[b[e\h]i_j!!'/
      data ca(15244:15245)/crlf/
      data ca(15246:15287)/'fqUl^UXj!a_UYj!a[UdUgVhWiZi^hbffdhbi^jUj!a'/
      data ca(15288:15309)/'dUfVgWhZh^gbefchai^j!!'/
      data ca(15310:15311)/crlf/
      data ca(15312:15353)/'gqUl^UXj!a_UYj!ac[ac!a[UjUi[iU!a\_b_!aUjdj'/
      data ca(15354:15359)/'fecj!!'/
      data ca(15360:15361)/crlf/
      data ca(15362:15403)/'hqUk^UXj!a_UYj!ac[ac!a[UjUi[iU!a\_b_!aUj\j'/
      data ca(15404:15405)/'!!'/
      data ca(15406:15407)/crlf/
      data ca(15408:15449)/'iqWmiWjWkUj[jYiWhVfUcU`V^X\[[^ZbZe[h\i_jaj'/
      data ca(15450:15491)/'difghc!acUaV_X][\^[b[e\h]i_j!aajcieggc!adc'/
      data ca(15492:15495)/'kc!!'/
      data ca(15496:15497)/crlf/
      data ca(15498:15539)/'jqTn]UWj!a^UXj!ajUdj!akUej!aZUaU!agUnU!a[_'/
      data ca(15540:15555)/'g_!aTj[j!aajhj!!'/
      data ca(15556:15557)/crlf/
      data ca(15558:15585)/'kq[hdU^j!aeU_j!aaUhU!a[jbj!!'/
      data ca(15586:15587)/crlf/
      data ca(15588:15629)/'lqXjgUbfah`i^j\jZiYgYeZd[eZf!afUaf`h^j!acU'/
      data ca(15630:15633)/'jU!!'/
      data ca(15634:15635)/crlf/
      data ca(15636:15677)/'mqUl^UXj!a_UYj!alU[b!ab^fj!aa^ej!a[UbU!ahU'/
      data ca(15678:15693)/'nU!aUj\j!abjhj!!'/
      data ca(15694:15695)/crlf/
      data ca(15696:15727)/'nqWk`UZj!aaU[j!a]UdU!aWjfjhdej!!'/
      data ca(15728:15729)/crlf/
      data ca(15730:15771)/'oqTo]UWj!a]U^j!a^U_h!akU^j!akUej!alUfj!aZU'/
      data ca(15772:15793)/'^U!akUoU!aTjZj!abjij!!'/
      data ca(15794:15795)/crlf/
      data ca(15796:15837)/'pqUn^UXj!a^Ueg!a^Xej!akUej!a[U^U!ahUnU!aUj'/
      data ca(15838:15841)/'[j!!'/
      data ca(15842:15843)/crlf/
      data ca(15844:15885)/'qqVlbU_V]X[[Z^YbYeZh[i]j`jcieggdhai]iZhWgV'/
      data ca(15886:15927)/'eUbU!abU`V^X\[[^ZbZe[h]j!a`jbidgfdgah]hZgW'/
      data ca(15928:15931)/'eU!!'/
      data ca(15932:15933)/crlf/
      data ca(15934:15975)/'rqUl^UXj!a_UYj!a[UgUjVkXkZj]h_d`\`!agUiVjX'/
      data ca(15976:15991)/'jZi]g_d`!aUj\j!!'/
      data ca(15992:15993)/crlf/
      data ca(15994:16035)/'sqVlbU_V]X[[Z^YbYeZh[i]j`jcieggdhai]iZhWgV'/
      data ca(16036:16077)/'eUbU!abU`V^X\[[^ZbZe[h]j!a`jbidgfdgah]hZgW'/
      data ca(16078:16119)/'eU!a[h[g\e^d_daebgbncoeofmfl!abgcmdnenfm!!'/
      data ca(16120:16121)/crlf/
      data ca(16122:16163)/'tqUm^UXj!a_UYj!a[UfUiVjXjZi]h^e_\_!afUhViX'/
      data ca(16164:16205)/'iZh]g^e_!aa_c`daeifjhjihig!adafhgihiih!aUj'/
      data ca(16206:16209)/'\j!!'/
      data ca(16210:16211)/crlf/
      data ca(16212:16253)/'uqVmiWjWkUj[jYiWhVeUaU^V\X\Z]\^]eagc!a\Z^\'/
      data ca(16254:16285)/'e`fagcgffheibj^j[iZhYfYdXjYhZh!!'/
      data ca(16286:16287)/crlf/
      data ca(16288:16323)/'vqWldU^j!aeU_j!a^U[[]UlUk[kU!a[jbj!!'/
      data ca(16324:16325)/crlf/
      data ca(16326:16367)/'wqUn]UZ`YdYgZi]jajdifggdkU!a^U[`ZdZg[i]j!a'/
      data ca(16368:16379)/'ZUaU!ahUnU!!'/
      data ca(16380:16381)/crlf/
      data ca(16382:16415)/'xqWk]U^j!a^U_h!akU^j!a[UaU!agUmU!!'/
      data ca(16416:16417)/crlf/
      data ca(16418:16459)/'yqTn\UZj!a]U[h!adUZj!adUbj!aeUch!alUbj!aYU'/
      data ca(16460:16469)/'`U!aiUoU!!'/
      data ca(16470:16471)/crlf/
      data ca(16472:16513)/'zqVl]Udj!a^Uej!akUWj!a[UaU!agUmU!aUj[j!aaj'/
      data ca(16514:16517)/'gj!!'/
      data ca(16518:16519)/crlf/
      data ca(16520:16561)/'{qWl]Ua_^j!a^Ub__j!alUb_!a[UaU!ahUnU!a[jbj'/
      data ca(16562:16563)/'!!'/
      data ca(16564:16565)/crlf/
      data ca(16566:16601)/'|qVljUWj!akUXj!a^U[[]UkU!aWjejgddj!!'/
      data ca(16602:16603)/crlf/
      data ca(16604:16631)/'èrZh^Q^q!a_Q_q!a^QeQ!a^qeq!!'/
      data ca(16632:16633)/crlf/
      data ca(16634:16643)/'ÑgZhZUhm!!'/
      data ca(16644:16645)/crlf/
      data ca(16646:16673)/'êrZhcQcq!adQdq!a]QdQ!a]qdq!!'/
      data ca(16674:16675)/crlf/
      data ca(16676:16695)/'ßrVlYca^ic!aYca_ic!!'/
      data ca(16696:16697)/crlf/
      data ca(16698:16707)/'bcWkXiji!!'/
      data ca(16708:16709)/crlf/
      data ca(16710:16729)/'Æv\geUcWbYbZc[dZcY!!'/
      data ca(16730:16731)/crlf/
      data ca(16732:16773)/'«qWlg\ecdgdiejhjjhkf!ah\fcegeifj!aece`d]b\'/
      data ca(16774:16813)/'`\]][`ZcZf[h\i^j`jbidfec!a`\^]\`[c[g\i!!'/
      data ca(16814:16815)/crlf/
      data ca(16816:16857)/'»qWj_U[b[e\h]i!a`U\b!a\b]__]a\c\e]f^g`gcff'/
      data ca(16858:16891)/'diaj_j]i\f\b!ae]f_fcefciaj!a\U`U!!'/
      data ca(16892:16893)/crlf/
      data ca(16894:16935)/'…qXjf_f`g`g_f]d\a\^]\`[c[f\h]i_jajdiff!aa\'/
      data ca(16936:16947)/'_]]`\c\g]i!!'/
      data ca(16948:16949)/crlf/
      data ca(16950:16991)/' qWliUecdgdiejhjjhkf!ajUfcegeifj!aece`d]b\'/
      data ca(16992:17033)/'`\]][`ZcZf[h\i^j`jbidfec!a`\^]\`[c[g\i!afU'/
      data ca(17034:17037)/'jU!!'/
      data ca(17038:17039)/crlf/
      data ca(17040:17081)/'ÀqXj\e`dccfag_f]d\a\^]\`[c[f\h]i_jajdifg!a'/
      data ca(17082:17095)/'a\_]]`\c\g]i!!'/
      data ca(17096:17097)/crlf/
      data ca(17098:17139)/'ÃqZiiVhWiXjWjViUgUeVdWcYb\_j^n]p!agUeWdYc]'/
      data ca(17140:17171)/'af`j_m^o]p[qYqXpXoYnZoYp!a^\h\!!'/
      data ca(17172:17173)/crlf/
      data ca(17174:17215)/'ÕqWkh\djcmap^q[qYpXoXnYmZnYo!ag\cjbm`p^q!a'/
      data ca(17216:17257)/'ece`d]b\`\]][`ZcZf[h\i^j`jbidfec!a`\^]\`[c'/
      data ca(17258:17263)/'[g\i!!'/
      data ca(17264:17265)/crlf/
      data ca(17266:17307)/'ŒqWl_UYj!a`UZj!a\c^_`]b\d\f]g^g`efeifj!ad\'/
      data ca(17308:17331)/'f^f`dfdiejhjjhkf!a\U`U!!'/
      data ca(17332:17333)/crlf/
      data ca(17334:17375)/'œq[hdUcVdWeVdU!a\`]^_\b\c]c`afaibj!aa\b]b`'/
      data ca(17376:17389)/'`f`iajdjfhgf!!'/
      data ca(17390:17391)/crlf/
      data ca(17392:17433)/'–q[heUdVeWfVeU!a]`^^`\c\d]d`aj`m_o^p\qZqYp'/
      data ca(17434:17459)/'YoZn[oZp!ab\c]c``j_m^o\q!!'/
      data ca(17460:17461)/crlf/
      data ca(17462:17503)/'—qWk_UYj!a`UZj!ag]f^g_h^h]g\f\d]`a^b\b!a^b'/
      data ca(17504:17533)/'`cbicj!a^b_caibjdjfihf!a\U`U!!'/
      data ca(17534:17535)/crlf/
      data ca(17536:17575)/'“q\hdU`c_g_i`jcjehff!aeUac`g`iaj!aaUeU!!'/
      data ca(17576:17577)/crlf/
      data ca(17578:17619)/'”qPqQ`R^T\W\X]X_WcUj!aV\W]W_VcTj!aWcY_[]]\'/
      data ca(17620:17661)/'_\a]b^b`_j!a_\a^a`^j!aacc_e]g\i\k]l^l`jfji'/
      data ca(17662:17685)/'kj!ai\k^k`ifiijjmjohpf!!'/
      data ca(17686:17687)/crlf/
      data ca(17688:17729)/'‘qUlV`W^Y\\\]]]_\cZj!a[\\]\_[cYj!a\c^_`]b\'/
      data ca(17730:17765)/'d\f]g^g`efeifj!ad\f^f`dfdiejhjjhkf!!'/
      data ca(17766:17767)/crlf/
      data ca(17768:17809)/'’qXja\^]\`[c[f\h]i_jajdiffgcg`f^e]c\a\!aa\'/
      data ca(17810:17835)/'_]]`\c\g]i!aajcieffcf_e]!!'/
      data ca(17836:17837)/crlf/
      data ca(17838:17879)/'÷qVkW`X^Z\]\^]^_]cYq!a\\]]]_\cXq!a]c^``]b\'/
      data ca(17880:17921)/'d\f]g^h`hcgfeibj`j^i]f]c!af]g_gcffdibj!aUq'/
      data ca(17922:17925)/'\q!!'/
      data ca(17926:17927)/crlf/
      data ca(17928:17969)/'◊qWkg\aq!ah\bq!aece`d]b\`\]][`ZcZf[h\i^j`j'/
      data ca(17970:17997)/'bidfec!a`\^]\`[c[g\i!a^qeq!!'/
      data ca(17998:17999)/crlf/
      data ca(18000:18041)/'ÿqXiY`Z^\\_\`]`__c]j!a^\_]__^c\j!a_ca_c]e\'/
      data ca(18042:18055)/'g\h]h^g_f^g]!!'/
      data ca(18056:18057)/crlf/
      data ca(18058:18099)/'ŸqYjg^g_h_h^g]d\a\^]]^]`^aeeff!a]_^`edfefh'/
      data ca(18100:18117)/'eibj_j\i[h[g\g\h!!'/
      data ca(18118:18119)/crlf/
      data ca(18120:18159)/'⁄qZhcU_c^g^i_jbjdhef!adU`c_g_i`j!a]\f\!!'/
      data ca(18160:18161)/crlf/
      data ca(18162:18203)/'€qUlV`W^Y\\\]]]`[f[h]j!a[\\]\`ZfZh[i]j_jai'/
      data ca(18204:18239)/'cgec!ag\ecdgdiejhjjhkf!ah\fcegeifj!!'/
      data ca(18240:18241)/crlf/
      data ca(18242:18283)/'‹qWkX`Y^[\^\_]_`]f]h_j!a]\^]^`\f\h]i_j`jci'/
      data ca(18284:18297)/'eggdh`h\g\h^!!'/
      data ca(18298:18299)/crlf/
      data ca(18300:18341)/'›qRoS`T^V\Y\Z]Z`XfXhZj!aX\Y]Y`WfWhXiZj\j^i'/
      data ca(18342:18383)/'`gae!ac\aeahbidjfjhijgkelal\k\l^!ad\bebhdj'/
      data ca(18384:18385)/'!!'/
      data ca(18386:18387)/crlf/
      data ca(18388:18429)/'ﬁqWkZ`\]^\a\b^ba!a`\a^aa`e_g]i[jZjYiYhZg[h'/
      data ca(18430:18471)/'Zi!a`e`hajdjfihf!ah]g^h_i^i]h\g\e]c_baaeah'/
      data ca(18472:18475)/'bj!!'/
      data ca(18476:18477)/crlf/
      data ca(18478:18519)/'ﬂqVkW`X^Z\]\^]^`\f\h^j!a\\]]]`[f[h\i^j`jbi'/
      data ca(18520:18561)/'dgfc!ai\ejdmbp_q\qZpYoYnZm[nZo!ah\djcmap_q'/
      data ca(18562:18563)/'!!'/
      data ca(18564:18565)/crlf/
      data ca(18566:18607)/'‡qWkh\g^e`]f[hZj!a[`\^^\a\e^!a\^^]a]e^g^!a'/
      data ca(18608:18631)/'[h]haidifh!a]hajdjfhgf!!'/
      data ca(18632:18633)/crlf/
      data ca(18634:18675)/'ërZhcQaR`S_U_W`YaZb\b^``!aaR`T`VaXbYc[c]b_'/
      data ca(18676:18717)/'^abccecgbiaj`l`nap!a`bbdbfah`i_k_m`oapcq!!'/
      data ca(18718:18719)/crlf/
      data ca(18720:18729)/'ïr]eaQaq!!'/
      data ca(18730:18731)/crlf/
      data ca(18732:18773)/'írZh_QaRbScUcWbYaZ`\`^b`!aaRbTbVaX`Y_[_]`_'/
      data ca(18774:18815)/'da`c_e_g`iajblbnap!abb`d`fahbickcmboap_q!!'/
      data ca(18816:18817)/crlf/
      data ca(18818:18859)/'¶rUmXdXbY_[^]^__cbecgcibj`!aXbY`[_]__`cced'/
      data ca(18860:18869)/'gdicj`j^!!'/
      data ca(18870:18871)/crlf/
      data ca(18872:18877)/'ôfaa!!'/
      data ca(18878:18879)/crlf/
       end
