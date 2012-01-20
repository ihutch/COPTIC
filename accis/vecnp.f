c*********************************************************************
c        Draw a vector to the plotter file, normalized coords   */
      subroutine vecnp(nx,ny,ud)
      real nx,ny
      integer ud
      include 'plotcom.h'
c      integer updown in common set to 99 by blockdata or color.
      integer pldx,pldy
      integer plx,ply
      common /plcoord/ plx,ply

      integer ipscale
      parameter (ipscale=10)
c Bigger scale is supposed to prevent rounding problems. See 1016 below.
      character*4 pu,pd, spc, endpair
      character*17 postlude
      character*300 prelude
      integer nr,nu,nd,ns,ne,no
      common /pltfc/pu,pd,spc,endpair,postlude,prelude
      common /pltfi/nr,nu,nd,ns,ne,no

c Difference form, requires M and D to be rmoveto, rdrawto: relative.
c Caused rounding erros in gs in some rare cases.
c      pldx=int(ipscale*(25+1000*nx)) - plx
c      pldy=int(ipscale*(25+1000*ny)) - ply
c      plx=plx+pldx
c      ply=ply+pldy
c Absolute form
      pldx=int(ipscale*(25+1000*nx))
      pldy=int(ipscale*(25+1000*ny))
      plx=pldx
      ply=pldy
c postscript section
      if(abs(pfsw).eq.2.or.abs(pfsw).eq.3)then
	 if(updown.ne.ud)then
	    if(ud.eq.0)call abufwrt(endpair,ne,12)
	    updown=ud
	 endif
	 call ibufwrt(pldx,12)
	 call abufwrt(spc,ns,12)
	 call ibufwrt(pldy,12)
	 call abufwrt(spc,ns,12)
	 if(ud.ne.0)then
	    call abufwrt(pd,nd,12)
c	    updown=ud
	 else
	    call abufwrt(pu,nd,12)
	 endif
      else
	 if(ud.ne.updown)then
	    updown=ud
	    if(ud.ne.0)then
	       call abufwrt(pd,nd,12)
	    else
	       call abufwrt(pu,nu,12)
	    endif
	 else 
	    call abufwrt(endpair,ne,12)
	 endif
	 call ibufwrt(pldx,12)
	 call abufwrt(spc,ns,12)
	 call ibufwrt(pldy,12)
      endif
      return
      end

c*********************************************************************/
      subroutine inib(iunit,filen)
      integer iunit
      character*(*) filen
      integer sblen
      include 'plotcom.h'

      character*80 sbuf
      common /wbuf/sblen,sbuf
      integer plx,ply
      common /plcoord/ plx,ply

      character*4 pu,pd, spc,endpair
      character*17 postlude
      character*300 prelude
      integer nr,nu,nd,ns,ne,no
      common /pltfc/pu,pd,spc,endpair,postlude,prelude
      common /pltfi/nr,nu,nd,ns,ne,no
      character*1745 ca
      character*2 crlf
      parameter (crlf=' \')
c See lnswrt below for the end-of line character definition here.
c g77 is incapable of handling non-printable characters properly.
      save ca

      data ca(00001:00023)/'%!PS-Adobe-3.0 EPSF-2.0'/
      data ca(00024:00025)/crlf/
      data ca(00026:00053)/'%%BoundingBox: 28 28 575 762'/
      data ca(00054:00055)/crlf/
      data ca(00056:00097)/'/SP0  {12 setlinewidth 1.000 1.000 1.000 s'/
      data ca(00098:00113)/'etrgbcolor } def'/
      data ca(00114:00115)/crlf/
      data ca(00116:00157)/'/SP1  {12 setlinewidth 0.000 0.000 0.545 s'/
      data ca(00158:00173)/'etrgbcolor } def'/
      data ca(00174:00175)/crlf/
      data ca(00176:00217)/'/SP2  {12 setlinewidth 0.180 0.545 0.341 s'/
      data ca(00218:00233)/'etrgbcolor } def'/
      data ca(00234:00235)/crlf/
      data ca(00236:00277)/'/SP3  {12 setlinewidth 0.125 0.698 0.667 s'/
      data ca(00278:00293)/'etrgbcolor } def'/
      data ca(00294:00295)/crlf/
      data ca(00296:00337)/'/SP4  {12 setlinewidth 0.698 0.133 0.133 s'/
      data ca(00338:00353)/'etrgbcolor } def'/
      data ca(00354:00355)/crlf/
      data ca(00356:00397)/'/SP5  {12 setlinewidth 0.580 0.000 0.827 s'/
      data ca(00398:00413)/'etrgbcolor } def'/
      data ca(00414:00415)/crlf/
      data ca(00416:00457)/'/SP6  {12 setlinewidth 0.957 0.643 0.376 s'/
      data ca(00458:00473)/'etrgbcolor } def'/
      data ca(00474:00475)/crlf/
      data ca(00476:00517)/'/SP7  {12 setlinewidth 0.439 0.502 0.565 s'/
      data ca(00518:00533)/'etrgbcolor } def'/
      data ca(00534:00535)/crlf/
      data ca(00536:00577)/'/SP8  {12 setlinewidth 0.745 0.745 0.745 s'/
      data ca(00578:00593)/'etrgbcolor } def'/
      data ca(00594:00595)/crlf/
      data ca(00596:00637)/'/SP9  {12 setlinewidth 0.000 0.000 1.000 s'/
      data ca(00638:00653)/'etrgbcolor } def'/
      data ca(00654:00655)/crlf/
      data ca(00656:00697)/'/SP10 {12 setlinewidth 0.000 1.000 0.000 s'/
      data ca(00698:00713)/'etrgbcolor } def'/
      data ca(00714:00715)/crlf/
      data ca(00716:00757)/'/SP11 {12 setlinewidth 0.000 1.000 1.000 s'/
      data ca(00758:00773)/'etrgbcolor } def'/
      data ca(00774:00775)/crlf/
      data ca(00776:00817)/'/SP12 {12 setlinewidth 1.000 0.000 0.000 s'/
      data ca(00818:00833)/'etrgbcolor } def'/
      data ca(00834:00835)/crlf/
      data ca(00836:00877)/'/SP13 {12 setlinewidth 1.000 0.000 1.000 s'/
      data ca(00878:00893)/'etrgbcolor } def'/
      data ca(00894:00895)/crlf/
      data ca(00896:00937)/'/SP14 {12 setlinewidth 1.000 1.000 0.000 s'/
      data ca(00938:00953)/'etrgbcolor } def'/
      data ca(00954:00955)/crlf/
      data ca(00956:00976)/'/D    {  lineto } def'/
      data ca(00977:00978)/crlf/
      data ca(00979:00999)/'/M    {  moveto } def'/
      data ca(01000:01001)/crlf/
      data ca(01002:01021)/'/MT   { moveto } def'/
      data ca(01022:01023)/crlf/
      data ca(01024:01062)/'/ST   { currentpoint stroke moveto} def'/
      data ca(01063:01064)/crlf/
      data ca(01065:01106)/'/LS   { exch dup 3 1 roll mul exch 8 div ['/
      data ca(01107:01122)/' exch dup } def '/
      data ca(01123:01124)/crlf/
      data ca(01125:01150)/'/LM   { mul exch dup } def'/
      data ca(01151:01152)/crlf/
      data ca(01153:01194)/'/LE   { mul exch pop ] exch setdash } def '/
      data ca(01195:01196)/crlf/
      data ca(01197:01224)/'/LT   { [ ] 0 setdash } def '/
      data ca(01225:01226)/crlf/
      data ca(01227:01252)/'/LT1  { LS 0 LM 8 LE } def'/
      data ca(01253:01254)/crlf/
      data ca(01255:01280)/'/LT2  { LS 4 LM 4 LE } def'/
      data ca(01281:01282)/crlf/
      data ca(01283:01308)/'/LT3  { LS 6 LM 2 LE } def'/
      data ca(01309:01310)/crlf/
      data ca(01311:01346)/'/LT4  { LS 6 LM 1 LM 0 LM 1 LE } def'/
      data ca(01347:01348)/crlf/
      data ca(01349:01384)/'/LT5  { LS 5 LM 1 LM 1 LM 1 LE } def'/
      data ca(01385:01386)/crlf/
      data ca(01387:01428)/'/LT6  { LS 3 LM 1 LM 1 LM 1 LM 1 LM 1 LE }'/
      data ca(01429:01432)/' def'/
      data ca(01433:01434)/crlf/
      data ca(01435:01476)/'/SP15 {12 setlinewidth 0.000 0.000 0.000 s'/
      data ca(01477:01492)/'etrgbcolor } def'/
      data ca(01493:01494)/crlf/
      data ca(01495:01533)/'/SF where { pop }  { /SF 1 def } ifelse'/
      data ca(01534:01535)/crlf/
      data ca(01536:01577)/'/sleep { SF mul cvi ST {gsave currentpoint'/
      data ca(01578:01579)/crlf/
      data ca(01580:01617)/'SP0 0 0 moveto 20 0 lineto 0 20 lineto'/
      data ca(01618:01654)/' fill moveto grestore }  repeat } def'/
      data ca(01655:01656)/crlf/
      data ca(01657:01698)/'gsave 1 setlinecap 1 setlinejoin 72 1016  '/
      data ca(01699:01729)/'div dup scale 426 426 translate'/
      data ca(01730:01731)/crlf/
      data ca(01732:01743)/'SP15 0 0 MT '/
      data ca(01744:01745)/crlf/

      
c The following ensures we can open the file for writing. If not then 
c an error condition exists.
      open(unit=iunit,FILE=filen,status='unknown',err=901)
      close(unit=iunit,status='delete',err=901)
      open(unit=iunit,FILE=filen,status='new',err=901)
      sblen=1
      if(abs(pfsw).eq.1)then
	 pu=';PU'
	 nu=3
	 pd=';PD'
	 nd=3
	 spc=','
	 ns=1
	 endpair=','
	 ne=1
	 postlude=';PG'//char(26)
	 no=4
c	 prelude= ' IN;SP1;PU0,0;PR;'   Caused rotation problem on hplj4.
	 prelude= '   ;SP1;PU0,0;PR;'
c        move plot pen to (0,0), switch to relative plotting.
	 nr=17
	 call abufwrt(prelude,nr,iunit)
         updown=99
      else
	 pu='M '
	 nu=2
	 pd='D '
	 nd=2
	 spc=' '
	 ns=1
	 endpair='ST '
	 ne=3
	 if(abs(pfsw).eq.3)then
	    ca(00026:00053)='%%BoundingBox: 30 30 620 480'
	    postlude='grestore'
	    no=8
	 elseif(abs(pfsw).eq.2)then
	    ca(00026:00053)='%%BoundingBox: 28 28 575 762'
	    postlude='showpage grestore'
	    no=17
	 endif
	 call lnswrt(iunit,ca,1745,'\',2)
	 if(abs(pfsw).eq.2)then
c Landscape, postscript.
	    write(iunit,*)' 8000  0 translate 90 rotate 0 0 MT'
            write(iunit,'(a)')'%%EndProlog'
            write(iunit,'(a)')'%%PageOrientation: Landscape'
	 elseif(abs(pfsw).eq.3)then
           write(iunit,*)' 0.8 dup scale'
         endif
         updown=0
      endif
      plx=0
      ply=0
      return
 901  write(*,*)'**** Accis Error! Cannot open file ',filen
      write(*,*)'Graphics files will not be written.'
      pfsw=0
      return
      end
c*********************************************************************/
      subroutine flushb(iunit)
      integer iunit
      integer sblen
      character*80 sbuf
      common /wbuf/sblen,sbuf

      character*4 pu,pd, spc, endpair
      character*17 postlude
      character*300 prelude
      integer nr,nu,nd,ns,ne,no
      common /pltfc/pu,pd,spc,endpair,postlude,prelude
      common /pltfi/nr,nu,nd,ns,ne,no

      if(sblen.gt.1)write(iunit,*)sbuf(1:sblen-1)
      write(iunit,*)postlude
      endfile(iunit)
      close(iunit)
      end
c********************************************************************
c Fill the path just drawn. A path is a set of vectors all drawn with
c pen down one after the other. The vecfill is dummy except for X11.
      subroutine pathfill()
      include 'plotcom.h'
c      write(*,*)'Called pathfill'
      if(pfsw.ge.0) call vecfill()
      if(abs(pfsw) .eq. 2 .or. abs(pfsw) .eq. 3)
     $     call abufwrt(' false upath ufill ',19,12)
      end

C********************************************************************
c*******************************************************************
      subroutine gradcolor(i)
      integer i
      include 'plotcom.h'
      character*4 pu,pd, spc, endpair
      character*17 postlude
      character*300 prelude
      integer nr,nu,nd,ns,ne,no
      common /pltfc/pu,pd,spc,endpair,postlude,prelude
      common /pltfi/nr,nu,nd,ns,ne,no

      character*60 string 
c      write(*,*)'gradcolor',i
      if(pfsw.ge.0) call acgradcolor(i)
      if(abs(pfsw).eq.3 .or. abs(pfsw).eq.2)then
c         write(*,*)'Calling getrgbcolor',i,ired,igreen,iblue
         call getrgbcolor(i,ired,igreen,iblue)
c         write(*,*)i,ired,igreen,iblue
         call abufwrt(spc,ns,12)
         red=ired/65535.
         green=igreen/65535.
         blue=iblue/65535.
c Stroke the previous path, if necessary.
         if(updown.ne.0) call abufwrt(endpair,ne,12)
c         if(red.gt.1.or.green.gt.1.or.blue.gt.1)then
c            write(*,*)i,red,green,blue,ired,igreen,iblue
c         endif
         write(string,'(3f8.4,'' setrgbcolor '')')
     $        red,green,blue
         call abufwrt(string,istlen(string,59)+1,12)
c Prevent the stroking of this path after the color change.         
         updown=0
      endif
      end
c*******************************************************************
c************************************************************************

      subroutine asleep(iusec)
c Generalized usleep routine including putting delay into ps files.
c Using repetitive drawing for timing.
      integer iusec
      include 'plotcom.h'
c Calibration
      isleep=.4*iusec
      if(abs(pfsw).eq.2.or.abs(pfsw).eq.3)then
         call ibufwrt(isleep,12)
         call abufwrt(' sleep ',7,12)
      endif
      call accisusleep(iusec)
      end

c***********************************************************************
      subroutine accisusleep(ius)
      integer ius
      include 'plotcom.h'
c Calibration 
      imax=ius*.25
      do i=1,imax
         call vec(0,0,0)
         call vec(1,0,1)
         call vec(0,1,1)
         call vecfill()
      enddo
      call vecn(crsrx,crsry,0)
      end
c*********************************************************************
c*********************************************************************
      integer function npgradtri(x,y,z,h,i3d)
c Draw triangle color gradient to the postscript file, 
c using level 3 command shfill.
c Returns 1 for success. 
c i3d indicates 3d x,y,z if 1, 2d x,y if 0. color in h. */
      real x(3),y(3),z(3),h(3)

      include 'plotcom.h'
      include 'world3.h'
      integer pldx,pldy
      integer plx,ply
      common /plcoord/ plx,ply

      integer ipscale
      parameter (ipscale=10)
c Bigger scale is supposed to prevent rounding problems. See 1016 below.
      character*4 pu,pd, spc, endpair
      character*17 postlude
      character*300 prelude
      integer nr,nu,nd,ns,ne,no
      common /pltfc/pu,pd,spc,endpair,postlude,prelude
      common /pltfi/nr,nu,nd,ns,ne,no

      character*60 string 
      real xn,yn,zn,xs,ys,zs

c If a dummy routine avoiding level-3 postscript, return 0.
c      npgradtri=0
c      return


      call nlwrt(12)
      call abufwrt('<< /ShadingType 4 /ColorSpace /DeviceRGB ',41,12)
      call abufwrt('/DataSource [ ',13,12)
c      call abufwrt('% (edge-flag x- y-positions',27,12)
c      call abufwrt(' color-triplet)x3',17,12)

      do i=1,3
         if(i3d.eq.1)then
            call wxyz2nxyz(x(i),y(i),z(i),xn,yn,zn)
            call trn32(xn,yn,zn,xs,ys,zs,0)
            xn=xs+xcbc2
            yn=ys+ycbc2
         else
            xn=wx2nx(x(i))
            yn=wy2ny(y(i))
            call tn2s(xn,yn,sx,sy)
         endif
         pldx=int(ipscale*(25+1000*xn))
         pldy=int(ipscale*(25+1000*yn))
         plx=pldx
         ply=pldy

         li=h(i)+1
         call getrgbcolor(li,ired,igreen,iblue)
         call abufwrt(spc,ns,12)
         red=ired/65535.
         green=igreen/65535.
         blue=iblue/65535.
         
         call nlwrt(12)
         call abufwrt('  0  ',5,12)
         call ibufwrt(pldx,12)
         call abufwrt(spc,ns,12)
         call ibufwrt(pldy,12) 
         call abufwrt(spc,ns,12)
         write(string,'(3f8.4,'' '')')red,green,blue
         call abufwrt(string,istlen(string,59)+1,12)

      enddo

      call abufwrt('] >> shfill ',12,12)
      call nlwrt(12)

      npgradtri=1

      end
