c*********************************************************************
c Routines for writing drawing commands to output 'plotter' files.
c The only format-dependent actions we need are
c 1. inib Initialize the file and some of the strings.
c 2. vecnp draw a line with pen up or down.
c 3. pathfill fill (non-destructively) the path just drawn.
c 4. gradcolor set a color from the gradient.
c 5. asleep sleep before executing next graphics commands
c 6. npgradtri attempt to draw a triangle gradient.
c 7. zerolinewidth set linewidth zero for color contouring.
c 8. npcolor write the 16-color change command.
c 9. flushb flush remaining writes, postlude, and close file.
c*********************************************************************
c Postscript (2,3) and PCL (1) are drawn to unit 12.
c pgf (4,5) to unit 13.
c*********************************************************************
c        Draw a vector to the plotter file, normalized coords   */
c On input ud=0 indicates moveto/pen-up; ud=1 lineto/pen-down; 
c ud=2 lift pen but don't draw; ud=3 stroke (don't move) and disable
c the stroke at the start of the next move by putting updown=99.
c ud=-1 draw point only: equivalent to a call with ud=0, then with ud=1.
      subroutine vecnp(nx,ny,ud)
      real nx,ny
      integer ud
      include 'plotcom.h'
c      integer updown in common set to 99 by blockdata or color.
      integer pldx,pldy
      integer plx,ply
      common /plcoord/ plx,ply

      integer ipsunit,ipgunit
      parameter (ipsunit=12,ipgunit=13)
      integer ipscale
      parameter (ipscale=10)
c Bigger scale is supposed to prevent rounding problems. See 1016 below.
      include 'npcom.f'
c Absolute form
      pldx=int(ipscale*(25+1000*nx))
      pldy=int(ipscale*(25+1000*ny))
      plx=pldx
      ply=pldy
      if(abs(pfsw).eq.2.or.abs(pfsw).eq.3)then
c Postscript section
         if(updown.ne.ud)then
c Not sure in PS.
            if(ud.ne.1.and.updown.ne.99)call abufwrt(endpair,ne,ipsunit)
c           if(ud.ne.1)call abufwrt(endpair,ne,ipsunit)
            updown=ud
         endif
c ud=2 signals pen-up but don't draw this position.
         if(ud.eq.2)return
         if(ud.eq.3)then
            updown=99
            return
         endif
         call ibufwrt(pldx,ipsunit)
         call abufwrt(spc,ns,ipsunit)
         call ibufwrt(pldy,ipsunit)
         call abufwrt(spc,ns,ipsunit)
         if(ud.eq.1)then
            call abufwrt(pd,nd,ipsunit)
         elseif(ud.eq.-1)then
c Point output: move to, then draw to same place, lift pen [and end].
            call abufwrt(pu,nd,ipsunit)
            call ibufwrt(pldx,ipsunit)
            call abufwrt(spc,ns,ipsunit)
            call ibufwrt(pldy,ipsunit)
            call abufwrt(spc,ns,ipsunit)
            call abufwrt(pd,nd,ipsunit)
            updown=0
c Apparently not needed:
c            call abufwrt(endpair,ne,ipsunit)
         else
            call abufwrt(pu,nd,ipsunit)
         endif
      elseif(abs(pfsw).eq.4.or.abs(pfsw).eq.5)then
c pgfgraphics
         if(updown.ne.ud)then
            if(ud.ne.1 .and. updown.ne.99)write(ipgunit,*)only(1:ns)
     $           ,'\pgfusepath{stroke}'
            updown=ud
         endif
         if(ud.eq.2)return
         if(ud.eq.3)updown=99
         if(ud.eq.1)then
            write(ipgunit,*)only(1:ns),'\pgfpathlineto{\pgfpointxy{'
     $           ,pldx,'}{',pldy ,'}}'
         elseif(ud.eq.0)then
            write(ipgunit,*)only(1:ns),'\pgfpathmoveto{\pgfpointxy{'
     $           ,pldx,'}{',pldy,'}}'
         endif
      else
c PCL
         if(ud.ne.updown)then
            updown=ud
            if(ud.eq.1)then
               call abufwrt(pd,nd,ipsunit)
            elseif(ud.eq.0)then
               call abufwrt(pu,nu,ipsunit)
            else
               return
            endif
         else 
            call abufwrt(endpair,ne,ipsunit)
         endif
         call ibufwrt(pldx,ipsunit)
         call abufwrt(spc,ns,ipsunit)
         call ibufwrt(pldy,ipsunit)
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

      include 'npcom.f'

      itunit=iunit+1
c      write(*,*)'Starting inib',iunit,'  ',filen
c The following ensures we can open the file for writing. If not then 
c an error condition exists.
      open(unit=iunit,FILE=filen,status='unknown',err=901)
      close(unit=iunit,status='delete',err=901)
      open(unit=iunit,FILE=filen,status='new',err=901)
      sblen=1
      if(abs(pfsw).eq.1)then
c PCL
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
         prelude= '   ;SP1;PU0,0;PR;'
c        move plot pen to (0,0), switch to relative plotting.
         nr=17
         call abufwrt(prelude,nr,iunit)
         updown=99
      elseif(abs(pfsw).eq.2 .or.abs(pfsw).eq.3)then
c PS
         pu='M '
         nu=2
         pd='D '
         nd=2
         spc=' '
         ns=1
         endpair='ST '
         ne=3
         if(abs(pfsw).eq.3)then
            postlude='grestore'
            no=8
         elseif(abs(pfsw).eq.2)then
            postlude='showpage grestore'
            no=17
         endif
         call psheader(iunit,pfsw)
         if(abs(pfsw).eq.2)then
c Landscape, PS'.
            write(iunit,*)' 8000  0 translate 90 rotate 0 0 MT'
            write(iunit,'(a)')'%%EndProlog'
            write(iunit,'(a)')'%%PageOrientation: Landscape'
         elseif(abs(pfsw).eq.3)then
            write(iunit,*)' 0.8 dup scale'
         endif
         updown=0
      elseif(abs(pfsw).eq.4)then
c PGF. 
c Define postlude etc.
         postlude='\end{pgfpicture}'
c Write out the header. Start with possible warning about \pgffillsave
 8       format(9a)
         write(iunit,8
     $        )'\ifx\undefined\pgffillsave\gdef\pgffillsave{\message{**'
     $        ,'**** You must'
         write(iunit,8
     $        )'    uncomment the pgfillsave code outside a beamer fra'
     $        ,'me ****** }}\fi'
         write(iunit,8)'%\makeatletter'
         write(iunit,8
     $        )'%\def\pgffillsave{\pgfsyssoftpath@getcurrentpath\mysof'
     $        ,'tpath\pgfusepath{fill}'
         write(iunit,8)'%\pgfsyssoftpath@setcurrentpath\mysoftpath}'
         write(iunit,8)'%\makeatother'
         write(iunit,*)'\begin{pgfpicture}','\pgfsetroundjoin'
c Set the pgfxy-scaling so that unit normal, which is 10000 plot units
c is 100 mm.
         write(iunit,*)'\pgfsetxvec{\pgfpoint{0.01mm}{0pt}}'
         write(iunit,*)'\pgfsetyvec{\pgfpoint{0pt}{.01mm}}'
c Use the integer ne  for frame counter
         ne=1
         only=''
         ns=lntrim(only)
      elseif(abs(pfsw).eq.5)then
c PGF animateinline with timeline
         ilast=lntrim(filen)
         filen(ilast:ilast)='t'
c Open timeline file for writing
         open(unit=itunit,FILE=filen,status='unknown',err=901)
         close(unit=itunit,status='delete',err=901)
         open(unit=itunit,FILE=filen,status='new',err=901)
         ne=0
         write(itunit,*)':: 0x0'
c Write header         
         write(iunit,8
     $        )'\ifx\undefined\pgffillsave\gdef\pgffillsave{\message{**'
     $        ,'**** You must'
         write(iunit,8
     $        )'    uncomment the pgfillsave code outside a beamer fra'
     $        ,'me ****** }}\fi'
         write(iunit,8)'%\makeatletter'
         write(iunit,8
     $        )'%\def\pgffillsave{\pgfsyssoftpath@getcurrentpath\mysof'
     $        ,'tpath\pgfusepath{fill}'
         write(iunit,8)'%\pgfsyssoftpath@setcurrentpath\mysoftpath}'
         write(iunit,8)'%\makeatother'
         write(iunit,8)'\begin{animateinline}['
         write(iunit,8)'  timeline=',filen(1:lntrim(filen)),','
         write(iunit,8)'  begin={'
         write(iunit,8)'   \begin{pgfpicture}'
         write(iunit,8)'    \pgfsetxvec{\pgfpoint{0.01mm}{0pt}}'
         write(iunit,8)'    \pgfsetyvec{\pgfpoint{0pt}{.01mm}}'
         write(iunit,8)'    \pgfpathrectangle{\pgfpointorigin}',
     $        '{\pgfpoint{10cm}{7.5cm}}'
         write(iunit,8)'    \pgfusepath{use as bounding box}'
     $        ,'\pgfsetroundjoin'
c An alternative might be \pgfsetmiterlimit{3}
         write(iunit,8)'  },'
         write(iunit,8)'  end={\end{pgfpicture}},'
         write(iunit,8)'  loop,autoplay'
         write(iunit,8)']{5}'
c Define postlude etc.
         postlude='\end{animateinline}'
         only=''
         ns=lntrim(only)
      endif
      plx=0
      ply=0
c      write(*,*)'Finished inib',iunit
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

      include 'npcom.f'

c      write(*,*)'Inside flushb',iunit
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
      integer ipsunit,ipgunit
      parameter (ipsunit=12,ipgunit=13)
      include 'npcom.f'

      if(pfsw.ge.0) call vecfill()
      if(abs(pfsw).eq.4.or. abs(pfsw).eq.5)then
c PGF
         write(ipgunit,*)only(1:ns),'\pgffillsave'
      elseif(abs(pfsw) .eq. 2 .or. abs(pfsw) .eq. 3)then
c This is level 2 postscript:
c     $     call abufwrt(' false upath ufill ',19,ipsunit)
c Level 1 is probably more compatible:
         call abufwrt(' gsave fill grestore ',21,ipsunit)
      endif
      end

C********************************************************************
c*******************************************************************
      subroutine gradcolor(i)
      integer i
      include 'plotcom.h'
      include 'npcom.f'
      integer ipsunit,ipgunit
      parameter (ipsunit=12,ipgunit=13)

      character*60 string 
c      write(*,*)'gradcolor',i
      ncolor=i+15
      if(pfsw.ge.0) call acgradcolor(i)
      if(abs(pfsw).ne.0)then
         call getrgbcolor(i,ired,igreen,iblue)
         call abufwrt(spc,ns,ipsunit)
         red=ired/65535.
         green=igreen/65535.
         blue=iblue/65535.
         if(abs(pfsw).eq.3 .or. abs(pfsw).eq.2)then
c Stroke the previous path, if necessary.
            if(updown.ne.0) call abufwrt(endpair,ne,ipsunit)
            write(string,'(3f8.4,'' setrgbcolor '')')
     $           red,green,blue
            call abufwrt(string,istlen(string,59)+1,ipsunit)
         elseif(abs(pfsw).eq.4.or.abs(pfsw).eq.5)then
c PGF
            if(updown.ne.0)write(ipgunit,*)only(1:ns)
     $           ,'\pgfusepath{stroke}'
 101        format(' \color[rgb]{',f5.3,',',f5.3,',',f5.3,'}')
            write(ipgunit,101)red,green,blue
         endif
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
      integer ipsunit,ipgunit
      parameter (ipsunit=12,ipgunit=13)
c common for ne integer
      include 'npcom.f'
c Lift the pen
      call vecn(crsrx,crsry,3)
c Calibration
      isleep=.4*iusec
      if(abs(pfsw).eq.2.or.abs(pfsw).eq.3)then
         call ibufwrt(isleep,ipsunit)
         call abufwrt(' sleep ',7,ipsunit)
      elseif(abs(pfsw).eq.4)then
c PGF. Increment the frame and update the onlystring.
         ne=ne+1
c         write(only,'(''\only<'',i,''->'')')ne
         write(only,*)'\only<',ne,'->'
         ns=lntrim(only)
      elseif(abs(pfsw).eq.5)then
         ne=ne+1
         write(ipgunit,*)'\newframe'
         write(ipgunit+1,*)'::',ne,'x0'
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

      integer ipsunit
      parameter (ipsunit=12)
      integer ipscale
      parameter (ipscale=10)
c Bigger scale is supposed to prevent rounding problems. See 1016 below.
      include 'npcom.f'
      character*60 string 
      real xn,yn,zn,xs,ys,zs

      if(abs(pfsw).eq.4.or.abs(pfsw).eq.5)then
c PGF not yet implemented.
c If a dummy routine avoiding level-3 postscript, return 0.
         npgradtri=0
         return
      endif

      call nlwrt(ipsunit)
      call abufwrt('<< /ShadingType 4 /ColorSpace /DeviceRGB '
     $     ,41,ipsunit)
      call abufwrt('/DataSource [ ',13,ipsunit)
c      call abufwrt('% (edge-flag x- y-positions',27,ipsunit)
c      call abufwrt(' color-triplet)x3',17,ipsunit)

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
         call abufwrt(spc,ns,ipsunit)
         red=ired/65535.
         green=igreen/65535.
         blue=iblue/65535.
         
         call nlwrt(ipsunit)
         call abufwrt('  0  ',5,ipsunit)
         call ibufwrt(pldx,ipsunit)
         call abufwrt(spc,ns,ipsunit)
         call ibufwrt(pldy,ipsunit) 
         call abufwrt(spc,ns,ipsunit)
         write(string,'(3f8.4,'' '')')red,green,blue
         call abufwrt(string,istlen(string,59)+1,ipsunit)

      enddo

      call abufwrt('] >> shfill ',12,ipsunit)
      call nlwrt(ipsunit)

      npgradtri=1

      end
c*********************************************************************
      subroutine zerolinewidth()
      include 'plotcom.h'
      integer ipsunit
      parameter (ipsunit=12)
c Make line width minimal for contouring
       if(abs(pfsw).eq.2 .or. abs(pfsw).eq.3)
     $     call abufwrt(' ST 0 setlinewidth ',19,ipsunit)
c Not yet implemented for PGF.
       end
c********************************************************************
c Write the color change commands to the plot file.
      subroutine npcolor(li)
      integer li
      integer wid
      character*6 spchr
      include 'plotcom.h'
      include 'npcom.f'
      integer ipsunit,ipgunit
      parameter (ipsunit=12,ipgunit=13)
      character*16 colorname(16)
      data colorname/'white','blue!50!black','green!50!black'
     $     ,'cyan!60!black','red!68!black','violet!90!white'
     $     ,'brown!90!white','gray','lightgray','blue','green','cyan'
     $     ,'red','magenta','yellow','black'/

c      write(*,*)'Color changing',li,pfsw,' ',colorname(li+1)
      if(pfsw.ne.0) then
         updown=99
         if(abs(pfsw).eq.2.or.abs(pfsw).eq.3)
     $        call abufwrt(' ST',3,ipsunit)
         if(abs(pfsw).eq.4.or.abs(pfsw).eq.5)then
c PGF
            write(ipgunit,*)only(1:ns),'\pgfusepath{stroke}'
            write(ipgunit,*)only(1:ns),'\pgfsetcolor{'
     $           ,colorname(li+1)(1:lntrim(colorname(li+1))),'}'
         else
c PS and PCL code combined.
            if(li.lt.15 )then
               spchr(1:3)=' SP'
               if(abs(pfsw).eq.2 .or. abs(pfsw).eq.3) then
                  call iwrite(mod(li,16),wid,spchr(4:6))
               else
                  call iwrite(mod(li,8)+1,wid,spchr(4:6))
               endif
               call abufwrt(spchr(1:3+wid)//' ',4+wid,ipsunit)
            else
               call abufwrt(' SP15 ',6,ipsunit)
            endif
         endif
      endif
      end

c*************************************************************************
      subroutine ibufwrt(i,iunit)
c Add a minimum length integer i to the line buffer for unit iunit
      integer i,iunit,iw
      character*80 str
         call iwrite(i,iw,str)
         call abufwrt(str,iw,iunit)
      end
c*************************************************************************
      subroutine fbufwrt(r,ip,iunit)
c Add a minimum length real (fx.ip) i to the line buffer for unit iunit
      integer iw,ip,iunit
      character*80 str
         call fwrite(r,iw,ip,str)
         call abufwrt(str,iw,iunit)
      end
c*************************************************************************
      subroutine abufwrt(str,la,iunit)
c Add a string str of length la to the line buffer for unit iunit
c If overflowing the line, write line.
      character*(*) str
      integer sblen,iunit,la
      character*80 sbuf
      common /wbuf/sblen,sbuf
c      write(*,*)'abufwrt:',str(1:la)
      if(sblen.le.0) write(*,*)'sblen error:',sblen,la,sbuf,iunit
      if(sblen+la.gt.78) then
         write(iunit,*)sbuf(1:sblen-1)
         sbuf=str(1:la)
         sblen=la+1
      else
         sbuf(sblen:sblen+la)=str(1:la)
         sblen=sblen+la
      endif
      end
c*********************************************************************/
      subroutine lnswrt(iunit,str,iln,ch,iomit)
c Write lines to unit iunit from the compact string str, length iln, using 
c character ch to define the end of line, omitting last iomit characters.
      integer iunit,iln,iomit
      character*(*) str
      character*1 ch
      integer iend,istpos,ilmin
      iend=0
    1 ilmin=iend+1
      iend=istpos(str,ilmin,iln,ch)
      if(iend.eq.0)iend=iln
      write(iunit,'(a)')str(ilmin:iend-iomit)
      if(iend.lt.iln) goto 1
      end
c********************************************************************
      subroutine nlwrt(iunit)
c Start a new line on the unit iunit.
      integer sblen,iunit,la
      character*80 sbuf
      common /wbuf/sblen,sbuf
      if(sblen.le.0) write(*,*)'sblen error:',sblen,la,sbuf,iunit
      write(iunit,*)sbuf(1:sblen-1)
      sblen=1
      end
c*********************************************************************
      subroutine psheader(iunit,pfsw)
      integer iunit,pfsw
c Write the psheader out to iunit.
      character*2 crlf
      parameter (crlf=' \')
c See lnswrt below for the end-of line character definition here.
c g77 is incapable of handling non-printable characters properly.'
c      save ca
      character*1745 ca
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

         if(abs(pfsw).eq.3)then
            ca(00026:00053)='%%BoundingBox: 30 30 620 480'
         elseif(abs(pfsw).eq.2)then
            ca(00026:00053)='%%BoundingBox: 28 28 575 762'
         endif
      call lnswrt(iunit,ca,1745,'\',2)
c ' Done
      end
c******************************************************************
c Obtain the length of a string omitting trailing blanks.
      function lntrim(string)
      character*(*) string
      do i=len(string),1,-1
         if(string(i:i).ne.' ') goto 101
      enddo
      i=0
 101  lntrim=i
      end
