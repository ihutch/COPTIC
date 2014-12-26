c*************************************************************************
      subroutine vec3n(x,y,z,iud)
c Plot a 3-vector in normal units. 
c Domain is normally (x,y,z)=scb(x,y,z)*(-1.,1.)
c World and scale factors for 3-D. World3.h:
c      real xcbc2,ycbc2,scbx3,scby3,scbx3,fixedn
c      integer ihiding
c      common/world3/xcbc2,ycbc2,scbx3,scby3,scbz3,
c     $ wx3min,wx3max,w3nx,wy3min,wy3max,w3ny,wz3min,wz3max,w3nz,
c     $    fixedn,ihiding
      real x,y,z
      integer iud
      real x2,y2,z2
      include 'world3.h'
      call trn32(x,y,z,x2,y2,z2,0)
      call optvecn(x2+xcbc2,y2+ycbc2,iud)
      end
c***********************************************************************
      subroutine scale3(x3min,x3max,y3min,y3max,z3min,z3max)
      real x3min,x3max,y3min,y3max,z3min,z3max
      include 'plotcom.h'
      include 'world3.h'
      wx3min=x3min
      wx3max=x3max
      wy3min=y3min
      wy3max=y3max
      wz3min=z3min
      wz3max=z3max
      call w2scl3
      end
c***********************************************************************
      subroutine getscale3(x3min,x3max,y3min,y3max,z3min,z3max)
      real x3min,x3max,y3min,y3max,z3min,z3max
      include 'plotcom.h'
      include 'world3.h'
      x3min=wx3min
      x3max=wx3max
      y3min=wy3min
      y3max=wy3max
      z3min=wz3min
      z3max=wz3max
      end
c***********************************************************************
      subroutine seti3trunc(itrunc)
c Set the value of the i3trunc switch that determines 3D truncation.
      include 'world3.h'
      i3trunc=itrunc
      end
c***********************************************************************
      subroutine togi3trunc()
c Toggle the value of the i3trunc switch that determines 3D truncation.
      include 'world3.h'
      if(i3trunc.ne.0)then
         i3trunc=0
      else
         i3trunc=1
      endif
      end
c***********************************************************************
      subroutine w2scl3
      include 'plotcom.h'
      include 'world3.h'
      if(wx3max.eq.wx3min)then
         wx3max=wx3min+1.
         write(*,*)'W2SCL3: wxmin=wxmax error.',wx3max
      endif
      w3nx=scbx3*2./(wx3max-wx3min)
      if(wy3max.eq.wy3min)then
         wy3max=wy3min+1.
         write(*,*)'W2SCL3: wymin=wymax error.',wy3max
      endif
      w3ny=scby3*2./(wy3max-wy3min)
      if(wz3max.eq.wz3min)then
         wz3max=wz3min+1.
         write(*,*)'W2SCL3: wzmin=wzmax error.',wx3max
      endif
      w3nz=scbz3*2./(wz3max-wz3min)
      end
c***********************************************************************
      subroutine setcube(cbx,cby,cbz,xcbc,ycbc)
c Set size of unit cell ("cube'), and 2D center position.
      include 'world3.h'
      real cbx,cby,cbz,xcbc,ycbc
      scbx3=cbx
      scby3=cby
      scbz3=cbz
      xcbc2=xcbc
      ycbc2=ycbc
c Probably should not do this here unless defaults are set, 
c which they are
      call w2scl3
      end
c***********************************************************************
      subroutine getcube(cbx,cby,cbz,xcbc,ycbc)
c Get size of unit cell ("cube'), and 2D center position.
      include 'world3.h'
      real cbx,cby,cbz,xcbc,ycbc
      cbx=scbx3
      cby=scby3
      cbz=scbz3
      xcbc=xcbc2
      ycbc=ycbc2
      end
c***********************************************************************
      subroutine vec3w(x,y,z,iud)
c Plot a vector in 3-D world coordinates.
      real x,y,z
      integer iud
      real xn,yn,zn
      include 'world3.h'
      include 'plotcom.h'
      if(i3trunc.eq.0)then
         xn=w3nx*(x-wx3min)-scbx3
         yn=w3ny*(y-wy3min)-scby3
         zn=w3nz*(z-wz3min)-scbz3
      else
         xn=w3nx*(max(min(x,wx3max),wx3min)-wx3min)-scbx3
         yn=w3ny*(max(min(y,wy3max),wy3min)-wy3min)-scby3
         zn=w3nz*(max(min(z,wz3max),wz3min)-wz3min)-scbz3
      endif
c      write(*,*)z,wz3max,w3zmin
c      write(*,*)'xn,yn,zn',xn,yn,zn
      call vec3n(xn,yn,zn,iud)
      end
c***********************************************************************
      subroutine nxyz2wxyz(xn,yn,zn,x,y,z)
c Transform from normal3 to world3
      real x,y,z
      real xn,yn,zn
      include 'world3.h'
      include 'plotcom.h'
      x=(xn+scbx3)/w3nx +wx3min
      y=(yn+scby3)/w3ny +wy3min
      z=(zn+scbz3)/w3nz +wz3min
      end
C********************************************************************
      subroutine wxyz2nxyz(x,y,z,xn,yn,zn)
c Transform from world3 to normal3
      real x,y,z
      real xn,yn,zn
      include 'world3.h'
      include 'plotcom.h'
      xn=w3nx*(x-wx3min)-scbx3
      yn=w3ny*(y-wy3min)-scby3
      zn=w3nz*(z-wz3min)-scbz3
      end
C********************************************************************
      subroutine xy2ncbc(x,y)
c Add centering offsets to x and y
      include 'world3.h'
      x=x+xcbc2
      y=y+ycbc2
      end
C********************************************************************
      function scbn(j)
      include 'world3.h'
      if(j.eq.1) then 
         scbn=scbx3
      elseif(j.eq.2)then 
         scbn=scby3
      elseif(j.eq.3)then 
         scbn=scbz3
      else
         scbn=0.
      endif
      end
C********************************************************************
      subroutine hdprset(isw,fxn)
c Set ihiding and the fixednormal coordinate.
c ihiding: Positive just hides. Negative also projects 3->2 as follows
c  -1,-2,-3 fixed coord, right handed; -4,-5,-6 (+3)= fixed, left handed.
c fixedn:  Value (normal) of the fixed coordinate. (Hence plane pos.)
      integer isw
      real fxn
      include 'world3.h'
      include 'plotcom.h'
      fixedn=fxn
      ihiding=isw
      ltlog=.false.
      if(isw.eq.0)then
         call axregion(0.31,0.91,0.1,0.7)
      else
         call axregion(-scbx3,scbx3,-scby3,scby3)
      endif
      end
C********************************************************************
      subroutine hdproject(hide,project,fixedcoord,fxw,rlsys)
c Set ihiding and the fixednormal coordinate.
c Like hdprset only more comprehensible.
c hide: do/ don't (0) set hiding.
c project: do/ don't (0) set projecting
c fixedcoord: tell which coordinate to fix for projection.
c fxw: World value of coordinate
c rlsys: right/ left (-1) handed system.
      integer hide,project,fixedcoord,rlsys
      real fxw
      include 'plotcom.h'
      include 'world3.h'
      if(project .eq. 0) then 
         ihiding=hide
      else
         ihiding=fixedcoord
         if(rlsys.eq.-1) ihiding=ihiding+3
         ihiding=-ihiding         
      endif
      ltlog=.false.
      if(fixedcoord.eq.1) fixedn=w3nx*(fxw-wx3min)-scbx3
      if(fixedcoord.eq.2) fixedn=w3ny*(fxw-wy3min)-scby3
      if(fixedcoord.eq.3) fixedn=w3nz*(fxw-wz3min)-scbz3
c      write(*,*)ihiding,fixedn
      call hdprset(ihiding,fixedn)
      end
C********************************************************************
      subroutine cubed(icin)
c Omit the lines to corner opposite that chosen.
c Code is 1,2,3,4 anticlockwise from x3min,y3min, +->top -->bottom.
c icorner=0 => omit none.
      integer icin
      include 'plotcom.h'
      include 'world3.h'
      integer ic,i,j,kx,ky,iu
c      write(*,*)'cubed call pfsw=',pfsw

c Chop the top bits off icorner
      icorner=icin-(icin/8)*8
      ic=-sign((mod(abs(icorner)+1,4)+1),icorner)
c      write(*,*)ic
      do 2 i=-1,1,2
         call hdprset(-3,i*scbz3)
         call vecn(-scbx3,-scby3,0)
         do 3 j=1,4
            kx=-1+mod((j+1)/2,2)*2
            ky=-1+mod((j/2),2)*2
            iu=1
            if(icorner.ne.0)then
            if(i*ic.gt.0)then
               if(j.eq.abs(ic).or.j.eq.(mod((abs(ic)+2),4)+1))then
                  iu=0
               endif
            endif
            endif
            call vecn(kx*scbx3,ky*scby3,iu)
    3    continue
    2 continue
      ic=abs(icorner)
      do 1 i=-1,1,2
         call hdprset(-1,i*scbx3)
         call vecn(scby3,-scbz3,0)
         if((ic.eq.1.and.i.eq.1).or.(ic.eq.2.and.i.eq.-1))then
c        call hidvecn(scby3,scbz3,1)
         else
            call vecn(scby3,scbz3,1)
         endif
         call vecn(-scby3,scbz3,0)
         if((ic.eq.3.and.i.eq.-1).or.(ic.eq.4.and.i.eq.1))then
c        call hidvecn(-scby3,-scbz3,1)
         else
            call vecn(-scby3,-scbz3,1)
         endif
    1 continue
      call hdprset(0,0.)
c      write(*,*)'cubed ending pfsw=',pfsw
      end
c***********************************************************************
c Version that finds the right corner to hide lines from.
c And returns it to caller.
      subroutine cubeproj(icorner)
      icorner=igetcubecorner()
      call cubed(icorner)
      end
c***********************************************************************
c Plan for a more systematic axproj. July 93
c ic indicates the corner and whether to flip.
c Vertical axis is always to the left.
c Flip means the vertical axis labels are in the plane including the
c corner, and the left horizontal axis has parallel labels.
c UnFlip means vertical labels in other plane, right horizontal parallel.
c Use upper bits of ic to determine horiz/vert decision:
c Bit0-2= corner (1-4)
c Bit4 (16) set => flip. 
c Bit5 (32) set => x-vertical else horizontal.
c Bit6 (64) set => y-vertical

      subroutine axproj(ic)
c Draw projected axes using the current projection according to ic.
      integer ic,ngpow
      include 'world3.h'
      real fixd
      integer ica,i1,i2
      integer ixc(0:5),iyc(0:5)
      logical flip,xhoriz,yhoriz,ltem
      data ixc/-1,-1,1,1,-1,-1/iyc/1,-1,-1,1,1,-1/
c
      
      yhoriz=(ic/64 - (ic/128)*2) .eq.0
      xhoriz=(ic/32 - (ic/64)*2) .eq.0
      flip=(ic/16 - (ic/32)*2) .ne.0
      ica=abs(ic - (ic/8)*8)
c      write(*,*)'flip,xhor,yhor,ica',flip,xhoriz,yhoriz,ica
      if(yhoriz)then
c Draw y-axis - to + if corner 1 unflipped, or 2, or 3 flipped.
c Else + to -. Parallel labels if 1fl,2un,3fl,4un. 
c Ticrev if perp and 2 or 4
c At -scbx3 if 1 or 4, else +
         call hdprset(-3,-scbz3)
         ltem=(flip.eqv.(mod(ica,2).eq.1))
         if(.not.ltem.and.(ica.eq.4.or.ica.eq.2))call ticrev()
         if(ica.eq.2 .or.(ica.eq.1.and..not.flip)
     $           .or.(ica.eq.3.and.flip))then
            call gaxis(wy3min,wy3max,ngpow,0.,0.,
     $   ixc(ica)*scbx3,ixc(ica)*scbx3,-scby3, scby3,
     $   ltem,.false.)
         else
            call gaxis(wy3max,wy3min,ngpow,0.,0.,
     $   ixc(ica)*scbx3,ixc(ica)*scbx3, scby3,-scby3,
     $   ltem,.false.)
         endif
         if(.not.ltem.and.(ica.eq.4.or.ica.eq.2))call ticrev()
      else
         if(ixc(ica).eq.1)then
c corner 2 or 3 
c Vertical. Draw y axis from - to + , as an x-axis
            call hdprset(-1,scbx3)
            call gaxis(wy3min,wy3max,ngpow,0.,0.,
     $     -scby3,scby3,-scbz3,-scbz3,.true.,.false.)
         else
c corner 1 or 4 Draw y axis from + to - , as a y-axis
         call hdprset(-4,-scbx3)
         call gaxis(wy3max,wy3min,ngpow,0.,0.,
     $     -scbz3,-scbz3,scby3,-scby3,.true.,.false.)
         endif
      endif
      if(xhoriz)then
c Draw x-axis - to + if corner 1, or 2fl, or 4 un.
c Else + to -. Parallel labels if 1un,2fl,3un,4fl. 
c Ticrev if perp and 1 or 3
c At -scby3 if 1 or 2, else +
         call hdprset(-3,-scbz3)
         ltem=(flip.neqv.(mod(ica,2).eq.1))
         if(.not.ltem.and.(ica.eq.3.or.ica.eq.1))call ticrev()
         if(ica.eq.1 .or.(ica.eq.4.and..not.flip)
     $           .or.(ica.eq.2.and.flip))then
            call gaxis(wx3min,wx3max,ngpow,0.,0.,
     $   -scbx3, scbx3,iyc(ica)*scby3,iyc(ica)*scby3,
     $   ltem,.false.)
         else
            call gaxis(wx3max,wx3min,ngpow,0.,0.,
     $    scbx3,-scbx3,iyc(ica)*scby3,iyc(ica)*scby3,
     $   ltem,.false.)
         endif
         if(.not.ltem.and.(ica.eq.3.or.ica.eq.1))call ticrev()
      else
         if(ica.le.2)then
c corner 1 or 2 draw x-axis - to + as an x-axis
         call hdprset(-5,-scby3)
         call gaxis(wx3min,wx3max,ngpow,0.,0.,
     $     -scbx3,scbx3,-scbz3,-scbz3,.true.,.false.)
         else
c corner 3 or 4 draw x-axis + to - as a y-axis
         call hdprset(-2,scby3)
         call gaxis(wx3max,wx3min,ngpow,0.,0.,
     $     -scbz3,-scbz3,scbx3,-scbx3,.true.,.false.)
         endif
      endif
c z- axis.
      i2=1
      if(ica.le.2)i2=-i2
      if(ica.eq.3.and..not.flip)i2=-i2
      if(ica.eq.1.and..not.flip)i2=-i2
      i1=ica
      if(flip)i1=i1-1
      if(i1.ne.1.and.i1.ne.3)then
         fixd=i2*scbx3
         if(i1.eq.2)then
            call hdprset(-1,fixd)
            call gaxis(wz3min,wz3max,ngpow,0.,0.
     $     ,-scby3,-scby3,-scbz3,scbz3,.false.,.false.)
         else
c        i1=0 or 4.
            call hdprset(-4,fixd)
            call gaxis(wz3min,wz3max,ngpow,0.,0.
     $     ,-scbz3,scbz3,scby3,scby3,.false.,.false.)
         endif
      else
         fixd=i2*scby3
         if(i1.eq.1)then
            call hdprset(-5,fixd)
            call gaxis(wz3min,wz3max,ngpow,0.,0.
     $     ,-scbx3,-scbx3,-scbz3,scbz3,.false.,.false.)
         else
c        i1=3
            call hdprset(-2,fixd)
            call gaxis(wz3min,wz3max,ngpow,0.,0.
     $     ,-scbz3,scbz3,scbx3,scbx3,.false.,.false.)
         endif
      endif
      call hdprset(0,0.)
      end
c****************************************************************************
      subroutine geteye(x2,y2,z2)
c Read eye.
      open(unit=2,file='eye.dat',status='old',err=99)
      read(2,*,err=99)x2,y2,z2
      close(2)
      if(x2**2+y2**2+z2**2.eq.0.)goto 99
      if(x2.eq.0. .and. y2.eq.0.)goto 99
      goto 98
 99   continue
c     If eye.dat has an error (e.g. does not exist) Default eye.
      write(*,*)'File eye.dat invalid. Using default'
      x2=4.
      y2=-10.
      z2=6.
 98   continue
      end
c****************************************************************************
      subroutine puteye(x2,y2,z2)
c Write eye.
      open(unit=2,file='eye.dat',status='unknown',err=99)
      write(2,*)x2,y2,z2
      close(2)
 99   continue
      end
c********************************************************************
c Return the nearest corner to eye in standard convention.
      function igetcorner()
c This was wrong. geteye reads the eyefile.
c      call geteye(x2,y2,z2)
      call trn32(xdum,ydum,zdum,x2,y2,z2,-1)
      if(y2.le.0.)then 
         if(x2.le.0.)then
            icorner=1
         else
            icorner=2
         endif
      else
         if(x2.le.0.)then
            icorner=4
         else
            icorner=3
         endif
      endif
c      write(*,*)'x2,y2,z2',x2,y2,z2,mod(icorner,2)
c If appropriate tell axproj to flip labels.
      if(abs(x2).gt.abs(y2).eqv.(mod(icorner,2).ne.0))
     $     icorner=icorner+16
c If appropriate use vertical labels
c        if(z2*z2.lt.x2*x2) icorner=icorner+32
c Trying for better results.
      xy2i=min(x2**2,y2**2)
      if(z2*z2.lt.xy2i+.2*y2**2) icorner=icorner+32
      if(z2*z2.lt.xy2i+.2*x2**2) icorner=icorner+64
      igetcorner=icorner
      end
c********************************************************************
c Return the nearest corner to eye for use with cubed
      function igetcubecorner()
      call trn32(xdum,ydum,zdum,x2,y2,z2,-1)
      if(y2.le.0.)then 
         if(x2.le.0.)then
            icorner=1
         else
            icorner=2
         endif
      else
         if(x2.le.0.)then
            icorner=4
         else
            icorner=3
         endif
      endif
      if(z2.lt.0)icorner=-icorner
      igetcubecorner=icorner
      end
