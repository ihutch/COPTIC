c*********************************************************************
      subroutine axlabels(xaxlab, yaxlab)
c    Axis labels.   */
      character*(*) xaxlab,yaxlab
      include 'plotcom.h'
      real ctemp,stemp,temp
      ctemp=chrscos
      stemp=chrssin
      chrscos=1.
      chrssin=0.
      call jdrwstr((naxmin+naxmax)*.5,
     $        naypt+xticoff-2.0*chrshght,xaxlab,0.)
      chrscos=0.
      chrssin=1.
      temp=(naxpt+5.0*yticoff)
      call jdrwstr(temp,(naymin+naymax)*.5,yaxlab,0.)
      chrscos=ctemp
      chrssin=stemp
      return
      end
c***********************************************************************
      subroutine axident3old()
      include 'world3.h'
c draw an axis orientation
      call trn32(-scbx3,-scby3,0.7*scbz3,x2,y2,z2,0)
      call vecn(x2+xcbc2,y2+ycbc2,0)
      call drcstr('!B z!@')
      call trn32(0.7*scbx3,-scby3,-scbz3,x2,y2,z2,0)
      call vecn(x2+xcbc2,y2+ycbc2,0)
      call drcstr('!B x!@')
      call trn32(-scbx3,0.7*scby3,-scbz3,x2,y2,z2,0)
      call vecn(x2+xcbc2,y2+ycbc2,0)
      call drcstr('!B y!@')
      end
c***********************************************************************
      subroutine axident3()
      call ax3labels('           x','           y','      z')
      end
c***********************************************************************
      subroutine boxtitle(title)
      character*(*) title
      include 'plotcom.h'
      real ctemp,stemp
      ctemp=chrscos
      stemp=chrssin
      chrscos=1.
      chrssin=0.
      call jdrwstr( (naxmin+naxmax)*0.5,(naymax-xticoff),title,0.)
      chrscos=ctemp
      chrssin=stemp
      end
c***********************************************************************
c Based on the same code decisions as axproj
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

      subroutine ax3labels(xaxlab,yaxlab,zaxlab)
      character*(*) xaxlab,yaxlab,zaxlab
c Draw projected axes using the current projection according to ic.
      integer ic
      include 'plotcom.h'
      include 'world3.h'
      real fixd
      integer ica,i1,i2
      integer ixc(0:5),iyc(0:5)
      logical flip,xhoriz,yhoriz,ltem
      data ixc/-1,-1,1,1,-1,-1/iyc/1,-1,-1,1,1,-1/
c Get the corner code.
      ic=igetcorner()
c
      ctemp=chrscos
      stemp=chrssin
      
      yhoriz=(ic/64 - (ic/128)*2) .eq.0
      xhoriz=(ic/32 - (ic/64)*2) .eq.0
      flip=(ic/16 - (ic/32)*2) .ne.0
      ica=abs(ic - (ic/8)*8)
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
c           call gaxis(wy3min,wy3max,ngpow,0.,0.,
c     $           ixc(ica)*scbx3,ixc(ica)*scbx3,-scby3, scby3,
c     $           ltem,.false.)
c y-axis as if x
            chrscos=0.
            chrssin=1.*ixc(ica)
            xlp=ixc(ica)*scbx3+3.*(xticoff)
            if(ltem)xlp=xlp-5.0*xticoff
c            xlp=ixc(ica)*scbx3-ixc(ica)*(xticoff-2.0*chrshght)
            call jdrwstr(xlp,0.,yaxlab,0.)
         else
c           call gaxis(wy3max,wy3min,ngpow,0.,0.,
c     $           ixc(ica)*scbx3,ixc(ica)*scbx3, scby3,-scby3,
c     $           ltem,.false.)
c y-axis as if x
            chrscos=0.
            chrssin=1.*ixc(ica)
            xlp=ixc(ica)*scbx3-3.*xticoff
            if(ltem) xlp=xlp+5.*xticoff
            call jdrwstr(xlp,0.,yaxlab,0.)
         endif
         if(.not.ltem.and.(ica.eq.4.or.ica.eq.2))call ticrev()
      else
         if(ixc(ica).eq.1)then
c corner 2 or 3 
c Vertical. Draw y axis from - to + , as an x-axis
            call hdprset(-1,scbx3)
c           call gaxis(wy3min,wy3max,ngpow,0.,0.,
c     $    -scby3,scby3,-scbz3,-scbz3,.true.,.false.)
c y-axis as if x
            chrscos=1.
            chrssin=0.
            call jdrwstr(0.,-scbz3+xticoff-2.0*chrshght,
     $           yaxlab,0.)
         else
c corner 1 or 4 Draw y axis from + to - , as a y-axis
         call hdprset(-4,-scbx3)
c        call gaxis(wy3max,wy3min,ngpow,0.,0.,
c     $    -scbz3,-scbz3,scby3,-scby3,.true.,.false.)
            chrscos=0.
            chrssin=1.*ixc(ica)
            call jdrwstr(-scbz3+xticoff-2.0*chrshght,0.,
     $           yaxlab,0.)
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
c           call gaxis(wx3min,wx3max,ngpow,0.,0.,
c     $  -scbx3, scbx3,iyc(ica)*scby3,iyc(ica)*scby3,
c     $   ltem,.false.)
            if(.not.ltem)then
               if(flip)then
                  chrscos=1.
                  ylp=iyc(ica)*scby3-3.*(xticoff)
               else
                  chrscos=-1.
                  ylp=iyc(ica)*scby3-xticoff+4.*chrshght
               endif
            else
               chrscos=1.
               ylp=iyc(ica)*scby3-(-xticoff+2.*chrshght)
            endif
            chrssin=0.
            call jdrwstr(0.,ylp,xaxlab,0.)

         else
c           call gaxis(wx3max,wx3min,ngpow,0.,0.,
c     $   scbx3,-scbx3,iyc(ica)*scby3,iyc(ica)*scby3,
c     $   ltem,.false.)
c               write(*,*)'ltem',ltem
            if(ltem)then
               chrscos=-1.
               ylp=iyc(ica)*scby3-xticoff+2.*chrshght
            else
               if(flip)then
                  chrscos=-1.
                  ylp=iyc(ica)*scby3+3.*xticoff
               else
                  chrscos=1.
                  ylp=iyc(ica)*scby3+xticoff-4.*chrshght
               endif
            endif
            chrssin=0.
            call jdrwstr(0.,ylp,xaxlab,0.)
         endif
         if(.not.ltem.and.(ica.eq.3.or.ica.eq.1))call ticrev()
      else
         if(ica.le.2)then
c corner 1 or 2 draw x-axis - to + as an x-axis
            call hdprset(-5,-scby3)
c            call gaxis(wx3min,wx3max,ngpow,0.,0.,
c     $           -scbx3,scbx3,-scbz3,-scbz3,.true.,.false.)
            chrscos=1.
            chrssin=0.
            call jdrwstr(0.,-scbz3+xticoff-2.*chrshght,xaxlab,0.)
         else
c corner 3 or 4 draw x-axis + to - as a y-axis
            call hdprset(-2,scby3)
c            call gaxis(wx3max,wx3min,ngpow,0.,0.,
c     $           -scbz3,-scbz3,scbx3,-scbx3,.true.,.false.)
            chrscos=0.
            chrssin=-1.
            call jdrwstr(-scbz3+xticoff-2.*chrshght,0.,xaxlab,0.)
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
c           call gaxis(wz3min,wz3max,ngpow,0.,0.
c     $    ,-scby3,-scby3,-scbz3,scbz3,.false.,.false.)
            chrscos=1.
            chrssin=0.
            call jdrwstr(-scby3,1.3*chrshght,zaxlab,-1.2)
         else
c        i1=0 or 4.
            call hdprset(-4,fixd)
c           call gaxis(wz3min,wz3max,ngpow,0.,0.
c     $    ,-scbz3,scbz3,scby3,scby3,.false.,.false.)
            chrscos=0.
            chrssin=-1.
            call jdrwstr(1.3*chrshght,scby3,zaxlab,-1.2)
         endif
      else
         fixd=i2*scby3
         if(i1.eq.1)then
            call hdprset(-5,fixd)
c           call gaxis(wz3min,wz3max,ngpow,0.,0.
c     $    ,-scbx3,-scbx3,-scbz3,scbz3,.false.,.false.)
            chrscos=1.
            chrssin=0.
            call jdrwstr(-scbx3,1.3*chrshght,zaxlab,-1.2)
         else
c        i1=3
            call hdprset(-2,fixd)
c           call gaxis(wz3min,wz3max,ngpow,0.,0.
c     $    ,-scbz3,scbz3,scbx3,scbx3,.false.,.false.)
            chrscos=0.
            chrssin=-1.
            call jdrwstr(1.3*chrshght,scbx3,zaxlab,-1.2)
         endif
      endif
      call hdprset(0,0.)
      chrscos=ctemp
      chrssin=stemp
      end
