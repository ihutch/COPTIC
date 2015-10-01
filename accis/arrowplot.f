c Vectorial Arrow-Plot
      subroutine arrowplot(E1,E2,Erange,Li,imax,jmax,x,y,iswin,is1,is2)
c Draw arrows representing the field magnitude and direction of
c E1/2(imax/of/Li,jmax) on a domain x(imax[,jmax]),y([imax,]jmax).
c On entry the domain should have been set up so that x, y span it.
c Erange is the scaling parameter for arrow length, relative to 
c 1 unit of normalized position.
c If on entry Erange=0, calculate it automatically from the arrays.
c isw: 0: x,y unused; just index. 1: x,y vectors, 2: x,y matrices
c isw: if 4 is added to isw, use is1, is2 as the steps in x,y position
c      if not, then is1, and is2 arguments may be omitted in call.
c isw: if 8 is added to isw, limit the arrow arrays to narm in each direction
c      i.e. calculate is1 and is2 algorithmically internally.


      include 'plotcom.h'
      integer Li,imax,jmax,isw
      real E1(Li,jmax),E2(Li,jmax),x(*),y(*)
      character*20 string
      integer narm
      parameter (narm=20)

      isw=iswin
      if(isw/4-2*(isw/8).ne.0)then
         isl1=is1
         isl2=is2
      else
         isl1=1
         isl2=1
      endif
      if(isw/8-2*(isw/16).ne.0)then
         isl1=(imax-1)/narm+1
         isl2=(jmax-1)/narm+1
      endif
      isw=isw-4*(isw/4)

      if(Erange.eq.0)then
c Find Erange by scaling.
         call minmax2(E1,Li,imax,jmax,e1min,e1max)
         call minmax2(E2,Li,imax,jmax,e2min,e2max)
         e2max=max(abs(e2min),abs(e2max))*(jmax/float(isl2))
         e1max=max(abs(e1min),abs(e1max))*(imax/float(isl1))
c         write(*,*)'e1max=',e1max,'  e2max=',e2max
         Erange=max(e1max,e2max)
         if(Erange.eq.0.)then
            write(*,*)'Arrowplot scaling error. Null vector field.'
            return
         endif
         call fitrange(0.,Erange*.15,2,ipow,fac10,delta,first,xlast) 
c         write(*,*)Erange,delta,first,xlast,fac10
         csize=delta/Erange
         call charangl(90.)
         call fwrite(delta,iw,1,string)
         call jdrwstr(naxmax+chrshght+.01,naymin+0.1,string(1:iw),1.)
         call charsize(csize,0.3*csize)
         call drcstr('!A_!@')
      endif

      do j=1,jmax,isl2
         do i=1,imax,isl1
            Ex=E1(i,j)
            Ey=E2(i,j)
            angle=atan2(Ey,Ex)
            Emag=sqrt(Ex**2 + Ey**2)
            csize=Emag/Erange+1.e-6
            if(isw.eq.0)then
               xc=i
               yc=j
            elseif(isw.eq.1)then
               xc=x(i)
               yc=y(j)
            elseif(isw.eq.2)then
               xc=x(i+Li*(j-1))
               yc=y(i+Li*(j-1))
            endif
            call charangl(180.*angle/3.1415926)
            call charsize(csize,0.3*csize)
            call jdrwstr(wx2nx(xc),wy2ny(yc),'!A_',0.)
c            call jdrwstr(wx2nx(xc),wy2ny(yc),'!A!!',0.)
c            write(*,*)Ex,Ey,Emag,angle
         enddo
      enddo

c Restore defaults.
      call charsize(0.,0.)
      call charangl(0.)
      call drcstr('!@')
      end
c*******************************************************************
      subroutine arg3arrow(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12)
c Convenience routine calls draw3arrow with explicit arguments
      real arrow(12)
      arrow(1)=a1
      arrow(2)=a2
      arrow(3)=a3
      arrow(4)=a4
      arrow(5)=a5
      arrow(6)=a6
      arrow(7)=a7
      arrow(8)=a8
      arrow(9)=a9
      arrow(10)=a10
      arrow(11)=a11
      arrow(12)=a12
      call draw3arrow(arrow)
      end
c*******************************************************************
      subroutine draw3arrow(arrow)
c Draw a 3-D arrow whose geometry is specified by data in the array
c arrow which are: Base position world coords (3), Point posn world (3),
c Barb parallel, perp size as fraction of length (2). 
c Feather size parallel, perp, frac of length (2).
c If barbperp is negative, read two additional parameters at the end.
c nangle (1) the number of barbs. isw a further switch.
c 
      integer ngp
      parameter (ibx=1,iby=2,ibz=3,ipx=4,ipy=5,ipz=6,ibba=7,ibbe=8,
     $     ifa=9,ife=10,ina=11,isw=12,ngp=isw)
      integer nanglemax
      parameter (nanglemax=32)
      real arrow(ngp)
      real v1(3)
      real avec(3),aperp(3),across(3),alamb(3)
      integer nangle
      real shaftfrac
      integer nalong
      parameter (nalong=4)
      real shaftcoords(nalong,nanglemax,3)
      real xshaft(nalong,nanglemax),yshaft(nalong,nanglemax)
     $     ,zshaft(nalong,nanglemax)
      real work(0:nalong+1,0:nanglemax)
      real d(5)
      equivalence (xshaft,shaftcoords(1,1,1))
     $     ,(yshaft,shaftcoords(1,1,2))
     $     ,(zshaft,shaftcoords(1,1,3))
      data nangle/4/shaftfrac/.4/

      iw=0
      ncolor=igetcolor()
      if(arrow(ibbe).lt.0)then
         nangle=int(arrow(ina))
         iw=int(arrow(isw))
      endif
      amag=0.
      do k=1,3
         avec(k)=(arrow(ipx-1+k)-arrow(ibx-1+k))
         amag=amag+avec(k)**2
      enddo
      amag=sqrt(amag)
c Decide a perpendicular direction.
      iperp=3
 1    do k=1,3
         aperp(k)=0.
         if(k.eq.iperp)aperp(k)=1.
      enddo
      call crossprod(aperp,avec,across,acmag)
      if(acmag.lt.1.e-5*amag)then
         iperp=iperp-1
         goto 1
      endif
c Now across is a vector perpendicular to both the arrow and coordinate
c direction iperp. Calculate the other perpendicular vector and 
c normalize all.
      call crossprod(avec,across,alamb,almag)
      do k=1,3
         avec(k)=avec(k)/amag
         across(k)=across(k)/acmag
         alamb(k)=alamb(k)/almag
      enddo
c Geometry completed.

c Draw shaft
      call vec3w(arrow(ibx),arrow(iby),arrow(ibz),0)
      call vec3w(arrow(ipx),arrow(ipy),arrow(ipz),1)
c Draw feather
      call vec3w(arrow(ibx),arrow(iby),arrow(ibz),0)
      do i=1,nangle+1
         angle=2.*3.141593*(i-1.)/(nangle)
         do k=1,3
            v1(k)=arrow(ibx-1+k)-arrow(ifa)*avec(k)
     $           +arrow(ife)*amag*
     $           (cos(angle)*across(k)+sin(angle)*alamb(k))
         enddo
c         call vec3w(v1(1),v1(2),v1(3),1)
         call vec3w(arrow(ibx),arrow(iby),arrow(ibz),1)
         call vec3w(v1(1),v1(2),v1(3),0)
      enddo
      call vec3w(arrow(ipx),arrow(ipy),arrow(ipz),0)
c Maybe draw a finite width shaft using face-orientation shading.
      if(iw/2-2*(iw/4).eq.1)then
         d(1)=-alamb(1)
         d(2)=-alamb(2)
         d(3)=-alamb(3)
c         d(1)=-across(1)
c         d(2)=-across(2)
c         d(3)=-across(3)
c 4,5 give distance limits.
c For absolute position we need this:
c         d(4)=-.5*abs(arrow(ibbe))
c         d(5)=.5*abs(arrow(ibbe))
c         do k=1,3
c            d(4)=d(4)+across(k)*arrow(k)
c            d(5)=d(5)+across(k)*arrow(k)
c         enddo
c Color limits of cosine angle. 
         d(4)=-1.2
         d(5)=1.2
         do i=1,nangle+1
            angle=2.*3.141593*(i-1.)/(nangle)
            do k=1,3
               transv=(cos(angle)*across(k)+sin(angle)*alamb(k))
               shaftcoords(2,i,k)=arrow(ipx-1+k)
     $              -arrow(ibba)*avec(k)
     $              +shaftfrac*abs(arrow(ibbe))*amag*transv
               shaftcoords(1,i,k)=arrow(ibx-1+k)
     $              -shaftfrac*arrow(ifa)*avec(k)
     $              +shaftfrac*abs(arrow(ife))*amag*transv
               shaftcoords(3,i,k)=arrow(ipx-1+k)
     $              -arrow(ibba)*avec(k)
     $              +abs(arrow(ibbe))*amag*transv
               shaftcoords(4,i,k)=arrow(ipx-1+k)
c Avoid degenerate point:
     $              +1.e-5*abs(arrow(ibbe))*amag*transv
            enddo
         enddo
         call color(0)
c Shade by orientation, periodic in j-index
         jsw=5+16
c         jsw=5
c         call accisgradinit(22000,-40000,-40000,64000,64000,64000)
         call surfdr3(xshaft,yshaft,zshaft,nalong,nalong,nangle+1,work
     $        ,jsw,d)
      endif
c Draw barb. Normally included in the finite thickness shaft.
      if(.true.)then
      call color(ncolor)
      call vec3w(arrow(ipx),arrow(ipy),arrow(ipz),0)
      do i=1,nangle+1
         angle=2.*3.141593*(i-1.)/(nangle)
         do k=1,3
            v1(k)=arrow(ipx-1+k)-arrow(ibba)*avec(k)
     $           +arrow(ibbe)*amag*
     $           (cos(angle)*across(k)+sin(angle)*alamb(k))
         enddo
         if(iw-2*(iw/2).eq.1)call vec3w(v1(1),v1(2),v1(3),1)
         call vec3w(arrow(ipx),arrow(ipy),arrow(ipz),1)
         call vecfill()
         call vec3w(v1(1),v1(2),v1(3),0)
      enddo
      endif
      call color(ncolor)
      end
c*****************************************************************
      subroutine crossprod(v1,v2,result,rmag)
c Form the cross-product between v1 and v2, returning also the magnitude
c of the result.
      real v1(3),v2(3),result(3),rmag
      rmag=0.
      do k=1,3
         k1=mod(k,3)+1
         k2=mod(k+1,3)+1
         result(k)=v1(k1)*v2(k2)-v1(k2)*v2(k1)
         rmag=rmag+result(k)**2
      enddo
      rmag=sqrt(rmag)
      end
c******************************************************************
c*******************************************************************
      subroutine arrow3path(path,nalen,nangle,isw,tr,bbl,bbw,nbb)
c Draw a 3-D path whose track is specified by data in the array 
c path(3,nalen)
c nangle is the number of angle positions of barbs or cylinders. 
c isw a further switch:
c isw bit0 .eq. 1 => draw feather for line drawing else no feather.
c isw bit1 .eq. 1 => barbfill for line drawing, else no fill, just line.
c isw bit3 .eq. 1 => Don't recalculate surface, just redraw.
c tr tube radius. Zero => line drawing, rather than surface.
c bbl, bbw, Barb parallel, perp size (absolute units)
c nbb indicates the number of path sections between barbs.
c The last path position (istep=nalen) always has a barb. 
c So, other barbs are at end of sections istep if mod(nalen-i,nbb).eq.0.
c Using nbb=1 => barb every section, nbb=2 every other section.
c There are approximately nalen/nbb barbs.  And barbs take 3 array
c positions instead of 1.  So the total length of the array needed is
c nalen*(1+2./nbb)
c 
      implicit none
      integer nalen,nbb,nangle,isw
      real path(3,nalen),bbl,bbw,tr

c Local Storage
      integer nanglemax,nalong,iaj
      parameter (nanglemax=32,nalong=1000)
      real v1(3)
      real avec(3),aperp(3),across(3),alamb(3)
      real shaftcoords(nalong,nanglemax,3)
      real xshaft(nalong,nanglemax),yshaft(nalong,nanglemax)
     $     ,zshaft(nalong,nanglemax)
      equivalence (xshaft,shaftcoords(1,1,1))
     $     ,(yshaft,shaftcoords(1,1,2))
     $     ,(zshaft,shaftcoords(1,1,3))
      real work(0:nalong+1,0:nanglemax)
      real d(5)
      real amag,acmag,almag,angle,transv
      integer ncolor,igetcolor,ipx,ibx,istep,i,k,iperp,jsw
      save shaftcoords,d,ipx

c      write(*,*)'nalen,nangle,isw,tr,bbl,bbw,nbb,path='
c     $     ,nalen,nangle,isw,tr,bbl,bbw,nbb,path

      ncolor=igetcolor()

c If we are just redrawing surface, skip to the end.
      if(isw/4-2*(isw/8).eq.1 .and. tr.gt.0.)goto 4

c Here we check whether the tube will fit in storage.
      if(.not.(int(nalen*(1+2./nbb)).le.nalong-2)
     $     .and.tr.gt.0.)then
         write(*,*)'Arrow3path parameters',nalen,nbb
     $        ,' exceed available storage',nalong
         stop
      endif
      iaj=-1
      if(isw-2*(isw/2).eq.0.and.tr.le.0.)then
c No feather for line drawing.         
         iaj=0
      endif
      ipx=2
c Start of path step iteration.
      do istep=2,nalen
         ibx=ipx-1
c Find avec, the path step.
         amag=0.
         do k=1,3
            avec(k)=path(k,istep)-path(k,istep-1)
            amag=amag+avec(k)**2
         enddo
         amag=sqrt(amag)
c Decide a perpendicular direction. aperp is along one coordinate 
c usually iperp=3, z, but a different one if avec is along z.
         iperp=3
 1       continue
         do k=1,3
            aperp(k)=0.
            if(k.eq.iperp)aperp(k)=1.
         enddo
         call crossprod(aperp,avec,across,acmag)
         if(acmag.lt.1.e-5*amag)then
            iperp=iperp-1
            goto 1
         endif
c Now find across, a vector perpendicular to both the path, avec, and
c coordinate direction iperp. Calculate the other perpendicular vector
c and normalize both.
         call crossprod(avec,across,alamb,almag)
         do k=1,3
            avec(k)=avec(k)/amag
            across(k)=across(k)/acmag
            alamb(k)=alamb(k)/almag
         enddo
c Geometry completed for this step.

         if(istep.eq.2)then
c First time, decide the direction of illumination (shading cosine).
            d(1)=-alamb(1)
            d(2)=-alamb(2)
            d(3)=-alamb(3)
c Color limits of cosine angle. 
            d(4)=-1.2
            d(5)=1.2
         endif

 2       continue
         if(tr.le.0)then
c Line-drawn arrow. 
c First time through draw the feather as indicated by non-zero iaj.
            ipx=ipx+iaj
c Draw line shaft for this step
            call vec3w(path(1,ibx),path(2,ibx),path(3,ibx),0)
            call vec3w(path(1,ipx),path(2,ipx),path(3,ipx),1)
            if(mod(nalen-istep,nbb).eq.0)then
c Draw line barb.
               call color(ncolor)
c               call vec3w(path(1,ipx),path(2,ipx),path(3,ipx),0)
c               write(*,*)'Drawing line barb',istep,ipx,bbl,bbw,isw,iaj
c     $              ,nangle
               do i=1,nangle+1
                  angle=2.*3.141593*(i-1.)/(nangle)
                  do k=1,3
                     transv=(cos(angle)*across(k)+sin(angle)*alamb(k))
                     v1(k)=path(k,ipx)-bbl*avec(k)+bbw*transv
                  enddo
                  if(isw/2-2*(isw/4).eq.1)call vec3w(v1(1),v1(2),v1(3)
     $                 ,1)
                  call vec3w(path(1,ipx),path(2,ipx),path(3,ipx),1)
                  call vecfill()
                  call vec3w(v1(1),v1(2),v1(3),0)
               enddo
            endif
            if(iaj.le.-1)then
c We have just done the feather; switch off and repeat leading point.
               ipx=ipx-iaj
               iaj=0
               goto 2
            endif
            ipx=ipx+1
         else
 3          continue
c Write into finite width shaft using face-orientation shading.
            if(mod(nalen-istep,nbb).eq.0 .and. iaj.ne.-1)then
c Barb at this step
c               write(*,*)'Barb write',ipx,iaj,ibx,istep,avec
c     $              ,(path(k,istep+iaj),k=1,3),bbl,bbw
c               write(*,*)'amag,avec=',amag,avec
               do i=1,nangle+1
                  angle=2.*3.141593*(i-1.)/(nangle)
                  do k=1,3
                     transv=(cos(angle)*across(k)+sin(angle)*alamb(k))
                     shaftcoords(ipx,i,k)=path(k,istep+iaj)-bbl*avec(k)
     $                    +tr*transv
                     shaftcoords(ipx+1,i,k)=path(k,istep+iaj)-bbl
     $                    *avec(k)+bbw*transv
                     shaftcoords(ipx+2,i,k)=path(k,istep+iaj)
     $                    +tr*transv
                     if(istep.eq.nalen)then
c Close the leading arrow avoiding degenerate point:
                        shaftcoords(ipx+2,i,k)=path(k,istep+iaj)
     $                    +1.e-5*tr*transv
                     endif
                  enddo
               enddo
               ipx=ipx+3
            else
c No barb at this step. Just add to the path:
               ipx=ipx+iaj
               do i=1,nangle+1
                  angle=2.*3.141593*(i-1.)/(nangle)
                  do k=1,3
                     transv=(cos(angle)*across(k)+sin(angle)*alamb(k))
                     shaftcoords(ipx,i,k)=path(k,istep+iaj)+tr*transv
                  enddo
               enddo
               if(iaj.le.-1)then
c Finish the extra first step.
                  ipx=ipx-iaj
                  iaj=0
                  goto 3
               else
                  ipx=ipx+1
               endif
            endif
         endif
c End of path step iteration
      enddo

 4    if(tr.gt.0.)then
         call color(0)
c Shade by orientation, periodic in j-index
         jsw=5+16
         call surfdr3(xshaft,yshaft,zshaft,nalong,ipx-1,nangle+1,work
     $        ,jsw,d)
      endif
      call color(ncolor)
      end
c*****************************************************************
