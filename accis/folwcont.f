c Modified Mar 94 to use parabolic interpolation.
c Aug 94 to use more compact character ppath. Less than 10% slower.
c Contouring routines using following algorithm.
c Summary: 
c       CONTREC(z,ppath,imax,jmax,cl,nc) - simple shell for elementary call.
c Use a rectangular mesh. If nc=0, fit contours, else use cl(nc).
c
c       CONTOURL(z,ppath,L,imax,jmax,cl,nc,x,y,consw) - Full call.
c Contour a function z(1:imax(of L),jmax), at nc levels, cl(nc), or fit
c if nc=0, on a mesh defined by x,y. If nc .lt. 0 use no labels.
c ppath is the character marker array and must have 
c dimension at least (imax,jmax)
c Switch consw determines the type of plot. Lower nibble
c   0   Equally spaced arrays. Arguments x and y are not used.
c   1   Vectors x and y determine the unequally spaced arrays.
c   2   Arrays x and y determine the arbitrary mesh.
c General call should be initialized by first.
c
c       CONSGEN - Searches for contour start and calls confol
c       CONFOL - Follows the contour.
c       MESH2W(xc,yc,ic,x,y,l,consw)
c         - Transform xc,yc(ic) to the mesh x,y, depending on switch consw.
c
c External calls: minmax2,PLTINIT, LABELINE, POLYLINE... 
c************************************************************************
      subroutine contrec(z,ppath,imax,jmax,cl,nc)
c Basic call: Contour a whole array on a rectangular equally spaced mesh.
c No axes.
      integer imax,jmax,nc
      real z(imax,jmax),cl(1)
      character ppath(imax,jmax)
      call pltinit(1.,float(imax),1.,float(jmax))
      call axbox
      call contourl(z,ppath,imax,imax,jmax,cl,nc,z,z,0)
      end
c************************************************************************
      subroutine autocolcont(z,ildx,imax,jmax)
      common/cont1stlast/c1st,clast
      integer ildx,imax,jmax
c Color (only) Contour of z on rectangular mesh.
      real z(ildx,jmax)
      real cldummy(1),xdummy(1),ydummy(1)
      character ppathdummy(1)

      call pltinit(1.,float(imax),1.,float(jmax))
      cldummy(1)=0.
      call contourl(z,ppathdummy,ildx,imax,jmax,
     $     cldummy,0,xdummy,ydummy,48)
      call color(15)
      call axbox
      call gradlegend(c1st,clast,-.2,0.,-.2,1.,.03,.false.) 
      end
c************************************************************************
      subroutine contourl(z,ppath,l,imax,jmax,cl,nc,x,y,consw)
c Contour a function on a mesh defined by x,y.
c 2 Aug 92 Labelled contouring. Feb 93.
c 1 Jan 2001 coloring by height.
c Plot should be initialized first.
      integer l,imax,jmax,nc
      real z(l,jmax),x(l),y(l),cl(*)
      character ppath(l,jmax)
      integer consw
c Switch consw determines the type of plot. Lower nibble
c   0   Equally spaced arrays. Arguments x and y are not used.
c   1   Vectors x and y determine the unequally spaced arrays.
c   2   Arrays x and y determine the arbitrary mesh.
c   Bit 5 (16) set, do coloring
c   Bit 6 (32) set, don't contour (only color).
c   Bit 7 (64) set, color using triangle gradients.
c     In which case second byte tells the level skip (default 1).
c If nc lt 0 no labels. If nc eq 0 then fit contours; use cl(1) for
c number of contours if non-zero. Else use array cl(nc).
c 
      include 'plotcom.h'
      real cv
      integer ic,ici,i
      integer point,width,cyc
      character str1*30,str2*30
      common/lablnc/width,str1
      external followdata
c Maximum length of a single contour. Increase if necessary.
      parameter (ici=4000)
      real xc(ici),yc(ici)
      real cw,ch
      logical labels
      integer inc,in,nxfac,theconsw,icfil,itri
      parameter (inc=7)
      real minz,maxz,xfac,xtic,x1st,xlast
      real xd(4),yd(4),zd(4)
      parameter (ngradcol=240)
      common/cont1stlast/c1st,clast

c      write(*,*)'Inside contourl'
c save charsize
      cw=chrswdth
      ch=chrshght
c Interpret switch.
      theconsw=consw- 16*(consw/16)
      icfil=consw/16 - 2*(consw/32)
      inocont=consw/32-2*(consw/64)
      itri=consw/64 -  2*(consw/128)
c Second byte
      istep=consw/256-256*(consw/(256*256))
c Third byte plus switch
      ifcol=(consw/65536-256*(consw/(256*65536)))*65536 + theconsw
c      write(*,*)'consw,theconsw,ifcol,icfil',consw,theconsw,ifcol,icfil
      if(nc.eq.0)then
c Contour level fitting.
         call minmax2(z,L,imax,jmax,minz,maxz)
         in=inc
         if(int(abs(cl(1))).ne.0) in=int(abs(cl(1)))
         if(.false.)then
c Obsolete very complicated scheme for label cycle calculation.
c Get a good spacing that spans the data range with no more than 8
c cycles of labels.  ctic is then a power of 10 times 1,2,4, or 5.
            call fitrange(minz,maxz,8,nxfac,xfac,ctic,c1st,clast)
c         write(*,*)ctic,c1st,clast,x1st,cyc,xcyc,in
c Determine a good cyc number xcyc that spans the desired number of
c contours in no more than 8 steps.
            call fitrange(0.,float(in),8,nxfac,xfac,xcyc,x1st,xlast)
            cyc=int(xcyc)
c Find good xtic number that spans ctic in no more than cyc+1 steps.
            call fitrange(c1st+ctic,c1st+2*ctic,cyc+1,
     $           nxfac,xfac,xtic,x1st,xlast)
c         write(*,*)'ctic,c1st,xtic,x1st,cyc,xcyc,in'
c         write(*,*)ctic,c1st,xtic,x1st,cyc,xcyc,in
c Round down, ensuring integers are ok, so that these are now commensurate.
            cyc=int(1.0001*ctic/xtic)
            xtic=ctic/cyc
            in=int((clast-c1st)/xtic)
            x1st=c1st+ctic-cyc*xtic
            cv=maxz-minz
c         write(*,*)ctic,c1st,xtic,x1st,cyc,xcyc,in
         else
c Current fitting of contours and labels code.
            call fitrange(minz,maxz,in,nxfac,xfac,xtic,x1st,xlast)
c Start assuming the label is every contour.
            cyc=1
 201        continue
            if(abs(xlast-x1st)/(cyc*xtic).gt.10)then
c If there would be more than 10 labels. Too many. Find the unit ixtic
               al10=alog10(xtic)
               ixtic=nint(10.**(al10-nint(al10-.5)))
c               write(*,*)xtic,xfac,al10,al10-nint(al10-.5),ixtic
               if(ixtic.eq.4)then
c If ixtic.eq.4, then increase xtic unit to 5.
                  xtic=1.25*xtic
                  goto 201
               endif
c Otherwise, increase the number of contours per label cycle.
               if(cyc.eq.4)then
                  cyc=5
               else
                  cyc=cyc*2
               endif
               goto 201
            endif
            in=int((xlast-x1st)/xtic)
c            write(*,*)xfac,x1st,xtic,cyc
         endif
      else
         in=abs(nc)
         cyc=1+in/8
         c1st=cl(1)
         if(in.gt.1)then 
            cv=abs(cl(in)-cl(1))
            clast=cl(in)
         else
            cv=abs(cl(1))
         endif
      endif
c At this point we need cyc, x1st, xtic, xlast
      if(lclog)then
         if(nc.eq.0)then
c Do logarithmic contour coloring fitting
            call minmax2(z,L,imax,jmax,minz,maxz)
            clast=10.**nint(alog10(maxz)+0.5)
            c1st=10.**nint(alog10(minz)-0.5)
         endif
         if(.not.c1st.gt.0.)then
            write(*,*)'Fixing Log contour fitting error',c1st,clast
            c1st=clast*1.e-4
c            stop
         endif
      endif
      c1stlog=log(c1st)
      cdiflog=1./(log(clast)-c1stlog)
c Coloring.
c      write(*,*)'minz,maxz,c1st,clast',minz,maxz,c1st,clast
      if(icfil.ne.0)then
      id=4
      incolor=igetcolor() 
c      write(*,*)'ncolor=',ncolor,'incolor=',incolor
c      write(*,*)'c1st,clast,ngradcol,itri',c1st,clast,ngradcol,itri
      do j=1,jmax
         y1=j-0.500001
         y2=j+0.500001
         if(j.eq.1)y1=1.001
         if(j.eq.jmax)y2=j-.001
         do i=1,imax
            if(itri.eq.0)then
c Block chunky.
               x1=i-0.500001
               x2=i+0.500001
               if(i.eq.1)x1=1.001
               if(i.eq.imax)x2=i-.001
               if(lclog)then
                  icolor=nint((log(max(z(i,j),c1st))-c1stlog)*(ngradcol
     $                 -1.)*cdiflog)
               else
                  icolor=nint((z(i,j)-c1st)*(ngradcol-1.)/(clast-c1st))
               endif
               if(icolor.gt.ngradcol-1)icolor=ngradcol-1
               if(icolor.lt.0)icolor=0
c               write(*,*)'icolor=',icolor
               call gradcolor(icolor) 
               xd(1)=x1
               yd(1)=y1
               xd(2)=x2
               yd(2)=y1
               xd(3)=x2
               yd(3)=y2
               xd(4)=x1
               yd(4)=y2 
               call mesh2w(xd,yd,id,x,y,l,theconsw) 
               call polyline(xd,yd,id)
               call pathfill() 
            else
c Triangle gradients.
               if(i.lt.imax .and. j.lt.jmax)then
                  xd(1)=i
                  yd(1)=j
                  xd(2)=i+1
                  yd(2)=j
                  xd(3)=i+1
                  yd(3)=j+1
                  xd(4)=i
                  yd(4)=j+1
                  call mesh2w(xd,yd,id,x,y,l,theconsw)
                  zd(1)=z(i,j)
                  zd(2)=z(i+1,j)
                  zd(3)=z(i+1,j+1)
                  zd(4)=z(i,j+1)
c                  write(*,*)c1st,clast,ngradcol
                  call gradquad(xd,yd,zd,zd,
     $                 c1st,clast,0,ngradcol-1,256*istep)
               endif
            endif
         enddo
      enddo
      call color(incolor)
      endif
      if(inocont.eq.0)then
      if(nc.lt.0)then
         labels=.false.
         cyc=1
      else
         labels=.true.
c    stored current size earlier
         call charsize(.01,.01)
c Decide on the format of the label based on the contour range.
         point=2-min(ifix(log10(cv)),2)
      endif
c Contour drawing
      do 1 i=1,in
         if(nc.ne.0)then
            cv=cl(i)
            idel=-1
            if(i.eq.1)idel=1
            cdel=abs(cl(i+idel)-cl(i))
         else
            cv=(x1st+i*xtic)
            cdel=xtic
         endif
         ic=ici
c This gave interesting random values. i1st not initialized.
c        if((mod(i+i1st,cyc).eq.0).and.labels)then
         if((mod(i,cyc).eq.0).and.labels)then
            ipoint=max(point,1-min(ifix(log10(cdel)+1.e-4),2))
            if(ipoint.gt.8)then
               write(str2,'(g12.6)')cv
c this circumlocution to work around an f2c/powerc bug. 
                str1=str2
               width=12
            else
               call fwrite(cv,width,ipoint,str1)
            endif
            if(cv.lt.0)then
               str2=str1
               str1=' '//str2(1:29)
               width=width+1
            endif
         else
            width=0
         endif
         call consgen(z,cv,l,imax,jmax,ppath,xc,yc,ic,x,y,ifcol)
    1 continue
      call charsize(cw,ch)
      endif
      end


c Contour-start searching for contourl.
c Contour labelling is controlled by the common:
c common/lablnc/width,str1(*30) giving the label length and string.
c Calls: Confol, polyline, labeline.
c*************************************************************************
      subroutine consgen(z,cv,l,ixmax,iymax,ppath,xc,yc,ic,x,y,theconsw)
c contour searching routine. 9 Aug 92
      real z(l,iymax),cv,x(1),y(1)
      integer l,ixmax,iymax
      character ppath(l,iymax)
      integer ic
      real xc(ic),yc(ic)
      integer consw,theconsw
c Switch determining the type of plot:
c       0   Equally spaced arrays. Arguments x and y are not used.
c       1   Vectors x and y determine the unequally spaced arrays.
c       2   Arrays x and y determine the arbitrary mesh.
c If the Third byte of theconsw is non-zero it is the gradcolor by
c which to fill a contour that is closed without encountering the bdy.
      integer i,j,k,kk,icc,ii,kk1
c labeling common
      integer width
      character str1*30
      common/lablnc/width,str1
      integer is(0:3),ioffs(0:3),iss
      integer idx(0:4), idy(0:4), id
      data is/0,1,0,1/ioffs/1,0,0,1/
      data idx/0,1,0,-1,0/
      data idy/-1,0,1,0,-1/

c consw is the lowest byte of the consw
      consw=theconsw-256*(theconsw/256)
c fillcolor is the next byte
      ifcolor=(theconsw/65536-256*(theconsw/(256*65536)))
c      write(*,*)'theconsw,ifcolor',theconsw,ifcolor

c Search for starting points and call confol.
      icc=ic
c Initialize the path record.
      do 1 j=1,iymax
         do 2 i=1,ixmax
c           do 21 k=1,4
               ppath(i,j)=char(0)
c   21      continue
    2    continue
    1 continue
c
c Search sides
      do 8 kk1=0,1
         do 6 kk=0,3
            i=is(kk)*(ioffs(kk))+(is(kk)-1)*(ioffs(kk)-1)*(ixmax-1)+1
            j=(1-is(kk))*(ioffs(kk))+is(kk)*ioffs(kk)*(iymax-1)+1
            iss=is(kk)
            do 7 ii=0,is(kk)*ixmax+(1-is(kk))*iymax-2
               if(z(i,j)-cv .ge. 0)then
                  k=kk+kk1
                  if(k.gt.3)k=k-4
                  if(z(i+idx(k),j+idy(k))-cv.lt.0)then
                     itest=ichar(ppath(i,j))/2**k
                     itest=itest-(itest/2)*2
                     if(itest.eq.0)then
c                    if(ichar(ppath(k+1,i,j)).eq.0)then
                        id=k+1
c Found a new starting point.
                        call confol(z,cv,l,ixmax,iymax,
     $                    i,j,id,ppath,xc,yc,ic)
                        call mesh2w(xc,yc,ic,x,y,l,consw)
                        call labeline(xc,yc,ic,str1,width)
                        ic=icc
                     endif
                  endif
               endif
               i=iss+i
               j=(1-iss)+j
    7       continue
    6    continue
    8 continue
c
c Search internal points.
      do 3 j=2,iymax-1
      do 4 i=2,ixmax-1
         if(z(i,j)-cv .ge. 0)then
c This might be a valid starting point. Check adjacent.
            do 5 k=0,3
               if(z(i+idx(k),j+idy(k))-cv.lt.0)then
                  itest=ichar(ppath(i,j))/2**k
                  itest=itest-(itest/2)*2
                  if(itest.eq.0)then
c                 if(ichar(ppath(k+1,i,j)).eq.0)then
                     id=k+1
c Found a new starting point.
                     call confol(z,cv,l,ixmax,iymax,
     $                 i,j,id,ppath,xc,yc,ic)
                     call mesh2w(xc,yc,ic,x,y,l,consw)
                     if(ifcolor.ne.0)then
c                        call getrgbcolor(icol,ired,igreen,iblue)
                        icol=igetcolor()
                        call acgradcolor(ifcolor)
                        call polyline(xc,yc,ic)
                        call pathfill()
                        call color(icol)
                     endif
                     call labeline(xc,yc,ic,str1,width)
                     ic=icc
                  endif
               endif
    5       continue
         endif
    4 continue
    3 continue
      end
c*************************************************************************
      subroutine mesh2w(xc,yc,ic,x,y,l,consw)
c Transform xc,yc(ic) to the mesh x,y, depending on switch consw.
c On Input: xc,yc are the cell (fractional) positions relative to the mesh. 
c x,y are the mesh itself.
c On Ouput: xc,yc are in the units of x and y.  
      integer ic,L,consw
      real xc(ic),yc(ic),x(*),y(*)
      real t1,t2,dx,dy
      integer ix,iy,iv,ipx,ipy,i

      if(consw.eq.0)return
c No scaling.
      if(consw.eq.1)then
c Use vector forms.
         do 1 i=1,abs(ic)
            ix=int(xc(i))
            xc(i)=x(ix)+(x(ix+1)-x(ix))*(xc(i)-ix)
            iy=int(yc(i))
            yc(i)=y(iy)+(y(iy+1)-y(iy))*(yc(i)-iy)
    1    continue
      elseif(consw.eq.2)then
c  Use array forms.
c  Mock up array behaviour: x(l,*),y(l,*)
         do 2 i=1,abs(ic)
            ix=int(xc(i))
            dx=xc(i)-ix
            iy=int(yc(i))
            dy=yc(i)-iy
            iv=ix+(iy-1)*L
            ipx=iv+1
            ipy=iv+L
            t1=x(iv)+(x(iv+1)-x(iv))*dx
            t2=x(ipy)+(x(ipy+1)-x(ipy))*dx
            xc(i)=(1.-dy)*t1+dy*t2
            t1=y(iv)+(y(iv+L)-y(iv))*dy
            t2=y(ipx)+(y(ipx+L)-y(ipx))*dy
            yc(i)=(1.-dx)*t1+dx*t2
    2    continue
      endif
      end
c************************************************************************
      subroutine confol(z,cv,l,ixmax,iymax,initx,inity,initd,
     $     ppath,xc,yc,i)
c Follow a contour of function z, at value cv. Path marking version.
c Inputs: function-array, contour value, dimensions: true,used,used.
      real z(l,iymax),cv
      integer l,ixmax,iymax
c Previous path record. Sets for each dir plotted for this point.
      character ppath(l,iymax)
c Initial point and direction:
      integer initx,inity,initd
c initd=[1,2,3,4] => direction of motion [E,N,W,S]
c Outputs:
      integer i,imax
      real xc(i),yc(i),z1,z2,di
c The interpolated contour and its length. i is equal to max length
c on input, and relevant length on output.
c
      integer idx(0:4), idy(0:4), id,ix,iy,ixn,iyn,ixn2,iyn2
c Direction delta indices idx,idy(id). The side we are crossing is 
c  (ix,iy)-(ix+idx(id-1),iy+idy(id-1)). Next point straight is
c  (ix+idx(id),iy+idy(id)).
      data idx/0,1,0,-1,0/
      data idy/-1,0,1,0,-1/
c Statement functions for interpolation.
c      xint(id)=(ix+idx(id)*(z(ix,iy)-cv)/
c     $      (z(ix,iy)-z(ix+idx(id),iy+idy(id))))
c      yint(id)=(iy+idy(id)*(z(ix,iy)-cv)/
c     $      (z(ix,iy)-z(ix+idx(id),iy+idy(id))))
c
      imax=i
      i=0
      ix=initx
      iy=inity
      id=initd
c Contour always has z1 (left) gt cv and z2 lt cv.
      z1=z(ix,iy)-cv
      z2=z(ix+idx(id-1),iy+idy(id-1))-cv
      if((z1 .lt. 0).or.(z2 .gt.0).or.(id.gt.4).or.(id.lt.1))then
         write(*,*)'CONFOL:Improper input values x,y,idir,z1-cv,z2-cv'
     $        ,ix,iy,id,z1,z2
         return
      endif
c Start of main loop.
    1 continue
c interpolate the next point.
      i=i+1
c      write(*,'('' new ix,iy,id,i'',4i4)')ix,iy,id,i
      call nextpt(z,l,ixmax,iymax,cv,
     $     ix,iy,idx(id-1),idy(id-1),di)
      xc(i)=ix+idx(id-1)*di
      yc(i)=iy+idy(id-1)*di
c This is presumably the place to put the function call to document the
c track followed. At the moment, all we are doing is storing the point
c in xc,yc. 
c Path documentation:
      call acpathdoc(z,cv,l,imax,xc(1),yc(1),i)

      itest=ichar(ppath(ix,iy))/2**(id-1)
      if(itest - (itest/2)*2 .eq. 0)
     $      ppath(ix,iy)=char(ichar(ppath(ix,iy))+2**(id-1))
      if(i.eq.imax)then
         write(*,*)'CONFOL: Contour length exhausted',imax
c Path document end
         call acpathdoc(z,cv,l,imax,0.,0.,1)
         return
      endif
      if((i.gt.1).and.(ix.eq.initx).and.(iy.eq.inity).and.
     $        (id.eq.initd))then
c        write(*,*)'Returned to initial point.',ix,iy
c Path document end
         call acpathdoc(z,cv,l,imax,0.,0.,1)
         return
      endif
c Decide which way to turn
      ixn=ix+idx(id)
      iyn=iy+idy(id)
      if((ixn.lt.1).or.(ixn.gt.ixmax).or.(iyn.lt.1).or.(iyn.gt.iymax))
     $     then
c        write(*,*)'Moved to edge.'
c Path document end
         call acpathdoc(z,cv,l,imax,0.,0.,1)
         return
      endif
      if(z(ixn,iyn)-cv .lt. 0)then
c Turn left
         id=id+1
         if(id.gt.4) id = 1
      else
         ixn2=ixn+idx(id-1)
         iyn2=iyn+idy(id-1)
         if(z(ixn2,iyn2)-cv.lt.0)then
c Go straight
            ix=ixn
            iy=iyn
         else
c Turn Right
            ix=ixn2
            iy=iyn2
            id=id-1
            if(id .lt. 1) id=4
         endif
      endif
      goto 1
      end
c**************************************************************************
c Dummy acpathdoc must be replaced by explicitly linked version if
c path documentation is actually desired.
      subroutine acpathdoc(z,cv,l,imax,xc,yc,i)
c Silence annoying warnings
      real z(*)
      r=z(1)
      r=cv
      j=l
      j=imax
      r=xc
      r=yc
      j=i
      end
c**************************************************************************
      subroutine  nextpt(z,l,ixm,iym,cv,ix,iy,idx,idy,di)
c Call the parabolic interpolation but fix boundaries to give zero
c second derivatives.
c Inputs: array z(ixm/l,iym); contour value cv; point ix,iy; direction idx,idy.
c Output: di the fraction of idx,idy to be moved to the point where z=cv.
c x+idx,y+idy and x,y must be in the mesh. Other points are checked. 
      integer l,ixm,iym,ix,iy,idx,idy
      real cv,di
      real z(l,iym)
      real ym1,y0,y1,y2
      integer i,j
      i=ix-idx
      j=iy-idy
      if(i.lt.1 .or. i.gt.ixm .or. j.lt.1 .or. j.gt.iym) then
         ym1=2.*z(ix,iy)-z(ix+idx,iy+idy)
      else
         ym1=z(i,j)
      endif
      y0=z(ix,iy)
      y1=z(ix+idx,iy+idy)
      i=ix+2*idx
      j=iy+2*idy
      if(i.lt.1 .or. i.gt.ixm .or. j.lt.1 .or. j.gt.iym) then
         y2=2.*z(ix+idx,iy+idy)-z(ix,iy)
      else
         y2=z(i,j)
      endif
      call intpara(ym1,y0,y1,y2,cv,di)
      end
c**************************************************************************
c Parabolic interpolation.
      subroutine intpara(ym1,y0,y1,y2,yv,di)
      real ym1,y0,y1,y2,yv,di
c  On entry:
c       ym1-y2  the four points to be interpolated
c       yv      value to be interpolated
c  On exit:
c       di      fractional index increment at which y(di)=yv.
      real g,yp,dy,dy1,gpyp
      dy=yv-y0
      dy1=y1-yv
      if(dy1.lt.0)then
         if(dy.le.0)goto 10
      else
         if(dy.gt.0)goto 10
      endif
      write(*,*)'INTPARA INPUT ERROR: Value not bracketed'
      write(*,'(a,g12.6,a,g12.6,a,i4,/)')
     $     ' y0=',y0,' y1=',y1
      return
   10 yp=y1-y0
      g=0.25*(-y2+y1+y0-ym1)
      gpyp=g+yp
      di=(abs(gpyp)+sqrt(gpyp*gpyp-4.*g*dy))
      if(di.ne.0.) di=abs(2*dy/di)
      end
c**************************************************************************
c Put an arbitrary label in the str1 contour labeling string(30 max).
      subroutine conlabel(str2)
      character*(*) str2
      integer width
      character str1*30
      common/lablnc/width,str1
      str1=str2
c cite the width without the closing 0.
      width=len(str2)-1
      end
c**************************************************************************
      logical function getconlog()
c Return the value of lclog
      include 'plotcom.h'
      getconlog=lclog
      end
c**************************************************************************
      subroutine setconlog(logval)
c Set the value of lclog that determines whether color contouring is done
c on logarithmic scale (.true.) or linear (.false.).
      logical logval
      include 'plotcom.h'
      lclog=logval
      end
c**************************************************************************
c Create a legend for the gradient color used in a contour plot
      subroutine gradlegend(c1,c2,x1,y1,x2,y2,colwidth,lpara)
      real c1,c2
      real x1,y1,x2,y2
      logical lpara
      real colwidth
c gradient goes from c1 to c2.
c Bar position in fractions of the plot box given by x1,y1;x2,y2
c Width of color bar colwidth times plot box width.
      logical laxlog
      integer ngradcol
c      parameter (ngradcol=256)
      parameter (ngradcol=240)
      real xd(4),yd(4)
      include 'plotcom.h'
c      data laxlog/.false./

      laxlog=lclog
      if(c1.eq.c2)then
         if(c1.eq.0.)then
            c2=1.
         else
            c2=1.001*c1
         endif
      endif

      th=atan2(y2-y1,x2-x1)
c      r=sqrt((y2-y1)**2+(x2-x1)**2)
      ct=cos(th)
      st=sin(th)
      if(colwidth.eq.0.)then
         cpb=0.01*(ct*(wxmax-wxmin)+st*(wymax-wymin))
      else
         cpb=colwidth*(ct*(wxmax-wxmin)+st*(wymax-wymin))
      endif
      incolor=igetcolor()

      xb=wxmin*(1.-x1) + wxmax*x1
      yb=wymin*(1.-y1) + wymax*y1
      xgmin=wx2nx(xb)
      ygmin=wy2ny(yb)

      do i=1,ngradcol
c         xc=x1+(x2-x1)*(i-1)/float(ngradcol-1)
         xc=x1+(x2-x1)*i/float(ngradcol-1)
         xc=wxmin*(1.-xc)+ wxmax*xc
         yc=y1+(y2-y1)*i/float(ngradcol-1)
         yc=wymin*(1.-yc)+ wymax*yc
c         write(*,*)xb,yb
         call gradcolor(i-1)
         xd(1)=xb
         xd(2)=xc
         xd(3)=xc+cpb*(-st)
         xd(4)=xb+cpb*(-st)

         yd(1)=yb
         yd(2)=yc
         yd(3)=yc+cpb*(ct)
         yd(4)=yb+cpb*(ct)

         call polyline(xd,yd,4)
         call pathfill()
         xb=xc
         yb=yc
      enddo
      call color(incolor)
      ngpow=0
      first=0.
      delta=0.
      xgmax=wx2nx(xb)
      ygmax=wy2ny(yb)
      call gaxis(c1,c2,ngpow,first,delta,
     $     xgmin,xgmax,ygmin,ygmax,lpara,laxlog)
      end








