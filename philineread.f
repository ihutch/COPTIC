      program philineread

      include 'examdecl.f'

      character*100 filenames(na_m)
      character*50 string
      parameter (nfx=20)
      integer ild,ilinechoice(ndims_mesh,nfx),ip(ndims_mesh)
      real philine(na_m),xline(na_m)
      real darray(nfx),pmax(nfx),punscale(nfx),rp(nfx)

c philineread commons
      logical lrange,lwrite
      integer iover
      character*(100)overfile
      common /linecom/xmin,xmax,ymin,ymax,lrange,lwrite
     $     ,iover,overfile

c xfig2trace parameters.
      parameter (np=200,nt=10)
      real xtraces(np,nt), ytraces(np,nt)
      integer nl(nt)
      real xyia(4)
c      character*100 xfigfile
c 
c silence warnings:
      zp(1,1,1)=0.
      fluxfilename=' '
      xmin=0.
      ymin=0.
      xmax=0.
      ymax=0.
      lwrite=.false.
      iover=0
      ild=3

c      write(*,*)nf,idl
      nf=nfx
      call lineargs(filenames,nf,ild,ilinechoice,rp)
      write(*,*)'Filenames',nf

      nplot=0
      do inm=1,nf
         phifilename=filenames(inm)
c         write(*,*)ifull,iuds
         call array3read(phifilename,ifull,iuds,u,ierr)
         if(ierr.eq.1)stop
         do i=1,ndims_mesh
            if(iuds(i).gt.na_m) then
               write(*,*)'Data too large for na_m',na_m,iuds(i)
               stop
            endif
         enddo

c      call sliceGweb(ifull,iuds,u,na_m,zp,
c     $     ixnp,xn,ifix,'potential:'//'!Ay!@')

c Select the lineout into the plotting arrays.      
         if(ild.ne.0)then
            do k=1,ndims_mesh-1
               ip(mod(ild+k-1,ndims_mesh)+1)=
     $              ilinechoice(mod(ild+k-1,ndims_mesh)+1,inm)
c               write(*,*)k,ilinechoice(mod(ild+k-1,ndims_mesh)+1,inm)
            enddo
c            write(*,*)'ild,ip(ild)',ild,ip(ild)
            do i=1,iuds(ild)
               ip(ild)=i
               xline(i)=xn(ixnp(ild)+i)/debyelen
               philine(i)=u(ip(1),ip(2),ip(3))
     $              /abs(rp(inm)*phip/debyelen)
               if(lwrite)write(*,'(2f10.5)')xline(i),philine(i)
            enddo
            call minmax(philine,iuds(ild),pmin,pa)
            pmax(inm)=pa
            punscale(inm)=pmax(inm)*abs(rp(inm)*phip/debyelen)
            darray(inm)=abs(rp(inm)*phip/debyelen)

            call winset(.true.)
            call pfset(3)
            if(inm.eq.1)then
               if(lrange)then
                  if(xmin-xmax.eq.0.)then
                     call minmax(xline,iuds(ild),xmin,xmax)
                  endif
                  if(ymin-ymax.eq.0.)then
                     call minmax(philine,iuds(ild),ymin,ymax)
                  endif
                  call pltinit(xmin,xmax,ymin,ymax)
                  call axis()
                  call axlabels('z/!Al!@',
     $                 '!Af!@/(|!Af!@!dp!d|r!dp!d/!Al!@)')
                  call winset(.true.)
                  call polyline(xline,philine,iuds(ild))
               else
                  call autoplot(xline,philine,iuds(ild))
               endif
               call axis2()
            else
               call color(inm)
c               call iwrite(inm,iwd,string)
c               call labeline(xline,philine,iuds(ild),string,iwd)
               call dashset(inm)
               call polyline(xline,philine,iuds(ild))
            endif
            string=' !Af!@!dp!d='
            call fwrite(phip,iwd,2,string(lentrim(string)+1:))
            string(lentrim(string):)='@'
            call fwrite(rp(inm)/debyelen,iwd,2,
     $           string(lentrim(string)+1:))
            call legendline(.5,(.01+inm*.05),0,
     $           string(1:lentrim(string)))
            nplot=nplot+1
         endif

         write(*,'(a,3i4,$)')'On grid',iuds
         write(*,*)(',',xn(ixnp(kk)+1),xn(ixnp(kk+1)),kk=1,3)
     $        ,ip
c         write(*,*)ild,ilinechoice
      enddo
      if(iover.gt.0)then
c Overplot traces from specified file.
         write(*,*)iover,' Overplot ',overfile(1:40)
         call xfig2trace(np,nt,xtraces,ytraces,il,nl,xyia,overfile)
         write(*,*)'Return from xfig2',il,(nl(k),k=1,il),xyia
            if(il.gt.0)then
               do k=1,il
                  call dashset(0)
                  call color(13)
                  call polyline(xtraces(1,k),ytraces(1,k),nl(k))
               enddo
            endif
         endif
      if(nplot.gt.0)then
         call pltend()
         call dashset(0)
         string=' r!dp!d/!Al!@='
         call lautomark(darray,pmax,nplot,.true.,.true.,0)
         imark=1
         do k=1,nplot
            call fwrite(rp(k)/debyelen,iwd,2,string(15:))
            call color(imark)
            if(k.eq.1)then
               call legendline(.1,.02+.05*imark,imark,string)
            elseif(rp(k).ne.rp(k-1))then
               imark=imark+1
               call color(imark)
               call legendline(.1,.02+.05*imark,imark,string)
            endif
            call polymark(darray(k),pmax(k),1,imark)
         enddo
         call color(15)
         call axlabels('|!Af!@!dp!d|r!dp!d/!Al!@'
     $        //'!A ~ !@Q/4!Ape!@!d0!d!Al!@',
     $        '!Af!@!dmax!d/(|!Af!@!dp!d|r!dp!d/!Al!@)')
c         call vecw(0.04,3.,0)
c         call vecw(1.,.12,1)
c         call vecw(0.01,2.,0)
c         call vecw(0.1,2.,1)
         call pltend()
         call lautomark(darray,punscale,nplot,.true.,.true.,0)
         call color(imark)
         imark=1
         do k=1,nplot
            call fwrite(rp(k)/debyelen,iwd,2,string(15:))
            call color(imark)
            if(k.eq.1)then
               call legendline(.1,.02+.05*imark,imark,string)
            elseif(rp(k).ne.rp(k-1))then
               imark=imark+1
               call color(imark)
               call legendline(.1,.02+.05*imark,imark,string)
            endif
            call polymark(darray(k),punscale(k),1,imark)
         enddo
         call color(15)
         write(*,*)'punscale',(punscale(k),k=1,nplot)
         call axlabels('|!Af!@!dp!d|r!dp!d/!Al!@'
     $        //'!A ~ !@Q/4!Ape!@!d0!d!Al!@',
     $        '!Af!@!dmax!d')
         call pltend()
      endif


      end

c*************************************************************
      subroutine lineargs(filenames,nf,ild,ilinechoice,rp)
c Deal with command line arguments.

c Non-switch arguments are the names of files to read from. For each of
c those, return the name in filename, the choice of line in ilinechoice,
c and the radius of object in rp. 
c On entry nf is array dimension. (IN)
c On exit nf is the number of files read. (OUT)
c Switch arguments set -x(min,max) -y(min,max) -l(linechoice for the
c subsequent files) -r radius (subsequent) -w writing to true. 
c -o overplot file (xfig with scaling). 

c ild is the dimension that is fixed, and must be the same for all files.
c the logic will break if it is changed by -l in the middle.

      include 'examdecl.f'
      character*100 filenames(na_m)
      real rp(nf)
      integer nf,ild,ilinechoice(ndims_mesh,nf)
      integer idj(ndims_mesh)

      logical lrange,lwrite
      integer iover
      character*(100) overfile
      common /linecom/xmin,xmax,ymin,ymax,lrange,lwrite
     $     ,iover,overfile

      ifull(1)=na_i
      ifull(2)=na_j
      ifull(3)=na_k
c Passed in array dimension
      nfx=nf

c Defaults and silence warnings.
      phifilename=' '
      fluxfilename=' '
      nf=1
      zp(1,1,1)=0.
      lrange=.false.
      do id=1,ndims_mesh
         idj(id)=0
         ilinechoice(id,1)=0
      enddo
      rread=1.
      iover=0

c Deal with arguments
      if(iargc().eq.0) goto 201
      do i=1,iargc()
         call getarg(i,argument)
         if(argument(1:1).eq.'-')then
            if(argument(1:2).eq.'-l')then
               read(argument(3:),*,end=11)ild,(idj(k),k=1,ndims_mesh-1)
 11            continue
            endif
            if(argument(1:2).eq.'-y')then
               read(argument(3:),*,end=12)ymin,ymax
               lrange=.true.
 12            continue
            endif
            if(argument(1:2).eq.'-x')then
               read(argument(3:),*,end=13)xmin,xmax
               lrange=.true.
 13            continue
            endif
            if(argument(1:2).eq.'-o')then
               read(argument(3:),*,end=14,err=14)overfile
               iover=1
 14            continue
            endif
            if(argument(1:2).eq.'-w')lwrite=.true.

            if(argument(1:13).eq.'--objfilename')
     $           read(argument(14:),'(a)',err=201)objfilename
            if(argument(1:2).eq.'-r')
     $           read(argument(3:),'(f8.4)',err=201)rread
            if(argument(1:2).eq.'-h')goto 203
            if(argument(1:2).eq.'-?')goto 203
            if(argument(1:2).eq.'-f')goto 204
         else
 204        read(argument(1:),'(a)',err=201)phifilename
c               write(*,*)ild, idj
            do k=1,ndims_mesh-1
               ilinechoice(mod(ild+k-1,ndims_mesh)+1,nf)=idj(k)
c                  write(*,*)'nf,k,idj(k)',nf,k,idj(k)
            enddo
            filenames(nf)(1:)=phifilename
            rp(nf)=rread
            if(nf.eq.nfx)then
               write(*,*)'Exhausted file dimension',nf
               return
            endif
            nf=nf+1
         endif
      enddo
      nf=nf-1

c      write(*,*)ild,
c     $     (ilinechoice(2,i),' ',rp(i),' ',filenames(i)(1:30),i=1,nf)
      goto 202
c------------------------------------------------------------
c Help text
 201  continue
      write(*,*)'=====Error reading command line argument'
 203  continue
 301  format(a,3i5)
 302  format(a,3f8.3)
      write(*,301)'Usage: philineread [switches] <phifile>'//
     $     ' [<phifile2> ...]'
c      write(*,301)' --objfile<filename>  set name of object data file.'
c     $     //' [ccpicgeom.dat'
      write(*,301)' -l<idim>,<irow1>,<irow2> set fixed dimension [',ild
      write(*,301)'    and row position in other dimensions. ['
     $     ,(ilinechoice(k,1),k=1,ndims_mesh)
      write(*,301)' -y<min>,<max>   -x<min>,<max>  set plot ranges'
      write(*,301)' -o<figfile> overplot traces using xfig2trace.'
      write(*,302)' -r<r> set radius [',rread
      write(*,301)' -h -?   Print usage.'
      call exit(0)
 202  continue
      if(lentrim(partfilename).lt.5)goto 203
      end
c*****************************************************************

c***********************************************************************
      subroutine xfig2trace(np,nt,xtraces,ytraces,il,nl,xyia,filename)

c Read up to nt traces each of length up to np (IN)
c from xfig format plot file filename (IN).
c Return il traces of lengths nl(il) in xtraces, ytraces (OUT).
c Use the optionally updated box xyia(4) (INOUT) to scale them.

c Read polylines in filename as traces. 
c Use the last rectangle in file as the corners of the plotted region,
c corresponding to xmin, xmax, ymin, ymax = xyia(4).

c If a comment exists in the xfig file starting "Box:" 
c (entered by doing an xfig edit of an object and typing it in)
c then read 4 values from this comment into xyia. 
c Otherwise use input values.

c For each comment starting "Scale:" read 1 value 
c and multiply the xyia y-values by it.
c For each comment starting "XScale:" read 1 value 
c and multiply the xyia x-values by it.
c (Must occur *after* any Box: comment to be effective. Can be 
c ensured by using the same xfig comment box with a new line.)

c Abnormal returns are indicated by il. il>=1000 (too many lines)
c il=-1 (no file) il=0 (no traces) il=-2 (box reading error).

      real xtraces(np,nt),ytraces(np,nt)
      integer il,nl(nt)
      real xyia(4)
      character*(*) filename

      real xrect(5),yrect(5)
      character*100 line
      real pvalues(16)

c      write(*,'(a,4f8.3)')'xmin,xmax,ymin,ymax'
c     $     ,xyia(1),xyia(2),xyia(3),xyia(4)
      il=1
      open(10,file=filename,status='old',err=304)
      do i=1,200
         read(10,'(a)',end=302,err=301)line
c         write(*,'(i4,a)')i,line
         if(line(1:3).eq.'2 1')then
            read(line,*)pvalues
c            write(*,'(a,16f4.0)')'Polyline',pvalues
            nl(il)=pvalues(16)
            if(nl(il).gt.np)then
               write(*,*)'Polyline ',il,' too long:',nl(il)
               nl(il)=np
            endif
            read(10,*)(xtraces(j,il),ytraces(j,il),j=1,nl(il))
c            write(*,'(10f7.1)')(xtraces(j,il),j=1,nl(il))
c            write(*,'(10f7.1)')(ytraces(j,il),j=1,nl(il))
            il=il+1
         elseif(line(1:3).eq.'2 2')then
            read(line,*)pvalues
c            write(*,'(a,16f4.0)')'Rectangle',pvalues
            nl(il)=pvalues(16)
            if(nl(il).ne.5)then
               write(*,*)'Error. Rectangle has not 5 points',nl(il)
               stop
            endif
            irect=il
            read(10,*)(xrect(j),yrect(j),j=1,nl(il))            
c            write(*,'(10f7.1)')(xrect(j),j=1,nl(il))
c            write(*,'(10f7.1)')(yrect(j),j=1,nl(il))
         elseif(line(1:6).eq.'# Box:')then
            read(line(7:),*,err=303,end=303)(xyia(k),k=1,4)
            write(*,*)'xfig2trace: xyia',(xyia(k),k=1,4)
         elseif(line(1:8).eq.'# Scale:')then
            read(line(9:),*,err=303,end=303)scale
            write(*,*)'xfig2trace: scale',scale
            do k=3,4
               xyia(k)=xyia(k)*scale
            enddo
         elseif(line(1:9).eq.'# XScale:')then
            read(line(10:),*,err=303,end=303)scale
            write(*,*)'xfig2trace: XScale',scale
            do k=1,2
               xyia(k)=xyia(k)*scale
            enddo
         endif
         if(il.gt.nt)goto 300
      enddo
 300  write(*,*)'Too many lines to read',nt
      il=il*1000
      return
 301  write(*,*)'xfig2trace Error reading line ',line
      il=il-1
      return
 303  write(*,*)'xfig2trace Error reading Box/Scale comment',line
      il=-2
      return
 302  continue
      write(*,*)'Completed reading file.'
      il=il-1
c      write(*,*)'il',il,(nl(k),k=1,il)
c Now we transform from xfig coordinates to plot world coordinates.
      xscale=(xyia(2)-xyia(1))/(xrect(2)-xrect(1))
      yscale=(xyia(4)-xyia(3))/(yrect(1)-yrect(4))
      xzero=-(xyia(2)*xrect(1)-xyia(1)*xrect(2))/(xrect(2)-xrect(1))
      yzero=-(xyia(4)*yrect(4)-xyia(3)*yrect(1))/(yrect(1)-yrect(4))

c      write(*,*)'Scalings ',xscale,xzero,yscale,yzero

c Transform to world coordinates.
      do i=1,il
         do j=1,nl(i)
            xtraces(j,i)=xtraces(j,i)*xscale+xzero
            ytraces(j,i)=ytraces(j,i)*yscale+yzero
         enddo
      enddo
      return
      
 304  write(*,*)'xfig2trace: Could not open: ',filename(1:50)
      il=-1

      end
c****************************************************************
