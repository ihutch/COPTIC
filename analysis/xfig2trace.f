c***********************************************************************
      subroutine xfig2trace(np,nt,xtraces,ytraces,il,nl,xyia,filename)

c Read up to nt traces each of length up to np (IN)
c from xfig format plot file filename (IN).
c Return il traces of lengths nl(il) in xtraces, ytraces (OUT).
c Use the optionally updated box xyia(4) (INOUT) to scale them.

c Read **polylines** in filename as traces. 

c Use the last **rectangle or closed quadrilateral** in file as the
c corners of the plotted region, corresponding to xmin, xmax, ymin, ymax
c = xyia(8).

c If a comment exists in the xfig file starting "Box:" 
c (entered by doing an xfig edit of an object and typing it in)
c then read 8 values x,y,x,y,x,y,x,y from this comment into xyia. 
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
      real xq(4),yq(4)
      character*(*) filename

      real xrect(5),yrect(5)
      character*100 line
      real pvalues(16)

c Defaults use xyia.
      xq(1)=0.
      xq(2)=0.

c      write(*,'(a,4f8.3)')'xmin,xmax,ymin,ymax'
c     $     ,xyia(1),xyia(2),xyia(3),xyia(4)
      il=1
      open(10,file=filename,status='old',err=304)
      do i=1,200
         read(10,'(a)',end=302,err=301)line
c         write(*,'(a)')line
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
         elseif(line(1:3).eq.'2 2'.or.line(1:3).eq.'2 3')then
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
c Old min max approach
            read(line(7:),*,err=303,end=303)(xyia(k),k=1,4)
            write(*,*)'xfig2trace: xyia',(xyia(k),k=1,4)
         elseif(line(1:7).eq.'# Quad:')then
c            read(line(8:),*,err=303,end=303)
            read(line(8:),*)(xq(k),yq(k),k=1,4)
            write(*,*)'xfig2trace: xq,yq',(xq(k),yq(k),k=1,4)
         elseif(line(1:8).eq.'# Scale:')then
            read(line(9:),*,err=303,end=303)scale
            write(*,*)'xfig2trace: scale',scale
            do k=3,4
               xyia(k)=xyia(k)*scale
            enddo
         elseif(line(1:9).eq.'# XScale:')then
            read(line(10:),*,err=303,end=303)scale
            write(*,*)'xfig2trace: Xscale',scale
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
c      write(*,*)'Completed reading file.'
      il=il-1
c      write(*,*)'il',il,(nl(k),k=1,il)
c---------------------------------------------------------------
      if(xq(2).eq.xq(1))then
c Old transform useful only for rectangle.
c Now we transform from xfig coordinates to plot world coordinates.
         write(*,*)'Using rectangular box transform', xyia,xrect,yrect
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
      else
c---------------------------------------------------------------
c Transform for general quadrilateral.
      
         do i=1,il
            do j=1,nl(i)
               x=xtraces(j,i)
               y=ytraces(j,i)
c Perpendicular distance from line 1-2
               p1=(-(x-xrect(1))*(yrect(2)-yrect(1))+(y-yrect(1))
     $              *(xrect(2)-xrect(1)))/sqrt((yrect(2)-yrect(1))**2
     $              +(xrect(2)-xrect(1))**2)
               p2=(-(x-xrect(2))*(yrect(3)-yrect(2))+(y-yrect(2))
     $              *(xrect(3)-xrect(2)))/sqrt((yrect(3)-yrect(2))**2
     $              +(xrect(3)-xrect(2))**2)
c Perpendicular distance from line 3-4
               p3=(-(x-xrect(3))*(yrect(4)-yrect(3))+(y-yrect(3))
     $              *(xrect(4)-xrect(3)))/sqrt((yrect(4)-yrect(3))**2
     $              +(xrect(4)-xrect(3))**2)
               p4=(-(x-xrect(4))*(yrect(1)-yrect(4))+(y-yrect(4))
     $              *(xrect(1)-xrect(4)))/sqrt((yrect(1)-yrect(4))**2
     $              +(xrect(1)-xrect(4))**2)
               f1=p1/(p1+p3)
               g1=1-f1
               f2=p2/(p2+p4)
               g2=1-f2
               xtraces(j,i)=g1*f2*xq(1)+g1*g2*xq(2)+f1*g2*xq(3)
     $              +f1*f2*xq(4)
               ytraces(j,i)=g1*f2*yq(1)+g1*g2*yq(2)+f1*g2*yq(3)
     $              +f1*f2*yq(4)
c               write(*,*)p1,p2,p3,p4
c               write(*,*)f1,g1,f2,g2
c               write(*,*)i,j,xtraces(j,i),ytraces(j,i)
            enddo
         enddo
c Kludge up the xyia
         xyia(1)=xq(1)
         xyia(2)=xq(2)
         xyia(3)=yq(1)
         xyia(4)=yq(3)

      endif
      return
      
 304  write(*,*)'xfig2trace: Could not open: ',filename(1:50)
      il=-1

      end
c****************************************************************
