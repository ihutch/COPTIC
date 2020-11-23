      program plottest
c
c Illustrative test program that exercises most low level routines.
c Various working variables.
      integer length
      parameter (length=100)
      real x(length),y(length),err(length),ym(length)
      integer i,irows,icolumns,itype,ipcolor(length)

c Make test arrays.
      do 2 i=1,length
       x(i)=float(i)*1.
       y(i)=1.*(1.4+sin(0.2*x(i)))
       err(i)=0.4*sin(x(i)*sin(x(i)))
       ym(i)=y(i)-0.5*err(i)
         ipcolor(i)=1+int((i-1.)/(length-1.)*240)
    2 continue

c Plot 1. Simplest one-call plot.
      call autoplot(x,y,length)
c If needed other calls can overwrite before terminating via:
      call pltend()

c Plot 2. A more complicated plot to illustrate different possibilities.
c Initialize the plot:
      call pltinit(1.,float(length),0.,2.)
c Draw the axes
      call axis()
c Label them if you like.
      call axlabels('X-axis','Y-axis')
c And put a title on the box
      call boxtitle('Test Plot')
c Example of a color set, 0:invisible, 1-7:dim, 8-15 bright:
      call color(12)
c Drawing a line exceeding the screen size, is not recommended but
c is safe because clipped automatically:
      call polyline(x,y,length)
c Return to default, bright black(!).
      call color(15)
c Set software truncation to a specified normal-units region;
      call truncf(.5,.85,.3,.5)
c Draw the line again. Only sections in the truncation region are overdrawn.
      call polyline(x,y,length)
c Switch off truncation.
      call truncf(0.,0.,0.,0.)
      call pltend

c Plot 3.
c Use the built in response facility by calling with negative switch.
c This will prompt for plotting to file.
c      call pfset(-3)
      call pfset(3)
c Set to dashed line plotting, only polylines are dashed:
      call dashset(2)
c Do a log autoplot of the arrays. x logarithmic, y linear.
      call lautoplot(x,y,length,.true.,.false.)
c Label the axes using different fonts.
      call axlabels('X-axis !Alabel!@ back','Y-!Baxis!A label')
c Overplot markers (6:stars) on the middle half of the points,
c not dashed because doesn't use polyline.
      call polymark(x(length/4),y(length/4),length/2,6)
c Set back to solid lines:
      call dashset(0)
c To illustrate a more general drawstring, draw a title, justified centered.
c (All strings drawn in normal units.)
c Make the characters .02 in size. (width,height)
      call charsize(.02,.02)
      call jdrwstr(0.65,0.72,' TITLE OF PLOT',0.)
c Restore default size.
      call charsize(0.,0.)
      call pltend()
c Switch off plotting to file.
c      call pfset(0)
c
c Plot 4 Simple log-log plot.
c But with reversed tics and axes
      call ticrev
c Make the axis point at the top right. Arguments are fractions of axregion.
      call axptset(1.,1.)
      call dashset(4)
      call lautoplot(x,y,length,.true.,.true.)
      call dashset(0)
      call axptset(0.,0.)
      call ticrev
c Draw additional axes in standard position, illustrating first,delta.
      call xaxis(5.,1.)
      call yaxis(3.,10.)
      call pltend
c
      call togminor()
c Plot 5 Simplest automatic scatter plot
c Illustrating the use of a general character as the marker.
      call automark(x,y,length,ichar('m'))
c Overplot some error bars.
      call polyerrs(x,y,err,length,.5,1.)
      call pltend
c
c Plot 6 Multiple Frame Plot
c (nrows,ncolumns,multype: x,y plot space from bits 0 and 1)
c      write(*,*)' Enter nrows and ncolumns,itype'
c      read(*,*)irows,icolumns,itype
        irows=3
        icolumns=4
        itype=3
      call multiframe(irows,icolumns,itype)
      do 10 i=1,irows*icolumns
         if((i/2)*2.ne.i)then
            call autoplot(x,y,length)
         else
            call autoplot(x,ym,length)
         endif
   10      continue
        call multiframe(0,0,0)
       call pltend

c Plot 7 Drawing a line whose color (gradient) encodes other array info
c blue-red-green gradient:
        call blueredgreenwhite
        call autoinit(x,y,length)
        call axis
        call polycolorline(x,y,length,ipcolor)
        call pltend
        call color(15)

      end

