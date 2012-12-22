/* Xwindow driver for accis plotting */
/* Fortran callable (f2c) routines */
/* ********************************************************************** */
/* Fix for Sun to put command line arguments in right place. */
int f__xargc;
char **f__xargv;
void f77_init(int *argc_ptr, char **argv_ptr, char **envp_ptr)
{
    f__xargc = *argc_ptr;
    f__xargv = argv_ptr;
}
/*
Refreshing version.
*/
#include <stdio.h>
#include  <X11/StringDefs.h> 
#include  <X11/Intrinsic.h> 
#include  <X11/Core.h> 

/* Globals for these routines*/
Widget accis_wshell;
Widget accis_drawing;
Display *accis_display;
Drawable accis_window;
Pixmap accis_pixmap;

GC accis_gc;
int accis_depth;
Colormap accis_colormap;

struct Screen_Size {  
    Dimension width; /* unsigned short */
    Dimension height;
  };
    static struct Screen_Size s_s;
/* Static maximum number of points in path */
#define accis_path_max 4000
XPoint accis_path[accis_path_max];
int accis_pathlen=0;

/* 256 color gradient globals */
/* Use only 240 colors so that 16 are left for the 16 colors and you can
   make gifs out of the resultant without a problem.*/
#define a_gradPixno 240
unsigned long a_gradPix[a_gradPixno];
int a_grad_inited=0;
int a_gradred[a_gradPixno];
int a_gradgreen[a_gradPixno];
int a_gradblue[a_gradPixno];

/* 16-color model globals. Pixels for truecolor display. */
#define a_maxPixels 16
unsigned long accis_pixels[a_maxPixels]={
  16777215,
       139,
   3050327,
   2142890,
  11674146,
   9699539,
  16032864,
   7372944,
  12500670,
       255,
     65280,
     65535,
  16711680,
  16711935,
  16776960,
         0,
};

char *accis_colornames[a_maxPixels]=
{
  "White",
  "DarkBlue",
  "SeaGreen",
  "LightSeaGreen",
  "Firebrick",
  "DarkViolet",
  "SandyBrown",
  "SlateGrey",
  "Grey",
  "Blue",
  "Green",
  "Cyan",
  "Red",
  "Magenta",
  "Yellow",
  "Black"
};



/* Subroutine */ 
int svga_(scrxpix, scrypix, vmode, ncolor)
short *scrxpix, *scrypix, *vmode, *ncolor;
{  
  static int second=0;
#ifndef IFC 
  extern int f__xargc;
  extern char **f__xargv;
#else
  int f__xargc=0;
  char **f__xargv;
#endif
  int svga_argc=0;
  char *svga_argv[1]; 
  static int n;
  static Arg wargs[10];
  int theDepth;
  Colormap theColormap;
  extern void config_handler();
  extern void accis_refresh();
  XWindowAttributes attributes_return;

  if(second == 0){
    accis_wshell = XtInitialize("accis","Accis", NULL, 0,
      &f__xargc, f__xargv);  /* This refers to the f2c command line args.
			It may only work with f2c, therefore. */
      /* &svga_argc, svga_argv); alternate */

    accis_drawing = XtCreateManagedWidget("drawing",coreWidgetClass,
		 accis_wshell, NULL, 0);

    *vmode=88;
    *ncolor=15;
    /* Set up a default size of the drawing window widget. 
       This is overruled by Accis*geometry resources setting wshell size.
       */
    n = 0;
    XtSetArg(wargs[n], XtNheight, 480); n++;
    XtSetArg(wargs[n], XtNwidth, 640); n++;
    XtSetValues(accis_drawing, wargs, n);
    XtRealizeWidget(accis_wshell);
    accis_display = XtDisplay(accis_drawing);
    accis_window = XtWindow(accis_drawing);
    accis_gc = XCreateGC(accis_display, accis_window, 0, NULL);
    accis_depth=DefaultDepth(accis_display,0);
    accis_colormap=DefaultColormap(accis_display,0);
    initDefaultColors();
    /* Leave setup for resizing. */
    n = 0;
    XtSetArg(wargs[n], XtNheight, &s_s.height); n++;
    XtSetArg(wargs[n], XtNwidth, &s_s.width); n++;
    /* Pixmap setup */
    XtGetValues(accis_wshell, wargs, n);
    accis_pixmap=XCreatePixmap(accis_display,accis_window,
			       s_s.width,s_s.height,accis_depth);
    XtAddEventHandler(accis_drawing,ExposureMask,FALSE,
		      accis_refresh,NULL);
    XSelectInput(accis_display,accis_window,
		 KeyPressMask | ExposureMask | ButtonPress
		 | FocusChangeMask | EnterWindowMask |
		 ButtonReleaseMask | ButtonMotionMask ); /* Tell events */
    second++;
  }else{
    /* Need to have this to make the get correct. */
    ManageEvents();
  } 
  /* This doesn't need a configuration event handler to work correctly. */
  XtGetValues(accis_wshell, wargs, n);
  *scrxpix=s_s.width;
  *scrypix=s_s.height;
  XFlush(accis_display);
  XClearWindow(accis_display,accis_window);
  XSetForeground(accis_display,accis_gc,accis_pixels[0]);  
  /* Clear the window to background color. VMS doesn't do correctly.
     But also this seems to fix the bad match error. */
  XFillRectangle(accis_display,accis_window,accis_gc,0,0,
		 s_s.width,s_s.height);
  /* Clear the pixmap  */
  XFillRectangle(accis_display,accis_pixmap,accis_gc,0,0,
		 s_s.width,s_s.height);
  XSetForeground(accis_display,accis_gc,accis_pixels[15]); 
  XRaiseWindow(accis_display, accis_window);
  XGetWindowAttributes (accis_display, accis_window,&attributes_return);
  /* Then set the focus. We should not get a bad match that way*/
  if(attributes_return.map_state == IsViewable){
    /* XSetInputFocus(accis_display, accis_window, RevertToPointerRoot,
       CurrentTime);
    Sometimes it is easier not to do this. */
  }
  return 0;
}

/* ******************************************************************** */
initDefaultColors()
{
  XColor theRGBColor;
  XColor theHardColor;
  int status,truecolor;
  int i;
  /*See if this is a sensible 24 bit display or not */
  theRGBColor.red=255*256;
  theRGBColor.green=127*256;
  theRGBColor.blue=2*256;
  XAllocColor(accis_display,accis_colormap,&theRGBColor);
  if(theRGBColor.pixel==((255*256+127)*256+2) ){
    /* If we decide it is a truecolor display. Then we don't have to 
       go through the expensive business of getting the pixels by lookups.*/
    truecolor=1;
/*      fprintf(stderr,"True color shortcut Default colors.\n"); */
    /*    for(i=0;i<a_maxPixels;i++){
      status=XParseColor(accis_display,accis_colormap,accis_colornames[i],
			 &theHardColor);
      accis_pixels[i]=
	((theHardColor.red/256)*256+(theHardColor.green/256))*256
	+(theHardColor.blue/256);
      fprintf(stderr,"%u,\n",accis_pixels[i]); 
    }
    This now not needed initialized in declaration */
  }else{
    truecolor=0;
    /*fprintf(stderr,"Status false: pixel:%d\n",theRGBColor.pixel);*/
/*      fprintf(stderr,"Looking up Default colors.\n"); */
    for(i=0;i<a_maxPixels;i++){
      status=XLookupColor(accis_display,accis_colormap,accis_colornames[i],
			  &theRGBColor,&theHardColor);
      if(status !=0){
	status=XAllocColor(accis_display,accis_colormap,&theHardColor);
	if(status !=0){
	  accis_pixels[i]=theHardColor.pixel;
	}else{
	accis_pixels[i]=BlackPixel(accis_display,0);
	}
      }else{
	accis_pixels[i]=BlackPixel(accis_display,0);
      }
/*        fprintf(stderr,"%10u,   %.3f %.3f %.3f\n",accis_pixels[i], */
/*  	      theHardColor.red/(256.*256.), */
/*  	      theHardColor.green/(256.*256.), */
/*  	      theHardColor.blue/(256.*256.));  */
    }
  }
}

/* ******************************************************************** */
/* End plotting and return to text editing. */
/* Subroutine */ 
/* #include <curses.h>*/
int txtmode_()
{
  XEvent event; 
  XFlush(accis_display);
  do{
    /*    printf("Executing XtNextEvent"); */
    XtNextEvent(&event);
    /* XNextEvent(accis_display,&event); is equivalent */
    XtDispatchEvent(&event);  
    /*    printf("The event type: %d\n",event); */
  }while(event.type != ButtonPress && event.type != KeyPress );
  /* Here we should give the focus back to parent, but I don't see how.
    XSetInputFocus(accis_display, ??, PointerRoot,
		   CurrentTime); */
}

/* ******************************************************************** */
/* Flush the plot buffer
/* Subroutine */ 
int accisflush_()
{
  XFlush(accis_display);
}

/* ********************************************************************* */
/* Subroutine */ int scolor_(li)
long *li;
{
  /* *ncolor=*li; */
  if((*li < a_maxPixels) && (*li >= 0)){
    XSetForeground(accis_display,accis_gc,accis_pixels[(int) *li]);
    return 1;
  }else{    
    return 0;
  }
} /* scolor_ */

/* ******************************************************************** */
/* Subroutine */ int vec_(px, py, ud)
long *px, *py, *ud;
{ /*  Draw vector on screen, with pen up or down. */
    static int px1=0,py1=0,px2=0,py2=0;
    extern XPoint accis_path[];
    extern int accis_pathlen;

    px1=px2;
    py1=py2;
    px2 = *px;
    py2 = *py;
    if( *ud != 0) {
      XDrawLine(XtDisplay(accis_drawing),XtWindow(accis_drawing), accis_gc,
		  px1,py1,px2,py2);
      XDrawLine(XtDisplay(accis_drawing),accis_pixmap, accis_gc,
		  px1,py1,px2,py2);
      if(accis_pathlen<accis_path_max){      /* Add point to path */
	accis_pathlen++;
      }
    }else{ /* Restart path */
      accis_pathlen=0;
    }
    accis_path[accis_pathlen].x=*px;
    accis_path[accis_pathlen].y=*py;
/*    XFlush(accis_display);
 Flush removed here. Now relies on txtmode to flush display.  */
    return 0;
} /* vec_ */
/* ******************************************************************** */
int vecfill_()
{
    extern XPoint accis_path[];
    extern int accis_pathlen;
    if(accis_pathlen>1){ /* If path is more than 2 points, fill. */
      XFillPolygon(accis_display,accis_window,accis_gc,
		   accis_path,accis_pathlen+1,Nonconvex,CoordModeOrigin);
      XFillPolygon(accis_display,accis_pixmap,accis_gc,
		   accis_path,accis_pathlen+1,Nonconvex,CoordModeOrigin);
    }
}
/* ******************************************************************** */
/* Not currently in use. */
void config_handler(w,cs_s,event)
Widget w;
struct Screen_Size cs_s;
XEvent *event;
{
cs_s.width=event->xconfigure.width;
cs_s.height=event->xconfigure.height;
printf("config_handler width= %d, height= %d\n",cs_s.width,cs_s.height);
}
/* ******************************************************************** */
/* ******************************************************************** */
void accis_refresh(w,data,event)
Widget w;
caddr_t data;
XEvent *event;
{
  XCopyArea(XtDisplay(w),accis_pixmap,accis_window,accis_gc,0,0,
	    s_s.width,s_s.height,0,0);
  XFlush(accis_display);
}
/* ******************************************************************** */
int accisrefresh_()
{
caddr_t data;
XEvent *event;
 accis_refresh(accis_wshell,data,event);
}
/* ******************************************************************** */

ManageEvents()
{
  XEvent event; 
  while(XtPending()){
    XtNextEvent(&event);
    XtDispatchEvent(&event);  
  }
}

float xeye,yeye,zeye;
float xeye0,yeye0,zeye0;
float accis_x0,accis_y0;

void accis_butdown(w,data,event)
Widget w;
caddr_t data;
XEvent *event;
{
  accis_x0=event->xbutton.x;
  accis_y0=event->xbutton.y;
  butdown_(&xeye0,&yeye0,&zeye0);
  xeye=xeye0; yeye=yeye0; zeye=zeye0;
}

void accis_butup(w,data,event)
Widget w;
caddr_t data;
XEvent *event;
{
  butup_(&xeye,&yeye,&zeye);
}


void accis_moved(w,data,event)
Widget w;
caddr_t data;
XEvent *event;
{
  float xmoved,ymoved;
  /* need screen to normal scaling not done here at present*/ 
  xmoved=event->xbutton.x-accis_x0;
  ymoved=event->xbutton.y-accis_y0;
/*    printf("eye: %f %f %f\n",xeye,yeye,zeye); */
  viewrot_(&xmoved,&ymoved,&xeye0,&yeye0,&zeye0,&xeye,&yeye,&zeye);
  cubeupd_(&xeye,&yeye,&zeye);
}
/* ******************************************************************** */
int eye3d_(value)
     int *value;
{
  extern void accis_butdown();
  extern void accis_butup();
  extern void accis_moved();
  XEvent event; 
  XFlush(accis_display);  
  XtAddEventHandler(accis_drawing,ButtonPressMask,FALSE,
		      accis_butdown,NULL);
  XtAddEventHandler(accis_drawing,ButtonReleaseMask,FALSE,
		      accis_butup,NULL);
  XtAddEventHandler(accis_drawing,ButtonMotionMask,FALSE,
		      accis_moved,NULL);
  do{
/*      printf("Executing XtNextEvent "); */
    XtNextEvent(&event);
    XtDispatchEvent(&event);  
/*      printf("The event type: %d\n",event); */
  }while(event.type != KeyPress && event.type != ButtonRelease  );
  if( event.type != ButtonRelease || 
      ( xeye==xeye0 && yeye==yeye0 && zeye==zeye0) ) *value=0; else *value=1;
}
/* ******************************************************************** */
/* Routines for using 240 color gradients in preference to 16 fixed colors.*/
/************** Setup The Gradient **********************/
int accisgradinit_(r1,g1,b1,r2,g2,b2)
     long *r1,*g1,*b1,*r2,*g2,*b2;
     /* RGB are specified in the range 0 to 65535 */
{
  int i,j,status;
  XColor theRGBcolor;
  /*See if this is a sensible 24 bit display or not */
  theRGBcolor.red=255*256;
  theRGBcolor.green=127*256;
  theRGBcolor.blue=2*256;
  XAllocColor(accis_display,accis_colormap,&theRGBcolor);
  if(theRGBcolor.pixel==((255*256+127)*256+2)){
    status=1;
    /*fprintf(stderr,"True color shortcut.\n");*/
  }else{
    status=0;
    /*fprintf(stderr,"Status false: pixel:%d\n",theRGBcolor.pixel);*/
  }
  for (i=0;i<a_gradPixno;i++){
    j=(i* *r2+(a_gradPixno-1-i)* *r1)/(a_gradPixno-1.) ;
    if(j<0)j=0; else if(j>65535)j=65535;
    a_gradred[i]=theRGBcolor.red=j;
    j=(i* *g2+(a_gradPixno-1-i)* *g1)/(a_gradPixno-1.) ;
    if(j<0)j=0; else if(j>65535)j=65535;
    a_gradgreen[i]=theRGBcolor.green=j;
    j=(i* *b2+(a_gradPixno-1-i)* *b1)/(a_gradPixno-1.) ;
    if(j<0)j=0; else if(j>65535)j=65535;
    a_gradblue[i]=theRGBcolor.blue=j;
    /* printf("theRGBcolor %d,%d,%d\n",theRGBcolor.red,theRGBcolor.green,theRGBcolor.blue);  */
    if(status){
      a_gradPix[i]=
	((theRGBcolor.red/256)*256+(theRGBcolor.green/256))*256
	+(theRGBcolor.blue/256);
      /* fprintf(stderr,"Allocated Color %d =%d\n",i,a_gradPix[i]);*/
    }else if(XAllocColor(accis_display,accis_colormap,&theRGBcolor)){
      /*      fprintf(stderr,"Allocated Color %d=%d\n",j,theRGBcolor.pixel);*/
      a_gradPix[i]=theRGBcolor.pixel;
    }else{
      a_gradPix[i]=BlackPixel(accis_display,0);
    }
  }
  a_grad_inited=2;
}
/************** Setup Default Gradient Grey-scale **********************/
int accisgraddef_()
{
  int i,j,status;
  unsigned long white;
  XColor theRGBcolor;
     /* RGB are specified in the range 0 to 65535 */
  theRGBcolor.red=65535;
  theRGBcolor.green=65535;
  theRGBcolor.blue=65535;
  /*Lookup numeric white*/
  if(XAllocColor(accis_display,accis_colormap,&theRGBcolor)){
    /*    fprintf(stderr,"Allocated Color %d\n",theRGBcolor.pixel);*/
    white=theRGBcolor.pixel;
  }else{
    white=0;
  }
  /* If we have a true or direct color display, just scale. */
  if(white>256){
    for (i=0;i<a_gradPixno;i++){
      j=(i* 65535)/(a_gradPixno-1.) ;
/*        a_gradPix[i]=white*i/(a_gradPixno-1.); */
      a_gradred[i]=theRGBcolor.red=j;
      a_gradgreen[i]=theRGBcolor.green=j;
      a_gradblue[i]=theRGBcolor.blue=j;
      a_gradPix[i]=
	((theRGBcolor.red/256)*256+(theRGBcolor.green/256))*256
	+(theRGBcolor.blue/256);
/*        fprintf(stderr,"a_gradPix[%d]=%d\n",i,a_gradPix[i]); */
    }
  }else{ /* Have to do the network-expensive lookups */
    for (i=0;i<a_gradPixno;i++){
      j=(i* 65535)/(a_gradPixno-1.) ;
      a_gradred[i]=theRGBcolor.red=j;
      a_gradgreen[i]=theRGBcolor.green=j;
      a_gradblue[i]=theRGBcolor.blue=j;
      if(XAllocColor(accis_display,accis_colormap,&theRGBcolor)){
	/*     fprintf(stderr,"Allocated Color %d=%d\n",j,theRGBcolor.pixel);*/
	a_gradPix[i]=theRGBcolor.pixel;
      }else{
	a_gradPix[i]=BlackPixel(accis_display,0);
      }
    }
  }
  a_grad_inited=1;
}
/********** Use a gradient color out of 240 *********************************/
/* Subroutine */ int acgradcolor_(li)
long *li;
{
    /*Default gradient is a gray-scale.*/
    /*If not inited already init by default.*/
  if(!a_grad_inited)accisgraddef_();
  /* *ncolor=*li; */
  if((*li < a_gradPixno) && (*li >= 0)){
    XSetForeground(accis_display,accis_gc,a_gradPix[(int) *li]);
    return 1;
  }else{    
    return 0;
  }
}
/********** Tell the current rgb color ********************************
/* Subroutine */ int getrgbcolor_(ipixel,red,green,blue)
long *ipixel,*red,*green,*blue;
{
    /*If not inited already init by default.*/
  if(!a_grad_inited)accisgraddef_();
  *red=a_gradred[(int) *ipixel];
  *green=a_gradgreen[(int) *ipixel];
  *blue=a_gradblue[(int) *ipixel];
  return 0;
}
/***********************************************************************
This is killingly expensive over the network, and gets called every time
the color is set when writing a print file. 
XQueryColor is painfully slow.
 So not using this.
int getrgbcolor_(ipixel,red,green,blue)
long *ipixel,*red,*blue,*green;
{
  XColor theRGBcolor;
  if(!a_grad_inited)accisgraddef_();
  theRGBcolor.pixel=a_gradPix[(int) *ipixel];
  XQueryColor(accis_display,accis_colormap,&theRGBcolor);
  *red=theRGBcolor.red;
  *green=theRGBcolor.green;
  *blue=theRGBcolor.blue;
  return 0;
}
************************************************************************/

