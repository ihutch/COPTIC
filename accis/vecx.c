/* Xwindow driver for accis plotting */
/* Fortran callable routines */
/* This version compatible with linux g77 and ifc and Sun Apr 2003*/
int accis_driver_(){return 1;}
/* ********************************************************************** */

#include <stdio.h>
#include  <X11/StringDefs.h> 
#include  <X11/Intrinsic.h> 
#include  <X11/Core.h> 

/* Argument-passing typdef */

typedef int FORT_INT;

/* Globals for these routines*/
Widget accis_wshell;
Widget accis_drawing;
Display *accis_display;
Drawable accis_window;
Pixmap accis_pixmap;
XWindowAttributes accis_attributes;
/* Drawing only to the back buffer (pixmap) default */
int accis_back=1;

GC accis_gc;
XGCValues accis_gcvalues;
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
/* Until svga is called, the default is that there is no display */
int accis_nodisplay=1;
static int accis_eye3d=9999;

/* 256 color gradient globals */
/* Use only 240 colors so that 16 are left for the 16 colors and you can
   make gifs out of the resultant without a problem.*/
#define a_gradPixno 240
/*unsigned long a_gradPix[a_gradPixno];*/
int a_gradPix[a_gradPixno];
int a_testedtrue=0;
int a_truecolor=0;
int a_grad_inited=0;
int a_gradred[a_gradPixno];
int a_gradgreen[a_gradPixno];
int a_gradblue[a_gradPixno];
int a_gradno=a_gradPixno;/* Publically available Pixno */

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

/* Cope gracefully with bad match errors occasionally thrown by XSetInputFocus
   Install our own non-fatal handler during the call, only once.*/
static XErrorHandler accis_old_handler = (XErrorHandler) 0 ;
static int accis_errorhandler(Display *display, XErrorEvent *theEvent) {
   fprintf(stderr, "Intercepted Xlib error: error code %d request code %d",
		theEvent->error_code, theEvent->request_code) ;
   if(theEvent->error_code == BadMatch){
     fprintf(stderr,"  BadMatch. Doing nothing.\n");
     return 0;
   }else{
     fprintf(stderr,"  Unfiltered error passed to Xlib.\n");
     accis_old_handler(display,theEvent);
     return 1;
   }
}
#define ACCIS_SET_FOCUS	  accis_set_focus()
void accis_set_focus(){
  if(accis_old_handler==0)
    accis_old_handler=XSetErrorHandler(accis_errorhandler) ;
  /*if(XGetWindowAttributes(accis_display, accis_window,&accis_attributes) \
    == Success && accis_attributes.map_state == IsViewable)/*DeliberateErrors*/
    XSetInputFocus(accis_display, accis_window, RevertToParent,CurrentTime);
    /* Immediately restoring the old one did not work. */	
    /*XSync(accis_display,True);XSetErrorHandler(accis_old_handler);*/
}

void badmatchtest_(){ /* Attempts to generate routine bad match events*/
  XSetWindowAttributes theAttributes;

  XGetWindowAttributes(accis_display, accis_window,&accis_attributes);
  printf("override_redirect=%d",accis_attributes.override_redirect);
  theAttributes.override_redirect=1; /*Make crazy to generate BadValue */
  XChangeWindowAttributes(accis_display,accis_window,CWOverrideRedirect,&theAttributes);
  XGetWindowAttributes(accis_display, accis_window,&accis_attributes);
  printf(" changed to override_redirect=%d\n",accis_attributes.override_redirect);
  /*XLowerWindow(accis_display, accis_window); Does not work */
  /*XIconifyWindow(accis_display, accis_window,0);*/
  XWithdrawWindow(accis_display, accis_window,0); /*Works crudely.*/
  /*  ACCIS_SET_FOCUS; subsequently generates bad match*/
}


#define ACCIS_NARGVS 30
char *accis_argv[ACCIS_NARGVS];
int accis_argc=0; 
/* Get the command line arguments from fortran main using getargs etc. */
void getcmdargs_()
{
  FORT_INT iargc=0; 
  static char argv[512];
  FORT_INT charlen=512;
  extern ifcargs_();
  int i,j;

  cmdlineargs_(&iargc,argv,charlen);
  accis_argc=iargc;
  /* Parse the returned command line string into standard C argv*/
  accis_argv[0]=argv+strspn(argv," ");
  for(i=1;i<=(iargc<ACCIS_NARGVS ? iargc : ACCIS_NARGVS-1);i++){
    accis_argv[i]=accis_argv[i-1]+strcspn(accis_argv[i-1]," ");
    j=strspn(accis_argv[i]," ");
    *(accis_argv[i])=0;
    accis_argv[i]=accis_argv[i]+j;
  }
}

int accisinit()
{
  FORT_INT xp, yp, vm, nc;
  svga_(&xp,&yp,&vm,&nc);
}

/* Subroutine */ 
int svga_(scrxpix, scrypix, vmode, ncolor)
FORT_INT *scrxpix, *scrypix, *vmode, *ncolor;
{  
  static int second=0;
  static int n;
  static Arg wargs[10];
  int theDepth;
  Colormap theColormap;
  extern void config_handler();
  extern void accis_refresh();
  int *svga_argc;
  char **svga_argv;
  int oldwidth,oldheight;

  accis_nodisplay=0;
  *ncolor=15;
  /* Call fortran routine to get arguments into accis_ globals*/
  if(second == 0) {
    getcmdargs_();
    svga_argc=&accis_argc;
    svga_argv=accis_argv;
/*      for(n=0;n<*svga_argc;n++){printf("%d:%s:\n",n,svga_argv[n]);} */
  }

  if(second == 0){
    accis_wshell = XtInitialize("accis","Accis", NULL, 0,
      svga_argc, svga_argv);

    accis_drawing = XtCreateManagedWidget("drawing",coreWidgetClass,
		 accis_wshell, NULL, 0);

    *vmode=88;
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
    accis_gc = XCreateGC(accis_display, accis_window, GCLineWidth, &accis_gcvalues);
    /* line_width 1 is needed to avoid rare pixel omission in fills */
    accis_gcvalues.line_width=0;
    XChangeGC(accis_display, accis_gc, GCLineWidth, &accis_gcvalues);
    accis_depth=DefaultDepth(accis_display,0);
    accis_colormap=DefaultColormap(accis_display,0);
    initDefaultColors();
    /* Set up arguments for resizing. */
    n = 0;
    XtSetArg(wargs[n], XtNheight, &s_s.height); n++;
    XtSetArg(wargs[n], XtNwidth, &s_s.width); n++;
    XtGetValues(accis_wshell, wargs, n);
    /* Pixmap setup */
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
  /* Get the possibly new values of wdth/height*/
  oldwidth=s_s.width;
  oldheight=s_s.height;
  XtGetValues(accis_wshell, wargs, n);
  *scrxpix=s_s.width;
  *scrypix=s_s.height;
  if(oldwidth!=s_s.width || oldheight!=s_s.height){
/*     printf("size changed %d,%d to %d,%d\n", */
/* 	   oldwidth,oldheight,s_s.width,s_s.height); */
    XFreePixmap(accis_display,accis_pixmap);
    accis_pixmap=XCreatePixmap(accis_display,accis_window,
			       s_s.width,s_s.height,accis_depth);
  }
  /*This spoils the display of rotating windows completely.
    It is unnecessary since we fill rectangle display and pixel.
    XClearWindow(accis_display,accis_window); */
  XSetForeground(accis_display,accis_gc,accis_pixels[0]);  
  if(accis_back==0)XFillRectangle(accis_display,accis_window,accis_gc,0,0,
		 s_s.width,s_s.height);
  /* Clear the pixmap  */
  XFillRectangle(accis_display,accis_pixmap,accis_gc,0,0,
		 s_s.width,s_s.height);
  XSetForeground(accis_display,accis_gc,accis_pixels[15]); 
  XRaiseWindow(accis_display, accis_window);
  return 0;
}

/************************************************************************/
int is_truecolor()
{  
  /*See if this is a sensible 24 bit display or not 
    If we decide it is a truecolor display. Then we don't have to 
    go through the expensive business of getting the pixels by lookups.*/
  XColor theRGBColor;
  if(accis_display==NULL)accisinit();
  if(!a_testedtrue){ /* We haven't already discovered and stored*/
    theRGBColor.red=255*256;
    theRGBColor.green=127*256;
    theRGBColor.blue=2*256;
    XAllocColor(accis_display,accis_colormap,&theRGBColor);
    if(theRGBColor.pixel==((255*256+127)*256+2) ){
      a_truecolor=1;
    }else{
      a_truecolor=0;
    }
    a_testedtrue=1;
  }
  return a_truecolor;
}
/* ******************************************************************** */
int initDefaultColors()
{
  XColor theRGBColor;
  XColor theHardColor;
  int status,truecolor;
  int i;
  if(accis_nodisplay){
    truecolor=1;
  }else{
    if(is_truecolor()){	
      truecolor=1;
/*      fprintf(stderr,"True color shortcut Default colors.\n"); */
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
}

/* ******************************************************************** */
/* End plotting and return to text editing. */
/* Subroutine */ 
/* #include <curses.h>*/
int txtmode_()
{
  XEvent event; 
  accisrefresh_();
  /*usleep(100000);/* Still struggling with bad match events */
  /* Try to get the focus into this window for keyboard control.*/
  ACCIS_SET_FOCUS;
  do{
    /*    printf("Executing XtNextEvent"); */
    XtNextEvent(&event);
    /* XNextEvent(accis_display,&event); is equivalent */
    XtDispatchEvent(&event);  
    /*    printf("The event type: %d\n",event); */
  }while(event.type != ButtonPress && event.type != KeyPress );
}

/* ******************************************************************** */
/* Flush the plot buffer
/* Subroutine */ 
int accisflush_()
{
  XFlush(accis_display);  
  accisrefresh_();
  ManageEvents();
}

/* ********************************************************************* */
/* Subroutine */ int scolor_(li)
FORT_INT *li;
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
FORT_INT *px, *py, *ud;
{ /*  Draw vector on screen, with pen up or down. Or point.*/
    static int px1=0,py1=0,px2=0,py2=0;
    extern XPoint accis_path[];
    extern int accis_pathlen;

    px1=px2;
    py1=py2;
    px2 = *px;
    py2 = *py;
    if( *ud > 0) {
      if(accis_back==0)XDrawLine(XtDisplay(accis_drawing),
				 XtWindow(accis_drawing),
				 accis_gc,px1,py1,px2,py2);
      XDrawLine(XtDisplay(accis_drawing),accis_pixmap, accis_gc,
		  px1,py1,px2,py2);
      if(accis_pathlen<accis_path_max){      /* Add point to path */
	accis_pathlen++;
      }
    }else{ 
      if(*ud == -1){ /* Write Point only */
	if(accis_back==0)XDrawPoint(XtDisplay(accis_drawing),
				    XtWindow(accis_drawing),
				    accis_gc,px2,py2);
	XDrawPoint(XtDisplay(accis_drawing),accis_pixmap, 
		  accis_gc,px2,py2);
      }
      /* Restart path for disconnected point */
      accis_pathlen=0;
    }
    accis_path[accis_pathlen].x=*px;
    accis_path[accis_pathlen].y=*py;
    return 0;
} /* vec_ */
/* ******************************************************************** */
int vecfill_()
{
    extern XPoint accis_path[];
    extern int accis_pathlen;
    if(accis_pathlen>1){ /* If path is more than 2 points, fill. */
      if(accis_back==0)XFillPolygon(accis_display,accis_window,accis_gc,
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
  /* Does not fix bad match errors XSync(accis_display,False);*/
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
/* ******************************************************************** */
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
  /* If we started moving without setting eye nonzero, set it first*/
    if(xeye0==0. && yeye0==0. && zeye0==0.) accis_butdown(w,data,event); 
  /* need screen to normal scaling not done here at present*/ 
  xmoved=event->xbutton.x-accis_x0;
  ymoved=event->xbutton.y-accis_y0;
/*    printf("eye: %f %f %f\n",xeye,yeye,zeye); */
  viewrot_(&xmoved,&ymoved,&xeye0,&yeye0,&zeye0,&xeye,&yeye,&zeye);
  cubeupd_(&xeye,&yeye,&zeye);
  accisrefresh_();
}

void accis_keypress(w,data,event)
Widget w;
caddr_t data;
XEvent *event;
{
  /* Testing only 
  printf("Key event keycode: %u keysym:%x",event->xkey.keycode,
	 (int)XLookupKeysym(&(event->xkey),0)); */
  /* This breaks the interface. Don't do it. ManageEvents();*/
  ACCIS_SET_FOCUS;
}

/* ******************************************************************** */
/* Externally callable routine to set noeye3d return value.
   Set to 9999 to disable, i.e. revert to waiting for input. */ 
int noeye3d_(value)
     int *value;
{
  if(*value>1000)accis_eye3d=9999; else accis_eye3d=*value;
}

/* ******************************************************************** */
int eye3d_(value)
     int *value;
{
  extern void accis_butdown();
  extern void accis_butup();
  extern void accis_moved();
  XEvent event; 

  if(accis_nodisplay){ *value=0; return 0; }
  if(accis_eye3d == 1){ *value=0; return 0; }
  accisrefresh_();
  /* Disabled waiting code: */
  if(accis_eye3d != 9999){
    if(XPending(accis_display)){
      XPeekEvent(accis_display,&event);
      accis_eye3d=9999;
      if(event.type==KeyPress){
	*value=(int)XLookupKeysym(&(event.xkey),0);
      }
    }else{
      *value=accis_eye3d; return 0; 
    }
  }
  /*Without this, the keyboard does not focus to this window.*/
  ACCIS_SET_FOCUS;
  XtAddEventHandler(accis_drawing,ButtonPressMask,FALSE,
		      accis_butdown,NULL);
  XtAddEventHandler(accis_drawing,ButtonReleaseMask,FALSE,
		      accis_butup,NULL);
  XtAddEventHandler(accis_drawing,ButtonMotionMask,FALSE,
		      accis_moved,NULL);
  /* A handler seems to be needed if we are to recognize keys. Even if
     it does not do anything much.*/
  XtAddEventHandler(accis_drawing,KeyPressMask,FALSE,
		      accis_keypress,NULL);
  XtAddEventHandler(accis_drawing,EnterWindowMask,FALSE,
		      accis_keypress,NULL);
  /* Wait for a key/button press */
  do{
    XtNextEvent(&event);
    XtDispatchEvent(&event);  
/*   printf("First loop:The event type: %d, %d\n",event.type,ButtonPress);  */
  }while(event.type != ButtonPress && event.type != KeyPress);
  /* Recognize KeyPress as sign to exit.*/
  if(event.type == KeyPress) {
    *value=(int)XLookupKeysym(&(event.xkey),0);
    return *value;
  /* Get rid of all the queued contiguous KeyPress/Release events so we don't 
   over-run the rotation when the key is lifted. */
    if(XPending(accis_display)) XPeekEvent(accis_display,&event);
    while(XPending(accis_display) &&
	  (event.type==KeyPress || event.type==KeyRelease) ){
      XNextEvent(accis_display,&event);
      if(XPending(accis_display)) XPeekEvent(accis_display,&event);	\
    }
  }
  do{
/*      printf("Executing XtNextEvent "); */
    XtNextEvent(&event);
    XtDispatchEvent(&event);  
/*        printf("The event type: %d\n",event);  */
  }while(event.type != KeyPress && event.type != ButtonRelease  );
  if( event.type != ButtonRelease || 
      ( xeye==xeye0 && yeye==yeye0 && zeye==zeye0) ) *value=0; else *value=1;
  return *value;
}
/* ******************************************************************** */
/* Routines for using 240 color gradients in preference to 16 fixed colors.*/
/************** Setup The Gradient *********************
  Linear gradients from value 1 to value 2. */
int accisgradinit_(r1,g1,b1,r2,g2,b2)
     FORT_INT *r1,*g1,*b1,*r2,*g2,*b2;
     /* On 64 bit machine this seems to need int *r1,*g1,*b1,*r2,*g2,*b2;*/
     /* RGB are specified in the range 0 to 65535 */
{
  int i,j,status;
  XColor theRGBcolor;
  if(accis_nodisplay){
    status=1;
  }else{
    if(status=is_truecolor() ==0){
      fprintf(stderr,"True Color 24 bit Status false: pixel:%lu\n",theRGBcolor.pixel);
    }
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
/*     printf("theRGBcolor %d,%d,%d\n",theRGBcolor.red,theRGBcolor.green,theRGBcolor.blue); */
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
/************** Setup The Gradient *********************
  Normalized Polynomial gradients from arrays r,g,b of length rl,gl,bl. */
/* Currenly unused. Use accisgradset instead */
int accisgradpoly_(rc,gc,bc,rl,gl,bl)
     FORT_INT *rc,*gc,*bc,*rl,*gl,*bl;
     /* RGB are specified in the range 0 to 65535 */
{
  int i,j,k,status;
  double color;
  XColor theRGBcolor;
  if(accis_nodisplay){
    status=1;
  }else{
    if(status=is_truecolor() ==0){
      fprintf(stderr,"True Color 24 bit Status false: pixel:%lu\n",theRGBcolor.pixel);
    }
  }
  for (i=0;i<a_gradPixno;i++){
    color=0.;
    for (k=*rl-1;k>=0;k--){
      color=(color*i)/(a_gradPixno-1.)+rc[k];
    }
    j=65535*color;
    if(j<0)j=0; else if(j>65535)j=65535;
    a_gradred[i]=theRGBcolor.red=j;
    color=0.;
    for (k=*gl-1;k>=0;k--){
      color=(color*i)/(a_gradPixno-1.)+gc[k];
    }
    j=65535*color;
    if(j<0)j=0; else if(j>65535)j=65535;
    a_gradgreen[i]=theRGBcolor.green=j;
    color=0.;
    for (k=*bl-1;k>=0;k--){
      color=(color*i)/(a_gradPixno-1.)+bc[k];
    }
    j=65535*color;
    if(j<0)j=0; else if(j>65535)j=65535;
    a_gradblue[i]=theRGBcolor.blue=j;
/*     printf("theRGBcolor %d,%d,%d\n",theRGBcolor.red,theRGBcolor.green,theRGBcolor.blue); */
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
  FORT_INT top=65535;
  FORT_INT bot=0;
  accisgradinit_(&bot,&bot,&bot,&top,&top,&top);
}
/**********************************************************************/
/********** Use a gradient color out of 240 *********************************/
/* Subroutine */ int acgradcolor_(li)
FORT_INT *li;
{
  int ii;
    /*Default gradient is a gray-scale.*/
    /*If not inited already init by default.*/
  if(!a_grad_inited)accisgraddef_();
  /* *ncolor=*li; */
  ii=( (*li < a_gradPixno) ? ( *li >= 0 ? *li : 0) : a_gradPixno-1 );
/*   if((*li < a_gradPixno) && (*li >= 0)){ */
  XSetForeground(accis_display,accis_gc,a_gradPix[(int) ii]);
  return 1;
/*   }else{     */
/*     return 0; */
/*   } */
}
/********** Tell the current rgb color ********************************
/* Subroutine */ int getrgbcolor_(ipixel,red,green,blue)
FORT_INT *ipixel,*red,*green,*blue;
{
  int ii;
    /*If not inited already init by default.*/
  if(!a_grad_inited)accisgraddef_();
  /* Note that for fortran calls ipixel goes from 1 to 240*/
  ii=( (*ipixel<a_gradPixno) ? ( *ipixel >= 0 ? *ipixel : 0) :a_gradPixno-1 );
  *red=a_gradred[ii];
  *green=a_gradgreen[ii];
  *blue=a_gradblue[ii];
  return 0;
}
/***********************************************************************/
/* Subroutine for dummy svga call when using no display */ 
int svganodisplay_(scrxpix, scrypix, vmode, ncolor)
FORT_INT *scrxpix, *scrypix, *vmode, *ncolor;
{
  /*Now nodisplay=1 until svga is called. Still, this setting allows us to
    call this routine after already having plotted something. */
  accis_nodisplay=1;
  *scrxpix=640;
  *scrypix=480;
  *ncolor=15;
  *vmode=0;
  return 0;
}
/***********************************************************************
This is killingly expensive over the network, and gets called every time
the color is set when writing a print file. 
XQueryColor is painfully slow.
 So not using this.
int getrgbcolor_(ipixel,red,green,blue)
FORT_INT *ipixel,*red,*blue,*green;
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
int accisgradnumber_() /* Fortran integer routine for gradcolor number 240*/
{  return a_gradPixno; }
/******************************************************************/
static int ilimit(b,i,t)
     int b,i,t;
{
  int j;
  j=(i<b ? 0 : i);
  j=(t<j ? t : j);
  return j;
}
/**********************************************************************/
/* Set the gradient colors from three fortran arrays.
   Limiting the range to 0-65535, warning if the pixel number is not right. */
int accisgradset_(red,green,blue,npixel)
     FORT_INT *red,*green,*blue,*npixel;
{
  int i,j;
  XColor theRGBcolor;
  if(*npixel!=a_gradPixno) 
    fprintf(stderr,"accisgradset ERROR: Incorrect array length:%d\n",*npixel);
  for (i=0;i<a_gradPixno;i++){
    a_gradred[i]=theRGBcolor.red=ilimit(0,*(red+i),65535);
    a_gradgreen[i]=theRGBcolor.green=ilimit(0,*(green+i),65535);
    a_gradblue[i]=theRGBcolor.blue=ilimit(0,*(blue+i),65535);
    /* printf("theRGBcolor %d,%d,%d\n",theRGBcolor.red,theRGBcolor.green,theRGBcolor.blue);  */
    if(is_truecolor()){
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
/*************************************************************************/
/* Dummy here, returning 0 */
int igradtri_() { return 0;}
/*************************************************************************/
/* G77 fortran callable usleep */
#include <unistd.h>
void usleep_(usecs)
     long *usecs;
{
  usleep(  (unsigned long) *usecs);
}
/************************************************************************/
/* Switch to drawing to the back buffer only */
void glback_()
{
  accis_back=1;
}
/************************************************************************/
/* Bring back buffer to front and return to writing there.*/
void glfront_()
{
  accis_back=0;
  accisrefresh_();
}

