/* Xlib driver for accis plotting */
/* Fortran callable routines */
/* This version based on Xlib only, no Xt, derived from glX version. */
int accis_driver_(){return 4;}
/* ********************************************************************** */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <X11/X.h>
#include <X11/Xlib.h>
#include    <X11/Intrinsic.h> 
#include    <X11/Core.h> 
/* Argument-passing typdef */
typedef int FORT_INT;

/* Globals for these routines*/
static Display *accis_display=NULL;
static Window accis_root;
static Window accis_window;
static Pixmap accis_pixmap=0;
static XVisualInfo             *accis_vi;
static Colormap                accis_cmap;
static XSetWindowAttributes    accis_swa;
static GC                      accis_gc;
static XGCValues               accis_gcvalues;
static XWindowAttributes       accis_gwa;
static XEvent                  accis_xev;
static Visual                  *accis_visual;
static Colormap accis_colormap;
static unsigned long accis_valuemask=0x0000;
static int accis_depth;
static int accis_screen;
static int accis_listing=0;
static int accis_eye3d=9999;
/* Drawing only to the back buffer (pixmap) default */
static int accis_back=1;


/* Static maximum number of points in path */
#define accis_path_max 4000
static XPoint accis_path[accis_path_max];
static int accis_pathlen=0;
/* Until svga is called, the default is that there is no display 
   Indicated by 99, not quite the same as nodisplay=1*/
static int accis_nodisplay=99;

/* 256 color gradient globals */
/* Use only 240 colors so that 16 are left for the 16 colors and you can
   make gifs out of the resultant without a problem.*/
#define a_gradPixno 240
static int a_gradPix[a_gradPixno];
static int a_testedtrue=0;
static int a_truecolor=0;
static int a_grad_inited=0;
static int a_gradred[a_gradPixno];
static int a_gradgreen[a_gradPixno];
static int a_gradblue[a_gradPixno];
static int a_gradno=a_gradPixno;/* Publically available Pixno */

/* 16-color model globals. Pixels for truecolor display. */
#define a_maxPixels 16
static unsigned long accis_pixels[a_maxPixels]={
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

static unsigned char accis_rgb[3*a_maxPixels];

static char *accis_colornames[a_maxPixels]=
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

/* ********************************************************************** */
#define ACCIS_NARGVS 30
static char *accis_argv[ACCIS_NARGVS];
static int accis_argc=0; 
static char *accis_geometry=NULL;
/* Get the command line arguments from fortran main using getargs etc. */
void getcmdargs_()
{
  FORT_INT iargc=0; 
  static char argv[512];
  FORT_INT charlen=512;
  extern void cmdlineargs_();
  int i,j;

  cmdlineargs_(&iargc,argv,charlen); 
 /* &charlen was wrong. It seems that length is passed by value */
  accis_argc=iargc;
  /* Parse the returned command line string into standard C argv*/
  accis_argv[0]=argv+strspn(argv," ");
  for(i=1;i<=(iargc<ACCIS_NARGVS ? iargc : ACCIS_NARGVS-1);i++){
    accis_argv[i]=accis_argv[i-1]+strcspn(accis_argv[i-1]," ");
    j=strspn(accis_argv[i]," ");
    *(accis_argv[i])=0;
    accis_argv[i]=accis_argv[i]+j;
    /* Set pointers to those we care about */
    if(strstr(accis_argv[i-1],"-geom")==accis_argv[i-1])
      {accis_geometry=accis_argv[i];}
  }
}

/* ********************************************************************** */
/* Cope gracefully with bad match errors occasionally thrown by XSetInputFocus
   Install our own non-fatal handler during the call, only once.*/
static int accis_errorcount = 0;
static XErrorHandler accis_old_handler = (XErrorHandler) 0 ;
int accis_errorhandler(Display *display, XErrorEvent *theEvent) {
   if(theEvent->error_code == BadMatch){
     if(++accis_errorcount>2){
       fprintf(stderr,"Too many badmatch errors, %d\n",accis_errorcount);
       accis_errorcount=0;
       return 2;
     }
     if(theEvent->request_code == 42){
       /* Fix focus errors quietly. Wait and try to set again.
	  This reentrant call could easily give an infinite loop. 
	  Hence the above errorcount trap. */
       usleep(10000);
       XSetInputFocus(accis_display, accis_window, RevertToParent,CurrentTime);
     }else{/*
       fprintf(stderr, "Intercepted Xlib error: error code %d request code %d",
		theEvent->error_code, theEvent->request_code) ;
		fprintf(stderr,"  BadMatch. Doing nothing.\n");*/
     }
     return 0;
   }else{
     fprintf(stderr, "Intercepted Xlib error: error code %d request code %d",
		theEvent->error_code, theEvent->request_code) ;
     fprintf(stderr,"  Unfiltered error passed to Xlib.\n");
     accis_old_handler(display,theEvent);
     return 1;
   }
}
/* ********************************************************************** */
#define ACCIS_SET_FOCUS	  accis_set_focus()
void accis_set_focus(){
  if(accis_old_handler==0)
    accis_old_handler=XSetErrorHandler(accis_errorhandler) ;
  XSetInputFocus(accis_display, accis_window, RevertToParent,CurrentTime);
}
/* ********************************************************************** */
/* Main setup subroutine */ 
int svga_(scrxpix, scrypix, vmode, ncolor)
FORT_INT *scrxpix, *scrypix, *vmode, *ncolor;
{  
  static int second=0;
  extern void initDefaultColors();
  XSizeHints hints;
  int x_size,y_size,x_off,y_off,gravity;
  int oldwidth=0,oldheight=0;

  accis_nodisplay=0;
  *ncolor=15;
  *vmode=88;
  if(second == 0) {
  /* Call fortran routine to get arguments into accis_ globals*/
    getcmdargs_();

/* Start of Xlib setup calls ******************/
    if( (accis_display = XOpenDisplay(NULL)) == NULL) {
        printf("\n\tcannot connect to X server\n\n");
        exit(0); 
    }
    accis_root=DefaultRootWindow(accis_display);
    accis_screen=DefaultScreen(accis_display);
    accis_visual=DefaultVisual(accis_display,accis_screen);
    accis_depth=DefaultDepth(accis_display,accis_screen);
    accis_cmap = XCreateColormap(accis_display, accis_root, accis_visual, 
				 AllocNone);
    accis_swa.colormap = accis_cmap;
    accis_swa.background_pixel=0;
    accis_swa.event_mask = ExposureMask | KeyPressMask 
      | ButtonPressMask | ButtonReleaseMask | ButtonMotionMask ;

/* Simple way to use the Xresources */
    if(accis_geometry==NULL)
      accis_geometry=XGetDefault(accis_display,"Accis","Geometry");
    /* Combine specified and program default size info; parse and return */
    hints.flags= USPosition | PPosition | USSize | PSize;
    XWMGeometry(accis_display,0,accis_geometry,"800x600+100+40",0,
		&hints,&x_off,&y_off,&x_size,&y_size,&gravity);

    accis_window=XCreateWindow(accis_display, accis_root, 
 			       x_off, y_off, x_size, y_size, 0,
			       accis_depth, 
			       InputOutput, 
			       accis_visual, 
			       CWBackPixel | CWColormap | CWEventMask, 
			       &accis_swa);
    XMapWindow(accis_display, accis_window);
    XStoreName(accis_display, accis_window, "Accis");
    /* printf("Finished initializing, and mapped window\n"); */

    accis_gc=XCreateGC(accis_display,accis_window,
		       accis_valuemask,&accis_gcvalues);
    accis_colormap=accis_cmap;
    initDefaultColors();
    second++;
    /* printf("Created GC\n"); */
  }
  /* Get the possibly new values of wdth/height*/
  XGetWindowAttributes(accis_display,accis_window, &accis_gwa);
#define pixscale 1
  *scrxpix=accis_gwa.width*pixscale;
  *scrypix=accis_gwa.height*pixscale;
  if(oldwidth!=accis_gwa.width || oldheight!=accis_gwa.height){
/*     printf("size changed %d,%d to %d,%d\n", */
/* 	   oldwidth,oldheight,s_s.width,s_s.height); */
    if(accis_pixmap!=0)XFreePixmap(accis_display,accis_pixmap);
    accis_pixmap=XCreatePixmap(accis_display,accis_window,
			       accis_gwa.width,accis_gwa.height,accis_depth);
    oldwidth=accis_gwa.width;oldheight=accis_gwa.height;
  }

  XSetForeground(accis_display,accis_gc,accis_pixels[0]);  
  if(accis_back==0)XFillRectangle(accis_display,accis_window,accis_gc,0,0,
		 accis_gwa.width,accis_gwa.height);
  /* Clear the pixmap  */
  XFillRectangle(accis_display,accis_pixmap,accis_gc,0,0,
		 accis_gwa.width,accis_gwa.height);
  /* This causes very bad redrawing behaviour on rotation etc.  */
  /* XCopyArea(accis_display,accis_pixmap,accis_window, */
  /* 			 accis_gc,0,0,accis_gwa.width,accis_gwa.height,0,0); */
  
  XSetForeground(accis_display,accis_gc,accis_pixels[15]); 

  /* This does not seem to be necessary and may be the cause of problems.
     XRaiseWindow(accis_display, accis_window);*/
  return 0;
}
/************************************************************************/
int accisinit()
{
  FORT_INT xp, yp, vm, nc;
  svga_(&xp,&yp,&vm,&nc);
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
void initDefaultColors()
{
  XColor theRGBColor;
  XColor theHardColor;
  int status,truecolor;
  int i,j,pix;
  if(accis_nodisplay){
    truecolor=1;
  }else{
    if(is_truecolor()){	
      truecolor=1;
/*      fprintf(stderr,"True color shortcut Default colors.\n"); */
    }else{
      truecolor=0;
      fprintf(stderr,"Looking up default colors the hard way.\n");
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
      }
    }
  }
  /* set up the rgb values */
  for(i=0;i<a_maxPixels;i++){
    pix=accis_pixels[i];
    j=i*3;
    accis_rgb[j+2]=(pix-(pix/256)*256);
    pix=pix/256;
    accis_rgb[j+1]=(pix-(pix/256)*256);
    pix=pix/256;
    accis_rgb[j]=(pix-(pix/256)*256);
  }
}

/* ******************************************************************** */
#define EXPOSE_ACTION      \
  /* Get all the queued contiguous expose events, before redrawing */ \
  if(XPending(accis_display)) XPeekEvent(accis_display,&event);	      \
      while(XPending(accis_display) && event.type==Expose){	      \
	XNextEvent(accis_display,&event);			      \
	if(XPending(accis_display)) XPeekEvent(accis_display,&event); \
      } /* Redraw everything.*/					      \
  XCopyArea(accis_display,accis_pixmap,accis_window,accis_gc,0,0,     \
	    accis_gwa.width,accis_gwa.height,0,0);                    \
  XFlush(accis_display);

/* ******************************************************************** */
/* End plotting Subroutine */ 
void txtmode_()
{
  XEvent event;
  ACCIS_SET_FOCUS;
  EXPOSE_ACTION;
  do{
    /*    printf("Executing XtNextEvent"); */
    XNextEvent(accis_display,&event);
/*     printf("Event: type=%d\n",event.type); */
    if(event.type == Expose){EXPOSE_ACTION;}
  }while(event.type != ButtonPress && event.type != KeyPress );
}
/* ******************************************************************** */
void accisclear_()  /* Simply clear to background fortran callable*/ 
{
  XSetForeground(accis_display,accis_gc,accis_pixels[0]);  
  if(accis_back==0)XFillRectangle(accis_display,accis_window,accis_gc,0,0,
		 accis_gwa.width,accis_gwa.height);
  /* Clear the pixmap  */
  XFillRectangle(accis_display,accis_pixmap,accis_gc,0,0,
		 accis_gwa.width,accis_gwa.height);
  XSetForeground(accis_display,accis_gc,accis_pixels[15]); 
}

/* ******************************************************************** */
/* Flush the plot buffer Subroutine */ 
void accisflush_()
{
  if(accis_nodisplay) return;
  if(accis_back==1)XCopyArea(accis_display,accis_pixmap,accis_window,
			 accis_gc,0,0,accis_gwa.width,accis_gwa.height,0,0);
  XFlush(accis_display);
}
/* ********************************************************************* */
/* Set the Color Subroutine */ int scolor_(li)
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
/* Draw a vector Subroutine */ int vec_(px, py, ud)
FORT_INT *px, *py, *ud;
{ /*  Draw vector on screen, with pen up or down. */
/*     static int px1=0,py1=0,px2=0,py2=0; */
  static float px1=0.,py1=0.,px2=0.,py2=0.;
    extern XPoint accis_path[];
    extern int accis_pathlen;

    px1=px2;
    py1=py2;
    px2 = *px;
    py2 = *py;
    if( *ud == 1) {
      if(accis_back==0)XDrawLine(accis_display,
				 accis_window,
				 accis_gc,px1,py1,px2,py2);
      XDrawLine(accis_display,accis_pixmap, accis_gc,
		  px1,py1,px2,py2);
      if(accis_pathlen<accis_path_max-1){      /* Add point to path */
	accis_pathlen++;
      }
    }else{ 
      if(*ud == -1){ /* Write Point only */
	if(accis_back==0)XDrawPoint(accis_display,
				    accis_window,
				    accis_gc,px2,py2);
	XDrawPoint(accis_display,accis_pixmap, 
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
/* Fill the current path */
void vecfill_()
{
  int i;
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
float xeye,yeye,zeye;
float xeye0,yeye0,zeye0;
float accis_x0,accis_y0;

/* ********************************* */
void accis_butdown(event)
XEvent *event;
{
  extern int butdown_();
/*   printf("entering accis_butdown\n"); */
  accis_x0=event->xbutton.x;
  accis_y0=event->xbutton.y;
  butdown_(&xeye0,&yeye0,&zeye0);
  xeye=xeye0; yeye=yeye0; zeye=zeye0;
/*   printf("leaving butdown\n"); */
}
/* ********************************* */
void accis_butup(event)
XEvent *event;
{
  extern int butup_();
/*   printf("entering butup\n"); */
  butup_(&xeye,&yeye,&zeye);
}
/* ********************************* */
void accis_moved(event)
XEvent *event;
{
  extern int viewrot_();
  extern int cubeupd_();
  float xmoved,ymoved;
  /* If we started moving without setting eye nonzero, set it first*/
    if(xeye0==0. && yeye0==0. && zeye0==0.) accis_butdown(event); 
  /* need screen to normal scaling not done here at present*/ 
  xmoved=event->xbutton.x-accis_x0;
  ymoved=event->xbutton.y-accis_y0;
/*    printf("eye: %f %f %f\n",xeye,yeye,zeye); */
  viewrot_(&xmoved,&ymoved,&xeye0,&yeye0,&zeye0,&xeye,&yeye,&zeye);
  cubeupd_(&xeye,&yeye,&zeye);
  /* This flush is necessary to see the black cube. */
  accisflush_();
}
/* ********************************* */
void accis_keypress(event)
XEvent *event;
{
  /* Testing only 
     printf("Key event keycode: %u keysym:%x",event->xkey.keycode,
	 (int)XLookupKeysym(&(event->xkey),0)); */
  ACCIS_SET_FOCUS;
}

/* ******************************************************************** */
/* Externally callable routine to set noeye3d return value.
   Set to 9999 to disable. */ 
int noeye3d_(value)
     int *value;
{
  if(*value>1000)accis_eye3d=9999;
  else accis_eye3d=*value;
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

  XCopyArea(accis_display,accis_pixmap,accis_window,accis_gc,0,0,     \
	    accis_gwa.width,accis_gwa.height,0,0);                    \
  XFlush(accis_display);
  /* Wait for a key press */
  if(accis_eye3d != 9999){
    if(XPending(accis_display)){
      XPeekEvent(accis_display,&event);
      if(event.type==KeyPress){ /* Halt continuous running on keypress */
	*value=(int)XLookupKeysym(&(event.xkey),0);
	accis_eye3d=9999;
      }else{ /* Discard other events */
	XNextEvent(accis_display,&event);
      }
    }
    if(accis_eye3d!=9999) {*value=accis_eye3d;return 0;}
  }
  ACCIS_SET_FOCUS;
  do{
    XNextEvent(accis_display,&event);
/*     printf("First loop:The event type: %d, %d\n",event.type,ButtonPress); */
    switch(event.type) {
    case Expose: EXPOSE_ACTION; break;
    case ButtonPress: accis_butdown(&event); break;
    case MotionNotify: accis_moved(&event); break;
    case ButtonRelease:
    case KeyPress: accis_keypress(&event); break;
    default: break;
    }
  }while(event.type != ButtonPress && event.type != KeyPress);
  /* Recognize KeyPress as sign to exit.*/
  if(event.type == KeyPress) {
    *value=(int)XLookupKeysym(&(event.xkey),0);
    if(XPending(accis_display)) XPeekEvent(accis_display,&event);
    /*    printf("Got value: %d, Event.type: %d, Pending? %d\n"
	   ,*value,event.type,XPending(accis_display));
  /* Get all the queued contiguous KeyPress events so we don't over-
     run the rotation when the key is lifted. The problem here is when
     any KeyPress events are buried in other event types. There seem to
     be such events under some circumstances. They are NoExpose.
     These come from XCopyArea calls in accisflush when doing animations.
  */
    if(XPending(accis_display)) XPeekEvent(accis_display,&event);
    while(XPending(accis_display) &&
	  (event.type==KeyPress || event.type==KeyRelease
	   || event.type==NoExpose 
	   )){
      XNextEvent(accis_display,&event);
      /*      printf("Next event: %d, %d ",event.type,
	      (int)XLookupKeysym(&(event.xkey),0));*/
      if(XPending(accis_display)) XPeekEvent(accis_display,&event);
    } 
     /* printf("Key value=%d\n",*value);  */
    return *value;
  }
  do{
    XNextEvent(accis_display,&event);
/*        printf("The event type: %d\n",event); */
    switch(event.type) {
    case Expose: EXPOSE_ACTION; break;
    case ButtonPress: accis_butdown(&event); break;
    case ButtonRelease: accis_butup(&event); break;
    case KeyPress: accis_keypress(&event); break;
    case MotionNotify: accis_moved(&event); break;
    case EnterNotify:
    default: break;
    }
  }while(event.type != KeyPress && event.type != ButtonRelease  );
  if( event.type != ButtonRelease ||
      ( xeye==xeye0 && yeye==yeye0 && zeye==zeye0) ) *value=0; else *value=1;
  return *value;
}
/* ******************************************************************** */
/* Routines for using 240 color gradients in preference to 16 fixed colors.*/
/************** Setup The Gradient **********************/
void accisgradinit_(r1,g1,b1,r2,g2,b2)
     FORT_INT *r1,*g1,*b1,*r2,*g2,*b2;
     /* On 64 bit machine this seems to need int *r1,*g1,*b1,*r2,*g2,*b2;*/
     /* RGB are specified in the range 0 to 65535 */
{
  int i,j,status;
  XColor theRGBcolor;
  if(accis_nodisplay){
    status=1;
  }else{
    if( (status=is_truecolor()) ==0){
      fprintf(stderr,"True Color 24 bit Status false: pixel:%lu\n",
	      theRGBcolor.pixel);
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
void accisgraddef_()
{ 
  FORT_INT top=65535;
  FORT_INT bot=0;
  accisgradinit_(&bot,&bot,&bot,&top,&top,&top);
}
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
}
/********** Tell the ipixel rgb color ********************************/
/* Subroutine */ int getrgbcolor_(ipixel,red,green,blue)
FORT_INT *ipixel,*red,*green,*blue;
{
  int ii;
    /*If not inited already, init by default.*/
  if(!a_grad_inited)accisgraddef_();
  /* Note that for fortran calls ipixel goes from 1 to 240*/
  /* No that seems to be incorrect 0 to 239 is right */
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
/************************************************************************/
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
void accisgradset_(red,green,blue,npixel)
     FORT_INT *red,*green,*blue,*npixel;
{
  int i;
  XColor theRGBcolor;
  if(*npixel!=a_gradPixno) 
    fprintf(stderr,"accisgradset ERROR: Incorrect array length:%d\n",*npixel);
  for (i=0;i<a_gradPixno;i++){
    a_gradred[i]=theRGBcolor.red=ilimit(0,*(red+i),65535);
    a_gradgreen[i]=theRGBcolor.green=ilimit(0,*(green+i),65535);
    a_gradblue[i]=theRGBcolor.blue=ilimit(0,*(blue+i),65535);
    /*printf("theRGBcolor %d,%d,%d\n",theRGBcolor.red,theRGBcolor.green,theRGBcolor.blue);*/
    if(accis_nodisplay!=1){
      if(is_truecolor()){
	a_gradPix[i]=
	  ((theRGBcolor.red/256)*256+(theRGBcolor.green/256))*256
	  +(theRGBcolor.blue/256);
	/*fprintf(stderr,"Allocated Color %d =%d\n",i,a_gradPix[i]);*/
      }else if(XAllocColor(accis_display,accis_colormap,&theRGBcolor)){
	/*fprintf(stderr,"Allocated Color %d=%ld\n",i,theRGBcolor.pixel);*/
	a_gradPix[i]=theRGBcolor.pixel;
      }else{
	a_gradPix[i]=BlackPixel(accis_display,0);
      }
    }
  }
  a_grad_inited=2;
}
/*************************************************************************/
/* Dummy, returning 0 defaults back to vec rendering. */
int igradtri_() { return 0;}
/*************************************************************************/
/*************************************************************************/
/* G77 fortran callable usleep */
void usleep_(usecs)
     long *usecs;
{
  usleep(  (unsigned long) *usecs);
}
/************************************************************************/
/* These are for compatibility with code written for OpenGL only */
/* Switch to drawing to the back buffer only */
void glback_()
{accis_back=1;
}
/************************************************************************/
/* Bring back buffer to front and return to writing there.*/
void glfront_()
{
  accis_back=0;
  XCopyArea(accis_display,accis_pixmap,accis_window,accis_gc,0,0,     \
	    accis_gwa.width,accis_gwa.height,0,0);                    \
  XFlush(accis_display);
}
