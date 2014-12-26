/* glX driver for accis plotting */
/* Fortran callable routines */
/* This version based on OpenGL and glx. */
int accis_driver_(){return 3;}
/* ********************************************************************** */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <X11/X.h>
#include <X11/Xlib.h>
#include    <X11/Intrinsic.h> 
#include    <X11/Core.h> 
#include <GL/gl.h>
#include <GL/glx.h>
#include <GL/glu.h>

/* Argument-passing typdef */
typedef int FORT_INT;

/* Globals for these routines*/
static Display *accis_display=NULL;
static Window accis_root;
static Window accis_window;
/*static Pixmap accis_pixmap;*/
/* Attributes to require of the visual chosen:*/
/* It proves advantageous to _avoid_ doublebuffering for incremental
   drawing. Speed seems best this way too. Yet there was a problem with
   single buffering from loki with glXSwapBuffers.*/
static GLint  accis_double[] = { GLX_RGBA, /* Truecolor and Directcolor */
			      GLX_DEPTH_SIZE, 24, /* Depth 24 */
			      GLX_DOUBLEBUFFER, /* */
			      None };
static GLint  accis_single[] = { GLX_RGBA, /* Truecolor and Directcolor */
			      GLX_DEPTH_SIZE, 24, /* Depth 24 */
				 /*GLX_DOUBLEBUFFER, */
			      None };
static XVisualInfo             *accis_vi;
static Colormap                accis_cmap;
static XSetWindowAttributes    accis_swa;
static GLXContext              accis_glc;
static XWindowAttributes       accis_gwa;
/*static XEvent                  accis_xev;*/
static Colormap accis_colormap;
static int accis_depth;
static int accis_listing=0;
static int accis_eye3d=9999;

/* Static maximum number of points in path */
#define accis_path_max 4000
static XPoint accis_path[accis_path_max];
static int accis_pathlen=0;
/* Until svga is called, the default is that there is no display */
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
/*static int a_gradno=a_gradPixno;*/ /* Publically available Pixno */

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
#define ACCIS_SET_FOCUS	  accis_set_focus()
void accis_set_focus(){
  if(accis_old_handler==0)
    accis_old_handler=XSetErrorHandler(accis_errorhandler) ;
    XSetInputFocus(accis_display, accis_window, RevertToParent,CurrentTime);
}
/* Main setup subroutine */ 
int svga_(scrxpix, scrypix, vmode, ncolor)
FORT_INT *scrxpix, *scrypix, *vmode, *ncolor;
{  
  static int second=0;
  extern void initDefaultColors();
  XSizeHints hints;
  int x_size,y_size,x_off,y_off,gravity;
  accis_nodisplay=0;
  *ncolor=15;
  *vmode=88;
  if(second == 0) {
  /* Call fortran routine to get arguments into accis_ globals*/
    getcmdargs_();
    if( (accis_display = XOpenDisplay(NULL)) == NULL) {
        printf("\n\tcannot connect to X server\n\n");
        exit(0); 
    }
    if((accis_vi=glXChooseVisual(accis_display, 0, accis_double)) == NULL) {
      printf("\n\tno appropriate visual found\n");
      if((accis_vi=glXChooseVisual(accis_display, 0, accis_single)) == NULL) {
        exit(0); 
      }else{
	printf("\tfell-back to single-buffering\n");
	printf("\tGLX visual %#x selected\n", (int)accis_vi->visualid); 
      }
    }else{
      printf("\tGLX visual %#x selected\n", (int)accis_vi->visualid); 
    }/* %p hexadecimal*/

/* Start of Xlib calls ******************/
    accis_root = DefaultRootWindow(accis_display);
    accis_cmap = XCreateColormap(accis_display, accis_root, accis_vi->visual, 
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
			       accis_vi->depth, 
			       InputOutput, accis_vi->visual, 
		CWBackPixel | CWColormap | CWEventMask, &accis_swa);
    XMapWindow(accis_display, accis_window);
    XStoreName(accis_display, accis_window, "Accis");

/*Start of OpenGL calls ******************/
    accis_glc = glXCreateContext(accis_display, accis_vi, NULL, GL_TRUE);
    glXMakeCurrent(accis_display, accis_window, accis_glc);
    /* All writes into both buffers at once  */
    glDrawBuffer(GL_FRONT_AND_BACK); 
    /* reads from the back buffer */
    glReadBuffer(GL_BACK);

    accis_depth=accis_vi->depth;
    accis_colormap=accis_cmap;
    initDefaultColors();
    second++;

  }else{
  } 
  XGetWindowAttributes(accis_display,accis_window, &accis_gwa);

#define pixscale 1
  *scrxpix=accis_gwa.width*pixscale;
  *scrypix=accis_gwa.height*pixscale;
  /*   printf("Called XGet. Height,Width=%d, %d\n",*scrxpix,*scrypix); */
  glViewport(0, 0, accis_gwa.width, accis_gwa.height);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  /*  printf("scrxpix=%d,scrypix=%d\n",*scrxpix,*scrypix);*/
  glOrtho(0.,(float)(*scrxpix-1),(float)(*scrypix-1),0.,1.,20.);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(0., 0., 10., 0., 0., 0., 0., 1., 0.);

  /* Start the main drawing list. Clear the frame. */
  /* Always clearing gives behavior equivalent to vecx. Otherwise, if
     pltend() is not called then drawing actions build up in the listing
     and are redrawn on refresh. There might be some cases where that 
     behavior is preferable, but such cases are not yet identified.
     if(accis_listing==0)*/
  {
    accis_listing=1;
    glNewList(1,GL_COMPILE_AND_EXECUTE);
    glClearColor(1.0, 1.0, 1.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  }

  glColor3f(0.,0.,0.);
  return 0;
}

/************************************************************************/
void accisinit()
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
  int status;
  int i,j,pix;
  if(accis_nodisplay){
  }else{
    if(is_truecolor()){	
/*      fprintf(stderr,"True color shortcut Default colors.\n"); */
    }else{
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
      glCallList(1);						      \
      /*glXSwapBuffers(accis_display,accis_window); Is this needed?*/ 	\
      glFlush(); /* Proves to be necessary for remote servers. */	\

/* ******************************************************************** */
/* End plotting Subroutine */ 
void txtmode_()
{
  XEvent event;
  glEndList(); /* Close the drawing list started in svga.*/
  accis_listing=0;
  ACCIS_SET_FOCUS;
  EXPOSE_ACTION;
  glFlush();
  do{
    /*    printf("Executing XtNextEvent"); */
    XNextEvent(accis_display,&event);
/*     printf("Event: type=%d\n",event.type); */
    if(event.type == Expose){EXPOSE_ACTION;}
  }while(event.type != ButtonPress && event.type != KeyPress );
    /* We don't do this; but here's how to terminate cleanly.
    glXMakeCurrent(accis_display, None, NULL);
    glXDestroyContext(accis_display, accis_glc);
    XDestroyWindow(accis_display, accis_window);
    XCloseDisplay(accis_display);  */
}
/* ******************************************************************** */
void accisclear_()  /* Simply clear to background fortran callable*/ 
{    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); }

/* ******************************************************************** */
/* Flush the plot buffer Subroutine */ 
void accisflush_()
{
  XFlush(accis_display);  
  glFlush();
}
/* ********************************************************************* */
/* Set the Color Subroutine */ int scolor_(li)
FORT_INT *li;
{
  /* *ncolor=*li; */
  if((*li < a_maxPixels) && (*li >= 0)){
    glColor3ubv(&accis_rgb[(int) *li*3]);
/*     printf("color3=%d,%d,%d\n", */
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
/*     printf("px %.0f, py %.0f\n",px2,py2);
    if( *ud > 0) { */
    if( *ud == 1) { 
      glBegin(GL_LINES);
      glVertex2f(px1,py1);
      glVertex2f(px2,py2);
      glEnd();
      if(accis_pathlen<accis_path_max){      /* Add point to path */
	accis_pathlen++;
      }
    }else{ 
      if(*ud==-1){
	glBegin(GL_POINTS);
	glVertex2f(px2,py2);
	glEnd();
      }
      /* Restart path */
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
      glBegin(GL_POLYGON);
      for(i=0;i<=accis_pathlen;i++){
	glVertex2i(accis_path[i].x,accis_path[i].y);
      }
      glEnd();
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
  glFlush();
  /* This is for the background case working around GL bugs. */
  /* glXSwapBuffers(accis_display, accis_window);*/
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
void noeye3d_(value)
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
  glEndList(); /* Close the drawing list started in svga.*/
  accis_listing=0;
  if(accis_eye3d == 1){ *value=0; return 0; }

  glXSwapBuffers(accis_display, accis_window);
  XFlush(accis_display);
  /* Wait for a key press */
  if(accis_eye3d != 9999){
    if(XPending(accis_display)){
      XPeekEvent(accis_display,&event);
      if(event.type==KeyPress){
	*value=(int)XLookupKeysym(&(event.xkey),0);
	accis_eye3d=9999;
      }
    }else{
      *value=accis_eye3d; return 0; 
    }
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
    /* Get all the queued contiguous KeyPress events so we don't over-
       run the rotation when the key is lifted. */
    if(XPending(accis_display)) XPeekEvent(accis_display,&event);
    while(XPending(accis_display) &&
	  (event.type==KeyPress || event.type==KeyRelease) )
      {
	XNextEvent(accis_display,&event);
	if(XPending(accis_display)) XPeekEvent(accis_display,&event);
      } 
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
GLushort ared,agreen,ablue;
    /*Default gradient is a gray-scale.*/
    /*If not inited already init by default.*/
  if(!a_grad_inited)accisgraddef_();
  /* *ncolor=*li; */
  ii=( (*li < a_gradPixno) ? ( *li >= 0 ? *li : 0) : a_gradPixno-1 );
  ared=a_gradred[ii];
  agreen=a_gradgreen[ii];
  ablue=a_gradblue[ii];
  /*  printf("ared=%d,agreen=%d,ablue=%d\n",ared,agreen,ablue);*/
  glColor3us(ared,agreen,ablue);
  return 1; 
}
/********** Tell the current rgb color ********************************/
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
    if(accis_nodisplay!=1){/*Unless screen display is off, initialize gradPix.*/
      if(is_truecolor()){
	a_gradPix[i]=
	  ((theRGBcolor.red/256)*256+(theRGBcolor.green/256))*256
	  +(theRGBcolor.blue/256);
	/*fprintf(stderr,"Allocated Color %d =%d\n",i,a_gradPix[i]);*/
      }else if(XAllocColor(accis_display,accis_colormap,&theRGBcolor)){
	fprintf(stderr,"Allocated Color %d=%ld\n",i,theRGBcolor.pixel);
	a_gradPix[i]=theRGBcolor.pixel;
      }else{
	a_gradPix[i]=BlackPixel(accis_display,0);
      }
    }
  }
  a_grad_inited=2;
}
/*************************************************************************/
/* Dummy, returning 0 defaults back to vec rendering.
int igradtri_() { return 0;} */
/*************************************************************************/
/* Fortran callable GL triangle drawing. Return 1 for success. 
i3d indicates 3d x,y,z if 1, 2d x,y if 0. color in h. 
This routine does not (yet) leave the normalized coordinates at the 
position of the last vertex. This is incompatible with the final postion
of the vecx driver that does not use igradtri.
*/
int igradtri_(x,y,z,h,i3d)
float *x,*y,*z,*h;
FORT_INT *i3d;
{
  int i,li;
  float xn,yn,zn,xw,yw,xs,ys;
  FORT_INT ixs,iys;
  float cbx,cby,cbz,xcbc,ycbc;
  float zs, zero=0.;
  /* Fortran functions called: */
  float extern wx2nx_(),wy2ny_(),getwx2nx_(),getwy2ny_();
  void extern tn2s_(),wxyz2nxyz_(),getcube_(),trn32_();

/*   if(h[1]>=a_gradPixno || h[3]>=a_gradPixno || h[3]>=a_gradPixno){ */
/*     return 1; */
/*   } */
  if(*i3d==0){/* 2-D case */
/*     printf(" 2D case\n"); */
    glBegin(GL_TRIANGLES);
    for (i=0;i<3;i++){
      xw=x[i]; yw=y[i];
      /* Real fortran _function_ calls don't work with g77 on AMD64
      unless one uses the extra flag -fno-f2. This is fraught with
      danger.  Probably better just to circumvent the whole thing by
      using fortran subroutines only. So instead of
      xn=wx2nx_(&xw);
      yn=wy2ny_(&yw); */ /*The following subroutines replace and work.*/
      getwx2nx_(&xw,&xn);
      getwy2ny_(&yw,&yn);
      tn2s_(&xn,&yn,&ixs,&iys);
      xs=ixs; ys=iys;
      li=(int) h[i] +1; /* +1 compensates for adj in gradtri */
      acgradcolor_(&li);
      glVertex2f(xs,ys);
    }
    glEnd();
    return 1;
  }else{
/*     printf(" 3D case\n"); */
    getcube_(&cbx,&cby,&cbz,&xcbc,&ycbc);
    glBegin(GL_TRIANGLES);
    for (i=0;i<3;i++){
      wxyz2nxyz_(x+i,y+i,z+i,&xn,&yn,&zn);
      trn32_(&xn,&yn,&zn,&xs,&ys,&zs,&zero);
      xn=xs+xcbc;
      yn=ys+ycbc;
      tn2s_(&xn,&yn,&ixs,&iys);
      xs=ixs; ys=iys;
      li=(int) h[i];
      acgradcolor_(&li);
      glVertex2f(xs,ys);	
    }
    glEnd();
    return 1;
  }
}
/*************************************************************************/
/* G77 fortran callable usleep */
void usleep_(usecs)
     long *usecs;
{
  usleep(  (unsigned long) *usecs);
}
/************************************************************************/
/* Switch to drawing to the back buffer only */
void glback_()
{glDrawBuffer(GL_BACK);
}
/************************************************************************/
/* Bring back buffer to front and return to writing there.*/
void glfront_()
{
  glFlush();
  /* The order of the next two statements is key to getting all drawing */
  /* but to avoid the X bugs we don't for now do this:
     glDrawBuffer(GL_FRONT_AND_BACK); */
  glXSwapBuffers(accis_display,accis_window); /* Should never be called*/
}
