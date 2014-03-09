/* Test glx visuals for compatibility */


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
static Pixmap accis_pixmap;
static XVisualInfo             *accis_vi;
static Colormap                accis_cmap;
static XSetWindowAttributes    accis_swa;
static GLXContext              accis_glc;
static XWindowAttributes       accis_gwa;
static XEvent                  accis_xev;
static Colormap accis_colormap;
static int accis_depth;
static int accis_listing=0;
static int accis_eye3d=9999;

/* Attributes to require of the visual chosen:*/
/* It proves advantageous to _avoid_ doublebuffering for incremental
   drawing. Speed seems best this way too. */
static GLint  accis_att[] = { GLX_RGBA, /* Truecolor and Directcolor */
			      /*GLX_ACCUM_RED_SIZE,8, /* Require accum */
			      GLX_DEPTH_SIZE, 24, /* Depth 24 */
			      /*GLX_DOUBLEBUFFER, /* */
			      None };
int main()
{
    if( (accis_display = XOpenDisplay(NULL)) == NULL) {
        printf("\n\tcannot connect to X server\n\n");
        exit(0); 
    }
    if((accis_vi=glXChooseVisual(accis_display, 0, accis_att)) == NULL) {
      printf("\n\tno appropriate visual found\n\n");
        exit(0); 
    }else{
      printf("\tGLX visual %#x selected\n", (int)accis_vi->visualid); 
    }/* %p hexadecimal*/



}
