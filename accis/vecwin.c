/* Now up to date except for the file-only code  Mar 05*/
/* Modified Sep 2017 to get lorentztrack animation working */

#include <windows.h>
#include <stdio.h>
#include <ctype.h>

#define accis_path_max 4000
POINT accis_path[accis_path_max];
int accis_pathlen=0;
/* Until svga is called, the default is that there is no display */
int accis_nodisplay=1;

#define accis_maxPixels 16
COLORREF accis_pixels[accis_maxPixels]=
{ 0x00FFFFFF,
  0x00770000,  
  0x00007700,  
  0x00777700,  
  0x000000AA,  
  0x00770077,  
  0x000077AA,
  0x00777777,  
  0x00AAAAAA,  
  0x00FF0000,  
  0x0000FF00,  
  0x00FFFF00,  
  0x000000FF,  
  0x00FF00FF,  
  0x0000FFFF,  
  0x00000000,  
};

/* Unused on windows.
char *accis_colornames[accis_maxPixels]=
{
  "White",
  "MediumBlue",
  "SeaGreen",
  "MediumTurquoise",
  "Firebrick",
  "Orchid",
  "Brown",
  "LightGrey",
  "SlateGrey",
  "Blue",
  "Green",
  "Cyan",
  "Red",
  "Magenta",
  "Yellow",
  "Black"
};
*/

#define ACCIS_DEFWIDTH 800
#define ACCIS_DEFHEIGHT 600

static char g_szClassName[] = "MyWindowClass";
static HINSTANCE g_hInst = NULL;
static int wWidth, wHeight; /* Window width and height */
static RECT wRect;
PAINTSTRUCT ps;
HDC hdcMemory, hdcWindow;
BITMAP bm;
HBITMAP hbm; 
int first=1;
static HBRUSH hbr;
static HBRUSH accis_brush;
static HPEN accis_pen;
  
static WNDCLASSEX WndClass;
static HWND hwnd;
static MSG Msg;
static int accis_eye3d=9999;

LRESULT CALLBACK WndProc(HWND hwnd, UINT Message, WPARAM wParam, LPARAM lParam)
{

   switch(Message){
   case WM_CREATE:
     break;
   case WM_SIZE:  
     /* This is called with 0,0 when iconized. That is bad. So don't resize
	now.*/
     wWidth = LOWORD(lParam);  // width of client area 
     wHeight = HIWORD(lParam); // height of client area 
/*       printf("Resized: %d %d\n",wWidth,wHeight); */
     /*     AccisRenewBitmap(hwnd);*/
     break;
   case WM_PAINT:
/*       printf("Paint BitBlt\n"); */
     /*Trying for private DC*/
/*       hdcWindow = BeginPaint(hwnd, &ps); */
     BitBlt(hdcWindow, 0, 0, wWidth, wHeight, hdcMemory, 0, 0, SRCCOPY);
/*         EndPaint(hwnd, &ps); */
     GetClientRect(hwnd,&wRect);
     ValidateRect(hwnd,&wRect);/* Needed without endpaint. Not perfect.
				Resize leads to continuous BitBlt.
			       Needs new/current window size.*/
     break;
   case WM_CLOSE:
       /* It is not safe to do the following from the right-hand X.
	  Makes the program hang. The Default is also bad so act the same.*/
/*       DestroyWindow(hwnd); */
/*       break; */
   case WM_CHAR: /* Recognize keyboard events as proceeding */
   case WM_DESTROY:
   case WM_LBUTTONDOWN: /* Click in window to exit from this message loop. */
     DeleteDC(hdcMemory);
     PostQuitMessage(0);
     break;
   default:
     return DefWindowProc(hwnd, Message, wParam, lParam);
   }
   return 0;
}


/**********************************************************************/
/* Delete the old bitmap and make a new one. */
void AccisRenewBitmap(HWND hwnd){
/*      WPARAM wParam; */
/*      LPARAM lParam; */
     if(!first){
/*         printf("Deleting hdcMemory and hbm\n"); */
       DeleteDC(hdcMemory);
       DeleteObject(hbm);
     }else{
       first=0;
     }
     GetClientRect(hwnd, &wRect);
     hdcWindow = GetDC(hwnd);
     /* Hoped this is not nec when we use GetDC instead. No worse.*/
/*         hdcWindow = BeginPaint(hwnd, &ps);   */
     hdcMemory = CreateCompatibleDC(hdcWindow);            
/*       printf("Creating bitmap, %d %d\n",wWidth,wHeight); */
     /* This needs to use the window device context to get color.*/
     hbm=CreateCompatibleBitmap(hdcWindow,wWidth,wHeight);
/*         EndPaint(hwnd, &ps);  */
     SelectObject(hdcMemory,hbm);
     GetObject(hbm,sizeof(bm),&bm);
     /* Have to set bitmap background. Not set automatically.*/
     hbr=CreateSolidBrush(0x00FFFFFF);
     /* Removed this to rely on WM_PAINT which ought to be called auto.*/
/*       BitBlt(hdcWindow, 0, 0, wWidth, wHeight,hdcMemory,0,0,SRCCOPY); */
     FillRect(hdcMemory,&wRect,hbr);
     InvalidateRect(hwnd,&wRect,TRUE); /* Force whole window draw */
     /*This does not seem to be necessary once we invalidate */
/*         WndProc(hwnd,WM_PAINT,wParam,lParam);  */
}
/**********************************************************************/
/* svga does the job of WinMain setting up the window etc.*/

int svga_(int *scrxpix, int *scrypix, int *vmode, int *ncolor)
{
  static int second=0;
  /*  extern int f__xargc;
      extern char **f__xargv;
      int svga_argc=0;
      char *svga_argv[1]; 
  */
  /* For some bizarre reason the following two definitions are important
     even though they do not seem to be used. Without them the window does
     not pop up.
  */
  HINSTANCE hInstance;
  HINSTANCE hPrevInstance; 
  LPSTR lpCmdLine; 
  /*nCmdShow tells how to open window */
  int nCmdShow=SW_SHOWNORMAL;

  accis_nodisplay=0;
  g_hInst=GetModuleHandle(NULL);
/*    printf("Entered svga\n"); */
  if(second == 0){
    WndClass.cbSize        = sizeof(WNDCLASSEX);
    WndClass.style         = CS_OWNDC;
    WndClass.lpfnWndProc   = WndProc;
    WndClass.cbClsExtra    = 0;
    WndClass.cbWndExtra    = 0;
    WndClass.hInstance     = g_hInst;
    WndClass.hIcon         = LoadIcon(NULL, IDI_APPLICATION);
    WndClass.hCursor       = LoadCursor(NULL, IDC_ARROW);
    WndClass.hbrBackground = (HBRUSH)(COLOR_BTNFACE+1);
    WndClass.lpszMenuName  = NULL;
    WndClass.lpszClassName = g_szClassName;
    WndClass.hIconSm       = LoadIcon(NULL, IDI_APPLICATION);
    
    if(!RegisterClassEx(&WndClass))
      {
	MessageBox(0, "Window Registration Failed!", "Error!",
		   MB_ICONEXCLAMATION | MB_OK | MB_SYSTEMMODAL);
	return 0;
      }
    
    hwnd = CreateWindowEx(
			  WS_EX_CLIENTEDGE,
			  g_szClassName,
			  "Accis Window",
			  WS_OVERLAPPEDWINDOW,
			  CW_USEDEFAULT, CW_USEDEFAULT,
			  ACCIS_DEFWIDTH+12, ACCIS_DEFHEIGHT+31,
			  NULL, NULL, g_hInst, NULL);
    
    if(hwnd == NULL)
      {
	MessageBox(0, "Accis Window Creation Failed!", "Error!",
		   MB_ICONEXCLAMATION | MB_OK | MB_SYSTEMMODAL);
	return 0;
      }
    
    /*      printf("About to show window\n"); */
    ShowWindow(hwnd, nCmdShow);
/*      printf("About to update window\n"); */
    UpdateWindow(hwnd);

    second++;
  }

  /* Originally this was implicit in resize. mingw did not like it.*/
  AccisRenewBitmap(hwnd);
  
/*    printf("Entering peek message loop. Second=%d\n",second); */
  while(PeekMessage(&Msg, hwnd, 0, 0, PM_REMOVE))
    {
      TranslateMessage(&Msg);
      DispatchMessage(&Msg);
    }
  hdcWindow = GetDC(hwnd); /*Ensure we have still selected the window.*/

  *vmode=88;
  *ncolor=15;
  *scrxpix=wWidth;
  *scrypix=wHeight;
  
  return Msg.wParam;
}
/*************************************************************************/
int txtmode_()
{
  while(GetMessage(&Msg, hwnd, 0, 0))
    {
      TranslateMessage(&Msg);
      DispatchMessage(&Msg);
    }
  return 0;
}
/* ******************************************************************** */
int vec_(long *px,long *py,long *ud)
{ /*  Draw vector on screen, with pen up or down. */
    static int px1=0,py1=0,px2=0,py2=0;
    extern POINT accis_path[];
    extern int accis_pathlen;

    px1=px2;
    py1=py2;
    px2 = *px;
    py2 = *py;
    if( *ud != 0) {
/*        MoveToEx(hdcWindow,px1,py1,NULL); */
/*        MoveToEx(hdcMemory,px1,py1,NULL); */
      LineTo(hdcWindow,px2,py2);
      LineTo(hdcMemory,px2,py2);
      if(accis_pathlen<accis_path_max){      /* Add point to path */
	accis_pathlen++;
      }
    }else{ /* Restart path */
      accis_pathlen=0;
      MoveToEx(hdcWindow,px2,py2,NULL);
      MoveToEx(hdcMemory,px2,py2,NULL);
    }
    accis_path[accis_pathlen].x=*px;
    accis_path[accis_pathlen].y=*py;

    return 0;
} /* vec_ */
/* ******************************************************************** */

int scolor_(long *li)
{
  if((*li < accis_maxPixels) && (*li >= 0)){
    DeleteObject(accis_pen);
    accis_pen=CreatePen(PS_SOLID,0,accis_pixels[(int)(*li)]);
    DeleteObject(accis_brush);
    accis_brush=CreateSolidBrush(accis_pixels[(int)(*li)]);
    SelectObject(hdcWindow,accis_pen);
    SelectObject(hdcMemory,accis_pen);
    SelectObject(hdcWindow,accis_brush);
    SelectObject(hdcMemory,accis_brush);
    return 1;
  }else{    
    return 0;
  }
}
/* ********************************************************************* */
int vecfill_()
{
  SetPolyFillMode(hdcWindow,WINDING);
  SetPolyFillMode(hdcMemory,WINDING);
  Polygon(hdcWindow,accis_path,accis_pathlen+1);
  Polygon(hdcMemory,accis_path,accis_pathlen+1);
  return 0;
}
/* ********************************************************************* */

/* Rotation code*/
float xeye,yeye,zeye;
float xeye0,yeye0,zeye0;
float accis_x0,accis_y0;

void extern butdown_();
void extern viewrot_();
void extern cubeupd_();
void extern butup_();

void accis_butdown(lMsg)
     MSG lMsg;
{ /* accis_butdown */
  accis_x0=LOWORD(lMsg.lParam);
  accis_y0=HIWORD(lMsg.lParam);
  butdown_(&xeye0,&yeye0,&zeye0);
  xeye=xeye0; yeye=yeye0; zeye=zeye0;
  }
/**********************************************************************/
int cmpmessage(UINT message)
{
   GetMessage(&Msg, hwnd, 0, 0);
   if(Msg.message == message){   return 1;   }else{  return 0;   }
}
/**********************************************************************/
int eye3d_(value)
     int *value;
{
  float xmoved,ymoved;
  if(accis_nodisplay){ *value=0; return 0; }
  if(accis_eye3d == 1){ *value=0; return 0; }  

  /* printf("accis_eye3d %d\n",accis_eye3d); */

  /*Disabled response route may not be working*/
  if(accis_eye3d != 9999){
    if(PeekMessage(&Msg,hwnd,0,0,PM_REMOVE)){
      if(Msg.message == WM_KEYDOWN){
	accis_eye3d=9999;
      }else{
	*value=accis_eye3d; return 0;
      }
    }else{
      *value=accis_eye3d; return 0;
    }
  }

  GetMessage(&Msg,hwnd,0,0);
  while( (Msg.message != WM_LBUTTONDOWN) && (Msg.message != WM_KEYDOWN)){
    /*  while(!cmpmessage(WM_LBUTTONDOWN)){*/
    TranslateMessage(&Msg);
    DispatchMessage(&Msg);
    GetMessage(&Msg,hwnd,0,0);
  }
/*      printf("Got first button press...");  */
  if(Msg.message == WM_LBUTTONDOWN){
    accis_butdown(Msg);
    while(!cmpmessage(WM_LBUTTONUP)){
      /* If we started moving without setting eye nonzero, set it first*/
      if(xeye0==0. && yeye0==0. && zeye0==0.) accis_butdown(Msg); 
      xmoved=LOWORD(Msg.lParam)-accis_x0;
      ymoved=HIWORD(Msg.lParam)-accis_y0;
/*          printf("eye: %f %f %f\n",xeye,yeye,zeye);  */
      viewrot_(&xmoved,&ymoved,&xeye0,&yeye0,&zeye0,&xeye,&yeye,&zeye);
      cubeupd_(&xeye,&yeye,&zeye);
      {/* Full redraw section not working.
	PostMessage(hwnd,WM_LBUTTONDOWN,0,0);
	break; */
      }
    }
/*      printf("Exited eye3d loop\n");  */
    butup_(&xeye,&yeye,&zeye);
    *value=0;
    if(  !( xeye==xeye0 && yeye==yeye0 && zeye==zeye0) ) *value=1;
    return 0;
  }else{/* KeyPress event.*/

    /*left up right down 37 38 39 40 -> 65361 -64*/
    if(Msg.wParam >36 && Msg.wParam < 41){
      *value=Msg.wParam+ 65324;
    }else{
    /* tolower equivalent */
    *value= (Msg.wParam > 96 ? Msg.wParam : Msg.wParam+32) ;    
    }
    /*printf("Msg.wParam=keycode=%d\n",Msg.wParam);*/
    return Msg.wParam;
  }
}
/* 10 May 2003 got rotation working.*/
/*******************************************************************/

/* 256 color gradient globals */
/* Use only 240 colors so that 16 are left for the 16 colors and you can
   make gifs out of the resultant without a problem.*/
#define a_gradPixno 240

COLORREF a_gradPix[a_gradPixno];

unsigned a_gradcurrent;
unsigned a_grad_inited=0;
unsigned a_gradred[a_gradPixno];
unsigned a_gradgreen[a_gradPixno];
unsigned a_gradblue[a_gradPixno];

/* ******************************************************************** */
/* Routines for using 240 color gradients in preference to 16 fixed colors.*/
/************** Setup The Gradient **********************/
int accisgradinit_(r1,g1,b1,r2,g2,b2)
     long *r1,*g1,*b1,*r2,*g2,*b2;
     /* On entry RGB are specified in the range 0 to 65535 */
     /* But we divide each by 256 to get down to a 256 color scale,
	which seems to be the easiest windows default. */
{
  int i,j;
  for (i=0;i<a_gradPixno;i++){
    j=(i* *r2+(a_gradPixno-1-i)* *r1)/(a_gradPixno-1.) ;
    if(j<0)j=0; else if(j>65535)j=65535;
    a_gradred[i]=j;
    j=(i* *g2+(a_gradPixno-1-i)* *g1)/(a_gradPixno-1.) ;
    if(j<0)j=0; else if(j>65535)j=65535;
    a_gradgreen[i]=j;
    j=(i* *b2+(a_gradPixno-1-i)* *b1)/(a_gradPixno-1.) ;
    if(j<0)j=0; else if(j>65535)j=65535;
    a_gradblue[i]=j;
    a_gradPix[i]=
      a_gradred[i]/256+256*(a_gradgreen[i]/256 +256*(a_gradblue[i]/256));
    /*fprintf(stderr,"Gradinit: a_gradPix(%d)= %d\n",i,a_gradPix[i]);*/
  }
  a_grad_inited=2;
  return 0;
}
/******************************************************************/
static int ilimit(b,i,t)
     int b,i,t;
{
  int j;
  j=(i<b ? 0 : i);
  j=(t<j ? t : j);
  return j;
}
/*****************Set the color gradient from arrays**************************/
int accisgradset_(red,green,blue,npixel)
   int *red,*green,*blue,*npixel;
{  unsigned long i,j;
  if(*npixel!=a_gradPixno) 
    fprintf(stderr,"accisgradset ERROR: Incorrect array length:%d\n",*npixel);
     /* RGB are specified in the range 0 to 65535 */
  for (i=0;i<a_gradPixno;i++){
     a_gradred[i]=ilimit(0,*(red+i),65535);
    a_gradgreen[i]=ilimit(0,*(green+i),65535);
    a_gradblue[i]=ilimit(0,*(blue+i),65535);
    a_gradPix[i]=
      a_gradred[i]/256+256*(a_gradgreen[i]/256 +256*(a_gradblue[i]/256));
    /*    printf("red,green,blue: %d, %d, %d\n",
	   a_gradred[i],a_gradgreen[i],a_gradblue[i]);
	   printf("a_gradPix[%d]=%6x\n",i,a_gradPix[i]);  */
  }
  a_grad_inited=1;
  return 0;
}
/************** Setup Default Gradient Grey-scale **********************/
void accisgraddef_(red,green,blue,npixel)
     int *red,*green,*blue,*npixel;
{
  unsigned long i,j;
     /* RGB are specified in the range 0 to 65535 */
  for (i=0;i<a_gradPixno;i++){
    j=(i* 65535)/(a_gradPixno-1.) ;
    /*        a_gradPix[i]=white*i/(a_gradPixno-1.); */
    a_gradred[i]=j;
    a_gradgreen[i]=j;
    a_gradblue[i]=j;
    a_gradPix[i]=
      a_gradred[i]/256+256*(a_gradgreen[i]/256 +256*(a_gradblue[i]/256));
    /*    printf("red,green,blue: %d, %d, %d\n",
	   a_gradred[i],a_gradgreen[i],a_gradblue[i]);
	   printf("a_gradPix[%d]=%6x\n",i,a_gradPix[i]);  */
  }
  a_grad_inited=1;
  return 0;
}

/**********************************************************************/
int acgradcolor_(long *li)
{
  if(!a_grad_inited)accisgraddef_();
  if((*li < a_gradPixno) && (*li >= 0)){
    DeleteObject(accis_pen);
    accis_pen=CreatePen(PS_SOLID,0,a_gradPix[(int)(*li)]);
    DeleteObject(accis_brush);
    accis_brush=CreateSolidBrush(a_gradPix[(int)(*li)]);
    SelectObject(hdcWindow,accis_pen);
    SelectObject(hdcMemory,accis_pen);
    SelectObject(hdcWindow,accis_brush);
    SelectObject(hdcMemory,accis_brush);
    a_gradcurrent=*li;
    return a_gradPix[(int)(*li)];
  }else{
    a_gradcurrent=a_gradPixno;
    return 0;
  }
}
/**********************************************************************/
int getrgbcolor_(long *i, long *ired, long *igreen, long *iblue)
{ 
  int ic;
  /*  ic=a_gradcurrent;
   *i=(long) ic; This would probably be more rational but ...*/
  ic=*i;
  if(ic >= a_gradPixno) ic=a_gradPixno-1;
  if(ic < 0) ic=0;
  *ired=a_gradred[ic];
  *igreen=a_gradgreen[ic];
  *iblue=a_gradblue[ic];
  return 0;
}
/**********************************************************************/
/* Subroutine for dummy svga call when using no display */ 
int svganodisplay_(int *scrxpix, int *scrypix, int *vmode, int *ncolor)
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
/***********************************************************************/
int accisclear_()
{
  AccisRenewBitmap(hwnd);
  return 0;
}
/***********************************************************************/
/* Windows version based on guess and experiment!*/
int accisflush_()
{
 if(PeekMessage(&Msg,hwnd,0,0,PM_REMOVE)){
    TranslateMessage(&Msg);
    DispatchMessage(&Msg);
}
  return 0;
}
/***********************************************************************/
#include <windows.h>	/* WinAPI */
/* Windows sleep in 100ns units from 
https://gist.github.com/Youka/4153f12cf2e17a77314c*/
BOOLEAN nanosleep(LONGLONG ns){
	/* Declarations */
	HANDLE timer;	/* Timer handle */
	LARGE_INTEGER li;	/* Time defintion */
	/* Create timer */
	if(!(timer = CreateWaitableTimer(NULL, TRUE, NULL)))
		return FALSE;
	/* Set timer properties */
	li.QuadPart = -ns;
	if(!SetWaitableTimer(timer, &li, 0, NULL, NULL, FALSE)){
		CloseHandle(timer);
		return FALSE;
	}
	/* Start & wait for timer */
	WaitForSingleObject(timer, INFINITE);
	/* Clean resources */
	CloseHandle(timer);
	/* Slept without problems */
	return TRUE;
}
/***********************************************************************/
/* G77 fortran callable usleep */
void usleep_(usecs)
     long *usecs;
{
  int msecs;
  int nsecs;
  msecs=*usecs/1000;
  nsecs=*usecs*10;

  nanosleep(nsecs);
  /*Sleep(msecs);*/
}
/* ******************************************************************** */
/* Externally callable routine to set noeye3d return value.
   Set to 9999 to disable. */ 
int noeye3d_(value)
     int *value;
{
  if(*value>1000)accis_eye3d=9999;
  accis_eye3d=*value;
  while (PeekMessage(&Msg, hwnd,  0, 0, PM_REMOVE)){};
  return 0;
}
/*************************************************************************/
/* Dummy here, returning 0 */
int igradtri_() { return 0;}
int accisrefresh_() { return 0;}
