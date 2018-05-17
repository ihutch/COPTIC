C********************************************************************
      block data initiald
c Generic parts of initializations.
      include 'plotcom.h'
      data nframe,nrows,ncolumns/0,0,0/
      data naxmin,naxmax,naymin,naymax,naxpt,naypt
     $  / 0.31,0.91,0.1,0.7,0.31,0.1/
      data xticlen,yticlen,xticoff,yticoff,nxlabw,nxlabp,nylabw,nylabp
     $  / 0.015,0.015,-0.03,-0.02,4,1,4,1 /
      data ticnum/6/
      data lxlog/.false./lylog/.false./lclog/.false./
      data lminor/.true./
c Now this is initialized by truncf call in pltinit.
c      data  trcxmi,trcxma,trcymi,trcyma,ltlog
c     $  / 0.,0.,1.,1.,.false. /
      data updown/99/
      data pfsw,pfilno/0,0/
      data pfPS/0/
      data vmode/111/
      end

C********************************************************************
      block data dashdata
      logical ldash
      real dashlen,dashdist
      integer MASKNO,dashmask,jmask
      parameter (MASKNO=4)
      dimension dashmask(MASKNO),dashlen(MASKNO)
      common/dashline/ldash,dashlen,dashdist,dashmask,jmask
      data dashmask/1,0,1,0/
      data ldash/.false./
      data dashlen/0.03,0.01,0.01,0.01/
      data dashdist/1.e-5/
      data jmask/1/
      end
c**************************************************************************
      block data followdata
      integer width
      character str1*30
      common/lablnc/width,str1
      data width/0/
      end
C********************************************************************
      block data tn2shidata
      include 'world3.h'
      data scbx3,scby3,scbz3/0.25,0.25,0.20/
      data xcbc2,ycbc2/0.5,0.40/
      data wx3min,wx3max,wy3min,wy3max,wz3min,wz3max
     $     /-1.,1.,-1.,1.,-1.,1./
      data w3nx/0.25/w3ny/0.25/w3nz/0.2/
      data ihiding/0/i3trunc/0/
      data z3sign/1./z3signch/1./
      data ax3chars/'axis-1','axis-2','axis-3'/
      end

c***************************************************************************
      block data hdcdata
      include 'hidcom.h'
      data ytop(1)/0./
      end

