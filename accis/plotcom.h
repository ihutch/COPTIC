c Worksurface common:
      integer nrows,ncolumns,nframe,multype
      real naxmin,naxmax,naymin,naymax,naxpt,naypt
      common/wksrfc/naxmin,naxmax,naymin,naymax,naxpt,naypt,
     $   nrows,ncolumns,nframe,multype

c Worlds and labels common:
      real xticlen,yticlen,xticoff,yticoff
      integer nxlabw,nxlabp,nylabw,nylabp,ticnum
      logical lxlog,lylog,lminor,lclog
      integer nxpow,nypow
      real xpow,ypow,wxmin,wxmax,wymin,wymax,w2nx,w2ny,n2sy
      logical ltlog
      common/worlds/
     $ xticlen,yticlen,xticoff,yticoff,
     $ nxlabw,nxlabp,nylabw,nylabp,ticnum,
     $ lxlog,lylog,nxpow,nypow,xpow,ypow,lclog,
     $ wxmin,wxmax,wymin,wymax,w2nx,w2ny,n2sy,ltlog,lminor
c
c  Cursor position
      real crsrx,crsry
      integer*2 updown
      common/crsrxy/crsrx,crsry,updown
c
c  Truncation monitoring
      real trcxmi,trcxma,trcymi,trcyma
      common/trncat/ trcxmi,trcxma,trcymi,trcyma
c
c  Characters
      integer BUFFER,NOCHARS
      parameter (BUFFER=110000,NOCHARS=2000)
      real chrscos,chrssin,chrswdth,chrshght,chrsslnt
      real pchrswdth,pchrshgt,pchrscos,pchrssin
      integer chrsaddr
      character*1 chrsfont
      common/chrcter/chrscos,chrssin,chrswdth,chrshght,chrsslnt,
     & pchrswdth,pchrshgt,pchrscos,pchrssin
      common/ac_chrfnt/chrsaddr(NOCHARS),chrsfont(BUFFER)
c
c  Plot-to-file control
      integer pfsw
      integer pfilno,pfnextsw,pfPS,psini
      common/pltfil/pfsw,pfilno,pfnextsw,pfPS,psini
c
c  Screen Parameters
      integer scrxpix,scrypix,ncolor,vmode
      real yoverx
      common/screenp/scrxpix,scrypix,ncolor,vmode,yoverx
