      character*4 pu,pd,spc,endpair
      character*15 only
      character*20 postlude
      character*300 prelude
      integer nr,nu,nd,ns,ne,no,nl
      common /pltfc/pu,pd,spc,endpair,postlude,prelude,only
      common /pltfi/nr,nu,nd,ns,ne,no,nl
c
c  Bounding Box in normal units.
      integer bxllx,bxlly,bxurx,bxury 
      common/bxcorners/bxllx,bxlly,bxurx,bxury 
