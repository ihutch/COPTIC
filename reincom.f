c********************************************************************
c Left-over requirements for orbitinject.
      real averein,adeficit
c Put them into a common.
      common /reinextra/averein,adeficit
c Parameters used by orbitinject.
      integer myid
      parameter (myid=0)
      integer nthsize
      parameter (nthsize=201)

c Reindiag parameters.
      integer ndth,ndpsi
      parameter (ndth=100,ndpsi=100)
      real reinpos(ndth,ndpsi),cthtot(ndth),psitot(ndpsi)
      real reincth(ndth),reinpsi(ndpsi)
      real fv(ndth,3),vfv(ndth),sv(ndth),vs(ndth)
      common /diagrein/reinpos,cthtot,psitot,reincth,reinpsi,
     $     fv,vfv,sv,vs
