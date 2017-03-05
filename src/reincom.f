!********************************************************************
! Random interpolate data.
      integer nvel,nQth
      parameter (nvel=50)
      parameter (nQth=200)
      real Qcom(nQth) 
      real Gcom(nvel,nQth)
      real Vcom(nvel)
      real pu1(nvel),pu2(nvel)
      real Pc(nQth,nvel)
! New BC
      integer bcphi,bcr
      logical infdbl
! Reinjection flux as a function of cos(theta) (line) and chi (column,
! from 0 to 9)
      common /ranQGcom/Gcom,Vcom,Qcom,pu1,pu2,Pc,infdbl,bcphi,bcr
!********************************************************************
! Left-over requirements for orbitinject.
      real averein,adeficit
! Put them into a common.
      common /reinextra/averein,adeficit
! Parameters used by orbitinject.
      integer myid
      parameter (myid=0)
      integer nthsize
      parameter (nthsize=201)

! Reindiag parameters.
      integer ndth,ndpsi
      parameter (ndth=100,ndpsi=100)
      real reinpos(ndth,ndpsi),cthtot(ndth),psitot(ndpsi)
      real reincth(ndth),reinpsi(ndpsi)
      real fv(ndth,3),vfv(ndth),sv(ndth),vs(ndth)
      common /diagrein/reinpos,cthtot,psitot,reincth,reinpsi,
     $     fv,vfv,sv,vs
