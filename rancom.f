c********************************************************************
c Random interpolate data.
      integer nvel,nQth
      parameter (nvel=50)
      parameter (nQth=200)
      real Qcom(nQth) 
      real Gcom(nvel,nQth)
      real Vcom(nvel)
      real pu1(nvel),pu2(nvel)
      real Pc(nQth,nvel)
c New BC
      integer bcphi,bcr
      logical infdbl
c Reinjection flux as a function of cos(theta) (line) and chi (column,
c from 0 to 9)
      common /rancom/Gcom,Vcom,Qcom,pu1,pu2,Pc,infdbl,bcphi,bcr
