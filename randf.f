c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___
c
c This code is copyright(c) (2003-5) Ian H Hutchinson hutch@psfc.mit.edu
c
c  It may be used freely with the stipulation that any scientific or
c scholarly publication concerning work that uses the code must give an
c acknowledgement referring to the papers I.H.Hutchinson, Plasma Physics
c and Controlled Fusion, vol 44, p 1953 (2002), vol 45, p 1477 (2003).
c  The code may not be redistributed except in its original package.
c
c No warranty, explicit or implied, is given. If you choose to build
c or run the code, you do so at your own risk.
c
c Version 2.6 Aug 2005.
c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___
c***********************************************************************
      FUNCTION GASDEV(IDUM)
      save
      DATA ISET/0/
      IF (ISET.EQ.0) THEN
 1       continue
         V1=2.*RANd()-1.
         V2=2.*RANd()-1.
c         V1=2.*RAN0(idum)-1.
c         V2=2.*RAN0(idum)-1.
        R=V1**2+V2**2
        IF(R.GE.1..OR.R.EQ.0.)GO TO 1
        FAC=SQRT(-2.*LOG(R)/R)
        GSET=V1*FAC
        GASDEV=V2*FAC
        ISET=1
      ELSE
        GASDEV=GSET
        ISET=0
      ENDIF
      RETURN
      END
c**********************************************************************
      FUNCTION RAN0(IDUM)
      save
c Version of July 06 that removes the argument dependence.
      DIMENSION V(97)
      DATA IFF /0/
      IF(IFF.EQ.0)THEN
        IFF=1
        DO 11 J=1,97
          DUM=RANd()
11      CONTINUE
        DO 12 J=1,97
          V(J)=RANd()
12      CONTINUE
        Y=RANd()
      ENDIF
c IHH hack to prevent errors when Y=1. Was 97.
      J=1+INT(96.9999*Y)
      IF(J.GT.97.OR.J.LT.1)then
         write(*,*)'RAN0 error: j=',j,'  y=',y
c         PAUSE 'RAN0 problem'
         j=1
      endif
      Y=V(J)
      RAN0=Y
      V(J)=RANd()
      RETURN
      END
c**********************************************************************
      FUNCTION RAN1(IDUM)
      save
      DIMENSION R(97)
      PARAMETER (M1=259200,IA1=7141,IC1=54773,RM1=3.8580247E-6)
      PARAMETER (M2=134456,IA2=8121,IC2=28411,RM2=7.4373773E-6)
      PARAMETER (M3=243000,IA3=4561,IC3=51349)
      DATA IFF /0/
      IF (IDUM.LT.0.OR.IFF.EQ.0) THEN
        IFF=1
        IX1=MOD(IC1-IDUM,M1)
        IX1=MOD(IA1*IX1+IC1,M1)
        IX2=MOD(IX1,M2)
        IX1=MOD(IA1*IX1+IC1,M1)
        IX3=MOD(IX1,M3)
        DO 11 J=1,97
          IX1=MOD(IA1*IX1+IC1,M1)
          IX2=MOD(IA2*IX2+IC2,M2)
          R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
11      CONTINUE
        IDUM=1
      ENDIF
      IX1=MOD(IA1*IX1+IC1,M1)
      IX2=MOD(IA2*IX2+IC2,M2)
      IX3=MOD(IA3*IX3+IC3,M3)
      J=1+(97*IX3)/M3
c      IF(J.GT.97.OR.J.LT.1)PAUSE
      RAN1=R(J)
      R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
      RETURN
      END
c

c********************************************************************
c Given a monotonic function Q(x) on a 1-D grid x=1..nq, 
c   solve Q(x)=y for x.
c That is, invert Q to give x=Q^-1(y).
      subroutine invtfunc(Q,nq,y,x)
      implicit none
      integer nq
      real Q(nq)
      real y,x
c
      integer iqr,iql,iqx
      real Qx,Qr,Ql
      Ql=Q(1)
      Qr=Q(nq)
      iql=1
      iqr=nq
      if((y-Ql)*(y-Qr).gt.0.) then
c Value is outside the range.
         x=0
         return
      endif
 200  if(iqr-iql.eq.1)goto 210
      iqx=(iqr+iql)/2
      Qx=Q(iqx)
c      write(*,*)y,Ql,Qx,Qr,iql,iqr
c Formerly .lt. which is an error.
      if((Qx-y)*(Qr-y).le.0.) then
         Ql=Qx
         iql=iqx
      else
         Qr=Qx
         iqr=iqx
      endif
      goto 200
 210  continue
c Now iql and iqr, Ql and Qr bracket Q
      x=(y-Ql)/(Qr-Ql)+iql
c      if(x.gt.float(nq))then
c         write(*,*)'invtfn error',x,nq,y,q(1),q(nq),
c     $     ql,qr,iql,iqr,qx,iqx
c         write(*,*)'q=',q
c      endif
      end
c********************************************************************
C complementary error function from NR.
      FUNCTION ERFCC(X)
      Z=ABS(X)      
      T=1./(1.+0.5*Z)
      ERFCC=T*EXP(-Z*Z-1.26551223+T*(1.00002368+T*(.37409196+
     *    T*(.09678418+T*(-.18628806+T*(.27886807+T*(-1.13520398+
     *    T*(1.48851587+T*(-.82215223+T*.17087277)))))))))
      IF (X.LT.0.) ERFCC=2.-ERFCC
      RETURN
      END

