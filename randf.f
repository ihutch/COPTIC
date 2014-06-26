c***********************************************************************
c Restartable gasdev based on NR.
      FUNCTION gasdev(ireset)
      integer ireset
      real vr(2),v1,v2
      equivalence (v1,vr(1)),(v2,vr(2))
      include 'rancom.f'
      if(ireset.lt.0) gd_iset=0
      if(gd_iset.ne.1)then
 1       continue
         call ranlux(vr,2)
         v1=2.*v1-1.
         v2=2.*v2-1.
         r=v1**2+v2**2
        if(r.ge.1..or.r.eq.0.)go to 1
        fac=sqrt(-2.*log(r)/r)
        gd_gset=v1*fac
        gasdev=v2*fac
        gd_iset=1
      else
        gasdev=gd_gset
        gd_iset=0
      endif
      end
c**********************************************************************
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

