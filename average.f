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
c****************************************************************
c This is exp(X^2)*erfc(X)
      FUNCTION expERFCC(X)
      Z=ABS(X)      
      T=1./(1.+0.5*Z)
      expERFCC=T*EXP(-1.26551223+T*(1.00002368+T*(.37409196+
     *    T*(.09678418+T*(-.18628806+T*(.27886807+T*(-1.13520398+
     *    T*(1.48851587+T*(-.82215223+T*.17087277)))))))))
      IF (X.LT.0.) expERFCC=2.*exp(z**2)-expERFCC
      END
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
c*********************************************************************
      subroutine average3d(q,qave,ifull,iuds,istepave)
c 3d only version of averaging. Obsolete.
      integer ifull(3),iuds(3)
      real q(ifull(1),ifull(2),ifull(3))
      real qave(ifull(1),ifull(2),ifull(3))

      do k=1,iuds(3)
         do j=1,iuds(2)
            do i=1,iuds(1)
               qave(i,j,k)=((istepave-1.)*qave(i,j,k)+q(i,j,k))
     $              /float(istepave)
            enddo
         enddo
      enddo

      end
c*********************************************************************
c General dimensional version of running average. See mditerate.f
      subroutine averagegd(q,qave,ifull,iuds,istepave)
      include 'ndimsdecl.f'
      ipin=0
      call mditerave(q,ndims,ifull,iuds,ipin,qave,istepave)
      end
