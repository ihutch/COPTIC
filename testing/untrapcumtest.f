!*********************************************************************
      program untrapcumtest
      integer nphi,nin,ntest
      parameter (nphi=501,nin=400,nt=32,ntot=nin+nt)
      integer np(0:nphi),nm(0:nphi)
!      real fpa(nin),fma(nin),ua(nin)
      real cump(-ntot:ntot,0:nphi),cumv(-ntot:ntot,0:nphi)
      idebug=2
      psi=1.
      um=0.
      holelen=4.
      tl=-1.
!      ntest=490
      ntest=20
      write(*,*)'ntest=',ntest
      call cumdis(nphi,um,psi,holelen,tl,nin,nt,np,nm,cump,cumv)
!      if(idebug.gt.1)then
!         call autoplot(ua,fma,nm)
!         call axlabels('u','fma')
!         call polyline(ua,fpa,np)
!         call pltend()
!      endif
      maxp=np(ntest)/8+nt
      maxm=nm(ntest)/8+nt
!      maxp=np(ntest)+nt
!      maxm=nm(ntest)+nt
!      write(*,*) cump(-nt:nt,ntest)!,cump(maxp,ntest)
      call autoplot(cumv(-maxm,ntest),cump(-maxm,ntest),maxp+maxm+1)
      call axis2()
      call axlabels('u','cumulative distribution')
      call polymark(cumv(-nt-1,ntest),cump(-nt-1,ntest),2*nt+3,1)
      call pltend()
      end
!*********************************************************************
