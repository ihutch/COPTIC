c Routine to take a discrete probability distribution p(x), at an
c array of length np, at ordinates x(np), calculate from it a 
c cumulative probability cum(np) by integrating by trapezoidal 
c summation. Or accept a pre-calculated such cumulative distribution.
c And, given a value r, between 0 and 1, interpolate the value of
c x such that cum(x)=r.
c On Entry
c   np is the array length
c   p  is the probability distribution array
c   x  is the ordinate array to be interpolated.
c   cum might be the cumulative probability array or not yet.
c   r  is the value of cum (0..1) to interpolate x to.
c   isw decides if cum is to be initialized (yes if isw=0).
c On Return
c   xofcum is the interpolated value of x (ordinate).
c   cum is the cumulatve probability array (0..1).
c   p   might have been renormalized so its integral is unity.
c   isw !=0 if successful.[Another call without resetting doesn't initialize]
c   im  is the index immediately below the interpolate.
c   y   is the fraction of the interval at which to interpolate.

      real function xofcum(np,p,x,cum,r,isw,im,y)
      integer np,isw,im
      real p(np),x(np),cum(np),r

      if(.not.r.ge.0..or.r.gt.1)then
         write(*,*)'Xofcum Improper interpolation value r=',r
         xofcum=98.
         return
      endif
      if(isw.eq.0)then ! Initialize cum
         cum(1)=0.
         do i=2,np
            cum(i)=cum(i-1)+0.5*(p(i-1)+p(i))*(x(i)-x(i-1))
         enddo
         do i=1,np     ! Ensure proper normalization.
            cum(i)=cum(i)/(0.999999*cum(np))
            p(i)=p(i)/cum(np)
         enddo
      endif
! Now cum is initialized (or ought to be).       
      il=1
      ir=np
      cl=0.
      cr=1.
!      do i=1,2**np+1   ! Bisect to find the active interval
      do i=1,20   ! Bisect to find the active interval
         im=(il+ir)/2
         cm=cum(im)
!         write(*,'(4i5,4f8.4)')i,il,im,ir,cl,cm,cr,r
         if(im.eq.il)goto 1   ! Converged
         if(cm.le.r)then
            il=im
            cl=cm
         else
            ir=im
            cr=cm
         endif
      enddo
      write(*,'(a,4i5,2f8.3)')'Xofcum ERROR. Interval not converged',il
     $     ,im,ir,i,cm,r
      xofcum=99.
      return
 1    continue    
! Now im, im+1 is the relevant interval.
      delta=x(im+1)-x(im)
      pim=p(im)
      pip=p(im+1)
      y=(sqrt(pim**2+(pip-pim)*2.*(r-cm)/delta)-pim)/(pip-pim)
      xofcum=x(im)+y*delta
      if(y.gt.1.or..not.y.ge.0)then
         write(*,*)'Xofcum ERROR. Exceeded fractional value',y
      else
         isw=1
      endif
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine testxofcum
      integer npmax,np
      parameter (npmax=100,ntest=20)
      real p(npmax),x(npmax),cum(npmax),r
      xmax=5.
      np=npmax
      do j=1,2
         if(j.eq.2)np=5
         write(*,*)'np=',np
         write(*,*)'                 r          analytic-x         cp-x'
         do i=1,np
            x(i)=(i-1.)*xmax/(np-1.)
            p(i)=exp(-x(i))
         enddo
         if(j.eq.1)then
            call autoplot(x,p,np)
            call axlabels('x','p(x), !AJ!@pdx')
         else
            call color(2)
            call polyline(x,p,np)
         endif
         call color(j)
         isw=0
         round=.999999
         do i=1,ntest
            r=(1.-exp(-(i-1.)*xmax/(ntest-1.)))/(1.-exp(-xmax))*round
            cp=xofcum(np,p,x,cum,r,isw,im,y)
            write(*,*)i,r,-alog(1.-r*(1.-exp(-xmax))*round),cp
            call polymark(cp,r,1,j)
         enddo
         call polyline(x,cum,np)
      enddo
      call pltend
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Optional main. Normally comment out.
c      call testxofcum
c      end
