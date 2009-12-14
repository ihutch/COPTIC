      program integratetest
      parameter (imax=200)
      real h(0:imax)
      external finttest
      xw=0.
      xc=0.
      call cumprob(finttest,xw,xc,imax,h,ginfty)
      call yautoplot(h(0),imax+1)
      call pltend()
      end
c**********************************************************************
      real function finttest(x)

c Inverse square
c      f=1./(abs(x)**2+1.)

c Gaussian
c      f=exp(-(x-1.)**2)

c Exponential
      f=exp(-abs(x-1.))

c Half exponential
c      if(x.gt.0)then
c          f=exp(-x)
c      else
c         f=0.
c      endif

c Triangle
c      f=1.-abs(x+2.)
c      if(f.lt.0.)f=0.

c Box
c      if(abs(x-2.).gt.1.)then
c         f=0.
c      else
c         f=1.
c      endif
      finttest=f

      end
c**********************************************************************
