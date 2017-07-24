c Main program to get and plot phasespace data from file[s] on command line.
      include 'phasecom.f'
      character*100 phasefilename
      real x(npsbuf),u(npsbuf)
      real phirange
      parameter (phirange=0.5)
      character*10 string

      n=npsbuf
      do i=1,iargc()
         call getarg(i,phasefilename)
         call phaseread(phasefilename,n,x,u,t)
         write(string,'(f10.2)')t
         call pfset(-3) ! Just output the files; don't plot.
         call multiframe(2,1,0)
         call pltinit(x(1),x(n),-phirange,phirange)
         call axis()
         call axlabels(' ','  !Af!@')
         call jdrwstr(wx2nx(x(n)),wy2ny(.9*phirange),string,-1.)
         call polyline(x,u,n)
         call phaseplot
         call pltend
      enddo

      end
