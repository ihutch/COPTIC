c Main program to get and plot phasespace data from file[s] on command line.
c If switch -t is given do triangular contouring.
c Trouble is, it gives too many colors for ps2gif, and pnmquant has
c problems with the colorscale legend. So use makepnganim
      include 'phasecom.f'
      character*100 phasefilename
      real x(npsbuf),u(npsbuf)
      real phirange,phirangeinit
      parameter (phirangeinit=0.5)
      character*10 string

      phirange=phirangeinit
      n=npsbuf
      do i=1,iargc()
         call getarg(i,phasefilename)
         write(*,*)phasefilename(1:50)
         if(phasefilename(1:2).eq.'-t'.and.ipsftri.eq.0)then
            call psftri
         else
            call phaseread(phasefilename,n,x,u,t)
            if(n.eq.0)goto 1
            call minmax(u,n,umin,umax)
            if(max(umax,abs(umin)).gt.phirange)phirange=phirange
     $           +phirangeinit
            write(string,'(f10.2)')t
            call pfset(-3)      ! Just output the files; don't plot.
            call multiframe(2,1,0)
            call pltinit(x(1),x(n),-phirange,phirange)
            call axis()
            call axlabels(' ','  !Af!@')
            call jdrwstr(wx2nx(x(n)),wy2ny(.9*phirange),string,-1.)
            call polyline(x,u,n)
            call phaseplot
            call pltend
         endif
 1       continue
      enddo
      if(i.lt.2)then
         write(*,*)'Usage phasereadplot file1 file2 ....'
      endif
      end
 
