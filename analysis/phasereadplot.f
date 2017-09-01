c Main program to get and plot phasespace data from file[s] on command line.
c If switch -t is given do triangular contouring.
c Trouble is, it gives too many colors for ps2gif, and pnmquant has
c problems with the colorscale legend. So use makepnganim
      implicit none
      include 'phasecom.f'
      include '../accis/plotcom.h'
      character*100 phasefilename
      real x(npsbuf),u(npsbuf),uave(npsbuf)
      real phirange,phirangeinit
      parameter (phirangeinit=0.5)
      character*10 string
      integer i,ii,n,Np,Nave,Nastep
      real t,umin,umax,wx2nx,wy2ny

      Nave=1
      Nastep=1
      phirange=phirangeinit
      do i=1,iargc()
         call getarg(i,phasefilename)
         if(phasefilename(1:2).eq.'-N')then
c Set the starting number of filewriting to be N
            read(phasefilename(3:),*,err=2,end=2)Np
            pfilno=Np
            write(*,*)'Plot file number offset:',Np
            goto 1
 2          write(*,*)'Garbled -N flag (needs number):',phasefilename
            goto 1
         endif
         write(*,*)phasefilename(1:50)
         if(phasefilename(1:2).eq.'-A')then
            read(phasefilename(3:),*,err=3,end=3)Nave
            write(*,*)'Averaging potential over steps:',Nave
            goto 1
 3          write(*,*)'Garbled -a flag (needs number):',phasefilename
            goto 1
         endif
         if(phasefilename(1:2).eq.'-t'.and.ipsftri.eq.0)then
            call psftri
         else
            n=npsbuf
            call phaseread(phasefilename,n,x,u,t)
            if(n.eq.0)goto 1
! Averaging process over up to Nave steps. (A decaying exponential)
            if(Nave.gt.1)then
               do ii=1,n
                  uave(ii)=(uave(ii)*(Nastep-1.)+u(ii))/float(Nastep)
               enddo
               if(Nastep.lt.Nave)Nastep=Nastep+1
            endif
            call minmax(uave,n,umin,umax)
            if(max(umax,abs(umin)).gt.phirange*1.2)phirange=phirange
     $           +phirangeinit
            write(string,'(f10.2)')
            call pfset(-3)      ! Just output the files; don't plot.
            call multiframe(2,1,0)
            call pltinit(x(1),x(n),-phirange,phirange)
            call axis()
            call axlabels(' ','  !Af!@')
            call jdrwstr(wx2nx(x(n)),wy2ny(.9*phirange),string,-1.)
            if(Nave.gt.1)then
               call dashset(1)
               call polyline(x,u,n)
               call dashset(0)
               call polyline(x,uave,n)
            else
               call polyline(x,u,n)
            endif
            call phaseplot
            call pltend
         endif
 1       continue
      enddo
      if(i.lt.2)then
         write(*,*)'Usage phasereadplot [Options] file1 [file2 ....]'
         write(*,*)'Options: -Annn average-number',' -N starting-number'
      endif
      end
 
