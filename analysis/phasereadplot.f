c Main program to get and plot phasespace data from file[s] on command line.
c If switch -t is given do triangular contouring.
c Trouble is, it gives too many colors for ps2gif, and pnmquant has
c problems with the colorscale legend. So use makepnganim
      implicit none
      include '../src/ndimsdecl.f'
      include '../src/phasecom.f'
      include '../src/partcom.f'
      include '../accis/plotcom.h'
      character*100 phasefilename
      real x(npsbuf),u(npsbuf),uave(npsbuf)
      real phirange,phirangeinit
      parameter (phirangeinit=0.5)
      character*10 string
      integer i,ii,n,Np,Nave,Nastep,thespecies,idone,irun
      integer ilab,ispecies
      real t,umin,umax,wx2nx,wy2ny
      character*12 nlabel(2)
      data nlabel/' !Bn!di!d!@',' !Bn!de!d!@'/
      data idone/0/irun/0/
      Nave=1
      Nastep=1
      phirange=phirangeinit
      call pfset(3)
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
         elseif(phasefilename(1:2).eq.'-t'.and.ipsftri.eq.0)then
            call psftri
         elseif(phasefilename(1:2).eq.'-r')then
            irun=1
         elseif(phasefilename(1:2).eq.'-q')then
            irun=1
            call pfset(-3)
         else
            if(idone.gt.0)call prtend(' ')
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
            call minmax(u,n,umin,umax)
            if(max(umax,abs(umin)).gt.phirange*1.2)phirange=phirange
     $           +phirangeinit
            write(string,'(f10.2)')t
            call multiframe(nspecies+1,1,1)
            call dcharsize(.018,.018)
            call pltinit(x(1),x(n),-phirange,phirange)
            call axis()
            call axlabels(' ','  !Af!@')
            call jdrwstr(wx2nx(x(n)),wy2ny(.9*phirange),string,-1.)
            if(Nave.gt.1)then
               call dashset(4)
               call polyline(x,u,n)
               call dashset(0)
               call polyline(x,uave,n)
            else
               call polyline(x,u,n)
            endif
!           Plot density in the same frame
            call scalewn(x(1),x(n),0.8,1.35,.false.,.false.)
            call axptset(1.,1.)
            call ticrev
            call axis
            call ticrev
            call axptset(0.,0.)
            call legendline(1.04,0.3,258,'!Bn!@')
            call color(15)
            do ispecies=1,nspecies
               call color(ispecies+4)
               call polyline(psx,psn(1,ispecies),npsx)
               if(nspecies.gt.1)then !Assume electrons are first species
                  ilab=1
                  if(ispecies.eq.1)ilab=2
                  call legendline(0.8,0.35-0.08*ispecies,0
     $                 ,nlabel(ilab))
               endif
            enddo
!            thespecies=1
            do thespecies=1,nspecies
               call color(7)
               call phaseplot(thespecies)
            enddo
            idone=idone+1
            call accisflush
         endif
 1       continue
      enddo
      if(irun.eq.0)then 
         call pltend
      else 
         call prtend(' ')
      endif
      if(i.lt.2)then
         write(*,*)'Usage phasereadplot [Options] file1 [file2 ....]'
         write(*,*)'Options: -Annn average-number',' -N starting-number'
         write(*,*)'-r run continuously','  -q no screen plots.'
      endif
      end
 
