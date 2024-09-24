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
      real phirange,phirangeinit,psnmin,psnmax,p1min,p1max,p2min,p2max
      real pbmax,pbmin
      parameter (phirangeinit=0.5)
      character*10 string
      integer i,ii,n,Np,Nave,Nastep,thespecies,idone,irun
      integer ilab,ispecies
      real t,umin,umax,wx2nx,wy2ny
      character*12 nlabel(2)
      data nlabel/' !Bn!di!d!@',' !Bn!de!d!@'/
      data idone/0/irun/0/
      data psnmin/0.9/psnmax/1.3/
      Nave=1
      Nastep=1
      phirange=phirangeinit
      lsideplot=.true. !default
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
         elseif(phasefilename(1:2).eq.'-h')then
            goto 4
         elseif(phasefilename(1:2).eq.'-q')then
            irun=1
            call pfset(-3)
         elseif(phasefilename(1:2).eq.'-s')then
            lsideplot=.not.lsideplot
         elseif(phasefilename(1:2).eq.'-l')then
            read(phasefilename(3:),*,err=5,end=5)ilogspec
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
 10         if(max(umax,abs(umin)).gt.phirange*1.2)then
               phirange=phirange+0.5*phirangeinit
               goto 10
            endif
 11         if(max(umax,abs(umin),phirangeinit).lt.phirange*0.6)then
               phirange=phirange-phirangeinit
               goto 11
            endif
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
            call minmax(psn(1,1),npsx,p1min,p1max)
            call minmax(psn(1,2),npsx,p2min,p2max)
            pbmin=min(p1min,p2min)
            pbmax=max(p1max,p2max)
 12         if(pbmax-pbmin.gt.(psnmax-pbmin)*1.3)then
               psnmax=(psnmax-pbmin)*1.25+pbmin
               goto 12
            endif
 13         if(pbmax-pbmin.lt.(psnmax-pbmin)*0.7.and.psnmax.gt.1.3)then
               psnmax=(psnmax-pbmin)*.75+pbmin
               goto 13
            endif
 14         if(pbmin.lt.psnmin)then
               psnmin=psnmin*.1
               goto 14
            endif
 15         if(pbmin.gt.psnmin*1.4)then
               psnmin=min(psnmin*1.1,.9)
               goto 15
            endif
 16         if(pbmax-pbmin.lt.(psnmax-psnmin)*0.2)then
               write(*,'(4f8.4)')psnmax,psnmin,pbmax,pbmin
               psnmin=pbmin-.3*(pbmax-pbmin)
               psnmax=pbmax+1.8*(pbmax-pbmin)
               write(*,'(4f8.4)')psnmax,psnmin
               goto 16
            endif
            call scalewn(x(1),x(n),psnmin,psnmax,.false.,.false.)
            call axptset(1.,1.)
            call ticrev
            call axis
            call ticrev
            call axptset(0.,0.)
            if(lsideplot)then
               call legendline(1.04,0.3,258,'    !Bn!@')
            else
               call legendline(1.04,0.3,258,'!Bn!@')
            endif
            call color(15)
            do ispecies=1,nspecies
               call color(ispecies+4)
               call polyline(psx,psn(1,ispecies),npsx)
               if(nspecies.gt.1)then !Assume electrons are first species
                  ilab=1
                  if(ispecies.eq.1)ilab=2
                  call legendline(0.05,0.95-0.08*ispecies,0
     $                 ,nlabel(ilab))
               endif
            enddo
            do thespecies=1,nspecies
               finfmax=maxval(finfofv(1:npsv,thespecies))
               call color(7)
               call phaseplot(thespecies)
               call color(15)
               call vecw(psx(1),0.,0)
               call vecw(psx(npsx),0.,1)
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
      if(i.ge.2)return
 5    write(*,*)'Could not read ilogspec'
 4    continue
      write(*,*)'Usage phasereadplot [Options] file1 [file2 ....]'
      write(*,*)'Options: -Annn average-number',' -N starting-number'
      write(*,*)'-r run continuously','  -q no screen plots',
     $     '  -s toggle sideways plot of f(v)'
      write(*,*)'-lnn set species for log phase contours (before files)'
      end
 
