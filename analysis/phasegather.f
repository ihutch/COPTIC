      program phasegather
! Analysis. Not to be confused with routines in src/phasespace.f.
! Read particle data and construct phase space plots.
      include 'examdecl.f' 
! (Examdecl itself includes meshcom.f plascom.f, objcom.f)
      parameter (nfilemax=999)
      include '../src/partcom.f'
      include '../src/ptaccom.f'
      include '../src/myidcom.f'
      include '../src/phasecom.f'
 
      character*10 chartemp
      character*100 name
!      character*100 string
      logical ldoc
!      real extra(nptdiagmax,ndimsmax)
! Spatial limits bottom-top, dimensions
      real xlimit(2,ndimsmax),xnewlim(2,ndimsmax)
      real elimit(2)
      real Bdirs(ndimsmax+1)
! Velocity limits
      real vlimit(2,ndimsmax)
      
      integer nx,nv
      parameter (nx=70,nv=70,ne=400)
      real hist(nx,nv),cworka(nx,ne),zclv(nx)
      real ehist(ne),ev(ne),ehistscaled(ne)
      real eh(0:ne)

      integer iuin(ndimsmax)
      integer ivtk
      integer ndfirst,ndlast
      integer iepeak(1)
      integer iaxis
!      character*200 ivarnames(2*ndimsmax)

! Maximum number of files to read for specification without extension.
      nfmax=nfilemax

! Defaults
      nprocs=1
      ndfirst=1
      ndlast=3
      iaxis=0
      do id=1,ndimsmax
         xlimit(1,id)=-500.
         xlimit(2,id)=500.
         xnewlim(1,id)=0.
         xnewlim(2,id)=0.
         vlimit(1,id)=-3.5
         vlimit(2,id)=3.5
      enddo

      elimit(2)=1.5
      elimit(1)=-1.
! Make box join at e=0. Decide how many boxes e<0
      ne0=int(ne*(-elimit(1))/(elimit(2)-elimit(1)))
      de=-elimit(1)/ne0
      elimit(2)=ne*de+elimit(1)

      nfvaccum=0
      nhist=0
      nehist=0
      iworking=0
! Iterate over possibly multiple file definitions.
      do ifile=1,nfilemax
         call partexamargs(xlimit,vlimit,iuin,cellvol,Bdirs,ldoc,ivtk
     $        ,ispecies,ndfirst,ndlast,iworking)
         if(iworking.eq.-1)goto 11 ! Exhausted arguments/files
         if(ifile.eq.1)then
! On first return, base filename is in partfilename.
! Do the file naming calculations.
            ip=lentrim(partfilename)-3
            if(partfilename(ip:ip).eq.'.')then
! If filename has a 3-character extension, we are reading one file at a time. 
               nfmax=0
               name=partfilename
               if(partfilename(ip:ip+3).eq.'.pex')then
                  write(*,*)'Using stored distribution file'
                  nfmax=-1
               endif
               phifilename=partfilename(1:ip)//'phi'
            elseif(partfilename(ip-4:ip-1).eq.'.pex')then
               nfmax=-1
               name=partfilename
               write(*,*)'Reading numbered pex file '
     $              ,name(1:lentrim(name))
               phifilename=partfilename(1:ip-4)//'phi'
            else ! Extensionless filename
               phifilename=partfilename(1:lentrim(partfilename))//'.phi'
            endif            
! Read in potential
            ied=1
            call array3read(phifilename,ifull,iuds,ied,u,ierr)
            if(ierr.eq.0.and.abs(phip).gt.1.)elimit(1)=phip
            
! Don't allow xlimits beyond mesh ends less one cell
            do id=1,ndimsmax
! Set meshstart/end for periodic particles
               xmeshstart(id)=(xn(ixnp(id)+1)+xn(ixnp(id)+2))/2
               xmeshend(id)=(xn(ixnp(id+1))+xn(ixnp(id+1)-1))/2
               xlimit(1,id)=max(xlimit(1,id),xmeshstart(id))
               xlimit(2,id)=min(xlimit(2,id),xmeshend(id))
!         write(*,*)'xlimit',xlimit(1,id),xlimit(2,id)
            enddo
         endif     ! End of first ifile initialization

! Multiple particles files, only if we haven't read one with extension
         do i=0,nfmax
            if(nfmax.ne.0)then
               write(chartemp,'(''.'',i3.3)')i
               name=partfilename(1:lentrim(partfilename))//chartemp
            else
               name=partfilename
            endif
            write(*,*)'Reading file ',name(1:lentrim(name))
            Bt=0.
            call partread(name,ierr)
            if(ispecies.gt.nspecies)then
               ispecies=nspecies
               write(*,*)'nspecies=',nspecies,'  Reset ispecies'
     $              ,ispecies
            endif
            if(ldoc)then
               write(*,*)'debyelen,Ti,vd,rs,phip=',debyelen,Ti,vd,rs
     $              ,phip
               write(*,*)'iregion_part,n_part,dt,ldiags,rhoinf,nrein,',
     $              'phirein,numprocs='
               write(*,*)iregion_part,n_part,dt,ldiags,rhoinf,nrein
     $              ,phirein,numprocs
               write(*,*)'eoverm,Bt,Bfield,vpar,vperp=',eoverm,Bt,Bfield
     $              ,vpar,vperp
               write(*,*)'nspecies',nspecies
            endif
            if(ierr-4*(ierr/4).ne.0)goto 11
            if(Bdirs(4).gt.0. .or. Bt.eq.0)then
! All directions were set by commandline. Or none were read from file.
               do k=1,ndims
                  Bfield(k)=Bdirs(k)
               enddo
            endif
            if(cellvol.eq.-1)write(*,*)'Bfield (projection)',Bfield
! Now act on the particle data we have read in:
            if(iaxis.eq.0)then
! Choose active axis based on b-settings
               do k=1,ndims
                  if(Bfield(k).ne.0)iaxis=k
               enddo
! When first determined, shave ends off iaxis direction, if needed:
               fr=1./nx
               if(xlimit(1,iaxis).eq. (xn(ixnp(iaxis)+1)+xn(ixnp(iaxis)
     $              +2))/2) xlimit(1,iaxis)=(1-fr)*xlimit(1,iaxis)+fr
     $              *xlimit(2,iaxis)
               if(xlimit(2,iaxis).eq. (xn(ixnp(iaxis+1))+xn(ixnp(iaxis
     $              +1)-1))/2) xlimit(2,iaxis)=fr*xlimit(1,iaxis)+(1-fr)
     $              *xlimit(2,iaxis)
            endif
! As we go, plot dots in phase-space
            call plotlocalphasepoints(xlimit,vlimit,iaxis,ispecies
     $           ,nfvaccum)
! Increment the phase-space histogram
            call fhistinc(xlimit,vlimit,iaxis,ispecies,nhist,hist,nx
     $           ,nv)
! Increment the total energy histogram. Use all axes if iexis=0
            iexis=1
            call ehistinc(xlimit,elimit,ispecies,nehist,ehist,ne,u,ifull
     $           ,iused,iexis,potentialmax)
         enddo                  ! of multiple file reading
      enddo  ! of argument reading
 11   continue
      call pltend()

      xl=xlimit(2,iaxis)-xlimit(1,iaxis)
      wp=8.   ! Estimate of hole length to be subtracted from box length
      xl=(xl-wp)    ! Untrapped length
! Find the maximum of the energy occurrence histogram epeak.
      iepeak=maxloc(ehist)
      i=iepeak(1)
! The histogram is assigned as ie=int(1+(energy-elimit(1))/de)
! Which makes box 1 contain 0 to de, etc: on average 0.5de.
! Therefore this is the center energy of the peak box:
!      epeak=(i-0.5)*(elimit(2)-elimit(1))/ne +elimit(1)
! Centroid of 3 peak boxes (not weighted appropriately)
!      epeak=elimit(1)+
!     $     (ehist(i-1)*(i-1.5)+ehist(i)*(i-0.5)+ehist(i+1)*(i+0.5))
!     $     /(ehist(i-1)+ehist(i)+ehist(i+1))*(elimit(2)-elimit(1))/ne
! Centroid of 3 peak boxes weighted by trap/untrap ratio.
      tup=sqrt(2.*3.1415926/xl) ! sqrt trapped to untrapped weight ratio.
      epeak=elimit(1)+
     $     (ehist(i-1)*(i-1.5)/tup+ehist(i)*(i-0.5)
     $              +ehist(i+1)*(i+0.5)*tup)
     $     /(ehist(i-1)/tup+ehist(i)+ehist(i+1)*tup)
     $     *(elimit(2)-elimit(1))/ne
!      epeak=epeak-.002 ! hack adjustment
      epeak=0.         ! no adjustment
!      write(*,'(a,i8,3f8.4)')'nehist,emin,emax,epeak='
!     $     ,nehist,elimit(1),elimit(2),epeak
!      write(*,'(a,10f8.4)')'xl,wp=',xl,wp
      eh(0)=elimit(1)
      twoqum1=sqrt(sqrt(2.))-1
      G1=twoqum1*sqrt(-eh(0)/potentialmax)
      sf1=4.*sqrt(2.)*3.141593*sqrt(-eh(0))*(1.+4.*G1/3.+G1**2/2.)
      do i=1,ne
         ev(i)=(i-0.5)*(elimit(2)-elimit(1))/ne +elimit(1)-epeak
         eh(i)=(i-0.)*(elimit(2)-elimit(1))/ne +elimit(1)-epeak
! The scaling to give f(W) needs to divide the occurrences by the 
! phasespace attributable to the dW.
         if(eh(i).le.0.)then  ! Trapped particle phasespace.
            G1=twoqum1*sqrt(-eh(i)/potentialmax) 
!            G1=0  ! To turn off the T_b correction for finite W
            sf2=4.*sqrt(2.)*3.141593*sqrt(-eh(i))*(1.+4.*G1/3.+G1**2/2.)
            ehistscaled(i)=ehist(i)/(sf1-sf2)
            sf1=sf2
         elseif(eh(i-1).ge.0)then ! Untrapped phasespace.
            ehistscaled(i)=ehist(i)*sqrt(2.)
     $           /(4.*xl*(sqrt(eh(i))-sqrt(eh(i-1))))
         else ! Mixed trap/untrap
            ehistscaled(i)=ehist(i)/(sf1+2.*sqrt(2.)*xl*sqrt(eh(i)))
         endif
      enddo
      if(.false.)then
! Plot the energy histogram.
      call pfset(3)
      call autoplot(ev,ehist,ne)
      call axis()
!      call axis2()
      call polybox(eh,ehist,ne)
      call axlabels('Total Energy: W','Occurences: Distribution')
      call color(3)
      call axptset(0.,1.)
      call ticrev
      call altxaxis(.05,.05)
      call axptset(0.,0.)
      call ticrev
      call legendline(.45,1.02,258,'W')
      call winset(.true.)
      call polybox(20.*eh,ehist,ne)
      call pltend()

! The scaled version that's supposed to give f(v)
      call autoplot(ev,ehistscaled,ne)
      call axis()
!      call axis2()
      call polybox(eh,ehistscaled,ne)
      call axlabels('W','f(v)=dN/d|W|!u1/2!u')
      call color(3)
      call axptset(0.,1.)
      call ticrev
      call altxaxis(.05,.05)
      call axptset(0.,0.)
      call ticrev
      call legendline(.45,1.02,258,'Energy')
      call winset(.true.)
      call polybox(20.*eh,ehistscaled,ne)
      call pltend()
      call pfset(0)
      else
         write(*,*)'Not plotting energy information'
      endif

! Color Contour of Phase space Histogram.
      call minmax2(hist,nx,nx,nv,hmin,hmax)
      ncont=10
      do i=1,ncont
         zclv(i)=hmin+(i-1)*(hmax-hmin)/(ncont-1)
      enddo
      icl=ncont
      icsw=0+16
      call pltinit(1.,float(nx),1.,float(nv))
      call blueredgreenwhite
      call contourl(hist,cworka,nx,nx,nv,zclv,icl,dummy,dummy,icsw) 
      call gradlegend(zclv(1),zclv(icl),-.2,.1,-.2,1.,.03,.false.) 
      call scalewn(xlimit(1,iaxis),xlimit(2,iaxis),vlimit(1,iaxis)
     $     ,vlimit(2,iaxis),.false.,.false.)
      call axis
      call pltend()
! Second example to illustrate autocolcont.
      if(.false.)then
      call pltinit(1.,float(nx),1.,float(nv))
      call autocolcont(hist,nx,nx,nv)
      call scalewn(xlimit(1,iaxis),xlimit(2,iaxis),vlimit(1,iaxis)
     $     ,vlimit(2,iaxis),.false.,.false.)
      call axis
      call pltend()
      endif
! Plot using the coptic code, needs xmeshstart/end set correctly.
      nstep=0
      call phasepscont(ifull,iuds,u,nstep,.true.,'',q)
      call pltend()

      end
!*************************************************************
      subroutine plotlocalphasepoints(xlimit,vlimit,iaxis,ispecies
     $     ,nfvaccum)
! Plot points in phase space in direction iaxis for particles that lie
! in the space and velocity ranges given by xlimit, vlimit.
      include '../src/ndimsdecl.f'
      include '../src/partcom.f'
      real xlimit(2,ndims),vlimit(2,ndims)
      character*1 axnames(3)
      data axnames/'x','y','z'/ 

      iwarned=0
      id=999
! Initialize the plot
      if(nfvaccum.eq.0)then
!         call pfset(3)
         call pltinit(xlimit(1,iaxis),xlimit(2,iaxis),vlimit(1,iaxis)
     $     ,vlimit(2,iaxis))
         call axis()
         call axlabels(axnames(iaxis),'v!d'//axnames(iaxis)//'!d')
         call color(ibrickred())
      endif
!      write(*,*)'x,vlimits'
!      write(*,'(6f8.3)')xlimit,vlimit
      do j=iicparta(ispecies),iocparta(ispecies)
! Only for filled slots
         if(x_part(iflag,j).ne.0)then
            do id=1,ndims
               x=x_part(id,j)
               v=x_part(id+ndims,j)
               if(x.lt.xlimit(1,id).or.x.gt.xlimit(2,id))then
                  if(iwarned.eq.0.and.id.ne.iaxis)then
                     iwarned=1
            write(*,*)'**** WARNING: Only partial transverse domain:'
     $                    ,id,x_part(id,j)

                  endif
                  goto 12                  
               endif
               if(v.lt.vlimit(1,id).or.v.gt.vlimit(2,id))goto 12
            enddo
            nfvaccum=nfvaccum+1
! Here for a particle lying within the limits.
! Plot a point in phase space
            call vecw(x_part(iaxis,j),x_part(iaxis+ndims,j),-1)
         endif
 12      continue
      enddo
!      write(*,*)'Plotted',nfvaccum,' particles'
      end
!*************************************************************
!*************************************************************
      subroutine fhistinc(xlimit,vlimit,iaxis,ispecies
     $     ,nhist,hist,nx,nv)
! Increment histogram of phase space in direction iaxis for particles
! that lie in the space and velocity ranges given by xlimit, vlimit.
      include '../src/ndimsdecl.f'
      include '../src/partcom.f'
      real xlimit(2,ndims),vlimit(2,ndims)
      real hist(nx,nv)
!      character*1 axnames(3)
!      data axnames/'x','y','z'/ 
      save x1,x2,xdiff,v1,vdiff

      if(nhist.eq.0)then
! Initialize
         x1=xlimit(1,iaxis)
         x2=xlimit(2,iaxis)
         xdiff=x2*1.00001-x1
         v1=vlimit(1,iaxis)
         v2=vlimit(2,iaxis)
         vdiff=v2*1.00001-v1
!         write(*,*)'x1,x2,xdiff,v1,v2,vdiff',x1,x2,xdiff,v1,v2,vdiff
         do j=1,nv
            do i=1,nx
               hist(i,j)=0.
            enddo
         enddo
      endif
      do j=iicparta(ispecies),iocparta(ispecies)
! Only for filled slots
         if(x_part(iflag,j).ne.0)then
            do id=1,ndims
!               x=x_part(id,j)
!Leapfrog correct:
               x=x_part(id,j)-x_part(id+ndims,j)*0.5*x_part(idtp,j) 
               v=x_part(id+ndims,j)
               if(x.lt.xlimit(1,id).or.x.gt.xlimit(2,id))goto 12
               if(v.lt.vlimit(1,id).or.v.gt.vlimit(2,id))goto 12
            enddo
! Here for a particle lying within the limits.
!            x=x_part(iaxis,j)
            x=x_part(iaxis,j)-x_part(iaxis+ndims,j)*0.5*x_part(idtp,j)
            v=x_part(iaxis+ndims,j)
            ix=int(nx*(x-x1)/xdiff+1)
            iv=int(nv*(v-v1)/vdiff+1)
            nhist=nhist+1
            hist(ix,iv)=hist(ix,iv)+1.
         endif
 12      continue
      enddo
      write(*,*)'Accumulated',nhist,' particles'
      end
!*************************************************************
      subroutine ehistinc(xlimit,elimit,ispecies
     $     ,nhist,ehist,ne,u,ifull,iused,iexis,potentialmax)
! Increment the histogram of total energy distribution using current
! particle data. 
! iexis if non-zero restricts the kinetic energy to that axis. 
! But then equal intervals of energy leads to a singularity at
! zero. We could perhaps scale them differently.
      include '../src/ndimsdecl.f'
      include '../src/partcom.f'
      include '../src/plascom.f'
      include '../src/meshcom.f'
      real u(*)
      integer ifull(ndims),iexis
      real xlimit(2,ndims),elimit(2)
      real ehist(ne)
      integer iLs(ndims+1)
! Dummies
      integer iregion
      real cij
      data de/0./  ! silence warning
      iused=0

      iLs(1)=1
      do i =1,ndims
         iLs(i+1)=iLs(i)*ifull(i)
      enddo
      ntrap=0
      nuntrap=0

      if(nhist.eq.0)then
! Initialize
         do j=1,ne
            ehist(j)=0.
         enddo
         de=(elimit(2)-elimit(1))/ne
      endif
!      write(*,*)'ehist axis',iexis,' dt=',dt
      potentialmax=0.
      do j=iicparta(ispecies),iocparta(ispecies)
! Only for filled slots
         if(x_part(iflag,j).ne.0)then
            v2=0.
            do id=1,ndims
! Leapfrog corrected position for getpotential
               x=x_part(id,j)-x_part(id+ndims,j)*0.5*x_part(idtp,j)
               if(x.lt.xlimit(1,id).or.x.gt.xlimit(2,id))goto 12
! Here for a particle lying within the space limits so far
! Locate the mesh full fractional postion of x. From potential at pt.
               ix=interp(xn(ixnp(id)+1),ixnp(id+1)-ixnp(id),x,xm)
               if(xm.lt.1.or.xm.gt.ixnp(id+1)-ixnp(id))then
                  write(*,*)'Skipping mesh error',j,xm,x,xmeshstart(id)
                  goto 12
               endif
               x_part(2*ndims+id,j)=xm  ! For getpotential
! Sum the kinetic energy for all directions if iexis=0.
               if(iexis.eq.0.or.iexis.eq.id)v2=v2+x_part(id+ndims,j)**2
            enddo            
! This is a wanted particle. Add its energy to histogram.
            potential=getpotential(u,cij,iLs,x_part(2*ndims+1,j),iregion
     $           ,0)
            if(potential.gt.potentialmax)potentialmax=potential
            energy=eoverms(ispecies)*potential+v2/2.
            ie=int(1+(energy-elimit(1))/de)
            if(.not.(ie.gt.0.and.ie.le.ne))goto 12
            ehist(ie)=ehist(ie)+1.
            nhist=nhist+1
            if(energy.lt.0)then
               ntrap=ntrap+1
            else
               nuntrap=nuntrap+1
            endif
         endif
 12      continue
      enddo
!      write(*,*)'ntrap/untrap=',ntrap,nuntrap,float(ntrap)/nuntrap
      end

!*************************************************************
      subroutine partexamargs(xlimit,vlimit
     $           ,iuin,cellvol,Bdirs,ldoc,ivtk,ispecies,ndfirst,ndlast
     $           ,iworking)
      include 'examdecl.f'
      include '../src/ptaccom.f'
      real xlimit(2,3),vlimit(2,3),Bdirs(4)
      integer iuin(3)
      logical ldoc
      integer ivtk
! I think unused here 26 May 12. But I'm not sure.
      ifull(1)=na_i
      ifull(2)=na_j
      ifull(3)=na_k
      ivtk=0
      ispecies=1

      do i=1,3
         iuin(i)=9
! Default field/projection direction:
         Bdirs(i)=0
         if(i.eq.3)Bdirs(i)=1
      enddo
! Use cellvol=0. by default.
      cellvol=0.
      ldoc=.false.

! silence warnings:
      fluxfilename=' '
      zp(1,1,1)=0.
! Defaults
      partfilename=' '

! Deal with arguments
      if(iargc().eq.0) goto 201
      istart=iworking+1
      do i=istart,iargc()
         call getarg(i,argument)
         iworking=iworking+1
         if(argument(1:1).eq.'-')then  ! Non-file switches
            if(argument(1:2).eq.'-x')then
               read(argument(3:),*,err=201) xlimit(1,1),xlimit(2,1)
            elseif(argument(1:2).eq.'-y')then
               read(argument(3:),*,err=201) xlimit(1,2),xlimit(2,2)
            elseif(argument(1:2).eq.'-z')then
               read(argument(3:),*,err=201) xlimit(1,3),xlimit(2,3)
            elseif(argument(1:2).eq.'-u')then
               read(argument(3:),*,err=201) vlimit(1,1),vlimit(2,1)
            elseif(argument(1:4).eq.'-vtk')then
               ivtk=1
            elseif(argument(1:2).eq.'-v')then
               read(argument(3:),*,err=201) vlimit(1,2),vlimit(2,2)
            elseif(argument(1:2).eq.'-w')then
               read(argument(3:),*,err=201) vlimit(1,3),vlimit(2,3)
            elseif(argument(1:3).eq.'-pu')then
               nptdiag=nsbins
            elseif(argument(1:2).eq.'-p')then
               cellvol=-1
               bread=0.
               read(argument(3:),*,err=102)(Bdirs(j),j=1,3)
               bread=1.
 102           continue
               Bdirs(4)=0.
               do j=1,3
                  Bdirs(4)=Bdirs(4)+Bdirs(j)**2
               enddo
               if(Bdirs(4).gt.0)then
                  Bdirs(4)=sqrt(Bdirs(4))
                  do j=1,3
                     Bdirs(j)=Bdirs(j)/Bdirs(4)
                  enddo
               endif
               Bdirs(4)=Bdirs(4)*bread
            elseif(argument(1:2).eq.'-b')then
               read(argument(3:),*,err=201)iuin
            elseif(argument(1:2).eq.'-q')then
               ldoc=.true.
            elseif(argument(1:3).eq.'-sp')then
               ispecies=ispecies+1
            elseif(argument(1:2).eq.'-d')then
               read(argument(3:),*,err=201)ndfirst,ndlast
            endif
            if(argument(1:13).eq.'--objfilename')
     $        read(argument(14:),'(a)',err=201)objfilename
            if(argument(1:2).eq.'-f')
     $           read(argument(3:),'(a)',err=201)partfilename
            if(argument(1:2).eq.'-h')goto 203
            if(argument(1:2).eq.'-?')goto 203
         else
            read(argument(1:),'(a)',err=201)partfilename
!            write(*,*)partfilename
            return  ! Return when read a filename.
         endif
      enddo
      iworking=-1 ! Exhausted all arguments
      return
!------------------------------------------------------------
! Help text
 201  continue
      write(*,*)'=====Error reading command line argument'
 203  continue
 301  format(a,i5,a)
 302  format(a,4f8.3)
 303  format(a,5i5)
      write(*,301)'Usage: phasespace [switches] <partfile> (no ext)'
      write(*,302)' -x -y -z<fff,fff>  set position range. [',
     $     xlimit(1,1),xlimit(2,1)
      write(*,302)' -u -v -w<fff,fff>  set velocity range. ['
     $     ,vlimit(1,1),vlimit(2,1)
      write(*,303)' -b<nx,ny,nz>  set spatial block range. [',iuin
      write(*,302)' -p[bx,by,bz]  project [in direction]   [',Bdirs
      write(*,'(a,$)')' -d<i,i> Set first,last v-index'
      write(*,303)' to plot [',ndfirst,ndlast
      write(*,301)' -sp     Increment species number to examine'
      write(*,301)' -pu     Set nptdiag for uniform bins   [',nptdiag
     $     ,'  (Hence 2-D f plots)'
      write(*,301)' -q      Output diagnostics of file reading'
      write(*,301)' --objfile<filename>  set name of object data file.'
     $     //' [copticgeom.dat'
      write(*,301)' -vtk    Output distribution function vtk files.'
      write(*,301)' -h -?   Print usage.'
      call exit(0)
      end
!************************************************************************
