!************************************************************************
!************************************************************************
      subroutine initializeparams(ifull,iuds,xlimit,vlimit,xnewlim
     $     ,boltzamp0,cellvol,ifix)
      implicit none
      include 'ndimsdecl.f'
      include 'griddecl.f'
      include 'partcom.f'
      include 'colncom.f'
      include 'ptaccom.f'
      include 'plascom.f'
      include 'ptchcom.f'
      include 'myidcom.f'
      integer ifull(ndimsmax),iuds(ndimsmax)
      real xlimit(2,ndimsmax),vlimit(2,ndimsmax)
      real xnewlim(2,ndimsmax)
      real boltzamp0,cellvol
      integer ifix

      integer id,i1
!-------------------------------------------------------------
! This circumlocution to silence warnings.
      i1=ndims+1
      do id=i1,ndimsmax
         ifull(id)=1
         iuds(id)=1
      enddo
      do id=1,ndimsmax
! Use very big xlimits by default to include whole domain
         xlimit(1,id)=-500.
         xlimit(2,id)=500.
         xnewlim(1,id)=0.
         xnewlim(2,id)=0.
         vlimit(1,id)=5.
         vlimit(2,id)=-5.
      enddo
      cellvol=0.
      rs=5.0
      caverein=0.
      chi=0.
      ifix=0
      rhoinf=0.
      ndropped=0
      boltzamp0=boltzamp
      boltzsign=sign(1.,eoverms(1))
      if(nptdiag.eq.0)nptdiag=nsbins
      if(n_part.ne.0)ripernode=0.

      call initdriftfield

! Set phip from the first object. Must be done before nameconstruct.
      call phipset(myid)

!-------------------------------------------------------------
      end
!********************************************************************

      subroutine restartnames(lrestart,partfilename,restartpath
     $     ,phifilename,phafilename,fluxfilename,nf_nsteps,nsteps)
      implicit none
      integer lrestart,nf_nsteps,nsteps
      character*100 partfilename,phifilename,fluxfilename,restartpath
      character*100 phafilename
      include 'myidcom.f'
      integer nb,nbcat,lentrim,nameappendint,iferr
      if(lrestart.ne.0)then
! Part of the restart code needs to be here to determine the required
! total number of steps for which the flux initialization is needed.
! Since some must be here, we construct the names here and not later.
         partfilename=restartpath
         if(lrestart/4-2*(lrestart/8).ne.0)then
            partfilename(lentrim(partfilename)+1:)='restartfile'
         else
            call nameconstruct(partfilename)
         endif
         phifilename=partfilename
         nb=nbcat(phifilename,'.phi')
         phafilename=partfilename
         nb=nbcat(phafilename,'.pha')
         fluxfilename=partfilename
         nb=nbcat(fluxfilename,'.flx')
         nb=nameappendint(partfilename,'.',myid,
     $     max(3,int(log10(float(nprocs))+1.)))
         if(lrestart/2-2*(lrestart/4).ne.0)then
            iferr=0
!            write(*,*)'Reading flux file:',fluxfilename
            call readfluxfile(fluxfilename,iferr)
! The total number of steps for fluxdatainit is the sum of what we
! just read out of fluxfile and the new nsteps:
            nf_nsteps=nf_nsteps+nsteps
         endif
      endif

      end
!***********************************************************************
      subroutine restartread(lrestart,fluxfilename,partfilename,nsteps
     $     ,nstep,nf_maxsteps,phifilename,ifull,iuds,ied,u,ierr
     $     ,lmyidhead,myid,phafilename,uave)
      implicit none
      integer lrestart,nsteps,nstep,nf_maxsteps,ied,ierr,myid
      logical lmyidhead
      character*100 partfilename,phifilename,fluxfilename,phafilename
      integer ifull(*),iuds(*)
      real u(*),uave(*)
      
      integer iferr
      integer lentrim
      external lentrim

      if(lrestart.ne.0)then
! names are constructed earlier.
         if(lrestart/2-2*(lrestart/4).ne.0)then
            iferr=0
            call readfluxfile(fluxfilename,iferr)
         else
            iferr=1
         endif
         if(lrestart-4*(lrestart/4).ne.0)then
            call partread(partfilename,ierr)
            if(ierr-4*(ierr/4).eq.0)then 
! We succeeded in reading the part-file. Relocate the particles.
               write(*,'(a,i4,a,a,i3)')' cpu',myid
     $              ,' Restart file read: '
     $              ,partfilename(1:lentrim(partfilename)+1),lrestart
               call locateinit()
               if(nsteps+nstep.gt.nf_maxsteps)then
                  if(lmyidhead)write(*,*)'Asked for',
     $                 nsteps,' in addition to',nstep,
     $                 ' Total',nsteps+nstep,
     $                 ' will exceed flux storage',nf_maxsteps
!     $                 ' too much; set to',nf_maxsteps-1
!                  nsteps=nf_maxsteps-nsteps-1
               endif
               ied=1
! Only read the phi-file if the flux file was present. Full restart.
               if(iferr.eq.0)then
!                  write(*,*)'Reading phifile',ierr,phifilename
                  ierr=0  ! to suppress commentary
                  call array3read(phifilename,ifull,iuds,ied,u,ierr)
                  call array3read(phafilename,ifull,iuds,ied,uave,ierr)
               endif
            else
               write(*,*)'Restartread failed to read',partfilename
               stop
            endif
         endif
! In case we have overwritten phip with the value from the restart file,
! try to set it again from the first object. But tell it we are not
! the head node, so it does not give out messages.
         call phipset(1)
      endif

      end
!**********************************************************************
      subroutine finaldiags(lmyidhead,linjplot,Ti,mf_obj,nstep,rinf
     $     ,ifobj,ifplot,rcij,rs,cv,iobpsw)
      implicit none
      logical lmyidhead,linjplot
      real Ti,rinf,rcij,rs,cv(*)
      integer mf_obj,nstep,ifobj,ifplot,iobpsw

! Check some flux diagnostics and writing.
      if(lmyidhead)then 
         if(linjplot)call plotinject(Ti)
         do ifobj=1,mf_obj
            call fluxave(nstep/2,nstep,ifobj,ifplot,rinf)
         enddo
!         write(*,*)'Calling objplot'
         if(ifplot.gt.0)then
            if(rcij.le.0)rcij=rs
            call objplot(1,rcij,cv,iobpsw,0)
         endif
      endif

      end
!***********************************************************************
      subroutine phasescatter(ifull,iuds,u)
      implicit none
      include 'ndimsdecl.f'
      include 'meshcom.f'
      include 'partcom.f'
      integer ifull(ndims),iuds(ndims)
      real u(ifull(1),ifull(2),ifull(3))
      integer id
      real vrange,phirange
      parameter (id=1,vrange=3.,phirange=.5)
! Only if this is a one-dimensional problem (for now)
      if(iuds(2).ge.4 .and. iuds(3).ge.4) return

      call multiframe(2,1,0)
      call pltinit(xmeshstart(id),xmeshend(id),-phirange,phirange)
      call axis()
      call axlabels(' ','  !Af!@')
      call polyline(xn(ixnp(1)+1),u(1,2,2),ixnp(2)-ixnp(1))
      call pltinit(xmeshstart(id),xmeshend(id),-vrange,vrange)
      call axis()
      call axlabels('x','  v')
      call winset(.true.)
      call color(1)
      call scatterxy(x_part(id,1),x_part(id+3,1),iocparta(1),idtp)
      call accisflush()
      call color(15)
      call multiframe(0,0,0)
      end
!***********************************************************************
      subroutine phasepscont(ifull,iuds,u,nstep,lplot,restartpath)
      implicit none
      include 'ndimsdecl.f'
      include 'meshcom.f'
      include 'partcom.f'
      include 'myidcom.f'
      include 'phasecom.f'
      integer ifull(ndims),iuds(ndims),nstep
      logical lplot,ldensac
      real u(ifull(1),ifull(2),ifull(3))
      character*(*) restartpath
      integer id,thespecies
      real vrange,phirange,umin,umax,psnmax
      real wx2nx,wy2ny
      parameter (id=1,vrange=3.)
      character*100 phasefilename
      character*10 string
      character*12 nlabel(2)
      integer lentrim
      external lentrim
      data phirange/0.5/thespecies/1/
      data nlabel/' !Bn!di!d!@',' !Bn!de!d!@'/
! Only if this is a one-dimensional problem (for now)
      if(iuds(2).ge.4 .and. iuds(3).ge.4) return
      ldensac=.true.
c psaccum must be asked for by all processes
      thespecies=mod(thespecies,nspecies)+1
      call psaccum(thespecies,1)
      write(string,'(f10.3)')nstep*dt
c but writing and plotting only by top process
      if(myid.eq.nprocs-1)then
         phasefilename=restartpath
         call nameconstruct(phasefilename)
         write(phasefilename(lentrim(phasefilename)+1:)
     $        ,'(''.pps'',i5.5)')nstep
         call phasewrite(phasefilename,
     $        ixnp(2)-ixnp(1),xn(ixnp(1)+1),u(1,2,2),nstep*dt)

         if(lplot)then
         call minmax(u(1,2,2),iuds(1),umin,umax)
         phirange=max(phirange,umax*.95)
         call pfset(3)
         call multiframe(2,1,1)
         call pltinit(xmeshstart(id),xmeshend(id),-phirange,phirange)
         call axis()
!         call axis2
         call axlabels(' ','  !Af!@')
         call polyline(xn(ixnp(1)+1),u(1,2,2),ixnp(2)-ixnp(1))
         call jdrwstr(wx2nx(xmeshend(id)),wy2ny(.9*phirange),string,-1.)
         if(ldensac)then
! Plot density in the same frame
            psnmax=nparta(thespecies)/npsx
            call scalewn(xmeshstart(id),xmeshend(id),
     $           0.8,1.35,.false.,.false.)
            call axptset(1.,1.)
            call ticrev
            call axis
            call ticrev
            call axptset(0.,0.)
            do thespecies=1,nspecies
               call psnaccum(thespecies,id)
               call color(thespecies+4)
               psn=psn/psnmax
               call polyline(psx,psn,npsx)
               call legendline(0.8,0.05+0.08*thespecies,
     $              0,nlabel(thespecies))
            enddo
            call color(15)
            call legendline(1.04,0.3,258,'!Bn!@')
         else
            call axis2
         endif
         call phaseplot
         write(string,'(''sp='',i1)')thespecies
         if(.not.ldensac)
     $        call jdrwstr(wx2nx(xmeshend(id)),wy2ny(0.),string,-1.)
         call color(15)
         call multiframe(0,0,0)
         call accisflush()
         call prtend(' ')
         endif
      endif
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Select the value of array arr that ranks k out of n.
! Median is quickselect(n/2,n,arr)
! The array is mixed up in the process.
! From http://www.stat.cmu.edu/~ryantibs/median/
! Apparently inspired by/copied from Numerical Recipes.
! An implementation of Hoare's algorithm TOMS64 I think. 
      real function quickselect(k,n,arr)
      integer k,n
      real arr(n)
      integer i,ir,j,l,mid
      real a,temp

      l = 1
      ir = n
 2    if (ir-l.le.1) then
         if (ir-1.eq.1) then
            if (arr(ir).lt.arr(l)) then
               temp = arr(l)
               arr(l) = arr(ir)
               arr(ir) = temp
            endif
         endif
         quickselect = arr(k)
         return
      else
         mid = (l+ir)/2
         temp = arr(mid)
         arr(mid) = arr(l+1)
         arr(l+1) = temp
         if (arr(l).gt.arr(ir)) then
            temp = arr(l)
            arr(l) = arr(ir)
            arr(ir) = temp
         endif
         if (arr(l+1).gt.arr(ir)) then
            temp = arr(l+1)
            arr(l+1) = arr(ir)
            arr(ir) = temp
         endif
         if (arr(l).gt.arr(l+1)) then
            temp = arr(l)
            arr(l) = arr(l+1)
            arr(l+1) = temp
         endif
         i = l+1
         j = ir
         a = arr(l+1)
 3       continue
         i = i+1
         if (arr(i).lt.a) goto 3
 4       continue
         j = j-1
         if (arr(j).gt.a) goto 4
         if (j.lt.i) goto 5
         temp = arr(i)
         arr(i) = arr(j)
         arr(j) = temp
         goto 3
 5       arr(l+1) = arr(j)
         arr(j) = a
         if (j.ge.k) ir = j-1
         if (j.le.k) l = i
      endif
      goto 2
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real function getmedian(ifull,iuds,u,scratch)
! Return the median of the ndims dimensional array u addressed linearly
! ifull and iuds are the full and used dimensions of u.
! Use scratch which must have size at least prod_i[iuds(i)].
      include 'ndimsdecl.f'
      integer ifull(ndims),iuds(ndims)
      real u(*),scratch(*)
      integer indi(ndims),iview(3,ndims),iscratch
! Copy the u array into scratch which is compact using iuds.
      iscratch=0
      icomplete=mditerator(ndims,iview,indi,4,iuds)
 1    ipointer=1+indexcontract(ndims,ifull,indi)
      iscratch=iscratch+1
      scratch(iscratch)=u(ipointer)
      if(mditerator(ndims,iview,indi,0,iuds).eq.0)goto 1
      getmedian=quickselect(iscratch/2,iscratch,scratch)
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine makemedianzero(ifull,iuds,u,scratch)
! Subtract its median from the array u.
      include 'ndimsdecl.f'
      integer ifull(ndims),iuds(ndims)
      real u(*),scratch(*)
      integer indi(ndims),iview(3,ndims)

      themedian=getmedian(ifull,iuds,u,scratch)
      if(.false.)write(*,*)'themedian=',themedian
      icomplete=mditerator(ndims,iview,indi,4,iuds)
 1    ipointer=1+indexcontract(ndims,ifull,indi)
      u(ipointer)=u(ipointer)-themedian
      if(mditerator(ndims,iview,indi,0,iuds).eq.0)goto 1

      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
      subroutine bckgdset(bckgd,bdt,bdtnow,dtf,maccel,boltzamp
     $     ,ispecies,pinjcomp0,ninjcomp0,voltotal,nf_step)
! Set the background subtraction in units of ninfinity.
      implicit none
      real bckgd,bdt,bdtnow,dtf,boltzamp,voltotal
      real pinjcomp0(nspecies),ninjcomp0(nspecies)
      integer maccel,ispecies,nf_step
      include 'ndimsdecl.f'
      include 'plascom.f'
      include 'partcom.f'
      real dum

      if(bdt.eq.0)then
         bckgd=(1.-boltzamp)*eoverms(1)
      elseif(bdt.gt.0)then
! Acceleration code.
         bdtnow=max(1.,(bdt-1.)*(maccel-nf_step+2)/(maccel+1.)+1.)
         dt=bdtnow*dtf
         bckgd=(1.-boltzamp)*eoverms(1)
      elseif(bdt.lt.0)then
! Density growth code. Negative -da switch instead says enhance the density
! of external plasma by increasing the injection rate bdt*t.
         rhoinf=0.
         call rhoinfcalc(dt)
! Subtract specified weight uniform background (for single-species running).
! Adjusted for prior step.
         if(nf_step.gt.1)then
            bckgd=numprocs*(1+bdt*dt)*(n_part)*(1.-boltzamp)
     $           *eoverms(1)/(voltotal*rhoinf)
            bdtnow=1.+abs(bdt)*nf_step*dt
            do ispecies=1,nspecies
               dum=bdtnow*(ninjcomp0(ispecies)+pinjcomp0(ispecies))
               ninjcompa(ispecies)=int(dum)
               pinjcompa(ispecies)=dum-ninjcompa(ispecies)
            enddo
         else
            bckgd=(1.-boltzamp)*eoverms(1)
         endif
      endif

      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
