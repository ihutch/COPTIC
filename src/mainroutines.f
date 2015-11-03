c************************************************************************
c************************************************************************
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
c-------------------------------------------------------------
c This circumlocution to silence warnings.
      i1=ndims+1
      do id=i1,ndimsmax
         ifull(id)=1
         iuds(id)=1
      enddo
      do id=1,ndimsmax
c Use very big xlimits by default to include whole domain
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

c Set phip from the first object. Must be done before nameconstruct.
      call phipset(myid)

c-------------------------------------------------------------
      end
c********************************************************************

      subroutine restartnames(lrestart,partfilename,restartpath
     $     ,phifilename,fluxfilename,nf_nsteps,nsteps)
      implicit none
      integer lrestart,nf_nsteps,nsteps
      character*100 partfilename,phifilename,fluxfilename,restartpath
      
      include 'myidcom.f'
      integer nb,nbcat,lentrim,nameappendint,iferr
      if(lrestart.ne.0)then
c Part of the restart code needs to be here to determine the required
c total number of steps for which the flux initialization is needed.
c Since some must be here, we construct the names here and not later.
         partfilename=restartpath
         if(lrestart/4-2*(lrestart/8).ne.0)then
            partfilename(lentrim(partfilename)+1:)='restartfile'
         else
            call nameconstruct(partfilename)
         endif
         phifilename=partfilename
         nb=nbcat(phifilename,'.phi')
         fluxfilename=partfilename
         nb=nbcat(fluxfilename,'.flx')
         nb=nameappendint(partfilename,'.',myid,3)
         if(lrestart/2-2*(lrestart/4).ne.0)then
            iferr=0
c            write(*,*)'Reading flux file:',fluxfilename
            call readfluxfile(fluxfilename,iferr)
c The total number of steps for fluxdatainit is the sum of what we
c just read out of fluxfile and the new nsteps:
            nf_nsteps=nf_nsteps+nsteps
         endif
      endif

      end
c***********************************************************************
      subroutine restartread(lrestart,fluxfilename,partfilename,nsteps
     $     ,nf_step,nf_maxsteps,phifilename,ifull,iuds,ied,u,ierr
     $     ,lmyidhead,myid)
      implicit none
      integer lrestart,nsteps,nf_step,nf_maxsteps,ied,ierr,myid
      logical lmyidhead
      character*100 partfilename,phifilename,fluxfilename
      integer ifull(*),iuds(*)
      real u(*)
      
      integer iferr
      integer lentrim
      external lentrim

      if(lrestart.ne.0)then
c names are constructed earlier.
         if(lrestart/2-2*(lrestart/4).ne.0)then
            iferr=0
            call readfluxfile(fluxfilename,iferr)
         else
            iferr=1
         endif
         if(lrestart-4*(lrestart/4).ne.0)then
            call partread(partfilename,ierr)
            if(ierr-4*(ierr/4).eq.0)then 
c We succeeded in reading the part-file. Relocate the particles.
               write(*,'(a,i4,a,a,i3)')' cpu',myid
     $              ,' Restart file read: '
     $              ,partfilename(1:lentrim(partfilename)+1),lrestart
               call locateinit()
               if(nsteps+nf_step.gt.nf_maxsteps)then
                  if(lmyidhead)write(*,*)'Asked for',
     $                 nsteps,' in addition to',nf_step,
     $                 ' Total',nsteps+nf_step,
     $                 ' too much; set to',nf_maxsteps-1
                  nsteps=nf_maxsteps-nsteps-1
               endif
               ied=1
c Only read the phi-file if the flux file was present. Full restart.
               if(iferr.eq.0)then
c                  write(*,*)'Reading phifile',ierr,phifilename
                  call array3read(phifilename,ifull,iuds,ied,u,ierr)
               endif
            endif
         endif
c In case we have overwritten phip with the value from the restart file,
c try to set it again from the first object. But tell it we are not
c the head node, so it does not give out messages.
         call phipset(1)
      endif

      end
c**********************************************************************
      subroutine finaldiags(lmyidhead,linjplot,Ti,mf_obj,nf_step,rinf
     $     ,ifobj,ifplot,rcij,rs,cv,iobpsw)
      implicit none
      logical lmyidhead,linjplot
      real Ti,rinf,rcij,rs,cv(*)
      integer mf_obj,nf_step,ifobj,ifplot,iobpsw

c Check some flux diagnostics and writing.
      if(lmyidhead)then 
         if(linjplot)call plotinject(Ti)
         do ifobj=1,mf_obj
            call fluxave(nf_step/2,nf_step,ifobj,ifplot,rinf)
         enddo
c         write(*,*)'Calling objplot'
         if(ifplot.gt.0)then
            if(rcij.le.0)rcij=rs
            call objplot(1,rcij,cv,iobpsw,0)
         endif
      endif

      end
