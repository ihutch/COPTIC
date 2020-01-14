      program coptic
! Main program of cartesian coordinate, oblique boundary, pic code.
      implicit none
      include 'ndimsdecl.f'
! Mesh spacing description structure includes griddecl too.
      include 'meshcom.f'
! Coptic here begins to assume it is 3-D. 
! However allocation is to ndimsmax to allow adjustment to ndims.
! Consequently it may be feasible to put na_k=1 and do 2-d.
! coptic runs correctly with unequal dimensions but phiexamine does not.
      integer Li1,Li2,Li3
      parameter (Li1=na_i,Li2=Li1*na_j,Li3=Li2*na_k)
      real u(na_i,na_j,na_k),q(na_i,na_j,na_k)
! Running averages.
      real qave(na_i,na_j,na_k),uave(na_i,na_j,na_k)
      real scratch(na_i,na_j,na_k)
      double complex cscratch(na_i+4,na_j+4,na_k+4)
      real psum(na_i,na_j,na_k),volumes(na_i,na_j,na_k),voltotal
      real cij(2*ndimsmax+1,na_i,na_j,na_k)
! Diagnostics (moments)
      integer ndiagmax
      parameter (ndiagmax=7)
      real diagsum(na_i,na_j,na_k,ndiagmax+1,nspeciesmax)
! Distribution function lengths
      integer ndistmax
      parameter (ndistmax=300)
! Used dimensions, Full dimensions. Used dims-2
      integer iuds(ndimsmax),ifull(ndimsmax),ium2(ndimsmax)
! Processor cartesian geometry can be set by default.
      integer nblksi,nblksj,nblksk
      parameter (nblksi=99,nblksj=99,nblksk=99)
      integer idims(ndimsmax)
! mpi process information.
      include 'myidcom.f'
! Common data containing the BC-object geometric information
      include '3dcom.f'
! Structure vector needed for finding adjacent u values.
      integer iLs(ndimsmax+1)
! Particle common data
      include 'partcom.f'
! Plasma common data
      include 'plascom.f'
! Collision common data
      include 'colncom.f'
! Point charge common data
      include 'ptchcom.f'
! Boundary setting common data
      include 'slpcom.f'
! Face boundary data
      include 'facebcom.f'
! Particle distribution accumulation data
      include 'ptaccom.f'
! Debugging
      include 'dbgcom.f'
      integer idebug
      external bdyshare,bdyset,cijroutine,cijedge,psumtoq
     $     ,quasineutral,fadcomp,qvary
      real fadcomp
      external volnode,linregion
      character*100 partfilename,phifilename,fluxfilename,objfilename
      character*100 restartpath,phafilename
      character*256 argline
!      common /ctl_sor/mi_sor,xjac_sor,eps_sor,del_sor,k_sor
      logical ltestplot,lcijplot,lsliceplot,lorbitplot,linjplot
      logical lmyidhead,lphiplot,ldenplot,lfftsucceeded
      integer ipstep,iwstep,idistp,idcount,icijcount,lrestart
! Diagnostics etc
      real zp(na_m,na_m)
      real xlimit(2,ndimsmax),vlimit(2,ndimsmax)
      real xnewlim(2,ndimsmax)
! Input for face boundary data:
      real CFin(3+ndimsmax,2*ndimsmax)
! Center of objplot final plot.
      real cv(ndimsmax)
! Particle initialization wave displacement specification
      integer nwspec
      parameter (nwspec=2*ndims+1)
      real wavespec(nwspec)
! Various local parameters
      real bckgd,bdt,bdtnow,boltzamp0,cellvol,dtf
      real dum,dum2,dum3,dum4,dum5
      real error,pinjcomp0(nspeciesmax),rc,rcij,rinf,rs1,thetain
      integer i,iavesteps,ibinit,iCFcount,ickst,ictl,id,idiag,idn
      integer ied,ierr,ierrsor,ifix,ifobj,ifplot,iobpl
      integer ipoint,ispecies,istat,istepave,iobpsw,j,k,maccel
      integer mbzero,ninjcomp0(nspeciesmax),nstep,nsteps,nth,ndiags
! And Functions
      integer lentrim,nbcat,nameappendint,oicijfunc
      external lentrim,nbcat,nameappendint,oicijfunc
      
! Set up the structure vector.
      data iLs/1,Li1,Li2,Li3/
! Mesh and mpi parameter defaults:
      data idims/nblksi,nblksj,nblksk/
      data ifull/na_i,na_j,na_k/
! Data for plotting etc.
      data iobpl/0/
      data idebug/0/
      data ltestplot,lcijplot,lsliceplot,lorbitplot,linjplot/
     $     .false.,.false.,.false.,.true.,.false./
      data lphiplot,ldenplot/.false.,.false./
      data lfftsucceeded/.false./
      data lrestart/0/cv/0.,0.,0./
      data ipstep/1/idistp/0/idcount/0/icijcount/0/
      data wavespec/nwspec*0/  ! Default no wave
! End of declarations      ###################################
!-------------------------------------------------------------
! Replace Block Data programs with this 
      call blockdatainit()
!-------------------------------------------------------------
! Initialize the fortran random number generator with a fixed number
! for solutions of volumes etc. Each node then does the same.
      call rluxgo(0,0,0,0)
! Determine what reinjection scheme we use. Sets rjscheme.
      include 'REINJECT.f'
!---------------------------------------------------------------------
! Find out early (here) the mpi structure and my mpi id number.
      call mpigetmyid(myid,nprocs,ierr)
      lmyidhead=myid.eq.0
! numprocs is the parameter in partcom kept separately for some reason
      numprocs=nprocs
! This barrier is to try to make any error in parameter setting stop
! all processes in such a way that the whole run stops.
      call mpibarrier(ierr)
!--------------------------------------------------------------
! Deal with command-line arguments and geometry/object file.
      call parametersetting
     $     (lmyidhead,ltestplot,iobpl,iobpsw,rcij
     $     ,lsliceplot,ipstep,ldenplot,lphiplot,linjplot,ifplot,norbits
     $     ,thetain,nth,iavesteps,nparta,ripernode,crelax,ickst
     $     ,colntime,dt,bdt,subcycle,dropaccel,eoverms,Bfield,Bt
     $     ,ninjcomp,nsteps,nf_maxsteps,vneutral,vds,ndiags,ndiagmax
     $     ,debyelen,Ts,iwstep,idistp,lrestart,restartpath,extfield
     $     ,objfilename,lextfield,vpars,vperps,ndims,islp,slpD,CFin
     $     ,iCFcount,LPF,ipartperiod,lnotallp,Tneutral,Enfrac,colpow
     $     ,idims,argline,vdrifts,ldistshow,gp0,gt,gtt,gn,gnt,nspecies
     $     ,nspeciesmax,numratioa,Tperps,boltzamp,nptdiag,nqblkmax
     $     ,holelen,holepsi,holeum,holeeta,holepow,holerad,hspecies
     $     ,holegfac,wavespec,LNPF,ifull,ierr)
! Hack to prevent incompatible particles
      if(nparta(1).ne.0 .and.
     $     (ipartperiod(1).ne.0.or.ipartperiod(2).ne.0
     $     .or.ipartperiod(3).ne.0)) stop '-ni not allowed with pp'
!-----------------------------------------------------------------
! Finalize initial parameters after switch and geometry reading.
      call initializeparams(ifull,iuds,xlimit,vlimit,xnewlim
     $     ,boltzamp0,cellvol,ifix)
!---------------------------------------------------------------
! Construct the mesh vector(s) and ium2
 250  call meshconstruct(ndims,iuds,ifull,ipartperiod,rs)
      if(bgn(1)+bgn(2)+bgn(3).gt.0)call bgaadj
      if(lmyidhead)write(*,'(a,3i4,6f8.3)')
     $     ' Constructed mesh',iuds
     $     ,(xmeshstart(k),xmeshend(k),k=1,ndims)
      rs1=0.5*abs(xmeshend(1)-xmeshstart(1))
      if(lmyidhead.and.rs.gt.rs1 .and. islp.eq.0.and.iCFcount.le.0)
     $     write(*,*)'UNUSUAL choice islp=',islp
     $     ,' when non-square domain in use.'
      voltotal=1.
      do i=1,ndims
         voltotal=voltotal*(xmeshend(i)-xmeshstart(i))
! If any dimension has more than one mesh step defeat fftw use.
         if(imeshstep(i,3).ne.0)LNPF=.true.
      enddo
! Initialize reinjection geometry if needed for particular case.
      call geominit(myid)
!-----------------------------------------------------------------
! Initialize the face phi boundary conditions if we are using them.
! Actually, always, sinc CFcount is unreliable in defaults.
      do idn=1,2*ndims
         call bdyfaceinit(idn,CFin(1,idn))
      enddo
! Now print out the result of the initialization.
      if(lmyidhead)call bdyfaceinit(-2,CFin(1,1))
!---------------------------------------------------------------
      call ninjcalc(dt)
!----------------------------------------------------------------
      nf_species=nspecies
      nf_nsteps=nsteps
! Part of the restart code needs to be here to determine the required
! total number of steps for which the flux initialization is needed.
! Since some must be here, we construct the names here and not later.
      call restartnames(lrestart,partfilename,restartpath,phifilename
     $     ,phafilename,fluxfilename,nf_nsteps,nsteps)
! Initialize the fluxdata storage and addressing before cijroutine
! That is necessary because ijbin addressing is used in cijroutine for
! subsequent access by cijdirect when calculating floating potential.
      call fluxdatainit(myid)
!-----------------------------------------------------------------
! Initialize cij:
      if(lmyidhead)write(*,*)'Initializing the stencil data cij.'
      error=0.
      ipoint=0
      do id=1,ndims
         ium2(id)=iuds(id)-2
         ipoint=ipoint+iLs(id)
      enddo
      call mditerarg(cijroutine,ndims,ifull,ium2,ipoint,
     $     cij(1,1,1,1),debyelen,error,dum4,dum5)
      if(error.ne.0.)then
         icijcount=icijcount+1
         if(icijcount.le.2)then
            if(lmyidhead)write(*,*)'cijroutine warnings',int(error)
     $        ,' Shifting mesh and recalculating.'
            call meshshift()
            goto 250
         else
            if(lmyidhead)write(*,*)'Failed to avoid cij warnings'
     $           ,int(error)
         endif
      endif
!---------------------------------------------
! Here we try to read the stored geometry volume data.
      istat=1
      call stored3geometry(volumes,iuds,ifull,istat,.true.)
      call mpibarrier(ierr)
! We need to wait till all nodes have tried to read before
! we go ahead and write, else we get inconsistent file status.
! An istat=1 return says we succeeded. If so skip the calculation. 
      if(istat.eq.0)then
! Calculate the nodal volumes for all non-edge points.
         ipoint=0
         do id=1,ndims
            ipoint=ipoint+iLs(id)
         enddo         
         if(lmyidhead)write(*,*)
     $        'Starting volume setting. Be patient this first time...'
         call mditerarg(volnode,ndims,ifull,ium2,ipoint,
     $        volumes,cij,dum3,dum4,dum5)
         if(lmyidhead)write(*,*)'Finished volume setting'
! If head, write the geometry data if we've had to calculate it.
! That means it isn't saved if writing to local slave disks.
         if(lmyidhead)call stored3geometry(volumes,iuds,ifull,istat
     $        ,.true.)
      endif
!---------------------------------------------
! Initialize the region flags in the object data
      call iregioninit(ifull)
      if(lmyidhead)call reportfieldmask()
! Set an object pointer region -1 for all the non periodic edges
      ipoint=0
      call mditerarg(cijedge,ndims,ifull,iuds,ipoint,cij,dum2,dum3,dum4
     $     ,dum5)
      if(lmyidhead)write(*,*) 'Used No of pointers:',oicijfunc(),' of'
     $     ,iuds(1)*iuds(2)*iuds(3) ,' points.'
!---------------------------------------------
      if(.not.lmyidhead)then
! Don't do plotting from any node except the master.
         ltestplot=.false.
         lcijplot=.false.
         lsliceplot=.false.
         lorbitplot=.false.
         norbits=0
         iobpl=0
      else
!---------------------------------------------
! Some simple graphics of cij, and volumes.
         if(ltestplot)call text3graphs(ndims,iuds,ifull,cij,volumes)
! More elaborate graphics of volumes. Was ltestplot.
!         if(.false.)call sliceGweb(ifull,iuds,volumes,na_m,zp,
!     $              ixnp,xn,ifix,'volumes:'//'!Ay!@'//char(0),dum,dum)
! Plot objects 0,1 and 2 (bits)
         if(iobpl.ne.0.and.lmyidhead)then
            if(rcij.eq.0.)rcij=rs
            write(*,*)'cijplot',rcij
            call cijplot(ifull,iuds,cij,rcij,iobpl)
          endif
      endif
!---------------------------------------------
! Initialize charge, potential (set q,u to zero over entire array).
      call mditerset(q,ndims,ifull,iuds,0,0.)
      call mditerset(qave,ndims,ifull,iuds,0,0.)
      call mditerset(u,ndims,ifull,iuds,0,0.)
      call mditerset(uave,ndims,ifull,iuds,0,0.)
! Initialize additional potential and charge.
      call setadfield(ifull,iuds,iptch_mask,ltestplot)
!---------------------------------------------------------------     
! Initial potential solution.
! Control. Bit 1, use my sor params (not here). Bit 2 use fadcomp (not)
      ictl=0
      do id=1,ndims
! Make dimensions periodic as necessary: Bit 4-6 periodicity.
         if(LPF(id))ictl=ictl+4*2**id
      enddo
! An initial solver call with zero density.
      ierrsor=0
      if(debyelen.ne.0.)call sormpi(ndims,ifull,iuds,cij,u,q,bdyshare
     $     ,bdyset,fadcomp,ictl,ierrsor,myid,idims)
!      write(*,*)'Monitor u:',u(iuds(1)/4,iuds(2)/4,iuds(3)/4)
      ictl=2+ictl
      if(ltestplot)call sliceGweb(ifull,iuds,u,na_m,zp,
     $              ixnp,xn,ifix,'potential:'//'!Ay!@'//char(0),dum,dum)
      if(lmyidhead)then
         call vaccheck(ifull,iuds,cij,u,thetain,nth,rs,ltestplot)
      endif
! End of plotting.
!------------------------------------------------------------------
! (Re)Initialize the fortran random number generator.
      call rluxgo(1,myid,0,0)
! Initialize particles  ! Only if we are not restarting.
      if(lrestart-4*(lrestart/4).eq.0)then
         if(nqblkmax.le.0.)then
! Standard noisy particle initialization
            if(holepsi.eq.0)then
               call pinit(subcycle)
            else
! Old Hole initialization 
               call trapinit(subcycle)
            endif
         else
! Quiet initialization. Use only with no objects.
            if(ngeomobj.gt.0)stop '!!Quiet init only with no objects!!'
            call qinit(wavespec)
         endif
      endif
!---------------------------------------------
! Initialize the force tracking.
      call forcetrackinit()
!---------------------------------------------
! Acceleration parameters etc.
      phirein=0.
      do i=1,nspecies
         ninjcomp0(i)=ninjcompa(i)
         pinjcomp0(i)=pinjcompa(i)
      enddo
      if(ninjcomp.ne.0.and.lmyidhead)
     $     write(*,'(a,i8,f8.5,i8,f8.5)')' Fixed injection count:'
     $     ,(ninjcompa(i),pinjcompa(i),i=1,nspecies)
      maccel=nsteps/3
      dtf=dt
      mbzero=maccel
!-----------------------------------------------
! Restart read in if we are restarting
      call restartread(lrestart,fluxfilename,partfilename,nsteps
     $        ,nf_step,nf_maxsteps,phifilename,ifull,iuds,ied,u,ierr
     $        ,lmyidhead,myid,phafilename,uave)
!-----------------------------------------------
      if(lmyidhead)then
         if(colntime.ne.0)then
            write(*,'(a,f8.2,a,f8.4)')' Collision time='
     $        ,colntime,', Eneutral=',Eneutral
         endif
         write(*,'(/,a)')'Step Iterations Flux:'
      endif
! This call is necessary for restart. It deposits the restarted particles
! into psum, giving the required status immediately after padvnc when
! padvnc includes the charge deposition.
      call mditerset(psum,ndims,ifull,iuds,0,0.)
      call chargetomesh(psum,iLs,diagsum,ndiags) 
!      write(*,*)'diagsum22221',diagsum(2,2,2,2,1)
! This writes out the initializing diagnostics. For debugging purposes
      if(iavesteps.eq.1)call periodicwrite(ifull,iuds,iLs
     $     ,diagsum,uave,lmyidhead,ndiags,ndiagmax,0,nsteps
     $     ,idistp,vlimit,xnewlim,cellvol,ibinit,idcount,restartpath
     $     ,iavesteps)
! Initialize diagsum here, not at potential initialization as before,
! because chargetomesh does deposition which duplicates that by padvnc.
      do ispecies=1,nspecies
         do idiag=1,ndiags+1
            call mditerset(diagsum(1,1,1,idiag,ispecies)
     $           ,ndims,ifull,iuds,0,0.)
         enddo
      enddo
      nstep=nf_step
!----------------------------------------------------------
! #### Main step iteration ##########################################
      do j=1,nsteps
         if(nspecies.gt.1.or.holepsi.ne.0.)
     &        boltzamp=max(0.,boltzamp0*(mbzero-nf_step+2)/(mbzero+1.))
         nstep=nstep+1
         nf_step=min(nf_maxsteps-1,nf_step+1) 
! Don't increment nf_step beyond its allowed maximum.
         istepave=min(nstep,iavesteps)
! Transfer deposits periodically or from ghost cells.
         call diagperiod(psum,ifull,iuds,iLs,1)
         call psumreduce(psum,nrein,phirein,ndims,ifull,iuds,iLs)
!------------- 
! Set background if needed.
         call bckgdset(bckgd,bdt,bdtnow,dtf,maccel,boltzamp
     $     ,ispecies,pinjcomp0,ninjcomp0,voltotal,nf_step)

! Calculate rhoinfinity, needed in psumtoq. Dependent on reinjection type.
! Does nothing if ninjcomp and rhoinf!=0
         call rhoinfcalc(dt)
!         if(j.eq.1)rhoinf=rhoinf*0.9999999
!-------------
! If more than one species, no background subtraction.
         if(nspecies.gt.1)bckgd=0.
! Convert psums to charge density, q. Remember external psumtoq!
         call mditerarg(psumtoq,ndims,ifull,ium2,0,psum(2,2,2)
     $        ,q(2,2,2),volumes(2,2,2),bckgd,rhoinf)
! Reset psum, after psumtoq.
         call mditerset(psum,ndims,ifull,iuds,0,0.)
! Possibly adjust q for variable background
         if(bckgd.ne.0.and.bgn(1)+bgn(2)+bgn(3).gt.0) call
     $        mditerarg(qvary,ndims,ifull,ium2,0,q(2,2,2),bckgd
     $        ,boltzamp)         

! -------------- Solve for the new potential:-------------------
         if(idebug.gt.0)write(*,*)'Calling solver',iuds,bckgd,boltzamp
         if(debyelen.eq.0)then
            call mditerarg(quasineutral,ndims,ifull,ium2,
     $        0,q(2,2,2),u(2,2,2),volumes(2,2,2),uc(2,2,2),dum5)
            call bdyslope0(ndims,ifull,iuds,cij,u,q)
         else
            if(boltzamp.ne.0)then
! Use some fraction of linearized Boltzmann electron response.
! To damp long wavelength oscillations.
               call sormpi(ndims,ifull,iuds,cij,u,q,bdyshare,bdyset
     $              ,fadcomp,ictl,ierr,myid,idims)
!               write(*,*)'StepMon u:',u(iuds(1)/4,iuds(2)/4,iuds(3)/4)
            else
               if(.not.LNPF)then ! All but x periodic. 
! Assume we want to use the FFT poisson solver when all periodic uniform.
                  if(LPF(1))then ! x also periodic
                     call fftphisolve3d(ifull,iuds,u,q,cscratch
     $                    ,xmeshend(1)-xmeshstart(1)
     $                    ,xmeshend(2)-xmeshstart(2)
     $                    ,xmeshend(3)-xmeshstart(3)
     $                    ,ierr)
                  else
                     call ffttrid(ifull,iuds,u,q,cscratch
     $                    ,xmeshend(1)-xmeshstart(1)
     $                    ,xmeshend(2)-xmeshstart(2)
     $                    ,xmeshend(3)-xmeshstart(3)
     $                    ,ierr)
                  endif
                  if(ierr.ne.0)then
                     lfftsucceeded=.false.
                     write(*,*)'FAILED ATTEMPT to use FFTW Solver'
                  else
                     if(.not.lfftsucceeded.and.lmyidhead)
     $                    write(*,*)'Using FFTW Solver'
                     lfftsucceeded=.true.
                  endif
               endif
               if(.not.lfftsucceeded)then
! SORMPI vanilla solve. Turn off use of fadcomp by ictl-2
                  call sormpi(ndims,ifull,iuds,cij,u,q,bdyshare,bdyset
     $              ,fadcomp,ictl-2,ierr,myid,idims)
               endif
            endif
         endif
         if(idebug.gt.0)write(*,*)'Returned from solver'
         if(LPF(1).and.LPF(2).and.LPF(3).and.ngeomobj.eq.0)
     $        call makemedianzero(ifull,iuds,u,scratch) ! remove offset.
! ------------------------------------------------
         call calculateforces(ndims,iLs,cij,u)
! ------------------------------------------------
         if(ipstep.eq.0.or.mod(j,ipstep).eq.0)then
! Slice plots
!            write(*,*)'lsliceplot,ldenplot,ldistshow,lphiplot'
!     $           ,lsliceplot,ldenplot,ldistshow,lphiplot
            if(lsliceplot.and.ldenplot)call sliceGweb(ifull,iuds,q,na_m
     $           ,zp,ixnp,xn,ifix,'density: n'//char(0),dum,dum)
            if((iuds(2).eq.3.and.iuds(3).eq.3))then
               if(.not.ldistshow)then
                  if(lphiplot)then
                     if(myid.eq.0)write(*,*
     $                 )'1-d does not plot potential; use both -gn -gp'
                     ldistshow=.true.
                  endif
               endif
            else
               if(lsliceplot.and.lphiplot)then
                  call sliceGweb(ifull,iuds,u,na_m,zp,ixnp,
     $                 xn,ifix,'potential:'//'!Af!@'//char(0),dum,dum)
               endif
            endif
! Phase space done by all processes, even though only one of them plots.
            if(ldistshow.and.iuds(2).eq.3.and.iuds(3).eq.3)then
               call phasepscont(ifull,iuds,u,nstep,lphiplot
     $              ,restartpath)
            endif
         endif

         if(nf_step.eq.ickst)call checkuqcij(ifull,u,q,psum,volumes,cij)
!----------- Particle advance:------------------------------
         if(idebug.gt.0)write(*,*)'Calling padvnc',ispecies
         do ispecies=1,nspecies
            call padvnc(iLs,cij,u,ndiags,psum,diagsum,ispecies,ndiagmax)
         enddo
!------------------------------------------------
! Optional checking step. Used for restart testing.
         if(nf_step.eq.ickst)call checkx

         call fluxreduce()
! Now do cij update
         call cijdirect(debyelen,error)
! Store the step's rhoinf, dt, npart.
         ff_rho(nf_step)=rhoinf
         ff_dt(nf_step)=dt
         nf_npart(nf_step)=n_part

! Report step values etc.
         if(lmyidhead)call reportprogress(nstep,nsteps,nsubc,ndropped
     $        ,ierr,iavesteps)

! These running and box averages do not include the updates for this step.
! Accumulate running q and u averages (replaced 3d routines):
         call averagegd(q,qave,ifull,iuds,istepave)
         call averagegd(u,uave,ifull,iuds,istepave)
! Every iavesteps, calculate the box average of the moments, and write it
! out, if we are doing diagnostics.
         if(mod(j,iavesteps).eq.0)then
            call periodicwrite(ifull,iuds,iLs
     $        ,diagsum,uave,lmyidhead,ndiags,ndiagmax,nstep,nsteps
     $        ,idistp,vlimit,xnewlim,cellvol,ibinit,idcount,restartpath
     $        ,iavesteps)
            if(bdt.lt.0.and.lmyidhead)then
               write(*,'(a,f8.2,3f8.4)')'injcomp,rhoinf,phirein,bckgd='
     $              ,dum,rhoinf,phirein,bckgd
            endif
            if(.not.ldistshow.and.idistp.ne.0.and.
     $           iuds(2).eq.3.and.iuds(3).eq.3)then
! If we have not already written phaseplot, but we are writing distribs
! and it is a 1-d calculation:
               call phasepscont(ifull,iuds,u,nstep,lphiplot
     $              ,restartpath)
            endif
         endif
! Particle distribution diagnostics.
! The serial cost for this call with 1M particles is about 1s in 17s.
! Or less than 2s if both partaccum and vaccum are called. Thus the
! total cost is roughly 10% of particle costs. Tell it that it must
! set dimensions the first time by appending 0.
! Update the particle distribution diagnostics for species 1:
         if(idcount.gt.0)call partdistup(xlimit,vlimit,xnewlim,
     $        cellvol,myid,0,1)
! Sometimes write them out:
         if(iwstep.gt.0.and.mod(nstep,iwstep).eq.0)call datawrite(myid
     $        ,partfilename,restartpath,ifull,iuds,u,uave,qave,nstep)

! This non-standard fortran call works with gfortran and g77 to flush stdout.
! Pathscale demands an argument number. So give it explicitly.
! Comment it out if it causes problems.
         if(lmyidhead)call flush(6)
      enddo
! #### End of Main Step Iteration ####################################
!----------------------------------------------------------
      if(norbits.ne.0)call cijplot(ifull,iuds,cij,rs,iobpl)
      if(lorbitplot.and.norbits.ne.0)call orbitplot(ifull,iuds,u,phip,rc
     $     ,rs)
!-------------------------------------------------------------------
! Everyone writes what they have to.
      if(iwstep.gt.0 .or. myid.eq.0)call datawrite(myid,partfilename
     $     ,restartpath,ifull,iuds,u,uave,qave,nstep)
!-------------------------------------------------------------------
      call mpifinalize(ierr)
! Check some flux diagnostics and writing. Really nf_step.
      if(nf_step.eq.nstep)call finaldiags(lmyidhead,linjplot,Ti,mf_obj
     $     ,nf_step,rinf,ifobj,ifplot,rcij,rs,cv,iobpsw)

      end
! End of Main Program.
!************************************************************************
