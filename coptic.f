      program coptic
c Main program of cartesian coordinate, oblique boundary, pic code.

c Object data storage.
      include 'objcom.f'
c Storage array spatial count size
      include 'griddecl.f'
c coptic runs correctly with unequal dimensions but phiexamine does not.
      parameter (Li1=na_i,Li2=Li1*na_j,Li3=Li2*na_k)
      real u(na_i,na_j,na_k),q(na_i,na_j,na_k)
c Running averages.
      real qave(na_i,na_j,na_k),uave(na_i,na_j,na_k)
c
      real psum(na_i,na_j,na_k),volumes(na_i,na_j,na_k)
      real cij(2*ndims_cij+1,na_i,na_j,na_k)
c Diagnostics (moments)
      integer ndiagmax
      parameter (ndiagmax=7)
      real diagsum(na_i,na_j,na_k,ndiagmax)
c Distribution functions
      integer ndistmax
      parameter (ndistmax=300)
c      real fv(ndistmax,na_i,na_j,na_k)

c Used dimensions, Full dimensions. Used dims-2
      integer iuds(ndims_cij),ifull(ndims_cij),ium2(ndims_cij)
c Mesh spacing description structure
      include 'meshcom.f'
c Processor cartesian geometry can be set by default.
      integer nblksi,nblksj,nblksk
      parameter (nblksi=1,nblksj=1,nblksk=1)
      integer idims(ndims_cij)
c mpi process information.
      include 'myidcom.f'
c Common data containing the BC-object geometric information
      include '3dcom.f'
c Structure vector needed for finding adjacent u values.
      integer iLs(ndims_cij+1)
c Particle common data
      include 'partcom.f'
c Plasma common data
      include 'plascom.f'
c Collision common data
      include 'colncom.f'
c Point charge common data
      include 'ptchcom.f'
c Boundary setting common data
      include 'slpcom.f'
c Face boundary data
      include 'facebcom.f'

      external bdyshare,bdyset,faddu,cijroutine,cijedge,psumtoq
     $     ,quasineutral
      real faddu
      external volnode,linregion
      character*100 partfilename,phifilename,fluxfilename,objfilename
      character*100 diagfilename,restartpath
      character*100 argument
      character*256 argline
c      common /ctl_sor/mi_sor,xjac_sor,eps_sor,del_sor,k_sor
      logical ltestplot,lcijplot,lsliceplot,lorbitplot,linjplot
      logical lmyidhead,lphiplot,ldenplot
      integer ipstep,iwstep,idistp,idcount,icijcount,lrestart
c Diagnostics etc
      real zp(na_m,na_m,ndims_mesh)
      real xlimit(2,ndims_mesh),vlimit(2,ndims_mesh)
      real xnewlim(2,ndims_mesh)
c Input for face boundary data:
      real CFin(3+ndims_mesh,2*ndims_mesh)
c Center of objplot final plot.
      real cv(ndims_mesh)

c Set up the structure vector.
      data iLs/1,Li1,Li2,Li3/
c Mesh and mpi parameter defaults:
      data idims/nblksi,nblksj,nblksk/
      data ifull/na_i,na_j,na_k/
c Data for plotting etc.
      data iobpl/0/
      data ltestplot,lcijplot,lsliceplot,lorbitplot,linjplot/
     $     .false.,.false.,.false.,.true.,.false./
      data lphiplot,ldenplot/.false.,.false./
c      data thetain,nth/.1,1/
      data lrestart/0/cv/0.,0.,0./
      data ipstep/1/idistp/0/idcount/0/icijcount/0/
c-------------------------------------------------------------
c Consistency checks
      if(ndims.ne.ndims_cij)then
         write(*,*)'Inconsistent ndims, ndims_cij',ndims,ndims_cij
         stop
      endif
c-------------------------------------------------------------
c Initialize the fortran random number generator with a fixed number
c for solutions of volumes etc. Each node then does the same.
      rs=ran1(-1)
c Defaults:
c Determine what reinjection scheme we use. Sets rjscheme.
      include 'REINJECT.f'

      do id=1,ndims_mesh
c Use very big xlimits by default to include whole domain
         xlimit(1,id)=-500.
         xlimit(2,id)=500.
         xnewlim(1,id)=0.
         xnewlim(2,id)=0.
         vlimit(1,id)=5.
         vlimit(2,id)=-5.
c Default zero field
         Bfield(id)=0.
      enddo
      phip=0.
      cellvol=0.
      rs=5.0
      rsmesh=rs
      caverein=0.
      chi=0.
      ierr=0
c---------------------------------------------------------------------
c This necessary here so one knows early the mpi structure.
c Otherwise could have been hidden in sormpi and pass back numprocs.
      myid=0
      call mpigetmyid(myid,nprocs,ierr)
      lmyidhead=.true.
      if(myid.ne.0) lmyidhead=.false.
      numprocs=nprocs
c--------------------------------------------------------------
c Deal with command-line arguments and geometry/object file.
c First time this routine just sets defaults and the object file name.
      call copticcmdline(lmyidhead,ltestplot,iobpl,iobpsw,rcij
     $     ,lsliceplot,ipstep,ldenplot,lphiplot,linjplot,ifplot,norbits
     $     ,thetain,nth,iavesteps,n_part,numprocs,ripernode,crelax,ickst
     $     ,colntime,dt,bdt,subcycle,dropaccel,rmtoz,Bfield,Bt,ninjcomp
     $     ,nsteps ,nf_maxsteps,vneutral,vd,ndiags,ndiagmax,debyelen,Ti
     $     ,iwstep ,idistp,lrestart,restartpath,extfield,objfilename
     $     ,lextfield ,vpar,vperp,ndims,islp,slpD,CFin,iCFcount,LPF
     $     ,ipartperiod,lnotallp,Tneutral,Enfrac,colpow,idims,argline
     $     ,vdrift,ldistshow,gp0,gt,gtt,gn,gnt
     $     )
c Read in object file information.
      call readgeom(objfilename,myid,ifull,CFin,iCFcount,LPF,ierr
     $     ,argline)
c Second time: deal with any other command line parameters.
      call copticcmdline(lmyidhead,ltestplot,iobpl,iobpsw,rcij
     $     ,lsliceplot,ipstep,ldenplot,lphiplot,linjplot,ifplot,norbits
     $     ,thetain,nth,iavesteps,n_part,numprocs,ripernode,crelax,ickst
     $     ,colntime,dt,bdt,subcycle,dropaccel,rmtoz,Bfield,Bt,ninjcomp
     $     ,nsteps ,nf_maxsteps,vneutral,vd,ndiags,ndiagmax,debyelen,Ti
     $     ,iwstep ,idistp,lrestart,restartpath,extfield,objfilename
     $     ,lextfield ,vpar,vperp,ndims,islp,slpD,CFin,iCFcount,LPF
     $     ,ipartperiod,lnotallp,Tneutral,Enfrac,colpow,idims,argline
     $     ,vdrift,ldistshow,gp0,gt,gtt,gn,gnt
     $     )
      if(ierr.ne.0)stop
c The double call enables cmdline switches to override objfile settings.
c-----------------------------------------------------------------
c Finalize parameters after switch reading.
      ndropped=0
c---------------------------------------------------------------
c Construct the mesh vector(s) and ium2
 250  call meshconstruct(ndims,iuds,ifull,ipartperiod)
      if(lmyidhead)write(*,'(a,3i4,6f8.3)')
     $     ' Constructed mesh',iuds
     $     ,(xmeshstart(k),xmeshend(k),k=1,ndims)
c Localizing use of rs as the size of domain.
      rs1=0.5*abs(xmeshend(1)-xmeshstart(1))
      do id=1,ndims
         rsi=0.5*abs(xmeshend(id)-xmeshstart(id))
         if(rsi.gt.rs)then 
            rs=rsi
         endif
      enddo
      if(lmyidhead.and.rs.gt.rs1 .and. islp.eq.0.and.iCFcount.le.0)
     $     write(*,*)'UNUSUAL choice islp=',islp
     $     ,' when non-square domain in use.'
c Initialize reinjection geometry if needed for particular case.
      call geominit(myid)
c-----------------------------------------------------------------
c Initialize the face phi boundary conditions if we are using them.
      if(iCFcount.ne.0)then
         do idn=1,2*ndims
            call bdyfaceinit(idn,CFin(1,idn))
         enddo
c Now print out the result of the initialization.
         if(lmyidhead)call bdyfaceinit(-2,CFin(1,1))
      endif
c-----------------------------------------------------------------
      do id=1,ndims
         ium2(id)=iuds(id)-2
      enddo         
c---------------------------------------------------------------
c      write(*,*)'Doing ninjcalc',n_part,ripernode,dt
      if(n_part.ne.0)ripernode=0.
c Set ninjcomp if we are using ripernode
c This does not work until after we've set mesh in cartesian.
      if(ripernode.ne.0)call ninjcalc(dt)
c----------------------------------------------------------------
c Initialize the fluxdata storage and addressing before cijroutine
      call fluxdatainit(myid)
      if(lmyidhead)write(*,*)'Initializing the stencil data cij.'
c Initialize cij:
      ipoint=iLs(1)+iLs(2)+iLs(3)
      error=0.
      call mditerarg(cijroutine,ndims,ifull,ium2,ipoint,
     $     cij(1,1,1,1),debyelen,error,dum4)
      if(error.ne.0.)then
         icijcount=icijcount+1
         if(icijcount.le.2)then
            if(lmyidhead)write(*,*)'cijroutine warnings',error
     $        ,' Shifting mesh and recalculating.'
            call meshshift()
            goto 250
         else
            if(lmyidhead)write(*,*)'Failed to avoid cij warnings',error
         endif
      endif
c      write(*,*)'Finished cijroutine iteration'
c---------------------------------------------
c Here we try to read the stored geometry volume data.
      istat=1
      call stored3geometry(volumes,iuds,ifull,istat,.true.)
c don't calculate volumes testing.      istat=1
      if(istat.eq.0)then
c Calculate the nodal volumes for all non-edge points.
         ipoint=iLs(1)+iLs(2)+iLs(3)
         if(lmyidhead)write(*,*)
     $        'Starting volume setting. Be patient this first time...'
         call mditerarg(volnode,ndims,ifull,ium2,ipoint,
     $        volumes,cij,dum3,dum4)
         if(lmyidhead)write(*,*)'Finished volume setting'
c If head, write the geometry data if we've had to calculate it.
         if(lmyidhead)call stored3geometry(volumes,iuds,ifull,istat
     $        ,.true.)
      endif
c---------------------------------------------
c Set an object pointer for all the edges so their regions get
c set by the iregioninit call
      call iregioninit(ndims,ifull)
      if(lmyidhead)call reportfieldmask()
      ipoint=0
      call mditerarg(cijedge,ndims,ifull,iuds,ipoint,cij,dum2,dum3,dum4)
c Initialize the region flags in the object data
c This old position overruled the new edge setting.
c      call iregioninit(ndims,ifull)
c---------------------------------------------
      if(.not.lmyidhead)then
c Don't do plotting from any node except the master.
         ltestplot=.false.
         lcijplot=.false.
         lsliceplot=.false.
         lorbitplot=.false.
         norbits=0
         iobpl=0
      else
c---------------------------------------------
c Some simple graphics of cij, and volumes.
         if(ltestplot)call text3graphs(ndims,iuds,ifull,cij,volumes)
c More elaborate graphics of volumes
         if(ltestplot)call sliceGweb(ifull,iuds,volumes,na_m,zp,
     $              ixnp,xn,ifix,'volumes:'//'!Ay!@'//char(0),dum,dum)
c---------------------------------------------
c The following requires include objcom.f
         if(lmyidhead)write(*,*)
     $      'Used No of pointers:',oi_cij,' of',iuds(1)*iuds(2)*iuds(3)
     $        ,' points.'
c Plot objects 0,1 and 2 (bits)
         if(iobpl.ne.0.and.lmyidhead)then
            if(rcij.eq.0.)rcij=rs
c            call cijplot(ndims,ifull,iuds,cij,rcij,iobpl,0)
            call cijplot(ndims,ifull,iuds,cij,rcij,iobpl)
          endif
      endif
c---------------------------------------------
c Initialize charge (set q to zero over entire array).
      call mditerset(q,ndims,ifull,iuds,0,0.)
c Initialize potential (set u to zero over entire array).
      call mditerset(u,ndims,ifull,iuds,0,0.)
c Initialize diagsum if necessary.
      do idiag=1,ndiags
         call mditerset(diagsum(1,1,1,idiag),ndims,ifull,iuds,0,0.)
      enddo
c Initialize additional potential and charge if needed.
      if(iptch_mask.ne.0 .or. gtt.ne.0.)
     $     call setadfield(ndims,ifull,iuds,iptch_mask,lsliceplot)
c      if(myid.eq.0)call sliceGweb(ifull,iuds,rhoci,na_m,zp,
c     $              ixnp,xn,ifix,'rhoci',dum,dum)

c---------------------------------------------------------------     
c Control. Bit 1, use my sor params (not here). Bit 2 use faddu (not)
c Bit 4-6 periodicity.
      ictl=0
c Make dimensions periodic:
      do id=1,ndims
         if(LPF(id))ictl=ictl+4*2**id
      enddo
c      write(*,*)'Calling sormpi'
c An initial solver call with zero density.
      ierrsor=0
      if(debyelen.ne.0.)call sormpi(ndims,ifull,iuds,cij,u,q,bdyshare
     $     ,bdyset,faddu,ictl,ierrsor,myid,idims)
      ictl=2+ictl
c      write(*,*)'Return from initial sormpi call.'
      if(ltestplot)call sliceGweb(ifull,iuds,u,na_m,zp,
     $              ixnp,xn,ifix,'potential:'//'!Ay!@'//char(0),dum,dum)
c
c-------------------------------------------------------------------
      if(lmyidhead)then
         call vaccheck(ifull,iuds,cij,u,thetain,nth,rs,ltestplot)
      endif
c End of plotting.
c------------------------------------------------------------------
c Set phip from the first object if it makes sense.
      call phipset(myid)
c------------------------------------------------------------------
c (Re)Initialize the fortran random number generator.
      idum=-myid-1
      rdum=ran1(idum)
c Initialize with a specified number of particles.
c      write(*,*)'ibool_part=',ibool_part
c      if(lmyidhead)write(*,*)'Initializing',n_part,' particles'
      call pinit(subcycle)
c      if(lmyidhead)write(*,*)'Return from pinit'
c---------------------------------------------
c Initialize the force tracking.
      call forcetrackinit()
c---------------------------------------------
      phirein=0.
      ninjcomp0=ninjcomp
      if(ninjcomp.ne.0.and.lmyidhead)
     $     write(*,*)'Fixed injection count:',ninjcomp
      maccel=nsteps/3
      dtf=dt
c-----------------------------------------------
c Restart code
      if(lrestart.ne.0)then
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
         else
            iferr=1
         endif
c         write(*,*)'iferr=',iferr
         if(lrestart-4*(lrestart/4).ne.0)then
            call partread(partfilename,ierr)
            if(ierr-4*(ierr/4).eq.0)then 
c We succeeded in reading the part-file. Relocate the particles.
               write(*,'(a,i4,a,a,i3)')' cpu',myid
     $              ,' Restart file read: '
     $              ,partfilename(1:lentrim(partfilename)+1),lrestart
               call locateinit()
               ied=1
c Only read the phi-file if the flux file was present. Full restart.
               if(iferr.eq.0)then
c                  write(*,*)'Reading phifile',ierr,phifilename
                  call array3read(phifilename,ifull,iuds,ied,u,ierr)
               endif
            endif
         endif
         if(nsteps+nf_step.gt.nf_maxsteps)then
            if(lmyidhead)write(*,*)'Asked for',
     $           nsteps,' in addition to',nf_step,
     $           ' Total',nsteps+nf_step,
     $           ' too much; set to',nf_maxsteps-1
            nsteps=nf_maxsteps-nsteps-1
         endif
c In case we have overwritten phip with the value from the restart file,
c try to set it again from the first object. But tell it we are not
c the head node, so it does not give out messages.
         call phipset(1)
c         if(lmyidhead)write(*,*)'nrein,n_part,ioc_part,rhoinf,dt=',
c     $        nrein,n_part,ioc_part,rhoinf,dt,(dtprec(k),k=1,4)
      endif
c-----------------------------------------------
      if(lmyidhead)then
         if(colntime.ne.0)then
            write(*,'(a,f8.2,a,f8.4)')' Collision time='
     $        ,colntime,', Eneutral=',Eneutral
         endif
         write(*,'(/,a)')'Step Iterations Flux:'
      endif
c Main step iteration -------------------------------------
      do j=1,nsteps
         nf_step=nf_step+1
c Acceleration code.
         bdtnow=max(1.,(bdt-1.)*(maccel-j+2)/(maccel+1.)+1.)
         dt=bdtnow*dtf
         ninjcomp=int(bdtnow*ninjcomp0)
         if(ninjcomp.ne.0)nrein=ninjcomp

         call mditerset(psum,ndims,ifull,iuds,0,0.)
c         write(*,*)'chargetomesh calling, ndiags',ndiags
         call chargetomesh(psum,ndims,iLs,diagsum,ndiags)
         call diagperiod(psum,ndims,ifull,iuds,iLs,1)
c The following call is marginally more efficient than diagperiod call.
c         call psumperiod(psum,ndims,ifull,iuds,iLs)
c But it is not necessary so simplify to one period routine.
c Psumreduce takes care of the reductions that were in rhoinfcalc 
c and explicit psum. It encapsulates the iaddtype iaddop generation.
c Because psumtoq internally compensates for faddu, we reduce here
         call psumreduce(psum,ndims,ifull,iuds,iLs)
c Calculate rhoinfinity, needed in psumtoq. Dependent on reinjection type.
         call rhoinfcalc(dt)
c Convert psums to charge density, q. Remember external psumtoq!
         call mditerarg(psumtoq,ndims,ifull,ium2,
     $        0,psum(2,2,2),q(2,2,2),volumes(2,2,2),u(2,2,2))
         istepave=min(nf_step,iavesteps)

         if(debyelen.eq.0)then
            call mditerarg(quasineutral,ndims,ifull,ium2,
     $        0,q(2,2,2),u(2,2,2),volumes(2,2,2),dum4)
            call bdyslope0(ndims,ifull,iuds,cij,u,q)
         else
            call sormpi(ndims,ifull,iuds,cij,u,q,bdyshare,bdyset,faddu
     $           ,ictl,ierr,myid,idims)
         endif
         call calculateforces(ndims,iLs,cij,u)

         if(lmyidhead)write(*,'(i4.4,i4,$)')nf_step,ierr
         if(lsliceplot)then
            if(ipstep.eq.0.or.mod(j,ipstep).eq.0)then
               if(ldenplot)call sliceGweb(ifull,iuds,q,na_m,zp,
     $              ixnp,xn,ifix,'density: n'//char(0),dum,dum)
               if(lphiplot)call sliceGweb(ifull,iuds,u,na_m,zp,
     $              ixnp,xn,ifix,'potential:'//'!Ay!@'//char(0),dum,dum)
            endif
         endif

         if(nf_step.eq.ickst) then
c      write(*,*)'Checking Step',nf_step
            call checkuqcij(ifull,u,q,psum,volumes,cij)
            call padvnc(ndims,iLs,cij,u)
            call checkx
         else
c The normal call:
            call padvnc(ndims,iLs,cij,u)
         endif
         call fluxreduce()
c Now do cij update
         call cijdirect(ndims,debyelen,error)
c Store the step's rhoinf, dt, npart.
         ff_rho(nf_step)=rhoinf
         ff_dt(nf_step)=dt
         nf_npart(nf_step)=n_part
         if(lmyidhead)then
c write out flux to object 1.
            write(*,'(f6.3,''| '',$)')fluxdiag()
c            write(*,'(f6.3,''| '',$)')n_part/(rhoinf*1000)
c            write(*,'(i6,''| '',$)')n_part
            if(mod(nf_step,5).eq.0)write(*,*)
            if(mod(nf_step,(nsteps/25+1)*5).eq.0)then
               write(*,
     $    '(''nrein,n_part,ioc_part,rhoinf,dt='',i6,i9,i9,2f10.3)')
     $              nrein,n_part,ioc_part,rhoinf,dt
               if(nsubc.ne.0)write(*,'(''Subcycled:'',i5,$)')nsubc
               if(ndropped.ne.0)then
c Report dropped ions because of excessive acceleration.
                  write(*,'(a,i5,a,f8.3)'
     $             )' dropped-ion period-total:',ndropped
     $                 ,'  per step average:'
     $                 ,float(ndropped/((nsteps/25+1)*5))
                  ndropped=0
               else
                  write(*,*)
               endif
            endif
         endif

c These running and box averages do not include the updates for this step.
c Accumulate running q and u averages:
         call average3d(q,qave,ifull,iuds,istepave)
         call average3d(u,uave,ifull,iuds,istepave)
c Every iavesteps, calculate the box average of the moments, and write it
c out, if we are doing diagnostics.
         if(mod(j,iavesteps).eq.0)then
            if(ndiags.gt.0)then
c Reduce the data
               call diagreduce(diagsum,ndims,ifull,iuds,iLs,ndiags)
               call diagperiod(diagsum,ndims,ifull,iuds,iLs,ndiags)
c Do any other processing? Here or later?
c               call diagstep(iLs,diagsum,ndiags)
               if(lmyidhead)then
c If I'm the head, write it.
                  write(*,'(a,i3,a)')'Diags',ndiags,' '
c                  write(argument,'(''.dia'',i4.4)')j
                  write(argument,'(''.dia'',i4.4)')nf_step
                  diagfilename=' '
                  call namewrite(diagfilename,ifull,iuds,ndiags,diagsum
     $                 ,argument)
               endif
c Now reinit diagsum
               do idiag=1,ndiags
                  call mditerset(diagsum(1,1,1,idiag),ndims,ifull,iuds,0
     $                 ,0.)
               enddo
            endif
            if(idistp.ne.0)then
c Particle distribution diagnostics
               if(idcount.ne.0)then
c Not for the first time, print or plot.
c Reduce the data from nodes.
                  call ptdiagreduce()
                  if(lmyidhead)then 
                     if(2*(idistp/2)-4*(idistp/4).ne.0)
     $                    call pltsubdist(5,9,9,vlimit,xnewlim,cellvol)
                     diagfilename=' '
                     call nameconstruct(diagfilename)
                     write(diagfilename(lentrim(diagfilename)+1:)
     $                    ,'(''.pex'',i4.4)')nf_step
                     if(idistp-2*(idistp/2).ne.0)
     $                    call distwrite(xlimit,vlimit,xnewlim,
     $                    diagfilename,cellvol)
                  endif
c (Re)initialize the accumulation
                  ibset=1
                  call fvxinit(xnewlim,cellvol,ibset)
                  call partacinit(vlimit)
c Unless we want fsv to accumulate all the particle counts, it also
c ought to be reinitialized here. (partexamine expects this.)
                  call fsvzero()
               endif
               idcount=idcount+1
            endif
         endif
c Particle distribution diagnostics.
c The serial cost for this call with 1M particles is about 1s in 17s.
c Or less than 2s if both partaccum and vaccum are called. Thus the
c total cost is roughly 10% of particle costs. Tell it that it must
c set dimensions the first time by appending 0.
         if(idcount.gt.0)call partdistup(xlimit,vlimit,xnewlim,
     $        cellvol,myid,0)
         if(iwstep.gt.0)then
            if(mod(nf_step,iwstep).eq.0)call datawrite(myid
     $           ,partfilename,restartpath,ifull,iuds,u,uave,qave)
         endif

c This non-standard fortran call works with gfortran and g77 to flush stdout.
c Comment it out if it causes problems. (E.g. pathscale gives segfaults.)
c Pathscale demands an argument number. So give it explicitly. Should fix
c the segfaults.
         if(lmyidhead)call flush(6)
      enddo
c-------- End of Main Step Iteration -------------------------------
c      write(*,*)iorbitlen(1),(xorbit(k,1),k=1,10)
      if(norbits.ne.0)
     $     call cijplot(ndims,ifull,iuds,cij,rs,iobpl)
c     $     call cijplot(ndims,ifull,iuds,cij,rs,iobpl,norbits)

      if(lorbitplot.and.norbits.ne.0)call orbitplot(ifull,iuds,u,phip,rc
     $     ,rs)

c Everyone writes what they have to.
      if(iwstep.gt.0 .or. myid.eq.0)call datawrite(myid,partfilename
     $     ,restartpath,ifull,iuds,u,uave,qave)
      
c-------------------------------------------------------------------
      call mpifinalize(ierr)
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
