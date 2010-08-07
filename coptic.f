      program coptic
c Main program of cartesian coordinate, oblique boundary, pic code.

      include 'mpif.h'
      include 'objcom.f'
c Storage array spatial count size
      include 'griddecl.f'
c coptic runs correctly with unequal dimensions but phiexamine does not.
      parameter (Li1=na_i,Li2=Li1*na_j,Li3=Li2*na_k)
      real u(na_i,na_j,na_k),q(na_i,na_j,na_k)
     $     ,cij(2*ndims_sor+1,na_i,na_j,na_k)
c Running averages.
      real qave(na_i,na_j,na_k),uave(na_i,na_j,na_k)
c
      real psum(na_i,na_j,na_k),volumes(na_i,na_j,na_k)
c Used dimensions, Full dimensions. Used dims-2
      integer iuds(ndims_sor),ifull(ndims_sor),ium2(ndims_sor)
c Mesh spacing description structure
      include 'meshcom.f'
c Processor cartesian geometry can be set by default.
      integer nblksi,nblksj,nblksk
      parameter (nblksi=1,nblksj=1,nblksk=1)
      integer idims(ndims_sor)
c mpi process information.
      include 'myidcom.f'
c Common data containing the BC-object geometric information
      include '3dcom.f'
c Structure vector needed for finding adjacent u values.
c Don't use the mditerate common. It might not be right.
      integer iLs(ndims_sor+1)
c Particle common data
      include 'partcom.f'
c Plasma common data
      include 'plascom.f'
c Collision common data
      include 'colncom.f'
c Point charge common data
      include 'ptchcom.f'
c      integer ndims
c      parameter (ndims=ndims_sor)
      external bdyset,faddu,cijroutine,cijedge,psumtoq
      external volnode,linregion
      character*100 partfilename,phifilename,fluxfilename,objfilename
      character*100 argument
c      common /ctl_sor/mi_sor,xjac_sor,eps_sor,del_sor,k_sor
      logical ltestplot,lcijplot,lsliceplot,lorbitplot,linjplot
      logical lrestart,lmyidhead,lphiplot,ldenplot
      integer ipstep
c Diagnostics
      real zp(na_m,na_m,ndims_mesh)
c Set up the structure vector.
      data iLs/1,Li1,Li2,Li3/
c Mesh and mpi parameter defaults:
      data idims/nblksi,nblksj,nblksk/
      data ifull/na_i,na_j,na_k/
c No longer set here. Done by meshconstruct. 
c      data iuds/ni,nj,nk/
c Data for plotting etc.
      data iobpl/0/
      data ltestplot,lcijplot,lsliceplot,lorbitplot,linjplot/
     $     .false.,.false.,.false.,.false.,.false./
      data lphiplot,ldenplot/.true.,.true./
c      data thetain,nth/.1,1/
      data lrestart/.false./
      data ipstep/1/
c-------------------------------------------------------------
c Consistency checks
      if(ndims.ne.ndims_sor)then
         write(*,*)'Inconsistent ndims, ndims_sor',ndims,ndims_sor
         stop
      endif
c-------------------------------------------------------------
c Initialize the fortran random number generator with a fixed number
c for solutions of volumes etc. Each node then does the same.
      rs=ran1(-1)
c Defaults:
c Determine what reinjection scheme we use. Sets rjscheme.
      include 'REINJECT.f'
c Fixed number of particles rather than fixed injections.
      ninjcomp=0
      n_part=0
c Default to constant ripernode not n_part.
      ripernode=100.
c Default edge-potential (chi) relaxation rate.     
      crelax=1.*Ti/(1.+Ti)
      averein=0.
      dt=.1
      objfilename='ccpicgeom.dat'
      nsteps=5
      debyelen=1.
      Ti=1.
      vd=0.
      subcycle=0.
      colntime=0.
      vneutral=0.
      rs=5.0
      rsmesh=rs
      numprocs=1
      bdt=1.
      thetain=.1
      nth=1
      norbits=0
      ickst=0
      iavesteps=100
      lmyidhead=.true.
      ifplot=-1
c---------------------------------------------------------------------
c This necessary here so one knows early the mpi structure.
c Otherwise could have been hidden in sormpi and pass back numprocs.
      call mpigetmyid(myid,nprocs,ierr)
      if(myid.ne.0) lmyidhead=.false.
      numprocs=nprocs
c--------------------------------------------------------------
c Deal with arguments
c      if(iargc().eq.0) goto "help"
      do i=1,iargc()
         call getarg(i,argument)
         if(argument(1:3).eq.'-gt')ltestplot=.true.
         if(argument(1:3).eq.'-gc')read(argument(4:),*,end=201)iobpl
         if(argument(1:3).eq.'-gs')then
            lsliceplot=.true.
            read(argument(4:),*,err=210,end=210)ipstep
            goto 211
 210        ipstep=1
 211        continue
            write(*,*)'ipstep=',ipstep
         endif
         if(argument(1:3).eq.'-gd')ldenplot=.false.
         if(argument(1:3).eq.'-gp')lphiplot=.false.
         if(argument(1:3).eq.'-gi')linjplot=.true.
         if(argument(1:3).eq.'-gf')read(argument(4:),*,err=201)ifplot
         if(argument(1:3).eq.'-go')read(argument(4:),*,err=201)norbits
         if(argument(1:3).eq.'-at')then
            read(argument(4:),*,err=201)thetain
         elseif(argument(1:3).eq.'-an')then
            read(argument(4:),*,err=201)nth
         elseif(argument(1:2).eq.'-a')then 
            read(argument(3:),*,err=201)iavesteps
         endif
         if(argument(1:3).eq.'-ni')read(argument(4:),*,err=201)n_part
         if(argument(1:3).eq.'-pn')read(argument(4:),*,err=201)numprocs
         if(argument(1:3).eq.'-ri')read(argument(4:),*,err=201)ripernode
         if(argument(1:3).eq.'-rx')read(argument(4:),*,err=201)crelax
         if(argument(1:3).eq.'-ck')read(argument(4:),*,err=201)ickst
         if(argument(1:3).eq.'-ct')read(argument(4:),*,err=201)colntime

         if(argument(1:3).eq.'-dt')read(argument(4:),*,err=201)dt
         if(argument(1:3).eq.'-da')read(argument(4:),*,err=201)bdt
         if(argument(1:3).eq.'-ds')read(argument(4:),*,err=201)subcycle
         if(argument(1:7).eq.'--reinj')
     $        read(argument(8:),*,err=201)ninjcomp
         if(argument(1:2).eq.'-s')then
            read(argument(3:),*,err=201)nsteps
            if(nsteps.gt.nf_maxsteps)then
               if(lmyidhead)write(*,*)'Asked for more steps',nsteps,
     $              ' than allowed. Limit ',nf_maxsteps-1
               nsteps=nf_maxsteps-1
            endif
         endif
         if(argument(1:3).eq.'-vn')then
            read(argument(4:),*,err=201)vneutral
         elseif(argument(1:2).eq.'-v')then
            read(argument(3:),*,err=201)vd
         endif
         if(argument(1:2).eq.'-l')read(argument(3:),*,err=201)debyelen
         if(argument(1:2).eq.'-t')read(argument(3:),*,err=201)Ti
         if(argument(1:9).eq.'--restart')lrestart=.true.
         if(argument(1:10).eq.'--extfield')then
            read(argument(11:),*,err=201)extfield
c            write(*,*)'||||||||||||||extfield',extfield
            lextfield=.true.
         endif
         if(argument(1:9).eq.'--objfile')
     $        read(argument(10:),'(a)',err=201)objfilename
         if(argument(1:3).eq.'-ho')then
            call geomdocument()
            call exit(0)
         endif
         if(argument(1:2).eq.'-h')goto 203
         if(argument(1:2).eq.'-?')goto 203
      enddo
      goto 202
c------------------------------------------------------------
c Help text
 201  continue
      if(lmyidhead)write(*,*)'=====Error reading command line argument '
     $     ,argument(:20)
 203  continue
      if(myid.ne.0)goto 202
 301  format(a,i5)
 302  format(a,f8.3)
      write(*,301)'Usage: ccpic [switches]'
      write(*,301)'Parameter switches.'
     $     //' Leave no gap before value. Defaults indicated [ddd'
      write(*,301)' -ni   set No of particles/node; zero => unset.    ['
     $     ,n_part
      write(*,301)' --reinj    set reinjection number at each step.   ['
     $     ,ninjcomp
      write(*,302)' -ri   set rhoinfinity/node => reinjection number. ['
     $     ,ripernode
      write(*,302)' -rx   set Edge-potl relax rate: 0=>off, 1=>immed. ['
     $     ,crelax
      write(*,302)' -dt   set Timestep.              [',dt,
     $     ' -da   set Initial dt accel-factor[',bdt
      write(*,302)' -ds   set subcycle fraction.     [',subcycle
      write(*,301)' -s    set No of steps.           [',nsteps
      write(*,302)' -v    set Drift velocity.        [',vd
      write(*,302)' -t    set Ion Temperature.       [',Ti
      write(*,302)' -l    set Debye Length.          [',debyelen
      write(*,301)' -a    set averaging steps.       [',iavesteps
      write(*,301)' -ct   set collision time.        [',colntime
      write(*,301)' -vn   set neutral drift velocity [',vneutral
c      write(*,301)' -xs<3reals>, -xe<3reals>  Set mesh start/end.'
      write(*,301)' --objfile<filename>  set name of object data file.'
     $     //' [ccpicgeom.dat'
      write(*,301)' --restart  Attempt to restart from saved state.'
      write(*,301)'Debugging switches for testing'
      write(*,301)' -gt   Plot regions and solution tests.'
      write(*,301)' -gi   Plot injection accumulated diagnostics.'
      write(*,301)' -gs[] Plot slices of solution potential, density. '
     $     //'[At step n]. [',ipstep
      write(*,301)' -gd -gp Turn off slicing of density, potential. '
      write(*,301)' -gf   set quantity plotted for flux evolution and'//
     $     ' final distribution. [',ifplot
      write(*,301)' -gc   set wireframe [& stencils(-)] mask.'//
     $     ' objects<->bits. [',iobpl
      write(*,301)' -go   set No of orbits'
     $     //'(to plot on objects set by -gc). [',norbits
      write(*,301)' -at   set test angle.'
     $     //' -an   set No of angles. '
      write(*,301)' -ck   set checking timestep No. [',ickst
      write(*,301)' -h -?   Print usage.'
      write(*,301)' -ho     Print geomobj file format description'
      call exit(0)
 202  continue
c-----------------------------------------------------------------
c Finalize parameters after switch reading.
c Geometry and boundary information. Read in.
      call readgeom(objfilename,myid)
c---------------------------------------------------------------
c Construct the mesh vector(s) and ium2
      call meshconstruct(ndims,iuds)
      if(lmyidhead)write(*,'(a,3i4,6f8.3)')
     $     ' Constructed mesh',iuds
     $     ,(xmeshstart(k),xmeshend(k),k=1,ndims)
c Initialize geometry if needed for particular case.
      call geominit(myid)
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
c Initializations
      if(lmyidhead)write(*,*)'Initializing the stencil data cij'
c Initialize cij:
      ipoint=iLs(1)+iLs(2)+iLs(3)
c Seems to work even though the cijroutine is
c not defined with sufficient arguments for all those in this call.
      call mditerarg(cijroutine,ndims,ifull,ium2,ipoint,
     $     cij(1,1,1,1),debyelen,dum3,dum4)

c---------------------------------------------
c Here we try to read the stored geometry volume data.
c This is probably a bottleneck for multi-process.
      istat=1
      call stored3geometry(volumes,iuds,ifull,istat)
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
         if(lmyidhead)call stored3geometry(volumes,iuds,ifull,istat)
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
      if(myid.ne.0)then
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
c---------------------------------------------
c The following requires include objcom.f
c         if(lmyidhead)write(*,*)'Finished mesh/stencil setup:',iuds
         if(lmyidhead)write(*,*)
     $      'Used No of pointers:',oi_sor,' of',iuds(1)*iuds(2)*iuds(3)
     $        ,' points.'
c Plot objects 0,1 and 2 (bits)
c      iobpl=-7
         if(iobpl.ne.0.and.lmyidhead)
     $        call cijplot(ndims,ifull,iuds,cij,rs,iobpl,0)
      endif

c Initialize charge (set q to zero over entire array).
      call mditerset(q,ndims,ifull,iuds,0,0.)
c Initialize potential (set u to zero over entire array).
      call mditerset(u,ndims,ifull,iuds,0,0.)
c Initialize additional potential and charge if needed.
      if(iptch_mask.ne.0)
     $     call setadfield(ndims,ifull,iuds,iptch_mask,lsliceplot)

c---------------------------------------------------------------     
c An inital vacuum solution with zero density. 
c Control. Bit 1, use my sor params (not here). Bit 2 use faddu (not)
      ictl=0
c      write(*,*)'Calling sormpi, ni,nj=',ni,nj
c An initial solver call.
      call sormpi(ndims,ifull,iuds,cij,u,q,bdyset,faddu,ictl,ierr
     $     ,myid,idims)
      ictl=2
c      write(*,*)'Return from initial sormpi call.'
c
c-------------------------------------------------------------------
      if(lmyidhead)then
         call vaccheck(ifull,iuds,cij,u,thetain,nth,rs,ltestplot)
      endif
c End of plotting.
c------------------------------------------------------------------
c Set phip from the first object if it makes sense.
      if(obj_geom(oabc,1).ne.0)then
         phip=-obj_geom(oabc+2,1)/obj_geom(oabc,1)
      elseif(obj_geom(oradius,1).ne.0.)then
         phip=obj_geom(omag,1)*obj_geom(oradius,1)
         if(myid.eq.0)write(*,*)'Potential from point charge'
     $        ,obj_geom(omag,1),' at radius ',obj_geom(oradius,1)
     $        ,' Charge:',phip
      else
         write(*,*)'Potential phip not set from objects.'
         phip=0.
      endif
c------------------------------------------------------------------
c Initialize with a specified number of particles.
c (Re)Initialize the fortran random number generator.
      idum=-myid-1
      rdum=ran1(idum)
c      if(lmyidhead)write(*,*)'Initializing',n_part,' particles'
      call pinit(subcycle)
c      if(lmyidhead)write(*,*)'Return from pinit'
c---------------------------------------------
c Initialize the fluxdata storage and addressing.
      call fluxdatainit(myid)
c Initialize the force tracking.
      call forcetrackinit()
      write(*,*)'mf_obj=',mf_obj
c---------------------------------------------
      phirein=0.
      ninjcomp0=ninjcomp
      if(ninjcomp.ne.0.and.lmyidhead)
     $     write(*,*)'Fixed injection count:',ninjcomp
      maccel=nsteps/3
      dtf=dt
c-----------------------------------------------
c Restart code
c 400  continue
      if(lrestart)then
         partfilename=' '
         call nameconstruct(partfilename)
         phifilename=partfilename
         nb=nbcat(phifilename,'.phi')
         fluxfilename=partfilename
         nb=nbcat(fluxfilename,'.flx')
         nb=nameappendint(partfilename,'.',myid,3)
         call readfluxfile(fluxfilename,ierr)
         if(ierr.ne.0)goto 401
         call partread(partfilename,ierr)
         if(ierr.ne.0)goto 401
c         call phiread(phifilename,ifull,iuds,u,ierr)
         call array3read(phifilename,ifull,iuds,u,ierr)
         if(ierr.ne.0)goto 401
         write(*,*)'Restart files read successfully.'
         if(nsteps+nf_step.gt.nf_maxsteps)then
            if(lmyidhead)write(*,*)'Asked for',
     $           nsteps,' in addition to',nf_step,
     $           ' Total',nsteps+nf_step,
     $           ' too much; set to',nf_maxsteps-1
            nsteps=nf_maxsteps-nsteps-1
         endif
c         if(lmyidhead)write(*,*)'nrein,n_part,ioc_part,rhoinf,dt=',
c     $        nrein,n_part,ioc_part,rhoinf,dt
         goto 402
 401     continue
         write(*,*)'Failed to read restart files',
     $        fluxfilename(1:lentrim(fluxfilename)-4)
         lrestart=.false.
 402     continue
      endif
c-----------------------------------------------
      if(lmyidhead)write(*,*)'Step Iterations Flux:'
c Main step iteration -------------------------------------
      do j=1,nsteps
         nf_step=nf_step+1
c Acceleration code.
         bdtnow=max(1.,(bdt-1.)*(maccel-j+2)/(maccel+1.)+1.)
         dt=bdtnow*dtf
         ninjcomp=int(bdtnow*ninjcomp0)
         if(ninjcomp.ne.0)nrein=ninjcomp

         call mditerset(psum,ndims,ifull,iuds,0,0.)
         call chargetomesh(psum,iLs,diags)
c Psumreduce takes care of the reductions that were in rhoinfcalc 
c and explicit psum. It encapsulates the iaddtype iaddop generation.
c Because psumtoq internally compensates for faddu, we reduce here
         call psumreduce(psum,ndims,ifull,iuds,iLs) 
c Calculate rhoinfinity, needed in psumtoq. Dependent on reinjection type.
         call rhoinfcalc(dt)
c Convert psums to charge, q. Remember external psumtoq!
         call mditerarg(psumtoq,ndims,ifull,ium2,
     $        0,psum(2,2,2),q(2,2,2),volumes(2,2,2),u(2,2,2))

         call sormpi(ndims,ifull,iuds,cij,u,q,bdyset,faddu,ictl,ierr
     $     ,myid,idims)

         call calculateforces(ndims,iLs,cij,u)

         if(lmyidhead)write(*,'(i4.4,i4,$)')nf_step,ierr
         if(lsliceplot)then
            if(ipstep.eq.0.or.mod(j,ipstep).eq.0)then
               if(ldenplot)call sliceGweb(ifull,iuds,q,na_m,zp,
     $              ixnp,xn,ifix,'density: n')
               if(lphiplot)call sliceGweb(ifull,iuds,u,na_m,zp,
     $              ixnp,xn,ifix,'potential:'//'!Ay!@')
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
c Store the step's rhoinf, dt.
         ff_rho(j)=rhoinf
         ff_dt(j)=dt
         if(lmyidhead)then
            call fluxdiag()
            if(mod(nf_step,5).eq.0)write(*,*)
         endif
         if(lmyidhead.and.mod(nf_step,(nsteps/25+1)*5).eq.0)then
            write(*,
     $    '(''nrein,n_part,ioc_part,rhoinf,dt='',i5,i9,i9,2f10.3)')
     $        nrein,n_part,ioc_part,rhoinf,dt
            if(nsubc.ne.0)write(*,'(''Subcycled:'',i6)')nsubc
         endif
         istepave=min(nf_step,iavesteps)
         call average3d(q,qave,ifull,iuds,istepave)
         call average3d(u,uave,ifull,iuds,istepave)

c This non-standard fortran call works with gfortran and g77 to flush stdout.
         if(lmyidhead)call flush()
c Comment it out if it causes problems.
      enddo
c End of Main Step Iteration -------------------------------
c      write(*,*)iorbitlen(1),(xorbit(k,1),k=1,10)
c      if(lorbitplot)call orbit3plot(ifull,iuds,u,phip,rc,rs)
      if(norbits.ne.0)
     $     call cijplot(ndims,ifull,iuds,cij,rs,iobpl,norbits)
c      write(*,*)'Finished orbitplot.'

      call partwrite(partfilename,myid)
      if(lmyidhead)then
         if(iptch_mask.ne.0)then
            call mditeradd(u,ndims,ifull,iuds,0,uci)
            call mditeradd(uave,ndims,ifull,iuds,0,uci)
         endif
         call namewrite(phifilename,ifull,iuds,uave,'.pha')
         call namewrite(phifilename,ifull,iuds,qave,'.den')
         call namewrite(phifilename,ifull,iuds,u,'.phi')
      endif
      
c-------------------------------------------------------------------
      call mpifinalize(ierr)
c Check some flux diagnostics and writing.
      if(lmyidhead)then 
         call writefluxfile(fluxfilename)
         if(linjplot)call plotinject(Ti)
         do ifobj=1,mf_obj
            call fluxave(nsteps/2,nsteps,ifobj,ifplot,rinf)
         enddo
      endif
c      call readfluxfile(fluxfilename)
      end
c**********************************************************************
