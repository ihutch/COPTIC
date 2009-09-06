      program ccpic
c Main program of cartesian coordinate pic code.

      include 'mpif.h'
      include 'objcom.f' 
c Storage array spatial count size
      integer Li,ni,nj,nk
c      parameter (Li=100,ni=40,nj=40,nk=20)
c      parameter (Li=100,ni=60,nj=60,nk=60)
      parameter (Li=100,ni=16,nj=16,nk=16)
c      parameter (Li=100,ni=32,nj=32,nk=32)
c      parameter (Li=100,ni=64,nj=64,nk=64)
c      parameter (Li=130,ni=128,nj=128,nk=128)
c      parameter (Li=6,ni=6,nj=6,nk=6)
c      parameter (Li=100,ni=26,nj=16,nk=20)
      parameter (Li2=Li*Li,Li3=Li2*Li)
      real u(Li,Li,Li),q(Li,Li,Li),cij(2*ndims_sor+1,Li,Li,Li)

c Additional variables for testing. Eventually should be removed.
      real u2(Li,Li,Li),q2(Li,Li,Li),cij2(2*ndims_sor+1,Li,Li,Li)
      real psum2(Li,Li,Li),volumes2(Li,Li,Li)
      real x_part2(9,1000000)
      integer n_part2,if_part2(1000000),iregion_part2,ioc_part2
      integer nrein2,numprocs2,ninjcomp2
      real dt2,rhoinf2,phirein2
      logical ldiags2
c End of extra variables.
c
      real psum(Li,Li,Li),volumes(Li,Li,Li)
c Used dimensions, Full dimensions. Used dims-2
      integer iuds(ndims_sor),ifull(ndims_sor),ium2(ndims_sor)
c Arrays defining mesh ends
      real xmstart(ndims_sor),xmend(ndims_sor)
c Mesh spacing description structure
      include 'meshcom.f'
c Processor cartesian geometry
      integer nblksi,nblksj,nblksk
      parameter (nblksi=2,nblksj=1,nblksk=1)
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
      integer ndims
      parameter (ndims=ndims_sor)

      external bdyset,faddu,cijroutine,cijedge,psumtoq
      external volnode,linregion
      character*100 partfilename
      character*100 phifilename
      character*100 fluxfilename
      character*100 objfilename
      character*100 argument
      common /ctl_sor/mi_sor,xjac_sor,eps_sor,del_sor,k_sor
      logical ltestplot,lcijplot,lsliceplot,lorbitplot,linjplot
      logical lrestart,lmyidhead,lfplot,lregion
      
c Diagnostics
c      real usave(Li,Li,Li),error(Li,Li,Li),cijp(2*ndims_sor+1,Li,Li)
      real zp(Li,Li)
c Point in the active region
      real xir(ndims)

c Set up the structure vector.
c      data iLs/1,Li,(Li*Li),(Li*Li*Li)/
      data iLs/1,Li,Li2,Li3/
c Mesh and mpi parameter defaults:
      data idims/nblksi,nblksj,nblksk/
      data ifull/Li,Li,Li/
      data iuds/ni,nj,nk/
c Point which lies in the plasma region:
      data xir/2.,2.,2./
c Data for plotting etc.
      data iobpl/0/
      data ltestplot,lcijplot,lsliceplot,lorbitplot,linjplot,lfplot/
     $     .false.,.false.,.false.,.false.,.false.,.false./
      data thetain,nth/.1,1/
      data lrestart/.false./

c-------------------------------------------------------------
c Initialize the fortran random number generator with a fixed number
c for solutions of volumes etc. Each node does the same.
      rs=ran1(-1)
c Defaults:
c Fixed number of particles rather than fixed injections.
      ninjcomp=0
      n_part=0
c Default to constant rhoinf not n_part.
      rhoinf=100.
      dt=.1
      objfilename='ccpicgeom.dat'
      nsteps=5
      debyelen=1.
      Ti=1.
      vd=0.
      rs=5.0
      rsmesh=rs*1.00001
      numprocs=1
      bdt=1.
      norbits=0
      ickst=0
      lmyidhead=.true.
      do id=1,ndims
         xmstart(id)=-rsmesh
         xmend(id)=rsmesh
      enddo
c---------------------------------------------------------------------
c mpi initialization only.
c Control bit 3 (=4) pure initialization, no communication (to get myid).
      call sormpi(ndims,ifull,iuds,cij,u,q,bdyset,faddu,4,ierr
     $     ,myid,idims)
      if(myid.ne.0) lmyidhead=.false.
c This necessary or could be hidden in sormpi and pass back numprocs.
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
c--------------------------------------------------------------
c Deal with arguments
c      if(iargc().eq.0) goto "help"
      do i=1,iargc()
         call getarg(i,argument)
         if(argument(1:3).eq.'-gt')ltestplot=.true.
         if(argument(1:3).eq.'-gc')read(argument(4:),*,end=201)iobpl
         if(argument(1:3).eq.'-gs')lsliceplot=.true.
         if(argument(1:3).eq.'-gi')linjplot=.true.
         if(argument(1:3).eq.'-gf')lfplot=.true.
         if(argument(1:3).eq.'-go')read(argument(4:),*,err=201)norbits
         if(argument(1:3).eq.'-at')read(argument(4:),*,err=201)thetain
         if(argument(1:3).eq.'-an')read(argument(4:),*,err=201)nth
         if(argument(1:3).eq.'-ni')read(argument(4:),*,err=201)n_part
         if(argument(1:3).eq.'-pn')read(argument(4:),*,err=201)numprocs
         if(argument(1:3).eq.'-ri')read(argument(4:),*,err=201)rhoinf
         if(argument(1:3).eq.'-ck')read(argument(4:),*,err=201)ickst
         if(argument(1:3).eq.'-dt')read(argument(4:),*,err=201)dt
         if(argument(1:3).eq.'-da')read(argument(4:),*,err=201)bdt
         if(argument(1:7).eq.'--reinj')
     $        read(argument(8:),*,err=201)ninjcomp
         if(argument(1:2).eq.'-s')then
            read(argument(3:),*,err=201)nsteps
            if(nsteps.gt.nf_maxsteps)then
               if(lmyidhead)write(*,*)'Asked for more steps',nsteps,
     $              ' than allowed. Limit ',nf_maxsteps
               nsteps=nf_maxsteps
            endif
         endif
         if(argument(1:9).eq.'--restart')lrestart=.true.
         if(argument(1:2).eq.'-l')read(argument(3:),*,err=201)debyelen
         if(argument(1:2).eq.'-v')read(argument(3:),*,err=201)vd
         if(argument(1:2).eq.'-t')read(argument(3:),*,err=201)Ti
         if(argument(1:9).eq.'--objfile')
     $        read(argument(10:),'(a)',err=201)objfilename
         if(argument(1:3).eq.'-ho')then
            call geomdocument()
            call exit(0)
         endif
         if(argument(1:2).eq.'-h')goto 203
         if(argument(1:2).eq.'-?')goto 203
      enddo
      if(n_part.ne.0)rhoinf=0.
c Set ninjcomp if we are using rhoinf
      if(rhoinf.ne.0)call nreincalc(dt)
      goto 202
c------------------------------------------------------------
c Help text
 201  continue
      if(lmyidhead)write(*,*)'=====Error reading command line argument'
 203  continue
      if(myid.ne.0)goto 202
 301  format(a,i5)
 302  format(a,f8.3)
      write(*,301)'Usage: ccpic [switches]'
      write(*,301)'Particle switches.'
     $     //' Leave no gap before value. Defaults indicated [ddd'
      write(*,301)' -ni   set No of particles/node; zero => unset. ['
     $     ,n_part
      write(*,301)' --reinj    set reinjection number at each step.['
     $     ,ninjcomp
      write(*,302)' -ri   set rhoinfinity/node => reinjection number. ['
     $     ,rhoinf
      write(*,302)' -dt   set Timestep.  [',dt,
     $     ' -da   set Initial dt accel-factor. [',bdt
      write(*,301)' -s    set No of steps. [',nsteps
      write(*,302)' -v    set Drift velocity. [',vd
      write(*,302)' -t    set Ion Temperature. [',Ti
      write(*,302)' -l    set Debye Length. [',debyelen
      write(*,301)' --objfile<filename>  set name of object data file.'
     $     //' [ccpicgeom.dat'
      write(*,301)' --restart  Attempt to restart from saved state.'
      write(*,301)'Debugging switches for testing'
      write(*,301)' -gt   Plot regions and solution tests.'
      write(*,301)' -gf   Plot flux evolution and final distribution.'
      write(*,301)' -gi   Plot injection accumulated diagnostics.'
      write(*,301)' -gs   Plot slices of solution potential, density. '
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
      do id=1,ndims
         ium2(id)=iuds(id)-2
      enddo         
c-----------------------------------------------------------------
c Geometry and boundary information. Read in and initialize:
      call readgeom(objfilename,myid)
c Used later. Ought to purge.
      rc=obj_geom(5,1)
      phip=-obj_geom(10,1)/obj_geom(8,1)
      if(lmyidhead)write(*,*)'rc=',rc,'  phip=',phip
c---------------------------------------------------------------
c Construct the mesh vector(s) and ium2
      call meshconstruct(ndims,iuds,xmstart,xmend)
c----------------------------------------------------------------
c Initializations
      if(lmyidhead)write(*,*)'Initializing the stencil data cij'
c Initialize cij:
      ipoint=iLs(1)+iLs(2)+iLs(3)

c Call to mditerarg instead: Seems to work even though the cijroutine is
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
c We are going to populate the region in which xir lies.
c         iregion=insideall(ndims,xir)
         iregion=insidemask(ndims,xir)
         region=iregion
         if(lmyidhead)write(*,*)
     $        'Starting volume setting. Be patient this first time...'
         call mditerarg(volnode,ndims,ifull,ium2,ipoint,
     $        volumes,region,cij,dum4)
         if(lmyidhead)write(*,*)'Finished volume setting'
c If head, write the geometry data if we've had to calculate it.
         if(lmyidhead)call stored3geometry(volumes,iuds,ifull,istat)
      endif
c---------------------------------------------
c Set an object pointer for all the edges so their regions get
c set by the iregioninit call
      ipoint=0
c      call mditerate(ndims,ifull,iuds,cijedge,cij,ipoint)
      call mditerarg(cijedge,ndims,ifull,iuds,ipoint,cij,dum2,dum3,dum4)
c Initialize the region flags in the object data
      call iregioninit(ndims,ifull)
c---------------------------------------------
c Old position of sormpi initialization.
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
         if(lmyidhead)write(*,*)'Finished mesh setup:',iuds
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
c---------------------------------------------------------------         
c Control. Bit 1, use my sor params (not here). Bit 2 use faddu (not)
      ictl=0
c      write(*,*)'Calling sormpi, ni,nj=',ni,nj
c An initial solver call.
      call sormpi(ndims,ifull,iuds,cij,u,q,bdyset,faddu,ictl,ierr
     $     ,myid,idims)
      ictl=2
c
      if(lmyidhead)then
c-------------------------------------------------------------------
c Do some analytic checking of the case with a fixed potential sphere
c inside a logarithmic derivative boundary condition. 1/r solution.
c Also write out some data for checking.
         call spherecheck(ifull,iuds,u,phip,rc)
         if(ltestplot)then
c Plot some of the initial-solver data.
            call solu3plot(ifull,iuds,u,cij,phip,rc,thetain,nth
     $        ,rs)
            if(lmyidhead)write(*,*)'Return from solu3plot.'
            stop
         endif
      endif
c End of plotting.
c------------------------------------------------------------------
c We are going to populate the region in which xir lies.
c This needs a more general definition to accommodate union regions.
c      iregion_part=insideall(ndims,xir)
      lregion=linregion(ibool_part,ndims,xir)
c With a specified number of particles.
c Initialize the fortran random number generator.
      idum=-myid-1
      blah=ran1(idum) 
      if(lmyidhead)write(*,*)'Initializing',n_part,' particles'
      call pinit()
c      if(lmyidhead)write(*,*)'Return from pinit'
c A special orbit.
c Pinit resets x_part. So set it for the special first particle.
      x_part(1,1)=2.
      x_part(2,1)=0.
      x_part(3,1)=0.
c Prior half-step radial velocity
      x_part(4,1)=0.5*dt*(abs(phip)/x_part(1,1)**2)
c Tangential velocity of circular orbit at r=4.
      x_part(5,1)=sqrt(abs(phip/x_part(1,1))-x_part(4,1)**2)
      x_part(6,1)=0.
c
c      if(lmyidhead)write(*,*)'dt=',dt,' vd=',vd
c ' dtheta=',dt*x_part(5,1)/x_part(1,1),
c     $     ' steps=',nsteps,' total theta/2pi='
c     $     ,nsteps*dt*x_part(5,1)/x_part(1,1)/2./3.1415927
c---------------------------------------------
c Initialize the fluxdata storage and addressing.
      call fluxdatainit()
c---------------------------------------------
      phirein=0.
      ninjcomp0=ninjcomp
      if(ninjcomp.ne.0)write(*,*)'Fixed injection count:',ninjcomp
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
         call phiread(phifilename,ifull,iuds,u,ierr)
         if(ierr.ne.0)goto 401
         write(*,*)'Restart files read successfully.'
         if(nsteps+nf_step.gt.nf_maxsteps)then
            if(lmyidhead)write(*,*)'Asked for',
     $           nsteps,' in addition to',nf_step,
     $           ' Total',nsteps+nf_step,
     $           ' too much; set to',nf_maxsteps
            nsteps=nf_maxsteps-nsteps
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
c Create addtype and operator for reduce sum. Transferred to psumreduce.f
c      call mpisubopcreate(ndims,ifull,iuds,addsubarray_mpi,
c     $     iaddtype,iaddop)
c-----------------------------------------------
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
c Because psumtoq internally compensates for faddu, we reduction here
         call psumreduce(psum,ndims,ifull,iuds,iLs) 
c Calculate rhoinfinity, needed in psumtoq. Dependent on reinjection type.
         call rhoinfcalc(dt)
c Convert psums to charge, q. Remember external psumtoq!
         call mditerarg(psumtoq,ndims,ifull,ium2,
     $        0,psum(2,2,2),q(2,2,2),volumes(2,2,2),u(2,2,2))

         call sormpi(ndims,ifull,iuds,cij,u,q,bdyset,faddu,ictl,ierr
     $     ,myid,idims)
c If there are processes in the main communicator that are absent from the
c cartesian communicator. They have to be sent the new potential here.

         if(lmyidhead)write(*,'(i4.4,i4,'' iterations.'',$)')
     $        nf_step,ierr
c         write(*,*)dt
         if(lsliceplot)then
            call sliceGweb(ifull,iuds,u,Li,zp,
     $        ixnp,xn,ifix,'potential:'//'!Ay!@')
            call sliceGweb(ifull,iuds,q,Li,zp,
     $        ixnp,xn,ifix,'density: n')
         endif

         if(nf_step.eq.ickst) then
c This test routine assumes 3 full dimensions all equal to Li are used.
            call checkuqcij(Li,u,q,psum,volumes,cij,
     $           u2,q2,psum2,volumes2,cij2)
            call padvnc(ndims,cij,u,iLs)
            call checkx(n_part2,x_part2,
     $           if_part2,iregion_part2,ioc_part2,dt2,
     $           ldiags2,rhoinf2,nrein2,phirein2,numprocs2,ninjcomp2)
         else
c The normal call:
            call padvnc(ndims,cij,u,iLs)
         endif
         call fluxreduce()
c Store the step's rhoinf, dt.
         ff_rho(j)=rhoinf
         ff_dt(j)=dt
         if(lmyidhead)call fluxdiag()
c         if(lmyidhead.and.mod(j,nsteps/5+1).eq.0)
c     $  write(*,
c       write(*,
c     $    '(''nrein,n_part,ioc_part,rhoinf,dt='',i5,i7,i7,2f10.3)')
c     $        nrein,n_part,ioc_part,rhoinf,dt
c
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
      if(lmyidhead)call phiwrite(phifilename,ifull,iuds,u)

c-------------------------------------------------------------------
      call MPI_FINALIZE(ierr)

c Check some flux diagnostics and writing.
      if(lmyidhead)then 
         call outputflux(fluxfilename)
         if(linjplot)call plotinject(Ti)
         call fluxave(nsteps/2,nsteps,lfplot)
      endif
c      call readfluxfile(fluxfilename)
      end
c**********************************************************************
