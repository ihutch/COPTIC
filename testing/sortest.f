      program sortest
      include '../ndimsdecl.f'
c Main program of cartesian coordinate solver 
      include '../objcom.f'
c Mesh spacing description structure
      include '../meshcom.f'
      real u(na_i,na_j,na_k),q(na_i,na_j,na_k)
     $     ,cij(2*ndims+1,na_i,na_j,na_k)
      real volumes(na_i,na_j,na_k)
c Used dimensions, Full dimensions. Used dims-2
      integer iuds(ndims),ifull(ndims),ium2(ndims)
c Processor cartesian geometry can be set by default.
      integer nblksi,nblksj,nblksk
      parameter (nblksi=1,nblksj=1,nblksk=1)
      integer idims(ndims) 
c mpi process information.
      include '../myidcom.f'
c Structure vector needed for finding adjacent u values.
      integer iLs(ndims+1)
      external bdyshare,bdyset,bdysetnull,faddu,cijroutine,cijedge
      character*100 objfilename
      character*100 argument
      character*256 argline
      logical ltestplot,lsliceplot,linjplot
      logical lmyidhead,lphiplot,lpgraph
      integer ipstep,idebug
      real CFin(3+ndims,2*ndims)
      integer ipartperiod(ndims)
      include '../facebcom.f'
c Either include plascom or define vperp and Bfield.
      include '../plascom.f'
c      real vperp(ndims),Bfield(ndims)
      real zp(na_m,na_m,ndims)
c sor control values
      common /ctl_sor/mi_sor,xjac_sor,eps_sor,del_sor,k_sor
c Set up the structure vector.
      parameter (Li1=na_i,Li2=Li1*na_j,Li3=Li2*na_k)
      data iLs/1,Li1,Li2,Li3/
c Mesh and mpi parameter defaults:
      data idims/nblksi,nblksj,nblksk/
      data ifull/na_i,na_j,na_k/ipartperiod/0,0,0/
c Data for plotting etc.
      data iobpl/0/idebug/1/
      data ltestplot,lsliceplot,linjplot/
     $     .true.,.false.,.false./
      data lphiplot,lpgraph/.true.,.false./
      data ipstep/1/
c-------------------------------------------------------------
c Consistency checks
      if(ndims.ne.ndims)then
         write(*,*)'Inconsistent ndims, ndims',ndims,ndims
         stop
      endif
c-------------------------------------------------------------
      rs=5.0
      rsmesh=rs
      lmyidhead=.true.
c---------------------------------------------------------------------
c This necessary here so one knows early the mpi structure.
c Otherwise could have been hidden in sormpi and pass back numprocs.
      call mpigetmyid(myid,nprocs,ierr)
      if(myid.ne.0) lmyidhead=.false.
      numprocs=nprocs
      if(idebug.gt.0)write(*,*)'numprocs,myid',numprocs,myid
c--------------------------------------------------------------
c Deal with command-line arguments; not all valid here.
      call copticcmdline(lmyidhead,ltestplot,iobpl,iobpsw,rcij
     $     ,lsliceplot,ipstep,ldenplot,lphiplot,linjplot,ifplot,norbits
     $     ,thetain,nth,iavesteps,n_part,numprocs,ripernode,crelax,ickst
     $     ,colntime,dt,bdt,subcycle,dropaccel,eoverm,Bfield,Bt,ninjcomp
     $     ,nsteps,nf_maxsteps,vneutral,vd,ndiags,ndiagmax,debyelen,Ti
     $     ,iwstep,idistp,lrestart,restartpath,extfield,objfilename
     $     ,lextfield ,vpar,vperp,ndims,islp,slpD,CFin,iCFcount,LPF
     $     ,ipartperiod,lnotallp,Tneutral,Enfrac,colpow,idims,argline
     $     ,vdrift,ldistshow,gp0,gt,gtt,gn,gnt,nspecies,nspeciesmax
     $     ,numratioa)
c-----------------------------------------------------------------
c Finalize parameters after switch reading.
c Geometry and boundary information. Read in.
      if(idebug.gt.0)write(*,*)'Calling readgeom',myid,ifull,lmyidhead
      call readgeom(objfilename,myid,ifull,CFin,iCFcount,LPF,ierr
     $     ,argline)
      if(idebug.gt.0)write(*,*)'Finished readgeom'
c---------------------------------------------------------------
c Construct the mesh vector(s) from the geometry info.
      call meshconstruct(ndims,iuds,ifull,ipartperiod)
      if(lmyidhead)write(*,'(a,3i4,6f8.3)')
     $     ' Constructed mesh',iuds
     $     ,(xmeshstart(k),xmeshend(k),k=1,ndims)
c-----------------------------------------------------------------
c Initialize the face phi boundary conditions if we are using them.
c      write(*,*)'ICFCOUNT',iCFcount
      if(iCFcount.ne.0)then
         do idn=1,2*ndims
            call bdyfaceinit(idn,CFin(1,idn))
         enddo
c Now print out the result of the initialization.
         if(lmyidhead)call bdyfaceinit(-2,CFin(1,1))
      endif
c----------------------------------------------------------------
c Initializations
      if(lmyidhead)write(*,*)'Initializing the stencil data cij'
c Initialize cij for just the inner part, not the edges.
      do id=1,ndims
         ium2(id)=iuds(id)-2
      enddo         
      ipoint=iLs(1)+iLs(2)+iLs(3)
c      write(*,*)ndims,ifull
      call mditerarg(cijroutine,ndims,ifull,ium2,ipoint,
     $     cij(1,1,1,1),debyelen,dum3,dum4)
c---------------------------------------------
c Initialize the region flags in the object data
      call iregioninit(ifull)
c      if(lmyidhead)call reportfieldmask()
      ipoint=0
      call mditerarg(cijedge,ndims,ifull,iuds,ipoint,cij,dum2,dum3,dum4)
c This old position overruled the new edge setting.
c      call iregioninit(ndims,ifull)
c---------------------------------------------
      if(myid.ne.0)then
c Don't do plotting from any node except the master.
         ltestplot=.false.
         lsliceplot=.false.
         iobpl=0
      else
c---------------------------------------------
c Some simple graphics of cij
         if(lpgraph)call pfset(3)
         if(ltestplot)call text3rgraph(ndims,iuds,ifull,cij,volumes)
c---------------------------------------------
c The following requires include objcom.f
         if(lmyidhead)write(*,*)
     $      'Used No of pointers:',oi_cij,' of',iuds(1)*iuds(2)*iuds(3)
     $        ,' points.'
c Plot objects
         if(iobpl.ne.0.and.lmyidhead)then
            if(rcij.eq.0.)rcij=rs
            call cijplot(ndims,ifull,iuds,cij,rcij,iobpl,0)
          endif
      endif
c---------------------------------------------
c Initialize charge (set q to zero over entire array).
      call mditerset(q,ndims,ifull,iuds,0,0.)
c Initialize potential (set u to zero over entire array).
      call mditerset(u,ndims,ifull,iuds,0,0.)
c---------------------------------------------------------------     
c Control. Bit 1, use my sor params (not here). Bit 2 use faddu (not)
c Bit 4-6 periodicity.
      ictl=0
c Make dimensions periodic; necessary for mpi:
      do id=1,ndims
         if(LPF(id))ictl=ictl+4*2**id
      enddo
      write(*,*)'Calling sormpi, ni,nj=',ni,nj
c Solver call.
      do k=1,2
         if(k.gt.1)ictl=1
         call sormpi(ndims,ifull,iuds,cij,u,q,bdyshare,bdyset,faddu,ictl
     $        ,ierr,myid,idims)
         if(idebug.gt.0)write(*,*)'Completed sormpi',ierr,del_sor
         if(ierr.gt.0)goto 2
         eps_sor=2.e-6
      enddo
 2    continue
c Putting the finalize here prevents end mpi-crashes.
      call mpifinalize(ierr)
c
c-------------------------------------------------------------------
      if(lmyidhead)then
         if(ltestplot)call sliceGweb(ifull,iuds,u,na_m,zp,
     $              ixnp,xn,ifix,'potential:'//'!Ay!@'//char(0),dum,dum)

c This only does anything if object 2 is an outer sphere.
         if(idebug.gt.0)write(*,*)'Calling vaccheck',rs,ltestplot
         call vaccheck(ifull,iuds,cij,u,thetain,nth,rs,ltestplot)
      endif
c End of plotting.
      call exit(0)
c------------------------------------------------------------
 201  continue
      if(lmyidhead)write(*,*)'=====Error reading command line argument '
     $     ,argument(:20)
 203  continue
      if(lmyidhead)then
         write(*,*)'Do ./coptic -h to get help on switches'
         write(*,*)'Not all switches apply to sortest'
c      if(lmyidhead)call helpusage2()
      endif
      
      end
c**********************************************************************
