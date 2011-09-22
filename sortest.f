      program sortest
c Main program of cartesian coordinate solver 

      include 'objcom.f'
c Storage array spatial count size
      include 'griddecl.f'
      real u(na_i,na_j,na_k),q(na_i,na_j,na_k)
     $     ,cij(2*ndims_sor+1,na_i,na_j,na_k)
      real volumes(na_i,na_j,na_k)
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
c Structure vector needed for finding adjacent u values.
      integer iLs(ndims_sor+1)
      external bdyset,faddu,cijroutine,cijedge
      character*100 objfilename
      character*100 argument
      logical ltestplot,lsliceplot,linjplot
      logical lmyidhead,lphiplot,lpgraph
      integer ipstep
c sor control values
      common /ctl_sor/mi_sor,xjac_sor,eps_sor,del_sor,k_sor
c Set up the structure vector.
      parameter (Li1=na_i,Li2=Li1*na_j,Li3=Li2*na_k)
      data iLs/1,Li1,Li2,Li3/
c Mesh and mpi parameter defaults:
      data idims/nblksi,nblksj,nblksk/
      data ifull/na_i,na_j,na_k/
c Data for plotting etc.
      data iobpl/0/
      data ltestplot,lsliceplot,linjplot/
     $     .true.,.false.,.false./
      data lphiplot,lpgraph/.true.,.false./
      data ipstep/1/
c-------------------------------------------------------------
c Consistency checks
      if(ndims.ne.ndims_sor)then
         write(*,*)'Inconsistent ndims, ndims_sor',ndims,ndims_sor
         stop
      endif
c-------------------------------------------------------------
c Defaults:
c Default edge-potential (chi) relaxation rate.     
      objfilename='copticgeom.dat'
      debyelen=1.
      Ti=1.
      crelax=1.*Ti/(1.+Ti)
      vd=0.
      rs=5.0
      rsmesh=rs
      numprocs=1
      thetain=.1
      nth=1
      norbits=0
      ickst=0
      lmyidhead=.true.
      ndiags=0
      ifplot=-1
      rcij=0
c---------------------------------------------------------------------
c This necessary here so one knows early the mpi structure.
c Otherwise could have been hidden in sormpi and pass back numprocs.
      call mpigetmyid(myid,nprocs,ierr)
      if(myid.ne.0) lmyidhead=.false.
      numprocs=nprocs
c--------------------------------------------------------------
c Deal with arguments
      do i=1,iargc()
         call getarg(i,argument)
         if(argument(1:3).eq.'-gt')ltestplot=.true.
         if(argument(1:3).eq.'-gc')read(argument(4:),*,end=201)iobpl
         if(argument(1:3).eq.'-gr')read(argument(4:),*,end=201)rcij
         if(argument(1:3).eq.'-gs')then
            lsliceplot=.true.
            read(argument(4:),*,err=210,end=210)ipstep
            goto 211
 210        ipstep=1
 211        continue
            write(*,*)'ipstep=',ipstep
         endif
         if(argument(1:3).eq.'-gw')lpgraph=.true.
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
         if(argument(1:3).eq.'-pn')read(argument(4:),*,err=201)numprocs
         if(argument(1:2).eq.'-l')read(argument(3:),*,err=201)debyelen
         if(argument(1:2).eq.'-t')read(argument(3:),*,err=201)Ti
         if(argument(1:3).eq.'-of')
     $        read(argument(4:),'(a)',err=201)objfilename
         if(argument(1:3).eq.'-ho')then
            call geomdocument()
            call exit(0)
         endif
         if(argument(1:2).eq.'-h')goto 203
         if(argument(1:2).eq.'-?')goto 203
      enddo
c-----------------------------------------------------------------
c Finalize parameters after switch reading.
c Geometry and boundary information. Read in.
      call readgeom(objfilename,myid,ifull)
c---------------------------------------------------------------
c Construct the mesh vector(s) from the geometry info.
      call meshconstruct(ndims,iuds,ifull)
      if(lmyidhead)write(*,'(a,3i4,6f8.3)')
     $     ' Constructed mesh',iuds
     $     ,(xmeshstart(k),xmeshend(k),k=1,ndims)
c----------------------------------------------------------------
c Initializations
      if(lmyidhead)write(*,*)'Initializing the stencil data cij'
c Initialize cij for just the inner part, not the edges.
      do id=1,ndims
         ium2(id)=iuds(id)-2
      enddo         
      ipoint=iLs(1)+iLs(2)+iLs(3)
      call mditerarg(cijroutine,ndims,ifull,ium2,ipoint,
     $     cij(1,1,1,1),debyelen,dum3,dum4)
c---------------------------------------------
c Initialize the region flags in the object data
      call iregioninit(ndims,ifull)
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
     $      'Used No of pointers:',oi_sor,' of',iuds(1)*iuds(2)*iuds(3)
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
c An inital vacuum solution with zero density. 
c Control. Bit 1, use my sor params (not here). Bit 2 use faddu (not)
      ictl=0
c      write(*,*)'Calling sormpi, ni,nj=',ni,nj
c Solver call.
      do k=1,2
         if(k.gt.1)ictl=1
         call sormpi(ndims,ifull,iuds,cij,u,q,bdyset,faddu,ictl,ierr
     $        ,myid,idims)
         write(*,*)'Completed sormpi',ierr,del_sor
         if(ierr.gt.0)goto 2
         eps_sor=2.e-6
      enddo
 2    continue
c Putting the finalize here prevents end mpi-crashes.
      call mpifinalize(ierr)
c
c-------------------------------------------------------------------
      if(lmyidhead)then
         call vaccheck(ifull,iuds,cij,u,thetain,nth,rs,ltestplot)
      endif
c End of plotting.
      call exit(0)
c------------------------------------------------------------
 201  continue
      if(lmyidhead)write(*,*)'=====Error reading command line argument '
     $     ,argument(:20)
 203  continue
      write(*,*)'Not all switches apply to sortest'
      if(lmyidhead)call helpusage2()
      end
c**********************************************************************
