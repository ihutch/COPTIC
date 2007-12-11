      program ccpic
c
      include 'objcom.f'
c Storage array spatial count size
      integer Li,ni,nj,nk
c      parameter (Li=100,ni=40,nj=40,nk=16)
      parameter (Li=100,ni=16,nj=16,nk=20)
      parameter (Li2=Li*Li,Li3=Li2*Li)
c      parameter (Li=100,ni=32,nj=32,nk=40)
      real u(Li,Li,Li),q(Li,Li,Li),cij(2*ndims_sor+1,Li,Li,Li)
      real psum(Li,Li,Li),volumes(Li,Li,Li)
c Used dimensions, Full dimensions. Used dims-2
      integer iuds(ndims_sor),ifull(ndims_sor),ium2(ndims_sor)
c Mesh spacing description structure
      include 'meshcom.f'
c Processor cartesian geometry
      integer nblksi,nblksj,nblksk
      parameter (nblksi=1,nblksj=1,nblksk=1)
      integer idims(ndims_sor)
c MPI information.
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
      integer ndims,nd2
      parameter (ndims=ndims_sor,nd2=ndims*2)

      external bdyset,faddu,cijroutine,cijedge,psumtoq
      external volnode
      character*100 fluxfilename
      character*100 objfilename
      character*100 argument
      common /ctl_sor/mi_sor,xjac_sor,eps_sor,del_sor,k_sor
      logical ltestplot,lcijplot,lsliceplot,lorbitplot

c Diagnostics
c      real error(Li,Li,Li)
      real zp(Li,Li),cijp(2*ndims_sor+1,Li,Li)
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
      data ltestplot,lcijplot,lsliceplot,lorbitplot/
     $     .false.,.false.,.false.,.false./
      data thetain,nth/.1,1/

c Defaults:
c Fixed number of particles rather than fixed injections.
      ninjcomp=0
      n_part=0
c Default to constant rhoinf not n_part.
      rhoinf=100.
      dt=.1
      objfilename='ccpicgeom.dat'
      nsteps=5.
      debyelen=1.
      Ti=1.
      vd=0.
      rs=5.
      numprocs=1
      bdt=1.
      norbits=5

c--------------------------------------------------------------
c Deal with arguments
c      if(iargc().eq.0) goto "help"
      do i=1,iargc()
         call getarg(i,argument)
         if(argument(1:3).eq.'-gt') ltestplot=.true.
         if(argument(1:3).eq.'-gc') lcijplot=.true.
         if(argument(1:3).eq.'-gs') lsliceplot=.true.
         if(argument(1:3).eq.'-go') lorbitplot=.true.
         if(argument(1:3).eq.'-at')read(argument(4:),*)thetain
         if(argument(1:3).eq.'-an')read(argument(4:),*)nth
         if(argument(1:7).eq.'--reinj')read(argument(8:),*)ninjcomp
         if(argument(1:3).eq.'-ni')read(argument(4:),*)n_part
         if(argument(1:3).eq.'-ri')read(argument(4:),*)rhoinf
         if(argument(1:3).eq.'-dt')read(argument(4:),*)dt
         if(argument(1:3).eq.'-da')read(argument(4:),*)bdt
         if(argument(1:2).eq.'-s')read(argument(3:),*)nsteps
         if(argument(1:2).eq.'-l')read(argument(3:),*)debyelen
         if(argument(1:13).eq.'--objfilename')
     $        read(argument(14:),'(a)')objfilename
         if(argument(1:2).eq.'-h')goto 201
         if(argument(1:2).eq.'-?')goto 201
      enddo
      if(n_part.ne.0)rhoinf=0.
c Set ninjcomp if we are using rhoinf
      if(rhoinf.ne.0)call nreincalc(dt)
      goto 202
c------------------------------------------------------------
c Help text
 201  continue
      write(*,'(a)')'Usage: ccpic [switches]'
     $     ,' -gt   Plot solution tests. -gc   Plot mesh and stencils. '
     $     ,' -gs   Plot slices of solution potential. '
     $     //' -go   Plot orbits. '
     $     ,' -atf.fff   set test angle. -anf.fff   set No of angles. '
     $     ,' --reinjnnn    set reinjection number at each step.'
     $     ,' -ni   set No of particles/node.'
     $     ,' -ri   set rhoinfinity instead of total particles.'
     $     ,' -dt   set timestep.    -da   set initial dt accel-factor'
     $     ,' -snnn  set No of steps.'
     $     //'   --objfile<filename>  set name of object data file.'
     $     ,' -h -?   Print usage.'
      call exit(0)
c-------------------------------------------------------------
c Geometry information read in.
 202  call readgeom(objfilename)
c First object is sphere of radius rc and potential phi.
      rc=obj_geom(5,1)
      phip=-obj_geom(10,1)/obj_geom(8,1)
      write(*,*)'rc=',rc,'  phip=',phip
c Second object is bounding sphere of radius rs. 
c But use a tad more for the mesh size
      rs=obj_geom(5,2)*1.00001
c---------------------------------------------------------------
c Construct the mesh vector(s) and ium2
      call meshconstruct(ndims,iuds,ium2,rs)
c----------------------------------------------------------------
c INITIALIZATIONS
      write(*,*)'Initializing the stencil data cij'
c Initialize cij:
      ipoint=iLs(1)+iLs(2)+iLs(3)
c      call mditerate(ndims,ifull,ium2,cijroutine,
c     $     cij(1,1,1,1),ipoint)

c Call to mditerarg instead: Seems to work even though the cijroutine is
c not defined with sufficient arguments for all those in this call.
      call mditerarg(cijroutine,ndims,ifull,ium2,ipoint,
     $     cij(1,1,1,1),debyelen,dum3,dum4)

c---------------------------------------------
c Here we try to read the geometry data.
      istat=1
      call stored3geometry(volumes,iuds,ifull,istat)
      if(istat.eq.0)then
c Calculate the nodal volumes for all non-edge points.
         ipoint=iLs(1)+iLs(2)+iLs(3)
c We are going to populate the region in which xir lies.
         iregion=insideall(ndims,xir)
         region=iregion
         write(*,*)
     $        'Starting volume setting. Be patient this first time...'
         call mditerarg(volnode,ndims,ifull,ium2,ipoint,
     $        volumes,region,cij,dum4)
         write(*,*)'Finished volume setting'
c Here we write the geometry data if we've had to calculate it.
         call stored3geometry(volumes,iuds,ifull,istat)
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
c Control bit 3 (=4) pure initialization, no communication (to get myid).
      call sormpi(ndims,ifull,iuds,cij,u,q,bdyset,faddu,4,ierr
     $     ,myid,idims)
      if(myid.ne.0)then
c Don't do plotting from any node except the master.
         ltestplot=.false.
         lcijplot=.false.
         lsliceplot=.false.
         lorbitplot=.false.
      endif
c---------------------------------------------
c Some simple graphics of cij, and volumes.
      if(ltestplot)call text3graphs(ndims,iuds,ifull,cij,volumes)
c---------------------------------------------
c The following requires include objcom.f
      write(*,*)'Finished mesh setup:',iuds,
     $     ' Used No of pointers:',oi_sor,' of',iuds(1)*iuds(2)*iuds(3)
     $     ,' points.'
c      write(*,*)'Finished mesh setup.'
      if(lcijplot)call cijplot(ndims,ifull,iuds,cij,rs)

c Initialize charge (set q to zero over entire array).
      call mditerset(q,ndims,ifull,iuds,0,0.)
c---------------------------------------------------------------         
c Control. Bit 1, use my sor params (not here). Bit 2 use faddu (not)
      ictl=0
c      write(*,*)'Calling sormpi, ni,nj=',ni,nj
c An initial solver call.
      call sormpi(ndims,ifull,iuds,cij,u,q,bdyset,faddu,ictl,ierr
     $     ,myid,idims)
      ictl=2

      if(myid.eq.0)then
c------------------------------------------------------------------
c         if(ierr.lt.0)then
c            write(*,*)'Not converged.',ierr
c            write(*,*) 'mi_sor,k_sor,xjac_sor,del_sor',
c     $           mi_sor,k_sor,xjac_sor,del_sor,myid
c         endif
c-------------------------------------------------------------------
c-------------------------------------------------------------------
c Do some analytic checking of the case with a fixed potential sphere
c inside a logarithmic derivative boundary condition. 1/r solution.
c Also write out some data for checking.
         call spherecheck(ifull,iuds,u,phip,rc)
      endif
      if(ltestplot)then
c Plot some of the initial-solver data.
         call solu3plot(ifull,iuds,u,cij,phip,rc,thetain,nth
     $        ,rs)
         write(*,*)'Return from solu3plot.'
      endif
c End of plotting.
c------------------------------------------------------------------
c We are going to populate the region in which xir lies.
      iregion_part=insideall(ndims,xir)
c With a specified number of particles.
c      n_part=2
      call srand(myid)
      write(*,*)'Initializing',n_part,' particles'
      call pinit()
c      write(*,*)'Return from pinit'
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
      write(*,*)'dt=',dt,' dtheta=',dt*x_part(5,1)/x_part(1,1),
     $     ' steps=',nsteps,' total theta/2pi='
     $     ,nsteps*dt*x_part(5,1)/x_part(1,1)/2./3.1415927
c---------------------------------------------
c Initialize the fluxdata storage and addressing.
      call fluxdatainit()
c---------------------------------------------

      phirein=0.
      ninjcomp0=ninjcomp
      maccel=nsteps/3
      dtf=dt

      do j=1,nsteps
         nf_step=j
c Acceleration code.
         bdtnow=max(1.,(bdt-1.)*(maccel-j+2)/(maccel+1.)+1.)
         dt=bdtnow*dtf
         ninjcomp=bdtnow*ninjcomp0
         if(ninjcomp.ne.0)nrein=ninjcomp
c
         write(*,'(i4,i4''  x_p='',7f9.5)')j,ierr, (x_part(k,1),k=1,6)
     $        ,sqrt(x_part(1,1)**2+x_part(2,1)**2)
         call mditerset(psum,ndims,ifull,iuds,0,0.)
         call chargetomesh(psum,iLs,diags)
c Calculate rhoinfinity, needed in psumtoq.
         call rhoinfcalc(dt)
         write(*,*)'nrein,n_part,rhoinf,dt=',nrein,n_part,rhoinf,dt
c Convert psums to charge, q. Remember external psumtoq!
         call mditerarg(psumtoq,ndims,ifull,ium2,
     $        0,psum(2,2,2),q(2,2,2),volumes(2,2,2),u(2,2,2))
c
         call sormpi(ndims,ifull,iuds,cij,u,q,bdyset,faddu,ictl,ierr
     $     ,myid,idims)
         if(lsliceplot) call slice3web(ifull,iuds,u,cij,Li,zp,cijp,
     $        ixnp,xn,ifix,'potential:'//'!Ay!@',1)
         if(lsliceplot)call slice3web(ifull,iuds,q,cij,Li,zp,cijp,
     $        ixnp,xn,ifix,'density: n',0)
c
         call padvnc(ndims,cij,u,iLs)
         call fluxdiag()
c         call partreduce
      enddo
c      write(*,*)iorbitlen(1),(xorbit(k,1),k=1,10)
c      call slice3web(ifull,iuds,psum,cij,Li,zp,cijp,
c     $        ixnp,xn,ifix,'psum',0)

      if(lorbitplot)call orbit3plot(ifull,iuds,u,phip,rc,rs)
c-------------------------------------------------------------------
      call MPI_FINALIZE(ierr)

c Check some flux diagnostics and writing.
c      call fluxave()
      call outputflux(fluxfilename)
c      call readfluxfile(fluxfilename)
c      call fluxave()
      end
c**********************************************************************
      subroutine meshconstruct(ndims,iuds,ium2,rs)
      integer iuds(ndims),ium2(ndims)
      include 'meshcom.f'

      iof=0
      do id=1,ndims
c Pointer to start of vector.
         ixnp(id)=iof
c Mesh data. Make it extend from -rs to +rs
         do i=1,iuds(id)
            xn(i+iof)=rs*(2*(i-1.)/(iuds(id)-1.) - 1.)
            if(.false.)then
c Two region Non-uniform assuming iuds even.
            iu2=iuds(id)/2
c Uniform -1 to 1:
            xi=(i-iu2-0.5)/(iu2-0.5)
c Half of mesh at each spacing
            xmid=0.5
c First half of mesh extends to this fraction of the total rs.
            xs=0.25
            if(abs(xi).lt.xmid)then
               xn(i+iof)=rs*xs*xi/xmid
            else
               xn(i+iof)=rs*sign(((abs(xi)-xmid)+xs*(1.-abs(xi))),xi)
     $              /(1-xmid)
            endif
            endif
         enddo
         iof=iof+iuds(id)
c Mesh size with the faces removed:
         ium2(id)=iuds(id)-2
      enddo
      ixnp(ndims+1)=iof

      end
c*********************************************************************
      subroutine text3graphs(ndims,iuds,ifull,cij,volumes)
      integer iuds(ndims),ifull(ndims)
      real cij(*)
      real volumes(*)
      include 'meshcom.f'
      character*40 form1
c Standard volume for uniform mesh:
         vs=.99999
         do id=1,ndims
            vs=vs*(xn(ixnp(id)+2)-xn(ixnp(id)+1))
         enddo
         vs=vs
c Text graphic of slice through volumes
         write(*,*)'Volumes percentages:'
         write(form1,'(''('',i2,''i4)'')')iuds(1)
         write(*,form1)((nint(100.*
     $        volumes(j-1+(iuds(2)/2-1)*ifull(1)
     $                +(k-1)*ifull(1)*ifull(2)+1)/vs),
     $        j=1,iuds(1)),k=1,iuds(3))

c Text graphic of slice through cij
c      write(*,*)'iregion:'
         write(form1,'(''('',i2,''i1)'')')iuds(1)
         write(*,form1)((ireg3(j,iuds(2)/2,k,ifull,cij),
     $        j=1,iuds(1)),k=1,iuds(3))
      end
