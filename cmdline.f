c**********************************************************************
c Encapsulation of parameter setting.
      subroutine copticcmdline(lmyidhead,ltestplot,iobpl,iobpsw,rcij
     $     ,lsliceplot,ipstep,ldenplot,lphiplot,linjplot,ifplot,norbits
     $     ,thetain,nth,iavesteps,n_part,numprocs,ripernode,crelax,ickst
     $     ,colntime,dt,bdt,subcycle,rmtoz,Bfield,ninjcomp,nsteps
     $     ,nf_maxsteps,vneutral,vd,ndiags,ndiagmax,debyelen,Ti,iwstep
     $     ,idistp,lrestart,restartpath,extfield,objfilename,lextfield
     $     ,vpar,vperp,ndims,islp,slpD)
      implicit none

      integer iobpl,iobpsw,ipstep,ifplot,norbits,nth,iavesteps,n_part
     $     ,numprocs,ickst,ninjcomp,nsteps,nf_maxsteps,ndiags,ndiagmax
     $     ,iwstep,idistp,ndims,islp
      logical lmyidhead,ltestplot,lsliceplot,ldenplot,lphiplot,linjplot
     $     ,lrestart,lextfield
      real rcij,thetain,ripernode,crelax,colntime,dt,bdt,subcycle,rmtoz
     $     ,vneutral,vd,debyelen,Ti,extfield,vpar,slpD
      real Bfield(ndims),vperp(ndims)
      character*100 restartpath,objfilename

c Local variables:
      integer lentrim
      external lentrim
      integer i
      real Bt
      character*100 argument

c Set defaults first.

c Fixed number of particles rather than fixed injections.
      ninjcomp=0
      n_part=0
c Default to constant ripernode not n_part.
      ripernode=100.
      debyelen=1.
      vd=0.
      Ti=1.
c Default edge-potential (chi) relaxation rate.     
      crelax=1.*Ti/(1.+Ti)
      dt=.1
      objfilename='copticgeom.dat'
      nsteps=5
      subcycle=0.
      colntime=0.
      vneutral=0.
      numprocs=1
      bdt=1.
      thetain=.1
      nth=1
      norbits=0
      ickst=0
      iavesteps=100
      ndiags=0
      ifplot=-1
      iwstep=99999
      rcij=0
      rmtoz=1.
      iobpsw=1
c Boundary condition switch and value. 0=> logarithmic.
      islp=0
      slpD=0.


c Deal with arguments
c      if(iargc().eq.0) goto "help"
      do i=1,iargc()
         call getarg(i,argument)
         if(argument(1:3).eq.'-gt')ltestplot=.true.
         if(argument(1:3).eq.'-gc')read(argument(4:),*,end=201)iobpl
         if(argument(1:3).eq.'-gw')read(argument(4:),*,end=201)iobpsw
         if(argument(1:3).eq.'-gr')read(argument(4:),*,end=201)rcij
         if(argument(1:3).eq.'-gs')then
            lsliceplot=.true.
            read(argument(4:),*,err=210,end=210)ipstep
            goto 211
 210        ipstep=1
 211        continue
            if(lmyidhead)write(*,*)'Plotting ipstep=',ipstep
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
         if(argument(1:3).eq.'-mz')read(argument(4:),*,err=201)rmtoz
         if(argument(1:3).eq.'-BC')read(argument(4:),*,err=201)islp
         if(argument(1:3).eq.'-Bx')then
            read(argument(4:),*,err=201)Bfield(1)
         elseif(argument(1:3).eq.'-By')then
            read(argument(4:),*,err=201)Bfield(2)
         elseif(argument(1:3).eq.'-Bz')then
            read(argument(4:),*,err=201)Bfield(3)
         endif
         if(argument(1:7).eq.'--reinj')
     $        read(argument(8:),*,err=201)ninjcomp
         if(argument(1:3).eq.'-rn')
     $        read(argument(4:),*,err=201)ninjcomp
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
c For Mach bdy, set slpD equal to M.
            slpD=vd
         endif
         if(argument(1:2).eq.'-md')then 
            read(argument(4:),*,err=201)ndiags
            if(ndiags.gt.ndiagmax)then
               write(*,*)'Error: Too many diag-moments',ndiags
               stop
            endif
         endif
         if(argument(1:2).eq.'-l')read(argument(3:),*,err=201)debyelen
         if(argument(1:2).eq.'-t')read(argument(3:),*,err=201)Ti
         if(argument(1:2).eq.'-w')read(argument(3:),*,err=201)iwstep
         if(argument(1:3).eq.'-pd')then
            idistp=1
            read(argument(4:),*,err=240)idistp
         endif
         if(argument(1:3).eq.'-fs')then
            lrestart=.true.
            read(argument(4:),'(a)',err=201)restartpath
         endif
         if(argument(1:10).eq.'--extfield')then
            read(argument(11:),*,err=201)extfield
c            write(*,*)'||||||||||||||extfield',extfield
            lextfield=.true.
         endif
         if(argument(1:3).eq.'-of')
     $        read(argument(4:),'(a)',err=201)objfilename
         if(argument(1:1).ne.'-')
     $        read(argument(1:),'(a)',err=201)objfilename

         if(argument(1:3).eq.'-ho')then
            call geomdocument()
            call exit(0)
         endif
c Help text options
         if(argument(1:2).eq.'-h')goto 203
         if(argument(1:3).eq.'--h')goto 203
         if(argument(1:2).eq.'-?')goto 203
c Indicator that coptic arguments are ended.
         if(argument(1:3).eq.'-ea'.or.argument(1:2).eq.'--')then
            write(*,*)'Arguments terminated by: ',argument(1:3)
            goto 202
         endif
 240     continue
      enddo
 202  continue
      Bt=0.
      do i=1,ndims
         Bt=Bt+Bfield(i)**2
      enddo
      if(Bt.ne.0)then
c Normalize magnetic field.
         Bt=sqrt(Bt)
         do i=1,ndims
            Bfield(i)=Bfield(i)/Bt
         enddo
c Assume that vd is in the z-direction
         vpar=vd*Bfield(ndims)
         do i=1,ndims
            vperp(i)=-Bfield(i)*vpar
         enddo
         vperp(ndims)=vperp(ndims)+vd
         write(*,*)'vpar,vperp',vpar,vperp
      else
c Zero the vparallel and vperp. Probably not necessary; but tidy.
         vpar=0.
         do i=1,ndims
            vperp(i)=0.
         enddo
      endif
c      write(*,*)'Bfield',Bfield

      return
c------------------------------------------------------------
c Help text
 201  continue
      if(lmyidhead)write(*,*)'=====Error reading command line argument '
     $     ,argument(:20)
 203  continue
      if(lmyidhead)then
c         call helpusage()
 301  format(a,i5,a,i5)
 302  format(a,3f8.3)
      write(*,301)'Usage: coptic [objectfile] [-switches]'
      write(*,301)'Parameter switches.'
     $     //' Leave no gap before value. Defaults or set values [ddd'
      write(*,301)' -ni   set No of particles/node; zero => unset.    ['
     $     ,n_part
      write(*,301)' -rn   set reinjection number at each step.        ['
     $     ,ninjcomp
      write(*,302)' -ri   set rhoinfinity/node => reinjection number. ['
     $     ,ripernode
      write(*,302)' -rx   set Edge-potl relax rate: 0=>off, 1=>immed. ['
     $     ,crelax
      write(*,302)' -dt   set Timestep.              [',dt
      write(*,302)' -da   set Initial dt accel-factor[',bdt
      write(*,302)' -ds   set subcycle frac. Bad!    [',subcycle
      write(*,301)' -s    set No of steps.           [',nsteps
      write(*,302)' -v    set Drift (z-)velocity.    [',vd
      write(*,302)' -t    set Ion Temperature.       [',Ti
      write(*,302)' -l    set Debye Length.          [',debyelen
      write(*,302)' -mz   set mass/Z ratio           [',rmtoz
      write(*,302)' -ct   set collision time.        [',colntime
      write(*,302)' -vn   set neutral drift velocity [',vneutral
      write(*,302)' -Bx -By -Bz set mag field compts [',Bfield
      write(*,301)' -a    set averaging steps.       [',iavesteps
     $     ,'     Also period of diagnostic writes.'
      write(*,301)' -w    set write-step period.     [',iwstep
      write(*,301)' -pd   set distribution diags     [',idistp
     $     ,'     Bits:1 write, 2 plot.'
      write(*,301)' -md   set No of diag-moments(7). [',ndiags
      write(*,301)' -BC   set boundary condition type[',islp
     $     ,'     Bits:1 Mach/Log[=0]  2-4:Face1-3=0'
c      write(*,301)' -xs<3reals>, -xe<3reals>  Set mesh start/end.'
      write(*,301)' [-of]<filename>  set name of object data file.'
     $     //'   ['//objfilename(1:lentrim(objfilename))
      write(*,301)
     $     ' -fs[path]  Attempt to restart from state saved [in path].'
      write(*,301)' -ea --  end argument parsing. Skip succeeding.'
      write(*,301)'Debugging switches for testing'
      write(*,301)' -gt   Plot regions and solution tests.'
      write(*,301)' -gi   Plot injection accumulated diagnostics.'
      write(*,301)' -gs[] Plot slices of solution potential, density. '
     $     //'[At step n]. [',ipstep
      write(*,301)' -gd -gp Turn off slicing of density, potential. '
      write(*,301)' -gf   set quantity plotted for flux evolution and'//
     $     ' final distribution. [',ifplot
      write(*,301)' -gw   set objplot sw. [+256:intercepts]'//
     $     ' Shade by 1:flux 2:flux-density[',iobpsw
      write(*,301)' -gc   set wireframe [& stencils(-)] mask.'//
     $     ' objects<->bits. [',iobpl
      write(*,301)' -gr   set override view scale (box size)'//
     $     ' for -gc, -gf plots.  [',rcij
      write(*,301)' -go   set No of orbits'
     $     //'(to plot on objects set by -gc). [',norbits
      write(*,301)' -at   set test angle.'
     $     //' -an   set No of angles. '
      write(*,301)' -ck   set checking timestep No. [',ickst
      write(*,301)' -h -?   Print usage.'
      write(*,301)' -ho     Print geomobj file format description'
      call exit(0)
      endif
      end