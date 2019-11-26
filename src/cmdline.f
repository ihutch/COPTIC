! Combination of the two calls to copticcmdline. 
      subroutine parametersetting
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
      implicit none
      integer ifull,ierr
      include 'myidcom.f'

      integer iobpl,iobpsw,ipstep,ifplot,norbits,nth,iavesteps
     $     ,ickst,ninjcomp,nsteps,nf_maxsteps,ndiags,ndiagmax
     $     ,iwstep,idistp,ndims,islp,lrestart,nptdiag,nqblkmax
      logical lmyidhead,ltestplot,lsliceplot,ldenplot,lphiplot,linjplot
     $     ,lextfield,LPF(ndims),lnotallp,ldistshow,LNPF
      real rcij,thetain,ripernode,crelax,colntime,dt,bdt,subcycle
     $     ,dropaccel,vneutral,debyelen,extfield,slpD ,Tneutral,Enfrac
     $     ,colpow,boltzamp
      real holepsi,holelen,holeum,holeeta,holepow,holerad,holegfac
      integer hspecies

      real wavespec(2*ndims+1)
      real Bfield(ndims),Bt,CFin(3+ndims,6)
      integer iCFcount,ipartperiod(ndims),idims(ndims)
      character*100 restartpath,objfilename
      character*256 argline
      real gt(ndims),gp0(ndims),gtt,gn(ndims),gnt
      real Ts(*),vds(*),eoverms(*),Tperps(*)
      integer nparta(*),numratioa(*)
      real vpars(*)
      real vperps(ndims,*),vdrifts(ndims,*)
      integer nspecies,nspeciesmax

!----------------------------------------------------------------------
! Deal with command-line arguments and geometry/object file.
! First time this routine just sets defaults and the object file name.
      call copticcmdline
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
     $     ,holegfac,wavespec,LNPF)
! Read in object file information.
      call readgeom(objfilename,myid,ifull,CFin,iCFcount,LPF,ierr
     $     ,argline)
      if(ierr.ne.0)stop 'Error in readgeom call'
! Second time: deal with any other command line parameters.
      call copticcmdline
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
     $     ,holegfac,wavespec,LNPF)
! The double call enables cmdline switches to override objfile settings.
!----------------------------------------------------------------------
      end
*********************************************************************
! Encapsulation of parameter setting.
      subroutine copticcmdline
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
     $     ,holegfac,wavespec,LNPF)

      implicit none

      integer iobpl,iobpsw,ipstep,ifplot,norbits,nth,iavesteps
     $     ,ickst,ninjcomp,nsteps,nf_maxsteps,ndiags,ndiagmax
     $     ,iwstep,idistp,ndims,islp,lrestart,nptdiag,nqblkmax
      logical lmyidhead,ltestplot,lsliceplot,ldenplot,lphiplot,linjplot
     $     ,lextfield,LPF(ndims),lnotallp,ldistshow,LNPF
      real rcij,thetain,ripernode,crelax,colntime,dt,bdt,subcycle
     $     ,dropaccel,vneutral,debyelen,extfield,slpD ,Tneutral,Enfrac
     $     ,colpow,boltzamp
     $     ,holelen,holepsi,holeum,holeeta,holepow,holerad,holegfac
      real wavespec(2*ndims+1)
      real Bfield(ndims),Bt,CFin(3+ndims,6)
      integer iCFcount,ipartperiod(ndims),idims(ndims)
      character*100 restartpath,objfilename
      character*256 argline
      real gt(ndims),gp0(ndims),gtt,gn(ndims),gnt
! Multiple species can have these things different.
! Only the first species is possibly collisional.
      real Ts(*),vds(*),eoverms(*),Tperps(*)
      integer nparta(*),numratioa(*)
      real vpars(*)
      real vperps(ndims,*),vdrifts(ndims,*)
      integer nspecies,nspeciesmax,hspecies

! Local variables:
      integer lentrim,iargc,ipfset
      external lentrim
      integer i,id,idn,idcn,i0,i1,iargcount,iargpos,ispecies
      real vwork,bdotgn,vdotgn
      character*100 argument,message
      logical lfirst
      data lfirst/.true./
      save argument,message

! Set defaults and objfilename only the first time, subsequently skip.
      if(lfirst)then
! Default (initial) number of particle species.
         nspecies=1
! Default to constant ripernode not npart.
         ripernode=100.
! Fixed number of particles rather than fixed injections.
         ninjcomp=0
         nparta(nspecies)=0
         vds(nspecies)=0.
         Ts(nspecies)=1.
         Tperps(nspecies)=1.
         eoverms(nspecies)=1.
         numratioa(nspecies)=1
         Tneutral=1.
         do id=1,ndims
            vdrifts(id,nspecies)=0.
         enddo
         vdrifts(ndims,nspecies)=1.
! Default edge-potential (chi) relaxation rate.     
!         crelax=1.*Ts(nspecies)/(1.+Ts(nspecies))
         crelax=0.
         debyelen=1.
         dt=.1
         restartpath=' '
         objfilename='copticgeom.dat'
         nsteps=5
         subcycle=0.
         dropaccel=10
         colntime=0.
         vneutral=0.
         Enfrac=0.
         colpow=0.
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
         iobpsw=1
         ldistshow=.false.
         boltzamp=1.
         holepsi=0.
         nqblkmax=0  ! Default 0=>pinit. >0 qinit
! Boundary condition switch and value. 0=> logarithmic.
         islp=0
         slpD=0.
         iCFcount=0
         do idn=1,2*ndims
            do id=1,3+ndims
               CFin(id,idn)=0.
            enddo
         enddo
         do id=1,ndims
            LPF(id)=.false.
            ipartperiod(id)=0
! Default zero field
            Bfield(id)=0.
         enddo
         LNPF=.false.
         gtt=0.
         gnt=0.
         do i=1,iargc()
            call getarg(i,argument)
            if(argument(1:3).eq.'-of')
     $           read(argument(4:),'(a)',err=501)objfilename
            if(argument(1:1).ne.'-')
     $           read(argument(1:),'(a)',err=501)objfilename
 501        continue
         enddo
         argline=' '
         message=' '
         lfirst=.false.
         return
      endif
! End of first time through setting
! ---------------------------------------
! Deal with arguments
      iargcount=iargc()
      iargpos=1
      do i=0,iargcount
! Start of argline internal iteration
 502     continue
         if(i.eq.0)then
! First time through, deal with argline arguments (from the objfile). 
            if(lentrim(argline(iargpos:)).ne.0)then
! Write them.
               if(lmyidhead.and.iargpos.eq.1)write(*,'(a,a)'
     $              )'File Arguments:',argline(iargpos:lentrim(argline))
               call argextract(argline,iargpos,argument)
            else
               goto 241
            endif
         else
! Subsequent: i>0 getarg.
            call getarg(i,argument)
         endif
!         write(*,*)i,argument
         if(argument(1:3).eq.'-gt')ltestplot=.true.
         if(argument(1:3).eq.'-gn')ldistshow=.not.ldistshow
         if(argument(1:3).eq.'-gc')read(argument(4:),*,end=201)iobpl
         if(argument(1:3).eq.'-gw')read(argument(4:),*,end=201)iobpsw
         if(argument(1:3).eq.'-gr')read(argument(4:),*,end=201)rcij
         if(argument(1:3).eq.'-gd'.or.argument(1:3).eq.'-gp')then
            lsliceplot=.true.
            read(argument(4:),*,err=211,end=211)ipstep
 211        continue
            if(lmyidhead)write(*,*)'Plotting ipstep=',ipstep
         endif
         if(argument(1:3).eq.'-gd')ldenplot=.not.ldenplot
         if(argument(1:3).eq.'-gp')lphiplot=.not.lphiplot
         if(argument(1:3).eq.'-gi')linjplot=.true.
         if(argument(1:3).eq.'-gf')read(argument(4:),*,end=201)ifplot
         if(argument(1:3).eq.'-go')read(argument(4:),*,end=201)norbits
         if(argument(1:3).eq.'-gg')call noeye3d(0)
         if(argument(1:3).eq.'-gx')then
            read(argument(4:),*,err=201,end=213)ipfset
            call pfset(ipfset)
            goto 240
 213        call pfset(-3)
         endif
         if(argument(1:3).eq.'-at')then
            read(argument(4:),*,end=201)thetain
         elseif(argument(1:3).eq.'-an')then
            read(argument(4:),*,end=201)nth
         elseif(argument(1:2).eq.'-a')then 
            read(argument(3:),*,end=201)iavesteps
         endif
         if(argument(1:3).eq.'-ni')read(argument(4:),*,err
     $        =201)nparta(nspecies)
         if(argument(1:3).eq.'-ri')read(argument(4:),*,end=201)ripernode
         if(argument(1:3).eq.'-rx')read(argument(4:),*,end=201)crelax
         if(argument(1:3).eq.'-ck')read(argument(4:),*,end=201)ickst
         if(argument(1:3).eq.'-ct')read(argument(4:),*,end=201)colntime
         if(argument(1:3).eq.'-Ef')read(argument(4:),*,end=201)Enfrac
         if(argument(1:3).eq.'-cp')read(argument(4:),*,end=201)colpow
         if(argument(1:3).eq.'-dt')read(argument(4:),*,end=201)dt
         if(argument(1:3).eq.'-da')read(argument(4:),*,end=201)bdt
         if(argument(1:3).eq.'-ds')read(argument(4:),*,end=201)subcycle
         if(argument(1:3).eq.'-dd')read(argument(4:),*,end=201)dropaccel
         if(argument(1:3).eq.'-zm')then
            read(argument(4:),*,err =201)eoverms(nspecies)
            numratioa(nspecies)=1
         endif
         if(argument(1:3).eq.'-nr')read(argument(4:),*,end=201)
     $        numratioa(nspecies)
         if(argument(1:3).eq.'-bc')read(argument(4:),*,end=201)islp
         if(argument(1:3).eq.'-bp')then
            read(argument(4:),*,end=201)idn
            if(0.lt.idn.and.idn.lt.4)then
               LPF(idn)=.not.LPF(idn)
               do id=idn,idn+ndims,ndims
                  CFin(1,id)=0.
                  CFin(2,id)=1.
               enddo
            endif
            iCFcount=iCFcount+1
         endif
         if(argument(1:3).eq.'-bf')then
            idn=-1
! Make sure at least the first parameter is readable and sensible
            read(argument(4:),*,end=201)idn
            if(idn.eq.0)then
! Reset all.
               do idn=1,2*ndims
                  do id=1,6
                     CFin(id,idn)=0.
                  enddo
               enddo
               iCFcount=0
            elseif(0.lt.idn .and. idn.le.2*ndims+1)then
! Read coefficients
               if(idn.eq.2*ndims+1)then
                  i0=1
                  i1=2*ndims
               else
                  i0=idn
                  i1=idn
               endif
               do idcn=i0,i1
! Initialize first, in case we are given a short read line.
                  do id=1,6
                     CFin(id,idcn)=0.
                  enddo
                  read(argument(4:),*,err=201,end=212)idn
     $                 ,(CFin(id,idcn),id=1,6)
 212              continue
! Diagnostics:
!                  write(*,*)'Set face',idcn,(CFin(id,idcn),id=1,6)
                  LPF(mod(idn-1,ndims)+1)=.false.
                  iCFcount=iCFcount+1
               enddo
            else
               write(*,*)'Error in bdyface cmdline switch:'
     $              ,argument(1:30)
               stop
            endif
         endif
         if(argument(1:3).eq.'-bn')LNPF=.true.
         if(argument(1:3).eq.'-pp')then
            read(argument(4:),*,err=201,end=201)ipartperiod
         endif
         if(argument(1:3).eq.'-id')then
            read(argument(4:),*,err=201,end=201)idims
         endif
         if(argument(1:3).eq.'-Bx')then
            read(argument(4:),*,end=201)Bfield(1)
         elseif(argument(1:3).eq.'-By')then
            read(argument(4:),*,end=201)Bfield(2)
         elseif(argument(1:3).eq.'-Bz')then
            read(argument(4:),*,end=201)Bfield(3)
         endif
         if(argument(1:7).eq.'--reinj')
     $        read(argument(8:),*,end=201)ninjcomp
         if(argument(1:3).eq.'-rn')
     $        read(argument(4:),*,end=201)ninjcomp
         if(argument(1:3).eq.'-sb')then
! Boltzamp setting
            read(argument(4:),*,end=201)boltzamp
         elseif(argument(1:4).eq.'-spb')then
! Boltzamp setting obsolete switch
            read(argument(5:),*,end=201)boltzamp
         elseif(argument(1:3).eq.'-sp')then
            if(nspecies+1.gt.nspeciesmax)then
               write(*,*)'***Disallowed more species than available'
     $              ,nspecies+1,nspeciesmax
               stop
            else
               nspecies=nspecies+1
! Default to electron Temp/mass for subsequent species.
               Ts(nspecies)=1.
               Tperps(nspecies)=1.
               eoverms(nspecies)=-1836.
               numratioa(nspecies)=1
               boltzamp=0.
! Inherit the drifts of the previous species until explicitly changed.
               vds(nspecies)=vds(nspecies-1)
               do id=1,ndims
                  vdrifts(id,nspecies)=vdrifts(id,nspecies-1)
               enddo
            endif
! Drifts must be respecified for this species if different.
         elseif(argument(1:2).eq.'-s')then
            read(argument(3:),*,end=201)nsteps
            if(nsteps.gt.nf_maxsteps)then
               if(lmyidhead)write(*,*)'Asked for more steps',nsteps,
     $              ' than nf_maxsteps-1=',nf_maxsteps-1
            endif
         endif
         if(argument(1:3).eq.'-vn')then
            read(argument(4:),*,end=201)vneutral
         elseif(argument(1:3).eq.'-vx')then
            read(argument(4:),*,end=201)vdrifts(1,nspecies)
         elseif(argument(1:3).eq.'-vy')then
            read(argument(4:),*,end=201)vdrifts(2,nspecies)
         elseif(argument(1:3).eq.'-vz')then
            read(argument(4:),*,end=201)vdrifts(3,nspecies)
         elseif(argument(1:2).eq.'-v')then
            read(argument(3:),*,end=201)vds(nspecies)
! For Mach bdy, set slpD equal to M.
            slpD=vds(nspecies)
! By default put the vneutral the same as first species
            vneutral=vds(1)
         endif
         if(argument(1:3).eq.'-md')then 
            read(argument(4:),*,end=201)ndiags
            if(ndiags.gt.ndiagmax)then
               write(*,*)'Error: Too many diag-moments',ndiags
               stop
            endif
         endif
         if(argument(1:2).eq.'-l')read(argument(3:),*,end=201)debyelen
         if(argument(1:3).eq.'-ng')then
            read(argument(4:),*,end=201)gn
            gnt=sqrt(gn(1)**2+gn(2)**2+gn(3)**2)
         endif
         if(argument(1:4).eq.'-tge')then
! Electron temperature gradient parameters
            read(argument(5:),*,end=201)gp0,gt
            gtt=sqrt(gt(1)**2+gt(2)**2+gt(3)**2)
         elseif(argument(1:3).eq.'-tn')then
            write(*,*)'***Tneutral setting is obsolete. No can do.'
            stop
!            read(argument(4:),*,end=201)Tneutral
         elseif(argument(1:3).eq.'-tp')then
            read(argument(4:),*,end=201)Tperps(nspecies)
         elseif(argument(1:2).eq.'-t')then
            read(argument(3:),*,end=201)Ts(nspecies)
            Tperps(nspecies)=Ts(nspecies)
! Default Tneutral=Ts(nspecies)
!            Tneutral=Ts(nspecies)
         endif
         
         if(argument(1:2).eq.'-w')read(argument(3:),*,end=201)iwstep
         if(argument(1:3).eq.'-pd')then
            idistp=1
            read(argument(4:),*,end=240,err=240)idistp
         endif
         if(argument(1:3).eq.'-pu')nptdiag=0
         if(argument(1:3).eq.'-pi')then
            if(holepsi.ne.0)then
               if(holelen.eq.4 .and. holepow.eq.0)then
                  if(lmyidhead)write(*,*)
     $                'Setting quieting after setting hole restores'
     $                ,' holelen and holepow defaults.' 
                  holelen=4.+holepsi/2.
                  holepow=-1/holepsi
               endif
            endif
            nqblkmax=30        ! If no argument still use qinit
            read(argument(4:),*,err=201,end=201)nqblkmax
         endif
         if(argument(1:3).eq.'-pw')then ! Wavedisplace initialization
            read(argument(4:),*,end=201,err=201)wavespec
         endif
         if(argument(1:3).eq.'-fs')then
            read(argument(4:),*,end=201)lrestart
         endif
         if(argument(1:3).eq.'-fn')then
            read(argument(4:),'(a)',end=201)restartpath
         endif
         if(argument(1:3).eq.'-fp')then
         endif         
         if(argument(1:10).eq.'--extfield')then
            read(argument(11:),*,end=201)extfield
!            write(*,*)'||||||||||||||extfield',extfield
            lextfield=.true.
         endif
         if(argument(1:3).eq.'-of')
     $        read(argument(4:),'(a)',end=201)objfilename
         if(argument(1:1).ne.'-')
     $        read(argument(1:),'(a)',end=201)objfilename
         if(argument(1:3).eq.'-ih')then
            ! Ensure we override prior settings.
            holepsi=0.
            holelen=999
            holeeta=2.   ! Needed default for trapinit.
            holeum=0.
            holepow=999
            holerad=0.
            holegfac=0.
            hspecies=nspecies
            read(argument(4:),*,err=201,end=231)holepsi,holeum,holelen
     $           ,holeeta,holepow,holerad,holegfac
 231        continue
            if(holepsi.ne.0)then ! We are setting 
               if(nqblkmax.gt.0)then ! qinit
                  if(holegfac.eq.0)then
                     if(holepsi.gt..6)then
                        if(holelen.eq.999)holelen=4.+holepsi/2.
                        if(holepow.eq.999)holepow=-1/holepsi !->toplen
                        if(lmyidhead)write(*,*)
     $                       'Deep Non-Debye Hole: wing scale',holelen
                     else
                        if(holelen.eq.999)holelen=4.   ! Exact Debye.
                        if(holepow.eq.999)holepow=-10. !->toplen
                     endif
                  else  ! Flattened distribution
                     if(holelen.eq.999)holelen=4. !+holepsi/10. Exact.
                     if(holepow.eq.999)holepow=
     $                    max(holepsi-1./holepsi,-10.)
                     if(lmyidhead)write(*,*)'Hole flattening gfac='
     $                    ,holegfac,' Normal value 0.6'
                  endif
               else               ! trapinit
                  if(holelen.eq.999)holelen=4.
                  if(holepow.eq.999)holepow=0.
               endif
            endif
         endif
         if(argument(1:3).eq.'-ho')then
            call geomdocument()
            call exit(0)
         endif
         if(argument(1:3).eq.'-hg')goto 402
! Help text options
         if(argument(1:2).eq.'-h')goto 203
         if(argument(1:3).eq.'--h')goto 203
         if(argument(1:2).eq.'-?')goto 203
! Indicator that coptic arguments are ended.
         if(argument(1:3).eq.'-ea'.or.argument(1:2).eq.'--')then
            write(*,*)'Arguments terminated by: ',argument(1:3)
            goto 202
         endif
! Finished parsing this argument.
 240     continue
         if(i.eq.0)goto 502
! End of internal argline iteration
 241     continue
      enddo
! End of command line parameter parsing.
!-------------------------------------------------------

 202  continue
      do ispecies=1,nspecies
         if(vdrifts(1,ispecies).ne.0. .or. vdrifts(2,ispecies).ne.0)then
! --- Drift in direction other than z. Normalize vdrifts cosines.
            vwork=0.
            do i=1,ndims
! Make all the vdrift components non-zero so that we can use that fact
! as an indicator of non-z drift.
!               if(vdrifts(i,ispecies).eq.0.)vdrifts(i,ispecies)=1.e-25
               vwork=vwork+vdrifts(i,ispecies)**2
            enddo
            vwork=sqrt(vwork)
            do i=1,ndims
               vdrifts(i,ispecies)=vdrifts(i,ispecies)/vwork
            enddo
         endif
      enddo
! --- Deal with B-field
      Bt=0.
      Bdotgn=0.
      vdotgn=0.
      do i=1,ndims
         Bt=Bt+Bfield(i)**2
         Bdotgn=Bdotgn+Bfield(i)*gn(i)
         vdotgn=vdotgn+vdrifts(i,1)*gn(i)
      enddo
      Bt=sqrt(Bt)
      if(gnt.ne.0.)then
         if(Bt.ne.0.)then
            if(Bdotgn.gt.1.e-5*Bt*gnt)then
! Complain if there's a B-component in the gradient direction.
               message='*** Density gradient parallel to Bt not allowed'
               goto 203
            endif
         else
            message='*** Zero B-field with density gradient not allowed'
            goto 203
         endif
         if(.not.vdotgn.eq.0)then
            message='*** V-drift component along density gradient'//
     $           ' not allowed.'
            goto 203
         endif
      endif
      if(Bt.ne.0)then
! Normalize magnetic field cosines.
         do i=1,ndims
            Bfield(i)=Bfield(i)/Bt
         enddo
         if(Bt.lt.1.e3 .and. Bt*dt.gt.1.)then
! Flag an inappropriate field and dt combination
            if(lmyidhead)write(message,'(a,f8.2,a,f8.5,a,a)')
     $           'UNWISE field',Bt,' and dt',dt
     $           ,' Use B.dt less than 1; else inaccurate.'
         endif
         do ispecies=1,nspecies
! vds(ispecies):
            vpars(ispecies)=0.
            do i=1,ndims
               vpars(ispecies)=vpars(ispecies)+vds(ispecies)*vdrifts(i
     $              ,ispecies)*Bfield(i)
            enddo
            do i=1,ndims
               vperps(i,ispecies)=-Bfield(i)*vpars(ispecies)
            enddo
            vwork=0.
            do i=1,ndims
               vperps(i,ispecies)=vperps(i,ispecies)+vds(ispecies)
     $              *vdrifts(i,ispecies)
               vwork=vwork+(vperps(i,ispecies)+vpars(ispecies)
     $              *Bfield(i))**2
            enddo
            vwork=sqrt(vwork)
            if(lmyidhead)write(*,'(a,f10.6,a,3f10.6,a,f10.6)'
     $           )'vpars(ispecies),vperp,vtot',vpars(ispecies),','
     $           ,(vperps(id,ispecies),id=1,3),',',vwork
         enddo
      else
! Zero the vparallel and vperp. Probably not necessary; but tidy.
         do ispecies=1,nspecies
            vpars(ispecies)=0.
            do i=1,ndims
               vperps(i,ispecies)=0.
            enddo
         enddo
      endif
! Set and check particle and potential periodicity logicals.
      lnotallp=.false.
      do i=1,ndims
         if(ipartperiod(i).ne.4)then
            lnotallp=.true.
         else
            if(crelax.ne.0. .and. bdt.ge.0.)then
! Don't allow unwise operation with periodic particles.
               write(*,*)'**** UNWISE operation with periodic particles'
     $              ,' and non-zero crelax=',crelax
               write(*,*)' Add switch -rx0.'
               stop
            endif
         endif
! Potential face non-periodic face existence logical.
! Removed now that ffttrid is available. 
!         if(i.ne.1)LNPF=LNPF.or..not.LPF(i) 
      enddo
! Consistency checks for holes
      if(holepsi.ne.0)then
         do i=1,ndims
            if(gn(i).ne.0)then ! bgn is not yet functional.
               if(lmyidhead)
     $        write(*,*)'Holes do not work with nonuniform background'
               stop
            endif
         enddo
      endif
      return
!------------------------------------------------------------
! Help text
 201  continue
      if(lmyidhead)write(*,*)'=====Error reading command line argument '
     $     ,argument(:20)
 203  continue
      if(.not.lmyidhead)return
 301  format(a,i5,a,i5)
 302  format(a,6f8.3)
 303  format(a,3L3,a)
 304  format(a,f8.3,a)
 305  format(a,3i3,a,3i5)
 306  format(a,6f7.3)
 307  format(a,6f8.1)
 308  format(a,6i8)
 309  format(a,7f6.2)
 310  format(a,i8,a,i8)
 311  format(a,7f6.2)
 312  format(a,L3,a)
      write(*,301)'Usage: coptic [objectfile] [-switches]'
      write(*,301)'Parameter switches.'
     $     //' Leave no gap before value. Defaults or set values [ddd'
      write(*,301)'[-of]<filename>  set name of object data file.'
     $     //'   ['//objfilename(1:30)
      write(*,310)' -ni   set No of particles/node   ['
     $     ,(nparta(ispecies),ispecies=1,1),'     zero => unset.'
      write(*,301)' -rn   set Reinjection number     [',ninjcomp
     $     ,'     fixed/step => parts/node unset.'
      write(*,304)' -ri   set rhoinfinity/node       [',ripernode
     $     ,'  => reinjection number.'
      write(*,304)' -rx   set Edge-potl relax rate   [',crelax
     $     ,'  0=>off, 1=>immed.'
      write(*,302)' -dt   set Timestep.              [',dt
      write(*,304)' -da   set Initial dt accel-factor[',bdt
     $     ,'  first 33% of timesteps larger.'
      write(*,'(3a)')'       Negative=> fractional density increase'
     $     ,' per unit time'
      write(*,304)' -ds   set Subcycle impulse/step  [',subcycle
     $     ,'  subcycling invoked above nonzero.'
      write(*,304)' -dd   set Drop-ion impulse/step  [',dropaccel
     $     ,'  greater impulse => drop this ion.'
      write(*,301)' -s    set No of steps.           [',nsteps
      write(*,302)' -t    set (Ion) Temperature.     [',(Ts(ispecies)
     $     ,ispecies=1,nspecies)
      write(*,302)' -tp   set Perp. Temperature.     ['
     $     ,(Tperps(ispecies),ispecies=1,nspecies)
      write(*,306)' -tge  set Elec Temp Center&Grad  [',gp0,gt
      write(*,306)' -ng   set Density gradient       [',gn
      write(*,302)' -l    set Debye Length.          [',debyelen
      write(*,302)' -v    set Drift speed.           [',(vds(ispecies)
     $     ,ispecies=1,nspecies)
      write(*,306)' -vx -vy -vz set velocity cosines [',((vdrifts(id
     $     ,ispecies),id=1,3),ispecies=1,min(2,nspecies))
      write(*,302)' -ct   set Collision time.        [',colntime
      write(*,302)' -vn   set Neutral drift velocity [',vneutral
      write(*,302)' -Ef   set Ext v-drive fraction   [',Enfrac
      write(*,302)' -cp   set v-power coln freq      [',colpow
      write(*,308)' -sp   add a Particle species     [',nspecies
      write(*,302)' -sb   Boltzmann fraction         [',boltzamp
      write(*,307)' -zm   set Z/mass ratio           ['
     $     ,(eoverms(ispecies),ispecies=1,nspecies)
      write(*,308)' -nr   set Species Number ratio   ['
     $     ,(numratioa(ispecies),ispecies=1,nspecies)
      write(*,302)' -Bx -By -Bz set Mag Field compts [',Bfield
      write(*,301)' -w    set Write-step period.     [',iwstep
     $     ,'     If <1, only myid writes final.'
      write(*,301)' -a    set Averaging steps.       [',iavesteps
     $     ,'     Also period of diagnostic writes.'
      write(*,301)' -pd   set Distribution Diags     [',idistp
     $     ,'     Bits:1 write, 2 plot.'
      write(*,301)' -pu   set Uniform p-distrib bins [',nptdiag
     $     ,'     nptdiag value. 0=>nsbins'
      write(*,301)' -md   set No of diag-moments(7). [',ndiags
      write(*,301)' -bc   set Boundary condition type[',islp
     $     ,'     E.g. 8201=1(1)+4(8)+13(8192)'
      write(*,301)
     $     '         Bits:1 Mach/Log[=0]  2-7:Face1-6=0  8-13:u''''=k2u'
!      write(*,301)' -xs<3reals>, -xe<3reals>  Set mesh start/end.'
      write(*,*)'-bf   set Face Bndry conditions  [ off'
     $     ,' Values:  idn,Ain,Bin,C0in[,Cx,Cy,Cz]'
      if(iCFcount.gt.0)write(*,'(10x,l3,i2,6f8.4)')(LPF(1+mod(idn-1,3))
     $     ,idn,(CFin(id,idn),id=1,6),idn=1,6)
      write(*,*)'        idn face-number (7=>all).'
     $     ,' ABC Robin coefs. Cxyz gradients.'
      write(*,303)' -bp<i>  toggle bndry Periodicity [',LPF
     $     ,'    in dimension <i>.'
      write(*,312)' -bn     defeat bndry nonPeriodcty[',LNPF
     $     ,'       When true, use no fftw solver'
      write(*,301)' -pi<i>  set Quiet part-init level[',nqblkmax
     $     ,'     1:no quieting, >>1:quiet'
      write(*,311)' -pw[..] apply initial Wave 1,n,xi[',wavespec
      write(*,305)' -pp<i,j,k>  Partcl bcs/periodcty [',ipartperiod
     $     ,'    Use sum of:'
      write(*,*)'     0 open; 1 lower, 2 upper, 3 both absorbing;'
     $     ,' 4 periodic, 5 v-reset'
      write(*,*)'     +Domain end between nodes (upper byte)',
     $     ' 64:lower, 128:upper, 192:both.'
      write(*,305)' -id<i,j,k> set MPI block dims    [',idims
      write(*,309)' -ih<P>[..] Hole Psi[u,l,h,p,r,g  [',holepsi,holeum
     $     ,holelen,holeeta,holepow,holerad,holegfac
      write(*,301)' -fs<i>  set Restart switch:      [',lrestart
     $     ,' bit1:partls+potl, bit2:flux, bit3:name'

      write(*,'(a,a)') ' -fn[path]  set particle reading/writing Path: '
     $     ,restartpath(1:lentrim(restartpath))
      write(*,301)' -ea --  end argument parsing. Skip succeeding.'
      goto 401
 402  continue
      write(*,301)'Debugging switches for testing and sequential plots'
      write(*,301)' -gt   Plot regions and solution tests.'
      write(*,301)' -gi   Plot injection accumulated diagnostics.'
      write(*,*)'-gn   Plot collisional reinjection distribution.'
     $     ,' Or 1-d phase space (clnless)'
      write(*,301)' -gg   Run continuously without pausing.'
      write(*,301)' -gx[i]Graphics filing [alone if i<0 or absent].'
      write(*,301)
     $      ' -gp -gd[] Plot slices of setup; plus potential, '
     $     //'density. [At step n]. [',ipstep
      write(*,301)' -gf   set quantity plotted for flux evolution and'//
     $     ' final distribution. [',ifplot
      write(*,301)' -gw   set objplot sw. [+256:intercepts]'//
     $     ' Shade by 1:flux 2:flux-density[',iobpsw
      write(*,301)' To plot intersections with objects use -gf and -gw'
      write(*,301)' -gc   set wireframe [& stencils(-)] mask.'//
     $     ' objects<->bits. [',iobpl
      write(*,302)' -gr   set override view scale (box size)'//
     $     ' for -gc, -gf plots.  [',rcij
      write(*,301)' -go   set No of orbits'
     $     //'(to plot on objects set by -gc). [',norbits
      write(*,301)' -at   set test angle.'
     $     //' -an   set No of angles. '
      write(*,301)' -ck   set checking timestep No. [',ickst
      write(*,'(a)')
     $ 'Examples of 1-d output controls.'
     $ ,'  -gn outputs just a pps file every step, no display.'
     $ ,'  -gp outputs a pps file and displays every -a steps.'
     $ ,'  -gn -gp outputs a pps file and displays every step.'
     $ ,'  -gn -gp5 outputs a pps file and displays every 5 steps.'
     $ ,'  -gn -gp -gx outputs a pps file and a ps file each step'
     $     //' but does not display.'
     $ ,'  -gn -gp2 -gx3 outputs pps and ps files'
     $     //' and displays every 2 steps.'
 401  write(*,301)' -h -?   Print usage.'
      write(*,301)' -hg     Print debugging/plotting switch usage.'
      write(*,301)' -ho     Print geomobj file format description'
      if(lentrim(message).gt.1)write(*,'(a)')message(1:lentrim(message))
      call exit(0)
      end
!*********************************************************************
      subroutine argextract(argline,ipos,argstring)
! Extract an argument from the string starting at ipos and return it in
! argstring, leaving ipos pointed at the terminating blank.
      character*(*) argline,argstring

      ias=0
 2    continue
! Skip blanks
      if(argline(ipos:ipos).eq.' ')then
         ipos=ipos+1
         goto 2
      endif
 3    continue
! Read until we reach a blank.
      ias=ias+1
      argstring(ias:)=argline(ipos:ipos)
      if(argline(ipos:ipos).ne.' ')then
         ipos=ipos+1
         goto 3
      endif
      end
