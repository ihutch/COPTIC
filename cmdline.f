*********************************************************************
c Encapsulation of parameter setting.
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
     $     ,nspeciesmax,numratioa,Tperps,boltzamp,nptdiag)

      implicit none

      integer iobpl,iobpsw,ipstep,ifplot,norbits,nth,iavesteps
     $     ,ickst,ninjcomp,nsteps,nf_maxsteps,ndiags,ndiagmax
     $     ,iwstep,idistp,ndims,islp,lrestart,nptdiag
      logical lmyidhead,ltestplot,lsliceplot,ldenplot,lphiplot,linjplot
     $     ,lextfield,LPF(ndims),lnotallp,ldistshow
      real rcij,thetain,ripernode,crelax,colntime,dt,bdt,subcycle
     $     ,dropaccel,vneutral,debyelen,extfield,slpD
     $     ,Tneutral,Enfrac,colpow,boltzamp
      real Bfield(ndims),Bt,CFin(3+ndims,6)
      integer iCFcount,ipartperiod(ndims),idims(ndims)
      character*100 restartpath,objfilename
      character*256 argline
      real gt(ndims),gp0(ndims),gtt,gn(ndims),gnt
c Multiple species can have these things different.
c Only the first species is possibly collisional.
      real Ts(*),vds(*),eoverms(*),Tperps(*)
      integer nparta(*),numratioa(*)
      real vpars(*)
      real vperps(ndims,*),vdrifts(ndims,*)
      integer nspecies,nspeciesmax

c Local variables:
      integer lentrim,iargc
      external lentrim
      integer i,id,idn,idcn,i0,i1,iargcount,iargpos,ispecies
      real vwork,bdotgn
      character*100 argument,message
      logical lfirst
      data lfirst/.true./
      save argument,message

c Set defaults and objfilename only the first time, subsequently skip.
      if(lfirst)then
c Default (initial) number of particle species.
         nspecies=1
c Default to constant ripernode not npart.
         ripernode=100.
c Fixed number of particles rather than fixed injections.
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
c Default edge-potential (chi) relaxation rate.     
c         crelax=1.*Ts(nspecies)/(1.+Ts(nspecies))
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
         boltzamp=0.
c Boundary condition switch and value. 0=> logarithmic.
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
         enddo
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
c End of first time through setting
c ---------------------------------------
c Deal with arguments
      iargcount=iargc()
      iargpos=1
      do i=0,iargcount
c Start of argline internal iteration
 502     continue
         if(i.eq.0)then
c First time through, deal with argline arguments (from the objfile). 
            if(lentrim(argline(iargpos:)).ne.0)then
c Write them.
               if(lmyidhead.and.iargpos.eq.1)write(*,'(a,a)'
     $              )'File Arguments:',argline(iargpos:lentrim(argline))
               call argextract(argline,iargpos,argument)
            else
               goto 241
            endif
         else
c Subsequent: i>0 getarg.
            call getarg(i,argument)
         endif
c         write(*,*)i,argument
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
         if(argument(1:3).eq.'-gf')read(argument(4:),*,err=201)ifplot
         if(argument(1:3).eq.'-go')read(argument(4:),*,err=201)norbits
         if(argument(1:3).eq.'-at')then
            read(argument(4:),*,err=201)thetain
         elseif(argument(1:3).eq.'-an')then
            read(argument(4:),*,err=201)nth
         elseif(argument(1:2).eq.'-a')then 
            read(argument(3:),*,err=201)iavesteps
         endif
         if(argument(1:3).eq.'-ni')read(argument(4:),*,err
     $        =201)nparta(nspecies)
         if(argument(1:3).eq.'-ri')read(argument(4:),*,err=201)ripernode
         if(argument(1:3).eq.'-rx')read(argument(4:),*,err=201)crelax
         if(argument(1:3).eq.'-ck')read(argument(4:),*,err=201)ickst
         if(argument(1:3).eq.'-ct')read(argument(4:),*,err=201)colntime
         if(argument(1:3).eq.'-Ef')read(argument(4:),*,err=201)Enfrac
         if(argument(1:3).eq.'-cp')read(argument(4:),*,err=201)colpow
         if(argument(1:3).eq.'-dt')read(argument(4:),*,err=201)dt
         if(argument(1:3).eq.'-da')read(argument(4:),*,err=201)bdt
         if(argument(1:3).eq.'-ds')read(argument(4:),*,err=201)subcycle
         if(argument(1:3).eq.'-dd')read(argument(4:),*,err=201)dropaccel
         if(argument(1:3).eq.'-zm')then
            read(argument(4:),*,err =201)eoverms(nspecies)
            numratioa(nspecies)=1
c            numratioa(nspecies)=sqrt(abs(eoverms(nspecies)/eoverms(1)))
         endif
         if(argument(1:3).eq.'-nr')read(argument(4:),*,err=201)
     $        numratioa(nspecies)
         if(argument(1:3).eq.'-bc')read(argument(4:),*,err=201)islp
         if(argument(1:3).eq.'-bp')then
            read(argument(4:),*,err=201)idn
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
c Make sure at least the first parameter is readable and sensible
            read(argument(4:),*,err=201)idn
            if(idn.eq.0)then
c Reset all.
               do idn=1,2*ndims
                  do id=1,6
                     CFin(id,idn)=0.
                  enddo
               enddo
               iCFcount=0
            elseif(0.lt.idn .and. idn.le.2*ndims+1)then
c Read coefficients
               if(idn.eq.2*ndims+1)then
                  i0=1
                  i1=2*ndims
               else
                  i0=idn
                  i1=idn
               endif
               do idcn=i0,i1
c Initialize first, in case we are given a short read line.
                  do id=1,6
                     CFin(id,idcn)=0.
                  enddo
                  read(argument(4:),*,err=201,end=212)idn
     $                 ,(CFin(id,idcn),id=1,6)
 212              continue
c Diagnostics:
c                  write(*,*)'Set face',idcn,(CFin(id,idcn),id=1,6)
                  LPF(mod(idn-1,ndims)+1)=.false.
                  iCFcount=iCFcount+1
               enddo
            else
               write(*,*)'Error in bdyface cmdline switch:'
     $              ,argument(1:30)
               stop
            endif
         endif
         if(argument(1:3).eq.'-pp')then
            read(argument(4:),*,err=201,end=201)ipartperiod
         endif
         if(argument(1:3).eq.'-id')then
            read(argument(4:),*,err=201,end=201)idims
         endif
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
         if(argument(1:4).eq.'-spb')then
c Boltzamp setting
            read(argument(5:),*,err=201)boltzamp
         elseif(argument(1:3).eq.'-sp')then
            if(nspecies+1.gt.nspeciesmax)then
               write(*,*)'***Disallowed more species than available'
     $              ,nspecies+1,nspeciesmax
               stop
            else
               nspecies=nspecies+1
c Default to electron Temp/mass for subsequent species.
               Ts(nspecies)=1.
               Tperps(nspecies)=1.
               eoverms(nspecies)=-1836.
               numratioa(nspecies)=1
c     $              sqrt(abs(eoverms(nspecies)/eoverms(1)))
c Inherit the drifts of the previous species until explicitly changed.
               vds(nspecies)=vds(nspecies-1)
               do id=1,ndims
                  vdrifts(id,nspecies)=vdrifts(id,nspecies-1)
               enddo
            endif
c Drifts must be respecified for this species if different.
         elseif(argument(1:2).eq.'-s')then
            read(argument(3:),*,err=201)nsteps
            if(nsteps.gt.nf_maxsteps)then
               if(lmyidhead)write(*,*)'Asked for more steps',nsteps,
     $              ' than allowed. Limit ',nf_maxsteps-1
               nsteps=nf_maxsteps-1
            endif
         endif
         if(argument(1:3).eq.'-vn')then
            read(argument(4:),*,err=201)vneutral
         elseif(argument(1:3).eq.'-vx')then
            read(argument(4:),*,err=201)vdrifts(1,nspecies)
         elseif(argument(1:3).eq.'-vy')then
            read(argument(4:),*,err=201)vdrifts(2,nspecies)
         elseif(argument(1:3).eq.'-vz')then
            read(argument(4:),*,err=201)vdrifts(3,nspecies)
         elseif(argument(1:2).eq.'-v')then
            read(argument(3:),*,err=201)vds(nspecies)
c For Mach bdy, set slpD equal to M.
            slpD=vds(nspecies)
c By default put the vneutral the same as first species
            vneutral=vds(1)
         endif
         if(argument(1:3).eq.'-md')then 
            read(argument(4:),*,err=201)ndiags
            if(ndiags.gt.ndiagmax)then
               write(*,*)'Error: Too many diag-moments',ndiags
               stop
            endif
         endif
         if(argument(1:2).eq.'-l')read(argument(3:),*,err=201)debyelen
         if(argument(1:3).eq.'-ng')then
            read(argument(4:),*,err=201)gn
            gnt=sqrt(gn(1)**2+gn(2)**2+gn(3)**2)
         endif
         if(argument(1:4).eq.'-tge')then
c Electron temperature gradient parameters
            read(argument(5:),*,err=201)gp0,gt
            gtt=sqrt(gt(1)**2+gt(2)**2+gt(3)**2)
         elseif(argument(1:3).eq.'-tn')then
            write(*,*)'***Tneutral setting is obsolete. No can do.'
            stop
c            read(argument(4:),*,err=201)Tneutral
         elseif(argument(1:3).eq.'-tp')then
            read(argument(4:),*,err=201)Tperps(nspecies)
         elseif(argument(1:2).eq.'-t')then
            read(argument(3:),*,err=201)Ts(nspecies)
            Tperps(nspecies)=Ts(nspecies)
c Default Tneutral=Ts(nspecies)
c            Tneutral=Ts(nspecies)
         endif
         
         if(argument(1:2).eq.'-w')read(argument(3:),*,err=201)iwstep
         if(argument(1:3).eq.'-pd')then
            idistp=1
            read(argument(4:),*,err=240)idistp
         endif
         if(argument(1:3).eq.'-pu')nptdiag=0
         if(argument(1:3).eq.'-fs')then
            read(argument(4:),*,err=201)lrestart
         endif
         if(argument(1:3).eq.'-fn')then
            read(argument(4:),'(a)',err=201)restartpath
         endif
         if(argument(1:3).eq.'-fp')then
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
         if(argument(1:3).eq.'-hg')goto 402
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
         if(i.eq.0)goto 502
c End of internal argline iteration
 241     continue
      enddo
c End of command line parameter parsing.
c-------------------------------------------------------

 202  continue
      do ispecies=1,nspecies
         if(vdrifts(1,ispecies).ne.0. .or. vdrifts(2,ispecies).ne.0)then
c --- Drift in direction other than z. Normalize vdrifts cosines.
            vwork=0.
            do i=1,ndims
c Make all the vdrift components non-zero so that we can use that fact
c as an indicator of non-z drift.
c               if(vdrifts(i,ispecies).eq.0.)vdrifts(i,ispecies)=1.e-25
               vwork=vwork+vdrifts(i,ispecies)**2
            enddo
            vwork=sqrt(vwork)
            do i=1,ndims
               vdrifts(i,ispecies)=vdrifts(i,ispecies)/vwork
            enddo
         endif
      enddo
c ---
      if(colntime.eq.0.)then
c Collisionless, also set vneutral to vd, else things are inconsistent:
c But vneutral is only the first species. This is redundant.
c         vneutral=vds(1)
      endif
c --- Deal with B-field
      Bt=0.
      Bdotgn=0.
      do i=1,ndims
         Bt=Bt+Bfield(i)**2
         Bdotgn=Bdotgn+Bfield(i)*gn(i)
      enddo
      Bt=sqrt(Bt)
      if(gnt.ne.0.)then
         if(Bt.ne.0.)then
            if(Bdotgn.gt.1.e-5*Bt*gnt)then
c Complain if there's a B-component in the gradient direction.
               message='**Density gradient parallel to Bt not allowed'
               goto 203
            endif
         else
            message='**Zero B-field with density gradient not allowed'
            goto 203
         endif
      endif
      if(Bt.ne.0)then
c Normalize magnetic field cosines.
         do i=1,ndims
            Bfield(i)=Bfield(i)/Bt
         enddo
         if(Bt.lt.1.e3 .and. Bt*dt.gt.1.)then
c Flag an inappropriate field and dt combination
            if(lmyidhead)write(message,'(a,f8.2,a,f8.5,a,a)')
     $           'UNWISE field',Bt,' and dt',dt
     $           ,' Use B.dt less than 1; else inaccurate.'
         endif
         do ispecies=1,nspecies
c This test is redundant now. Should use the general (else) case.
c            if(vdrifts(1,ispecies).eq.0)then
            if(.false.)then
c Assume that vds(ispecies) is in the z-direction
               vpars(ispecies)=vds(ispecies)*Bfield(ndims)
               do i=1,ndims
                  vperps(i,ispecies)=-Bfield(i)*vpars(ispecies)
               enddo
            vperps(ndims,ispecies)=vperps(ndims,ispecies)+vds(ispecies)
            else
c vds(ispecies) non-z:
               vpars(ispecies)=0.
               do i=1,ndims
                  vpars(ispecies)=vds(ispecies)*vdrifts(i,ispecies)
     $                 *Bfield(i)
               enddo
               do i=1,ndims
                  vperps(i,ispecies)=-Bfield(i)*vpars(ispecies)
               enddo
               vwork=0.
               do i=1,ndims
                  vperps(i,ispecies)=vperps(i,ispecies)+vds(ispecies)
     $                 *vdrifts(i,ispecies)
                  vwork=vwork+(vperps(i,ispecies)+vpars(ispecies)
     $                 *Bfield(i))**2
               enddo
               vwork=sqrt(vwork)
               if(lmyidhead)write(*,'(a,f10.6,a,3f10.6,a,f10.6)'
     $              )'vpars(ispecies),vperp,vtot',vpars(ispecies),','
     $              ,(vperps(id,ispecies),id=1,3),',',vwork
            endif
         enddo
      else
c Zero the vparallel and vperp. Probably not necessary; but tidy.
         do ispecies=1,nspecies
            vpars(ispecies)=0.
            do i=1,ndims
               vperps(i,ispecies)=0.
            enddo
         enddo
      endif
c Set and check particle periodicity logical.
      lnotallp=.false.
      do i=1,ndims
         if(ipartperiod(i).ne.4)then
            lnotallp=.true.
         else
            if(crelax.ne.0.)then
c Don't allow unwise operation with periodic particles.
               write(*,*)'**** UNWISE operation with periodic particles'
     $              ,' and non-zero crelax=',crelax
               write(*,*)' Add switch -rx0.'
               stop
            endif
         endif
      enddo

      return
c------------------------------------------------------------
c Help text
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
      write(*,301)'Usage: coptic [objectfile] [-switches]'
      write(*,301)'Parameter switches.'
     $     //' Leave no gap before value. Defaults or set values [ddd'
      write(*,301)'[-of]<filename>  set name of object data file.'
     $     //'   ['//objfilename(1:30)
      write(*,301)' -ni   set No of particles/node   ['
     $     ,(nparta(ispecies),ispecies=1,1),'     zero => unset.'
      write(*,301)' -rn   set reinjection number     [',ninjcomp
     $     ,'     fixed/step => parts/node unset.'
      write(*,304)' -ri   set rhoinfinity/node       [',ripernode
     $     ,'  => reinjection number.'
      write(*,304)' -rx   set Edge-potl relax rate   [',crelax
     $     ,'  0=>off, 1=>immed.'
      write(*,302)' -dt   set Timestep.              [',dt
      write(*,304)' -da   set Initial dt accel-factor[',bdt
     $     ,'  first 33% of timesteps larger.'
      write(*,304)' -ds   set Subcycle impulse/step  [',subcycle
     $     ,'  subcycling invoked above nonzero.'
      write(*,304)' -dd   set Drop-ion impulse/step  [',dropaccel
     $     ,'  greater impulse => drop this ion.'
      write(*,301)' -s    set No of steps.           [',nsteps
      write(*,302)' -t    set (Ion) Temperature.     [',(Ts(ispecies)
     $     ,ispecies=1,nspecies)
      write(*,302)' -tp   set perp. Temperature.     ['
     $     ,(Tperps(ispecies),ispecies=1,nspecies)
      write(*,306)' -tge  set Elec Temp Center&Grad  [',gp0,gt
      write(*,306)' -ng   set Density gradient       [',gn
      write(*,302)' -l    set Debye Length.          [',debyelen
      write(*,302)' -v    set drift speed.           [',(vds(ispecies)
     $     ,ispecies=1,nspecies)
      write(*,306)' -vx -vy -vz set velocity cosines [',((vdrifts(id
     $     ,ispecies),id=1,3),ispecies=1,min(2,nspecies))
      write(*,302)' -ct   set collision time.        [',colntime
      write(*,302)' -vn   set neutral drift velocity [',vneutral
      write(*,302)' -Ef   set Ext v-drive fraction   [',Enfrac
      write(*,302)' -cp   set v-power coln freq      [',colpow
      write(*,308)' -sp   add a particle species     [',nspecies
      write(*,302)' -spb  extra Boltzmann fraction   [',boltzamp
      write(*,307)' -zm   set Z/mass ratio           ['
     $     ,(eoverms(ispecies),ispecies=1,nspecies)
      write(*,308)' -nr   set species number ratio   ['
     $     ,(numratioa(ispecies),ispecies=1,nspecies)
      write(*,302)' -Bx -By -Bz set mag field compts [',Bfield
      write(*,301)' -w    set write-step period.     [',iwstep
     $     ,'     If <1, only myid writes final.'
      write(*,301)' -a    set averaging steps.       [',iavesteps
     $     ,'     Also period of diagnostic writes.'
      write(*,301)' -pd   set distribution diags     [',idistp
     $     ,'     Bits:1 write, 2 plot.'
      write(*,301)' -pu   set uniform p-distrib bins [',nptdiag
     $     ,'     nptdiag value. 0=>nsbins'
      write(*,301)' -md   set No of diag-moments(7). [',ndiags
      write(*,301)' -bc   set boundary condition type[',islp
     $     ,'     E.g. 8201=1(1)+4(8)+13(8192)'
      write(*,301)
     $     '         Bits:1 Mach/Log[=0]  2-7:Face1-6=0  8-13:u''''=k2u'
c      write(*,301)' -xs<3reals>, -xe<3reals>  Set mesh start/end.'
      write(*,*)'-bf   set face bndry conditions  [ off'
     $     ,' Values:  idn,Ain,Bin,C0in[,Cx,Cy,Cz]'
      if(iCFcount.gt.0)write(*,'(10x,l3,i2,6f8.4)')(LPF(1+mod(idn-1,3))
     $     ,idn,(CFin(id,idn),id=1,6),idn=1,6)
      write(*,*)'        idn face-number (7=>all).'
     $     ,' ABC Robin coefs. Cxyz gradients.'
      write(*,303)' -bp<i>  toggle bndry periodicity [',LPF
     $     ,'    in dimension <i>.'
      write(*,305)' -pp<i,j,k>  partcl bcs/periodcty [',ipartperiod 
      write(*,*)'     0 open; 1 lower absorbing; 2 upper absorbing;'
     $     ,' 3 both absorb; 4 periodic'
      write(*,305)' -id<i,j,k> Set MPI block dims    [',idims
      write(*,301)' -fs<i>  set restart switch:      [',lrestart
     $     ,' bit1:partls+potl, bit2:flux, bit3:name'

      write(*,'(a,a)') ' -fn[path]  set particle reading/writing path: '
     $     ,restartpath(1:lentrim(restartpath))
      write(*,301)' -ea --  end argument parsing. Skip succeeding.'
      goto 401
 402  continue
      write(*,301)'Debugging switches for testing'
      write(*,301)' -gt   Plot regions and solution tests.'
      write(*,301)' -gi   Plot injection accumulated diagnostics.'
      write(*,301)' -gn   Plot collisional reinjection distribution.'
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
 401  write(*,301)' -h -?   Print usage.'
      write(*,301)' -hg     Print debugging/plotting switch usage.'
      write(*,301)' -ho     Print geomobj file format description'
      if(lentrim(message).gt.1)write(*,'(a)')message(1:lentrim(message))
      call exit(0)
      end
c*********************************************************************
      subroutine argextract(argline,ipos,argstring)
c Extract an argument from the string starting at ipos and return it in
c argstring, leaving ipos pointed at the terminating blank.
      character*(*) argline,argstring

      ias=0
 2    continue
c Skip blanks
      if(argline(ipos:ipos).eq.' ')then
         ipos=ipos+1
         goto 2
      endif
 3    continue
c Read until we reach a blank.
      ias=ias+1
      argstring(ias:)=argline(ipos:ipos)
      if(argline(ipos:ipos).ne.' ')then
         ipos=ipos+1
         goto 3
      endif
      end
