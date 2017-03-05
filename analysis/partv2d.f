      program part2d
! Examine the particle data, showing distribution function(s)
! in 2d. A hacked program at present mostly to check the magnetized
! collisional distributions. Much came from partexamine.
      include 'examdecl.f'
! (Examdecl itself includes meshcom.f plascom.f, objcom.f)
      parameter (nfilemax=999)
      include '../src/partcom.f'
      include '../src/ptaccom.f'

 
      character*10 chartemp
      character*100 name
      character*100 string

      integer nx,ny
      parameter (nxbin=50,nybin=50)
      real fvxvy(nxbin,nybin),vxa(nxbin,nybin),vya(nxbin,nybin)
      real work(0:nxbin+1,0:nybin+1)
      integer nl
      parameter (nl=10)
      real cl(nl),ht
      character pp(4,nxbin,nybin)

      logical ldoc
!      real extra(nptdiag,mdims),diff(nptdiag,mdims)
! Spatial limits bottom-top, dimensions
      real xlimit(2,mdims),xnewlim(2,mdims)
      real Bdirs(mdims+1)
! Velocity limits
      real vlimit(2,mdims)
      integer iuin(mdims)

      nfmax=nfilemax
! silence warnings:
      zp(1,1,1)=0.

! Defaults
      do id=1,mdims
! Use very big xlimits by default to include whole domain
! They are then reset by the accumulation itself.
         xlimit(1,id)=-500.
         xlimit(2,id)=500.
         xnewlim(1,id)=0.
         xnewlim(2,id)=0.
         vlimit(1,id)=3.
         vlimit(2,id)=-2.
      enddo

      call partexamargs(xlimit,vlimit,iuin,cellvol,Bdirs,ldoc)
      do id=1,mdims
! Needed initialization removed from partacinit.
         xmeshstart(id)=min(-5.,xlimit(1,id))
         xmeshend(id)=max(5.,xlimit(2,id))
         isuds(id)=iuin(id)
      enddo
      write(*,*)' isuds:',isuds
      if(isuds(1)*isuds(2)*isuds(3).gt.nsub_tot)then
         write(*,*)'Too many blocks requested. Reduce or recompile.'
         stop
      endif
      write(*,*)'nptdiag,nsbins',nptdiag,nsbins
      if(nptdiag.lt.nsbins)then
         write(*,*)'WARNING nptdiags=',nptdiag,' smaller than nsbins='
     $        ,nsbins,' Incorrect array size choice.'
      endif
! Now the base filename is in partfilename.
      ip=lentrim(partfilename)-3
      if(partfilename(ip:ip).eq.'.')then
! If filename has a 3-character extension or is a numbered pex
! file. Assume it is complete and that we are reading just one
! file. Should do this only the first time.
         nfmax=0
         name=partfilename
         write(*,*)'Reading single file ',name(1:lentrim(name))
         if(partfilename(ip:ip+3).eq.'.pex')then
            write(*,*)'Using stored distribution file'
            nfmax=-1
         endif
      elseif(partfilename(ip-4:ip-1).eq.'.pex')then
         nfmax=-1
         name=partfilename
         write(*,*)'Reading numbered pex file ',name(1:lentrim(name))
      endif

! Possible multiple files.
      do i=0,nfmax
         if(nfmax.ne.0)then
            write(chartemp,'(''.'',i3.3)')i
            name=partfilename(1:lentrim(partfilename))//chartemp
            write(*,*)'Reading file ',name(1:lentrim(name))
         endif
         Bt=0.
         call partread(name,ierr)
         if(ldoc)then
            write(*,*)'debyelen,Ti,vd,rs,phip=',debyelen,Ti,vd,rs,phip
            write(*,*)'iregion_part,n_part,dt,ldiags,rhoinf,nrein,',
     $           'phirein,numprocs='
            write(*,*)iregion_part,n_part,dt,ldiags,rhoinf,nrein,phirein
     $           ,numprocs
            write(*,*)'eoverm,Bt,Bfield,vpar,vperp=',eoverm,Bt,Bfield
     $           ,vpar,vperp
            stop
         endif
         if(ierr-4*(ierr/4).ne.0)goto 11
         if(Bdirs(4).gt.0. .or. Bt.eq.0)then
! All directions were set by commandline. Or none were read from file.
            do k=1,ndims
               Bfield(k)=Bdirs(k)
            enddo
         endif
         if(cellvol.eq.-1)write(*,*)'Bfield (projection)',Bfield
! Accumulate the fvxvy distribution
         do k=1,nybin
            do j=1,nxbin
               fvxvy(j,k)=0.
            enddo
         enddo
         do j=1,ioc_part
            if(x_part(iflag,j).ne.0)then
               vx=x_part(1+npdim,j)
               vy=x_part(3+npdim,j)
               nx=nint(nxbin*(vx-vlimit(1,1))/(vlimit(2,1)-vlimit(1,1)))
               ny=nint(nybin*(vy-vlimit(1,2))/(vlimit(2,2)-vlimit(1,2)))
               if(nx.ge.1.and.nx.le.nxbin.and.ny.ge.1.and.ny.le.nybin)
     $              then
                  fvxvy(nx,ny)=fvxvy(nx,ny)+1.
               endif
            endif
         enddo

      enddo
      do i=1,nybin
         write(*,'(80i1)')(int(2.*alog10(fvxvy(j,i))),j=1,nxbin)
      enddo


      do i=1,nxbin
      do j=1,nybin
         vxa(i,j)=vlimit(1,1)+(vlimit(2,1)-vlimit(1,1))*(i-1.)/(nxbin
     $        -1.)
         vya(i,j)=vlimit(1,2)+(vlimit(2,2)-vlimit(1,2))*(j-1.)/(nybin
     $        -1.)
      enddo
      enddo
!      write(*,*)vxa
98     call pltinit(0.,1.,0.,1.)
!       Plot the surface. With axes (2-D). Web color 10, axis color 7.
        j=2 + 256*10 + 256*256*7
        call surf3d(vxa,vya,fvxvy,nxbin,nxbin,nybin,j,work)
!       Draw a contour plot in perspective. Need to reset color anyway.
        call color(4)
        call axregion(-scbx3,scbx3,-scby3,scby3)
        call scalewn(vxa(1,1),vxa(nxbin,1),vya(1,1),vya(1,nybin),.false.
     $       ,.false.)
        call hdprset(-3,scbz3)
!       write(*,*) 'Done hdprset'
!       call axis   ! if desired.
        call scalewn(1.,float(nxbin),1.,float(nybin),.false.,.false.)
!       Contour without labels, direct on mesh.
!       write(*,*) 'Done scalewn'
! Not yet.
!        call contourl(fvxvy,pp,nxbin,nxbin,nybin,cl,nl,vxa,vya,16)
        call color(15)
        if(ieye3d().ne.0) goto 98


 11   continue



      end
!*************************************************************
      subroutine partexamargs(xlimit,vlimit,iuin,cellvol,Bdirs,ldoc)
      include 'examdecl.f'
      real xlimit(2,3),vlimit(2,3),Bdirs(4)
      integer iuin(3)
      logical ldoc

! I think unused here 26 May 12. But I'm not sure.
      ifull(1)=na_i
      ifull(2)=na_j
      ifull(3)=na_k

      do i=1,3
         iuin(i)=9
! convenient default field:
         Bdirs(i)=2-i
      enddo
! Use cellvol=0. by default.
      cellvol=0.
      ldoc=.false.

! silence warnings:
      fluxfilename=' '
      zp(1,1,1)=0.
! Defaults
      partfilename=' '

! Deal with arguments
      if(iargc().eq.0) goto 201
      do i=1,iargc()
         call getarg(i,argument)
         if(argument(1:1).eq.'-')then
            if(argument(1:2).eq.'-x')then
               read(argument(3:),*,err=201) xlimit(1,1),xlimit(2,1)
            elseif(argument(1:2).eq.'-y')then
               read(argument(3:),*,err=201) xlimit(1,2),xlimit(2,2)
            elseif(argument(1:2).eq.'-z')then
               read(argument(3:),*,err=201) xlimit(1,3),xlimit(2,3)
            elseif(argument(1:2).eq.'-u')then
               read(argument(3:),*,err=201) vlimit(1,1),vlimit(2,1)
            elseif(argument(1:2).eq.'-v')then
               read(argument(3:),*,err=201) vlimit(1,2),vlimit(2,2)
            elseif(argument(1:2).eq.'-w')then
               read(argument(3:),*,err=201) vlimit(1,3),vlimit(2,3)
            elseif(argument(1:2).eq.'-p')then
               cellvol=-1
               bread=0.
               read(argument(3:),*,err=102)(Bdirs(j),j=1,3)
               bread=1.
 102           continue
               Bdirs(4)=0.
               do j=1,3
                  Bdirs(4)=Bdirs(4)+Bdirs(j)**2
               enddo
               if(Bdirs(4).gt.0)then
                  Bdirs(4)=sqrt(Bdirs(4))
                  do j=1,3
                     Bdirs(j)=Bdirs(j)/Bdirs(4)
                  enddo
               endif
               Bdirs(4)=Bdirs(4)*bread
            elseif(argument(1:2).eq.'-b')then
               read(argument(3:),*,err=201)iuin
            elseif(argument(1:2).eq.'-q')then
               ldoc=.true.
            endif
            if(argument(1:13).eq.'--objfilename')
     $        read(argument(14:),'(a)',err=201)objfilename
            if(argument(1:2).eq.'-f')
     $           read(argument(3:),'(a)',err=201)partfilename
            if(argument(1:2).eq.'-h')goto 203
            if(argument(1:2).eq.'-?')goto 203
         else
            read(argument(1:),'(a)',err=201)partfilename
!            write(*,*)partfilename
         endif
         
      enddo
      goto 202
!------------------------------------------------------------
! Help text
 201  continue
      write(*,*)'=====Error reading command line argument'
 203  continue
 301  format(a,i5)
 302  format(a,4f8.3)
      write(*,301)'Usage: partexamine [switches] <partfile> (no ext)'
      write(*,301)' --objfile<filename>  set name of object data file.'
     $     //' [ccpicgeom.dat'
      write(*,301)' -x -y -z<fff,fff>  set position range. [ -5,5'
      write(*,301)' -u -v -w<fff,fff>  set velocity range. [ -5,5'
      write(*,301)' -b<nx,ny,nz>  set spatial block range. [',iuin(1)
      write(*,302)' -p[bx,by,bz]  project [in direction]   [',Bdirs
      write(*,301)' -f   set name of partfile.'
      write(*,301)' -h -?   Print usage.'
      call exit(0)
 202  continue
      if(lentrim(partfilename).lt.5)then
         write(*,*)'Short filename, length<5 not allowed'
         goto 203
      endif
      end
