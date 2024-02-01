      program partexamine
! Examine the particle data, showing distribution function(s) for 
! selected cell ranges.
      include 'examdecl.f' 
! (Examdecl itself includes meshcom.f plascom.f, objcom.f)
      parameter (nfilemax=999)
      include '../src/partcom.f'
      include '../src/ptaccom.f'
      include '../src/myidcom.f'
 
      character*10 chartemp
      character*20 theformat
      character*100 name
      character*100 string
      logical ldoc
      real extra(nptdiagmax,ndimsmax),diff(nptdiagmax,ndimsmax)
! Spatial limits bottom-top, dimensions
      real xlimit(2,ndimsmax),xnewlim(2,ndimsmax)
      real Bdirs(ndimsmax+1)
! Velocity limits
      real vlimit(2,ndimsmax),wicell,wjcell,wkcell
      character*1 axnames(3)
      real vmean(ndimsmax)
      real vtkxn(nsub_i+1),vtkyn(nsub_j+1),vtkzn(nsub_k+1)
      integer iuin(ndimsmax),vtkifull(ndimsmax),vtkiuds(ndimsmax)
      integer ibinary
      integer ivtk,ip3index,ip,icentering(2*ndimsmax)
      integer ivardims(2*ndimsmax),ivardims_alloc(2*ndimsmax)
      integer ndfirst,ndlast
      external ip3index
      data axnames/'x','y','z'/ 
      character*200 ivarnames(2*ndimsmax)

      nfmax=nfilemax
! silence warnings:
      zp(1,1,1)=0.

! Defaults
      nprocs=1
      ndfirst=1
      ndlast=3
      do id=1,ndimsmax
! Use very big xlimits by default to include whole domain
! They are then reset by the accumulation itself.
         xlimit(1,id)=-500.
         xlimit(2,id)=500.
         xnewlim(1,id)=0.
         xnewlim(2,id)=0.
! Overlapping vlimits make limitdeterm the usual setting method.
         vlimit(1,id)=5.
         vlimit(2,id)=-5.
      enddo

      call partexamargs(xlimit,vlimit,iuin,cellvol,Bdirs,ldoc,ivtk
     $     ,ispecies,ndfirst,ndlast)
      do id=1,ndimsmax
! Needed initialization removed from partacinit.
         xmeshstart(id)=min(-5.,xlimit(1,id))
         xmeshend(id)=max(5.,xlimit(2,id))
         isuds(id)=iuin(id)
      enddo
!      write(*,*)' isuds:',isuds
      if(isuds(1)*isuds(2)*isuds(3).gt.nsub_tot)then
         write(*,*)'Too many blocks requested. Reduce or recompile.'
         stop
      endif
!      write(*,*)'nptdiag,nsbins',nptdiag,nsbins
      if(nptdiag.lt.nsbins)then
         write(*,*)'WARNING nptdiags=',nptdiag,' smaller than nsbins='
     $        ,nsbins,' Incorrect array size choice.'
      endif
! Now the base filename is in partfilename.

!      ip=lentrim(partfilename)-3
      ip=istrstr(partfilename,'.')
! If filename has an character extension or is a (numbered) pex
! file. Assume it is complete and that we are reading just one
! file. Should do this only the first time.
      if(ip.gt.0)then
      if(partfilename(ip:ip+3).eq.'.pex')then
         nfmax=-1
         name=partfilename
         write(*,*)'Reading pex file ',name(1:lentrim(name))
      elseif(partfilename(ip:ip).eq.'.')then
         nfmax=0
         name=partfilename
         write(*,*)'Reading single file ',name(1:lentrim(name))
         if(partfilename(ip:ip+3).eq.'.pex')then
            write(*,*)'Using stored distribution file'
            nfmax=-1
         endif
      endif
      endif
! Possible multiple files.
      il=3
 12   do i=0,nfmax
         if(nfmax.ne.0)then
            write(theformat,'(a,i1,a,i1,a)')'(''.'',i',il,'.',il,')'
            write(chartemp,theformat)i
            name=partfilename(1:lentrim(partfilename))//chartemp
            write(*,'(3a,$)')' Reading ',name(1:lentrim(name)),' '
         endif
         Bt=0.
         call partread(name,ierr)
         if(ierr.ne.0.and.i.eq.0)
     $        write(*,*)'partread returns error:',ierr
         if(ierr-4*(ierr/4).ne.0)then
            goto 11
         endif
         if(ispecies.gt.nspecies)then
            ispecies=nspecies
            write(*,*)'nspecies=',nspecies,'  Reset ispecies',ispecies
         endif
         if(ldoc)then
            write(*,'(a,5f9.4)')' debyelen,Ti,vd,rs,phip=',debyelen,Ti
     $           ,vd,rs,phip
!            write(*,*)'iregion_part,n_part,dt,ldiags,rhoinf,nrein,',
!     $           'phirein,numprocs='
!            write(*,*)iregion_part,n_part,dt,ldiags,rhoinf,nrein,phirein
!     $           ,numprocs
!            write(*,*)'eoverm,Bt,Bfield,vpar,vperp=',eoverm,Bt,Bfield
!     $           ,vpar,vperp
!            write(*,*)'nspecies',nspecies
!            stop
         endif
         if(Bdirs(4).gt.0. .or. Bt.eq.0)then
! All directions were set by commandline. Or none were read from file.
            do k=1,ndims
               Bfield(k)=Bdirs(k)
            enddo
         endif
         if(cellvol.eq.-1)write(*,*)'Bfield (projection)',Bfield
! The cellvol==-1 call will set ivproj=1 in ptaccom.
         call partdistup(xlimit,vlimit,xnewlim,cellvol,0,isuds,ispecies)
!         write(*,*)'partdistup completed',i,fsv(10,1),ifsv(10,1)
      enddo
 11   continue
      il=il+1
      if(il.le.4)goto 12

      if(nfmax.eq.-1)then
         call distread(xlimit,vlimit,xnewlim,name,cellvol)
! Calculate the binned data.
         nfvaccum=0
         do i=1,nptdiag
            nfvaccum=nfvaccum+int(px(i,1))
         enddo
!         write(*,*)'nfvaccum,cellvol',nfvaccum,cellvol
! Don't call bincalc if we read data back.         
!         call bincalc()
      else
         ip=istrstr(name,'.')
         name(ip+1:)='pex'
         call distwrite(xlimit,vlimit,xnewlim,name,cellvol)
      endif
      
!      write(*,*)'isfull'
!      write(*,*)'isuds'
!      write(*,*)'vsbin',vsbin
!      write(*,*)'fvx',fvx
      if (ivtk.eq.1)then
         wicell=(xnewlim(2,1)-xnewlim(1,1))/isuds(1)
         wjcell=(xnewlim(2,2)-xnewlim(1,2))/isuds(2)
         wkcell=(xnewlim(2,3)-xnewlim(1,3))/isuds(3)
         do i=1,isuds(1)+1
            vtkxn(i)=xnewlim(1,1)+(i-1)*wicell
         enddo
         do j=1,isuds(2)+1
            vtkyn(j)=xnewlim(1,2)+(j-1)*wjcell
         enddo
         do k=1,isuds(3)+1
            vtkzn(k)=xnewlim(1,3)+(k-1)*wkcell
         enddo
         do i=1,ndimsmax
            vtkifull(i)=isfull(i)+1
         enddo
         do j=1,ndimsmax
            vtkiuds(j)=isuds(j)+1
         enddo
         ibinary=0

         do i=1,isuds(1)
            do j=1,isuds(2)
               do k=1,isuds(3)
                  do l=0,nsbins
                     do m=1,ndimsmax
                        vtkudata(i,j,k,l,m)=vhbin(l,m)
                     enddo
                  enddo
                  do l=1,nsbins
                     do m=ndimsmax+1,2*ndimsmax
                        ip=ip3index(isuds,i,j,k)+1
                        vtkudata(i,j,k,l-1,m)=fvx(l,m-ndimsmax,ip)
                     enddo
                  enddo
               enddo
            enddo
         enddo
         ivarnames(1)='vx'//char(0)
         ivarnames(2)='vy'//char(0)
         ivarnames(3)='vz'//char(0)
         ivarnames(4)='fvx'//char(0)
         ivarnames(5)='fvy'//char(0)
         ivarnames(6)='fvz'//char(0)
         do i=1,2*ndimsmax
            icentering(i)=0
         enddo         
         do i=1,ndimsmax
            ivardims(i)=nsbins+1
            ivardims_alloc(i)=nsbins+1
         enddo
         do i=ndimsmax+1,2*ndimsmax
            ivardims(i)=nsbins
            ivardims_alloc(i)=nsbins+1
         enddo
         call vtkwritevars(vtkifull,vtkiuds
     $        ,vtkudata
     $        ,vtkxn,vtkyn
     $        ,vtkzn
     $        ,ibinary,ivardims
     $        ,partfilename(1:lentrim(partfilename))//char(0)
     $        ,ivarnames,200
     $        ,2*ndimsmax,ivardims_alloc
     $        ,icentering) 
      else
             
      call multiframe(2,1,2)
      do id=1,ndimsmax
         do k=1,nsbins
            fk=fsv(k,id)
               fsv(k,id)=fk*csbin(k,id)
     $           *nptdiag/(vlimit(2,id)-vlimit(1,id))
!            write(*,*)k,id,fsv(k,id),csbin(k,id)
            extra(nsbins+1-k,id)=fsv(k,id)
         enddo
         do k=1,nptdiag
               fv(k,id)=fv(k,id)
     $           *nptdiag/(vlimit(2,id)-vlimit(1,id))
         enddo
         do k=1,nsbins
            diff(k,id)=fsv(k,id)-fsv(nsbins+1-k,id)
         enddo
!         write(*,*)nsbins
         call ticnumset(10)
         call autoplot(vdiag(1,id),fv(1,id),nptdiag)
         call boxtitle('Total particles per unit velocity')
         if(ivproj.eq.0)then
            write(string,'(a,i3)')'Distribution dimension',id
         else
            write(string,'(a,i3)')'Distribution projection',id
         endif
         call axlabels('velocity',string(1:lentrim(string)))
         call color(12)
         call polymark(vsbin(1,id),fsv(1,id),nsbins,1)
         call polybox(vhbin(0,id),fsv(1,id),nsbins)
         call color(13)
!         call polybox(vhbin(0,id),extra(1,id),nsbins)
!         call polybox(vhbin(0,id),diff(1,id),nsbins)
!         write(*,*) 'vdiag'
!         write(*,*)(vdiag(kk,id),kk=1,nptdiag)
!         write(*,*) 'fv'
!         write(*,*) (fv(kk,id),kk=1,nptdiag)
!         write(*,*)(vsbin(kk,id),kk=1,nsbins)
!         write(*,*)(fsv(kk,id),kk=1,nsbins)
         call color(15)
!         call pltend()
         call autoplot(xdiag(1,id),px(1,id),nptdiag)
         if(nptdiag.eq.nsbins)then
            write(*,*)'Density as a function of position dimension',id
            do kk=1,nptdiag
               write(*,*)xdiag(kk,id),px(kk,id)
            enddo
         endif
         call axlabels('position',string(1:lentrim(string)))
         call pltend()
      enddo

      if(nfmax.eq.-1)then
! Print to file the distribution.
         il=lentrim(name)
         name(il-2:il-1)='fv'
         open(25,file=name,status='unknown',err=101)
         close(25,status='delete')
         open(25,file=name,status='new',err=101)
         do id=1,ndimsmax
!            name(il:il)=char(48+id)
            write(25,*)'Dimension ',id
            write(25,'(a,a,a)')'legend: f(v!d',axnames(id),'!d)'
! The following is for ndimsmax>3.
!            write(25,'(a,i1,a)')'legend: f(v!d',id,'!d)'
            write(25,*)nptdiag
            write(25,'(2g14.6)')(vdiag(i,id),fv(i,id),i=1,nptdiag)
         enddo
         close(25)
         goto 102
 101     write(*,*)'Error opening file:',name
         close(25,status='delete')
         call exit(1)
 102     continue
      endif
      do id=1,ndimsmax
         vmean(id)=0.
         ftot=0.
         do i=1,nptdiag
            vmean(id)=vmean(id)+fv(i,id)*vdiag(i,id)
            ftot=ftot+fv(i,id)
         enddo
         vmean(id)=vmean(id)/ftot
      enddo
      write(*,'(''Mean velocity:'',3f10.4)')(vmean(i),i=1,ndimsmax)

! Plot the subdistributions at a particular cell.
      icell=isuds(1)/2+1
      jcell=isuds(2)/2+1
      kcell=isuds(3)/2+1
      call pltsubdist(icell,jcell,kcell,vlimit,xnewlim,cellvol
     $     ,ndfirst,ndlast)
      endif                   
      end
!*************************************************************
      subroutine partexamargs(xlimit,vlimit
     $           ,iuin,cellvol,Bdirs,ldoc,ivtk,ispecies,ndfirst,ndlast)
      include 'examdecl.f'
      include '../src/ptaccom.f'
      real xlimit(2,3),vlimit(2,3),Bdirs(4)
      integer iuin(3)
      logical ldoc
      integer ivtk
! I think unused here 26 May 12. But I'm not sure.
      ifull(1)=na_i
      ifull(2)=na_j
      ifull(3)=na_k
      ivtk=0
      ispecies=1

      do i=1,3
         iuin(i)=9
! Default field/projection direction:
         Bdirs(i)=0
         if(i.eq.3)Bdirs(i)=1
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
            elseif(argument(1:4).eq.'-vtk')then
               ivtk=1
            elseif(argument(1:2).eq.'-v')then
               read(argument(3:),*,err=201) vlimit(1,2),vlimit(2,2)
            elseif(argument(1:2).eq.'-w')then
               read(argument(3:),*,err=201) vlimit(1,3),vlimit(2,3)
            elseif(argument(1:3).eq.'-pf')then
               call pfset(3)
            elseif(argument(1:3).eq.'-pu')then
               nptdiag=nsbins
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
            elseif(argument(1:3).eq.'-sp')then
               ispecies=ispecies+1
            elseif(argument(1:2).eq.'-d')then
               read(argument(3:),*,err=201)ndfirst,ndlast
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
 301  format(a,i5,a)
 302  format(a,4f8.3)
 303  format(a,5i5)
      write(*,301)'Usage: partexamine [switches] <partfile> (no ext)'
      write(*,302)' -x -y -z<fff,fff>  set position range. [',
     $     xlimit(1,1),xlimit(2,1)
      write(*,302)' -u -v -w<fff,fff>  set velocity range. ['
     $     ,vlimit(1,1),vlimit(2,1)
      write(*,303)' -b<nx,ny,nz>  set spatial block range. [',iuin
      write(*,302)' -p[bx,by,bz]  project [in direction]   [',Bdirs
      write(*,'(a,$)')' -d<i,i> Set first,last v-index'
      write(*,303)' to plot [',ndfirst,ndlast
      write(*,301)' -sp     Increment species number to examine'
      write(*,301)' -pu     Set nptdiag for uniform bins   [',nptdiag
     $     ,'  (Hence 2-D f plots)'
      write(*,301)' -pf     Output ps files of plots'
      write(*,301)' -q      Output diagnostics of file reading'
      write(*,301)' --objfile<filename>  set name of object data file.'
     $     //' [copticgeom.dat'
      write(*,301)' -vtk    Output distribution function vtk files.'
      write(*,301)' -h -?   Print usage.'
      call exit(0)
 202  continue
      if(lentrim(partfilename).lt.5)then
         write(*,*)'Short filename, length<5 not allowed'
         goto 203
      endif
      end
