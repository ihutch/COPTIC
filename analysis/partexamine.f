      program partexamine
c Examine the particle data, showing distribution function(s) for 
c selected cell ranges.
      include 'examdecl.f' 
c (Examdecl itself includes meshcom.f plascom.f, objcom.f)
      parameter (nfilemax=999)
      include '../partcom.f'
      include '../ptaccom.f'
 
      character*10 chartemp
      character*100 name
      character*100 string
      logical ldoc
      real extra(nptdiag,mdims),diff(nptdiag,mdims)
c Spatial limits bottom-top, dimensions
      real xlimit(2,mdims),xnewlim(2,mdims)
      real Bdirs(mdims+1)
c Velocity limits
      real vlimit(2,mdims),wicell,wjcell,wkcell
      character*1 axnames(3)
      real vmean(mdims)
      real vtkxn(nsub_i+1),vtkyn(nsub_j+1),vtkzn(nsub_k+1)
      integer iuin(mdims),vtkifull(mdims),vtkiuds(mdims),ibinary
      integer ivtk,ip3index,ip,icentering(2*mdims)
      integer ivardims(2*mdims),ivardims_alloc(2*mdims)
      external ip3index
      data axnames/'x','y','z'/ 
      character*200 ivarnames(2*mdims)

      nfmax=nfilemax
c silence warnings:
      zp(1,1,1)=0.

c Defaults
      do id=1,mdims
c Use very big xlimits by default to include whole domain
c They are then reset by the accumulation itself.
         xlimit(1,id)=-500.
         xlimit(2,id)=500.
         xnewlim(1,id)=0.
         xnewlim(2,id)=0.
c Overlapping vlimits make limitdeterm the usual setting method.
         vlimit(1,id)=5.
         vlimit(2,id)=-5.
      enddo

      call partexamargs(xlimit,vlimit,iuin,cellvol,Bdirs,ldoc,ivtk)
      do id=1,mdims
c Needed initialization removed from partacinit.
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
c Now the base filename is in partfilename.

      ip=lentrim(partfilename)-3
      if(partfilename(ip:ip).eq.'.')then
c If filename has a 3-character extension or is a numbered pex
c file. Assume it is complete and that we are reading just one
c file. Should do this only the first time.
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
c Possible multiple files.
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
c All directions were set by commandline. Or none were read from file.
            do k=1,ndims
               Bfield(k)=Bdirs(k)
            enddo
         endif
         if(cellvol.eq.-1)write(*,*)'Bfield (projection)',Bfield
c The cellvol==-1 call will set ivproj=1 in ptaccom.
         call partdistup(xlimit,vlimit,xnewlim,cellvol,0,isuds)

      enddo
 11   continue

      if(nfmax.eq.-1)then
         call distread(xlimit,vlimit,xnewlim,name,cellvol)
c Calculate the binned data.
         nfvaccum=0
         do i=1,nptdiag
            nfvaccum=nfvaccum+px(i,1)
         enddo
c         write(*,*)'nfvaccum,cellvol',nfvaccum,cellvol
c Don't call bincalc if we read data back.         
c         call bincalc()
      else
         name(lentrim(name)-2:lentrim(name))='pex'
         call distwrite(xlimit,vlimit,xnewlim,name,cellvol)
      endif
      
c      write(*,*)'isfull'
c      write(*,*)'isuds'
c      write(*,*)'vsbin',vsbin
c      write(*,*)'fvx',fvx
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
         do i=1,mdims
            vtkifull(i)=isfull(i)+1
         enddo
         do j=1,mdims
            vtkiuds(j)=isuds(j)+1
         enddo
         ibinary=0

         do i=1,isuds(1)
            do j=1,isuds(2)
               do k=1,isuds(3)
                  do l=0,nsbins
                     do m=1,mdims
                        vtkudata(i,j,k,l,m)=vhbin(l,m)
                     enddo
                  enddo
                  do l=1,nsbins
                     do m=mdims+1,2*mdims
                        ip=ip3index(isuds,i,j,k)+1
                        vtkudata(i,j,k,l-1,m)=fvx(l,m-mdims,ip)
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
         do i=1,2*mdims
            icentering(i)=0
         enddo         
         do i=1,mdims
            ivardims(i)=nsbins+1
            ivardims_alloc(i)=nsbins+1
         enddo
         do i=mdims+1,2*mdims
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
     $        ,2*mdims,ivardims_alloc
     $        ,icentering) 
      else
             
      call multiframe(2,1,2)
      do id=1,mdims
         do k=1,nsbins
            fk=fsv(k,id)
               fsv(k,id)=fk*csbin(k,id)
     $           *nptdiag/(vlimit(2,id)-vlimit(1,id))
c            write(*,*)k,id,fsv(k,id),csbin(k,id)
            extra(nsbins+1-k,id)=fsv(k,id)
         enddo
         do k=1,nptdiag
               fv(k,id)=fv(k,id)
     $           *nptdiag/(vlimit(2,id)-vlimit(1,id))
         enddo
         do k=1,nsbins
            diff(k,id)=fsv(k,id)-fsv(nsbins+1-k,id)
         enddo
c         write(*,*)nsbins
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
c         call polybox(vhbin(0,id),extra(1,id),nsbins)
c         call polybox(vhbin(0,id),diff(1,id),nsbins)
c         write(*,*) 'vdiag'
c         write(*,*)(vdiag(kk,id),kk=1,nptdiag)
c         write(*,*) 'fv'
c         write(*,*) (fv(kk,id),kk=1,nptdiag)
c         write(*,*)(vsbin(kk,id),kk=1,nsbins)
c         write(*,*)(fsv(kk,id),kk=1,nsbins)
         call color(15)
c         call pltend()
         call autoplot(xdiag(1,id),px(1,id),nptdiag)
         call axlabels('position',string(1:lentrim(string)))
         call pltend()
      enddo

      if(nfmax.eq.-1)then
c Print to file the distribution.
         il=lentrim(name)
         name(il-2:il-1)='fv'
         open(25,file=name,status='unknown',err=101)
         close(25,status='delete')
         open(25,file=name,status='new',err=101)
         do id=1,mdims
c            name(il:il)=char(48+id)
            write(25,*)'Dimension ',id
            write(25,'(a,a,a)')'legend: f(v!d',axnames(id),'!d)'
c The following is for mdims>3.
c            write(25,'(a,i1,a)')'legend: f(v!d',id,'!d)'
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
      do id=1,mdims
         vmean(id)=0.
         ftot=0.
         do i=1,nptdiag
            vmean(id)=vmean(id)+fv(i,id)*vdiag(i,id)
            ftot=ftot+fv(i,id)
         enddo
         vmean(id)=vmean(id)/ftot
      enddo
      write(*,'(''Mean velocity:'',3f10.4)')(vmean(i),i=1,mdims)

c Plot the subdistributions at a particular cell.
      icell=isuds(1)/2+1
      jcell=isuds(2)/2+1
      kcell=isuds(3)/2+1
      call pltsubdist(icell,jcell,kcell,vlimit,xnewlim,cellvol)
      endif                   
      end
c*************************************************************
      subroutine partexamargs(xlimit,vlimit
     $           ,iuin,cellvol,Bdirs,ldoc,ivtk)
      include 'examdecl.f'
      include '../ptaccom.f'
      real xlimit(2,3),vlimit(2,3),Bdirs(4)
      integer iuin(3)
      logical ldoc
      integer ivtk
c I think unused here 26 May 12. But I'm not sure.
      ifull(1)=na_i
      ifull(2)=na_j
      ifull(3)=na_k
      ivtk=0

      do i=1,3
         iuin(i)=9
c convenient default field:
         Bdirs(i)=2-i
      enddo
c Use cellvol=0. by default.
      cellvol=0.
      ldoc=.false.

c silence warnings:
      fluxfilename=' '
      zp(1,1,1)=0.
c Defaults
      partfilename=' '

c Deal with arguments
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
c            write(*,*)partfilename
         endif
         
      enddo
      goto 202
c------------------------------------------------------------
c Help text
 201  continue
      write(*,*)'=====Error reading command line argument'
 203  continue
 301  format(a,i5)
 302  format(a,4f8.3)
 303  format(a,5i5)
      write(*,301)'Usage: partexamine [switches] <partfile> (no ext)'
      write(*,301)' --objfile<filename>  set name of object data file.'
     $     //' [copticgeom.dat'
      write(*,301)' -x -y -z<fff,fff>  set position range. [ -5,5'
      write(*,301)' -u -v -w<fff,fff>  set velocity range. [ -5,5'
      write(*,303)' -b<nx,ny,nz>  set spatial block range. [',iuin
      write(*,302)' -p[bx,by,bz]  project [in direction]   [',Bdirs
      write(*,301)' -f   set name of partfile.'
      write(*,301)' -h -?   Print usage.'
      write(*,301)' -vtk   Output distribution function vtk files.'
      call exit(0)
 202  continue
      if(lentrim(partfilename).lt.5)then
         write(*,*)'Short filename, length<5 not allowed'
         goto 203
      endif
      end
