      program partexamine
c Examine the particle data, showing distribution function(s) for 
c selected cell ranges.
      include 'examdecl.f'
c (Examdecl itself includes meshcom.f)
      parameter (nfilemax=999)
      include '../partcom.f'
      include '../ptaccom.f'
 
      character*10 chartemp
      character*100 name
      character*100 string

      real extra(nptdiag,mdims),diff(nptdiag,mdims)
c Spatial limits bottom-top, dimensions
      real xlimit(2,mdims),xnewlim(2,mdims)
c Velocity limits
      real vlimit(2,mdims)
      character*1 axnames(3)
      data axnames/'x','y','z'/ 

c      write(*,*)'na_i',na_i
      write(*,*)'nsub_i,nsub_j,nsub_k',nsub_i,nsub_j,nsub_k
      write(*,*)'nptdiag,nsbins',nptdiag,nsbins
      if(nptdiag.lt.nsbins)then
         write(*,*)'WARNING nptdiags=',nptdiag,' smaller than nsbins='
     $        ,nsbins,' Incorrect array size choice.'
      endif

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
c Needed initialization removed from partacinit.
         xmeshstart(id)=min(-5.,xlimit(1,id))
         xmeshend(id)=max(5.,xlimit(2,id))
      enddo
c Use cellvol=0. to indicate the first call. Nonzero subsequent.
      cellvol=0.

      call partexamargs(xlimit,vlimit)
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
         call partread(name,ierr)
         if(ierr.ne.0)goto 11

         call partdistup(xlimit,vlimit,xnewlim,cellvol,0)
c If we wish to accumulate to the uniform mesh for other than the first
c occasion (when cellvol is zero) we need to do:
c         call partsaccum

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
         write(string,'(a,i3)')'Distribution dimension',id
         call axlabels('velocity',string(1:lentrim(string)))
         call color(12)
         call polymark(vsbin(1,id),fsv(1,id),nsbins,1)
         call polybox(vhbin(0,id),fsv(1,id),nsbins)
         call color(13)
c         call polybox(vhbin(0,id),extra(1,id),nsbins)
c         call polybox(vhbin(0,id),diff(1,id),nsbins)
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

c Plot the subdistributions at a particular cell.
      icell=nsub_i/2+1
      jcell=nsub_j/2+1
      kcell=nsub_k/2+1
      call pltsubdist(icell,jcell,kcell,vlimit,xnewlim,cellvol)

      end
c*************************************************************
      subroutine partexamargs(xlimit,vlimit)
      include 'examdecl.f'
      real xlimit(2,3),vlimit(2,3)

      ifull(1)=na_i
      ifull(2)=na_j
      ifull(3)=na_k

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
            elseif(argument(1:2).eq.'-v')then
               read(argument(3:),*,err=201) vlimit(1,2),vlimit(2,2)
            elseif(argument(1:2).eq.'-w')then
               read(argument(3:),*,err=201) vlimit(1,3),vlimit(2,3)
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
c 302  format(a,f8.3)
      write(*,301)'Usage: partexamine [switches] <partfile> (no ext)'
      write(*,301)' --objfile<filename>  set name of object data file.'
     $     //' [ccpicgeom.dat'
      write(*,301)' -x -y -z<fff,fff>  set position range. [ -5,5'
      write(*,301)' -u -v -w<fff,fff>  set velocity range. [ -5,5'
      write(*,301)' -f   set name of partfile.'
      write(*,301)' -h -?   Print usage.'
      call exit(0)
 202  continue
      if(lentrim(partfilename).lt.5)then
         write(*,*)'Short filename, length<5 not allowed'
         goto 203
      endif
      end
