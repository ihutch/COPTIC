      program partexamine
c Examine the particle data, showing distribution function(s) for 
c selected cell ranges.
      include 'examdecl.f'
      parameter (nfilemax=999)
      include 'partcom.f'
 
      character*10 chartemp
      character*100 name
      character*50 string

      include 'ptaccom.f'
c Distributions in ptaccom.f

c Spatial limits bottom-top, dimensions
      real xlimit(2,mdims)
c Velocity limits
      real vlimit(2,mdims)

      nfmax=nfilemax
c silence warnings:
      zp(1,1,1)=0.
      xr(1)=0.

c Defaults
      do id=1,mdims
         xlimit(1,id)=-5.
         xlimit(2,id)=5.
c Overlapping vlimits make limitdeterm the usual setting method.
         vlimit(1,id)=5.
         vlimit(2,id)=-5.
      enddo

      call partexamargs(xlimit,vlimit)
c Now the base filename is in partfilename.
c Initialize partaccum.      Now later.
c      call partacinit(xlimit,vlimit)

      ip=lentrim(partfilename)-3
      if(partfilename(ip:ip).eq.'.')then
c If filename has a 3-character extension. Assume it is complete
c Should do this only the first time.
c and that we are reading just one file.
         nfmax=0
         name=partfilename
         write(*,*)'Reading single file ',name(1:lentrim(name))
         if(partfilename(ip:ip+3).eq.'.pex')then
            write(*,*)'Using stored distribution file'
            nfmax=-1
         endif
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
c Use the first file to establish the accumulation range.         
         if(i.eq.0)then
c This ought to be determined based on the number of samples.
c But xlimits mean that's problematic.
            nvlist=100
            call vlimitdeterm(npdim,x_part,if_part,ioc_part,xlimit
     $           ,vlimit,nvlist)
            call partacinit(xlimit,vlimit)
            write(*,'('' Velocity limits:'',6f7.3)') vlimit
         endif
c Do the accumulation for this file up to maximum relevant slot. 
         naccum=0
         call accumulate(npdim,x_part,if_part,ioc_part,naccum,xlimit
     $        ,vlimit)
         write(*,*)'Accumulated',naccum,' of',ioc_part,' total'
     $        ,' in',xlimit
c         write(*,*)'Calling bincalc'
c Should do this only the first time.
         if(i.eq.0)call bincalc(naccum)

      enddo
 11   continue

      if(nfmax.eq.-1)then
         open(25,file=name,status='old',form='unformatted',err=101)
         read(25)ndiagfile,mdimsfile
         read(25)(xlimit(1,j),xlimit(2,j),vlimit(1,j),vlimit(2,j),
     $        j=1,mdims)
         read(25)((xdiag(i,j),px(i,j),i=1,ndiag),j=1,mdims)
         read(25)((vdiag(i,j),fv(i,j),i=1,ndiag),j=1,mdims)
         close(25)
      else
         name(lentrim(name)-2:lentrim(name))='pex'
         open(25,file=name,status='unknown',err=101)
         close(25,status='delete')
         open(25,file=name,status='new',form='unformatted',err=101)
         write(25)ndiag,mdims
         write(25)(xlimit(1,j),xlimit(2,j),vlimit(1,j),vlimit(2,j),
     $        j=1,mdims)
         write(25)((xdiag(i,j),px(i,j),i=1,ndiag),j=1,mdims)
         write(25)((vdiag(i,j),fv(i,j),i=1,ndiag),j=1,mdims)
         close(25)
      endif
      goto 102
 101  write(*,*)'Error opening file:',name
      close(25,status='delete')
 102  continue


      call multiframe(2,1,2)
      do id=1,mdims
         do k=1,nsbins
            fk=fsv(k,id)
            fsv(k,id)=fk/csbin(k,id)
c            write(*,*)k,id,fsv(k,id),csbin(k,id)
         enddo
         call ticnumset(10)
         call autoplot(vdiag(1,id),fv(1,id),ndiag)
         write(string,'(a,i3)')'Distribution dimension',id
         call axlabels('velocity',string(1:lentrim(string)))
         call color(12)
         call polymark(vsbin(1,id),fsv(1,id),nsbins,1)
         call polybox(vhbin(0,id),fsv(1,id),nsbins)
         call color(15)
c         call pltend()
         call autoplot(xdiag(1,id),px(1,id),ndiag)
         call axlabels('position',string(1:lentrim(string)))
         call pltend()
      enddo

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
c****************************************************************
      subroutine partacinit(xlimit,vlimit)
c Accumulate the particles into bins.
      include 'ptaccom.f'
      include 'plascom.f'
      include 'meshcom.f'
      real xlimit(2,mdims),vlimit(2,mdims)

c Silence warning
      xr(1)=0.
c Indicate csbin not initialized yet:
      csbin(1,1)=0.
c Initialization.
      do id=1,mdims
         xmeshstart(id)=min(-5.,xlimit(1,id))
         xmeshend(id)=max(5.,xlimit(2,id))
         do i=1,ndiag
            vdiag(i,id)=vlimit(1,id)
     $           +(vlimit(2,id)-vlimit(1,id))*(i-0.5)/ndiag
            fv(i,id)=0.
            px(i,id)=0.
            xdiag(i,id)=xmeshstart(id)+(i-0.5)*
     $           (xmeshend(id)-xmeshstart(id))/(ndiag)
         enddo
c         write(*,*)'Position cell-center range',
c     $        id,xdiag(1,id),xdiag(ndiag,id)
c         write(*,*)'Velocity cell-center range',
c     $        id,vdiag(1,id),vdiag(ndiag,id)
      enddo
      return
      end
c*****************************************************************
      subroutine partaccum(xr,xlimit,vlimit)
c Accumulate a particle into bins.
      include 'ptaccom.f'
      include 'plascom.f'
      include 'meshcom.f'
      real xlimit(2,mdims),vlimit(2,mdims)

      do id=1,mdims
c Assign velocities to bins.
         v=xr(mdims+id)
         if(v.lt.vlimit(1,id))v=vlimit(1,id)
         if(v.gt.vlimit(2,id))v=vlimit(2,id)
         ibin=nint(ndiag*(.000005+0.99999*
     $        (v-vlimit(1,id))/(vlimit(2,id)-vlimit(1,id)))+0.5)
         if(ibin.lt.1.or.ibin.gt.ndiag)
     $        write(*,*)k,nin,id,' ibin',ibin,v
         fv(ibin,id)=fv(ibin,id)+1.
         if(csbin(1,1).ne.0.)then
c Doing summed bin accumulation
            ibs=ibinmap(ibin,id)
            fsv(ibs,id)=fsv(ibs,id)+1.
         endif
c Assign positions to bins
         x=(xr(id)-xmeshstart(id))/(xmeshend(id)-xmeshstart(id))
         ibin=nint(0.50000+x*(ndiag-.00000))
         if(ibin.lt.1 .or. ibin.gt.ndiag)then
            write(*,*)'ibin=',ibin,id,x,xr(id)
         else
            px(ibin,id)=px(ibin,id)+1.
         endif
      enddo

      end
c**********************************************************************
      subroutine accumulate(mdims,xpart,ifpart,iocpart,naccum,xlimit
     $     ,vlimit)
      real xpart(3*mdims,iocpart)
      integer ifpart(iocpart)
c Spatial limits bottom-top, dimensions
      real xlimit(2,mdims)
c Velocity limits
      real vlimit(2,mdims)

      do j=1,iocpart
c Only for filled slots
         if(ifpart(j).eq.1)then
            do id=1,mdims
               x=xpart(id,j)
               if(x.lt.xlimit(1,id).or.x.gt.xlimit(2,id))goto 12
            enddo
            naccum=naccum+1
            call partaccum(xpart(1,j),xlimit,vlimit)
         endif
 12      continue
      enddo
c         write(*,*)'Accumulated',naccum,' in',xlimit,' of',iocpart
c     $        ,' total'
             
      end
c**********************************************************************
      subroutine vlimitdeterm(mdims,xpart,ifpart,iocpart,xlimit,vlimit
     $     ,nvlist)
      real xpart(3*mdims,iocpart)
      integer ifpart(iocpart)
c Spatial limits bottom-top, dimensions
      real xlimit(2,mdims)
c Velocity limits
      real vlimit(2,mdims)
c Velocity Sorting arrays
      parameter (nvlistmax=200)
      real vtlist(nvlistmax),vblist(nvlistmax)
      if(nvlist.gt.nvlistmax)nvlist=nvlistmax
      do id=1,mdims
         do j=1,nvlist
            vblist(j)=vlimit(1,id)
            vtlist(j)=vlimit(2,id)
         enddo
         do j=1,iocpart
c Only for filled slots
            if(ifpart(j).eq.1)then
               v=xpart(id+mdims,j)
               call sorttoplimit(v,vtlist,nvlist)
               call sortbottomlimit(v,vblist,nvlist)
            endif
         enddo
c         write(*,*)vblist
c         write(*,*)vtlist
         vlimit(1,id)=vblist(nvlist)
         vlimit(2,id)=vtlist(nvlist)
      enddo
c      write(*,*)'vlimits',vlimit
      end
c========================================================================
      subroutine sortbottomlimit(v,vlist,nvlist)
      real vlist(nvlist)
c Insert the value v into its ordered place in vlist, retaining the bottom
c nvlist values.
c      write(*,*)'bot',v,nvlist,vlist
      if(v.ge.vlist(nvlist))return
      if(v.ge.vlist(1))then
         i1=1
         i2=nvlist
c Find my position by bisection.
 1       i=(i1+i2)/2
         if(v.lt.vlist(i))then
            i2=i
         else
            i1=i
         endif
         if(i2.gt.i1+1)goto 1
      else
         i2=1
      endif
c Here i2 is the position of list value just greater than v.
c      write(*,*)'botend',i1,i2,i,vlist(i2),v
      do i=nvlist,i2+1,-1
         vlist(i)=vlist(i-1)
      enddo
      vlist(i2)=v

      end
c*********************************************************************
      subroutine sorttoplimit(v,vlist,nvlist)
      real vlist(nvlist)
c Insert the value v into its reverse-ordered place in vlist, retaining
c the top nvlist values.
c      write(*,*)'topstart',nvlist,v,vlist
      if(v.le.vlist(nvlist))return
      if(v.le.vlist(1))then
         i1=1
         i2=nvlist
c Find my position by bisection.
 1       i=(i1+i2)/2
         if(v.gt.vlist(i))then
            i2=i
         else
            i1=i
         endif
         if(i2.gt.i1+1)goto 1
      else
         i2=1
      endif
c Here i2 is the position of list value just less than v.
c      write(*,*)'topend',i1,i2,i,vlist(i2),v
      do i=nvlist,i2+1,-1
         vlist(i)=vlist(i-1)
      enddo
      vlist(i2)=v

      end

c***********************************************************************
      subroutine bincalc(naccum)
      integer naccum
      include 'ptaccom.f'

c silence warning
      b=xr(1)
c 
      do id=1,mdims
         cumfv(0,id)=0.
         do j=1,nsbins
            vsbin(j,id)=0.
            csbin(j,id)=0.
            fsv(j,id)=0.
         enddo
         dv=(vdiag(ndiag,id)-vdiag(1,id))/(ndiag-1)
         vhbin(0,id)=vdiag(1,id)-dv*0.5
         ib=1
         do k=1,ndiag
            cumfv(k,id)=cumfv(k-1,id)+fv(k,id)/float(naccum)
c This linear mapping does not work well.
c            ib=1+ int(cumfv(k,id)*(nsbins)*(.99999))
            bx=float(ib)/nsbins
c cubic progression.
            cfn=1.0001*(3.*bx**2-2.*bx**3)
            if(cumfv(k,id).gt.cfn)then
c find the histogram bin-boundary.
               vhbin(ib,id)=vhbin(0,id)+(k-1)*dv
               ib=ib+1
            endif
            ibinmap(k,id)=ib
            if(ib.lt.1 .or. ib.gt.nsbins)then
               write(*,*)'ibinmap error',ib,cumfv(k,id),cfn
               stop
            endif
            vsbin(ib,id)=vsbin(ib,id)+vdiag(k,id)
            csbin(ib,id)=csbin(ib,id)+1.
c Also accumulate this data into the summed bins
            fsv(ib,id)=fsv(ib,id)+fv(k,id)
         enddo
         vhbin(nsbins,id)=vhbin(0,id)+ndiag*dv
c Now cumfv is the cumulative probability distribution (to 1.0) over the
c uniform bins, and ibinmap maps those bins to nonuniform bins.  The
c uniform bin centers are in vdiag. Each non-uniform bin is composed of
c the sum of the uniform bins that map to it. Its center is therefore at
c the centroid of those bins = vsbin. The number of uniform bins that
c they contain is csbin. 
         do k=1,nsbins
            vsbin(k,id)=vsbin(k,id)/csbin(k,id)
         enddo
      enddo
      write(*,*)'Bincalc has chosen',nsbins,' bin placement.'
c      write(*,*)' vsbin',vsbin
c      write(*,*)' csbin',csbin
c      write(*,*)' fsv  ',fsv
c      write(*,*)' ibinmap',ibinmap
      end
