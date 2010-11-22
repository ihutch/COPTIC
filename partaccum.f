
c********************************************************************
      subroutine partdistup(xlimit,vlimit,xnewlim,cellvol,myid)
c Routine to update the particle distribution diagnostics for using the
c current particle information. If this is called for the first time,
c then calculate the bins for accumulating the distribution. Otherwise
c just accumulate to them.
      implicit none
      include 'meshcom.f'
      include 'partcom.f'
      include 'ptaccom.f'
      real xlimit(2,mdims),vlimit(2,mdims),xnewlim(2,mdims)
      real cellvol
      integer myid
      integer nvlist
c Use the first occasion to establish the accumulation range.
c Indicated by zero cellvolume.
         if(cellvol.eq.0.)then
c This ought to be determined based on the number of samples.
c But xlimits mean that's problematic.
            nvlist=100
            call vlimitdeterm(npdim,x_part,if_part,ioc_part,xlimit
     $           ,vlimit,nvlist)
c            write(*,'('' Velocity limits:'',6f7.3)') vlimit
            call minmaxreduce(mdims,vlimit)
c            write(*,'('' Velocity reduced:'',6f7.3)') vlimit
c Indicate csbin not initialized and start initialization
            csbin(1,1)=-1.
            call partacinit(xlimit,vlimit)

c Do the accumulation for this file up to maximum relevant slot. 
            nfvaccum=0
            call partsaccum(npdim,x_part,if_part,ioc_part,xlimit
     $        ,vlimit,xnewlim,nfvaccum)
            write(*,*)'Accumulated',nfvaccum,' of',ioc_part,' total'
     $        ,' in',xlimit
c Reduce back the data for MPI cases.
            call ptdiagreduce()
            call minmaxreduce(mdims,xnewlim)
            if(myid.eq.0)write(*,'(a,i8,a,6f8.4)')'Reduced',nfvaccum
     $           ,' Xnewlim=',xnewlim
c Should do this only the first time.
            call bincalc()
            call fvxinit(xnewlim,cellvol)
         else
            call partsaccum(npdim,x_part,if_part,ioc_part,xlimit
     $        ,vlimit,xnewlim,nfvaccum)
         endif
c         write(*,*)'isfull',isfull,cellvol
c         write(*,*)'calling subaccum'
         call subaccum(mdims,x_part,if_part,ioc_part,
     $     isfull,vlimit,xnewlim)
         end
c****************************************************************
      subroutine partacinit(xlimit,vlimit)
c Initialize bins for accumulation of the particles.
      include 'meshcom.f'
      include 'ptaccom.f'
      include 'plascom.f'
      include 'griddecl.f'
      real xlimit(2,mdims),vlimit(2,mdims)

c Initialization.
      do id=1,mdims
         do i=1,nptdiag
            fv(i,id)=0.
            px(i,id)=0.
            vdiag(i,id)=vlimit(1,id)
     $           +(vlimit(2,id)-vlimit(1,id))*(i-0.5)/nptdiag
            xdiag(i,id)=xmeshstart(id)+(i-0.5)*
     $           (xmeshend(id)-xmeshstart(id))/(nptdiag)
         enddo
c         write(*,*)'Position cell-center range',
c     $        id,xdiag(1,id),xdiag(nptdiag,id)
c         write(*,*)'Velocity cell-center range',
c     $        id,vdiag(1,id),vdiag(nptdiag,id)
      enddo
      end
c*****************************************************************
      subroutine oneaccum(xr,vlimit)
c Accumulate a particle into velocity bins.
      include 'meshcom.f'
      include 'ptaccom.f'
      include 'plascom.f'
      include 'griddecl.f'
      real vlimit(2,mdims)
      real xr(3*mdims)

      do id=1,mdims
c Assign velocities to bins.
         v=xr(mdims+id)
         if(v.lt.vlimit(1,id))v=vlimit(1,id)
         if(v.gt.vlimit(2,id))v=vlimit(2,id)
         ibin=nint(nptdiag*(.000005+0.99999*
     $        (v-vlimit(1,id))/(vlimit(2,id)-vlimit(1,id)))+0.5)
         if(ibin.lt.1.or.ibin.gt.nptdiag)
     $        write(*,*)k,nin,id,' ibin',ibin,v
         fv(ibin,id)=fv(ibin,id)+1.
         if(csbin(1,1).ne.-1.)then
c Doing summed bin accumulation
            ibs=ibinmap(ibin,id)
            fsv(ibs,id)=fsv(ibs,id)+1.
         endif
c Assign positions to bins
         x=(xr(id)-xmeshstart(id))/(xmeshend(id)-xmeshstart(id))
         ibin=nint(0.50000+x*(nptdiag-.00000))
         if(ibin.lt.1 .or. ibin.gt.nptdiag)then
            write(*,*)'ibin=',ibin,id,x,xr(id)
         else
            px(ibin,id)=px(ibin,id)+1.
         endif
      enddo

      end
c**********************************************************************
      subroutine partsaccum(mdims,xpart,ifpart,iocpart,xlimit
     $     ,vlimit,xnewlim,nfvaccum)
c Accumulate all particles in the xlimit range into velocity bins.
c If on entry xnewlim(2).eq.xnewlim(1) then adjust those limits.
      real xpart(3*mdims,iocpart)
      integer ifpart(iocpart)
c Spatial limits bottom-top, dimensions
      real xlimit(2,mdims),xnewlim(2,mdims)
c Velocity limits
      real vlimit(2,mdims)

      if(xnewlim(1,1).eq.xnewlim(2,1))then
         limadj=1
      else
         limadj=0
      endif
      do j=1,iocpart
c Only for filled slots
         if(ifpart(j).eq.1)then
            do id=1,mdims
               x=xpart(id,j)
               if(x.lt.xlimit(1,id).or.x.gt.xlimit(2,id))goto 12
               if(limadj.eq.1)then
c Adjust the xnewlim.
                  if(x.lt.xnewlim(1,id))xnewlim(1,id)=x
                  if(x.gt.xnewlim(2,id))xnewlim(2,id)=x
               endif
            enddo
            nfvaccum=nfvaccum+1
            call oneaccum(xpart(1,j),vlimit)
         endif
 12      continue
      enddo
c      write(*,*)'Accumulated',nfvaccum,' in',xlimit,' of',iocpart
c     $        ,' total'
c      if(limadj.eq.1)write(*,'(a,6f10.5)')' xnewlim=',xnewlim

             
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
      subroutine bincalc()
      include 'meshcom.f'
      include 'ptaccom.f'
      include 'griddecl.f'

c 
      do id=1,mdims
         cumfv(0,id)=0.
         do j=1,nsbins
            vsbin(j,id)=0.
            csbin(j,id)=0.
            fsv(j,id)=0.
         enddo
         dv=(vdiag(nptdiag,id)-vdiag(1,id))/(nptdiag-1)
         vhbin(0,id)=vdiag(1,id)-dv*0.5
         ib=1
         do k=1,nptdiag
            cumfv(k,id)=cumfv(k-1,id)+fv(k,id)/float(nfvaccum)
c            write(*,*)id,k,fv(k,id),cumfv(k,id),ib
c This linear mapping does not work well.
c            ib=1+ int(cumfv(k,id)*(nsbins)*(.99999))
            bx=float(ib)/nsbins
c cubic progression.
            cfn=1.0001*(3.*bx**2-2.*bx**3)
c            write(*,*)'cfn',cfn,cumfv(k,id)
            if(cumfv(k,id).gt.cfn)then
c find the histogram bin-boundary.
               vhbin(ib,id)=vhbin(0,id)+(k-1)*dv
               ib=ib+1
            endif
            ibinmap(k,id)=ib
            if(ib.lt.1 .or. ib.gt.nsbins)then
               write(*,*)'ibinmap error',ib,cumfv(k,id),cfn,k,id
               stop
            endif
            vsbin(ib,id)=vsbin(ib,id)+vdiag(k,id)
            csbin(ib,id)=csbin(ib,id)+1.
c Also accumulate this data into the summed bins
            fsv(ib,id)=fsv(ib,id)+fv(k,id)
         enddo
c Fix any summed bins that are off the end, and hence zero width.
         do ibk=ib,nsbins
            vhbin(ibk,id)=vhbin(0,id)+nptdiag*dv
         enddo
c Now cumfv is the cumulative probability distribution (to 1.0) over the
c uniform bins, and ibinmap maps those bins to nonuniform bins.  The
c uniform bin centers are in vdiag. Each non-uniform bin is composed of
c the sum of the uniform bins that map to it. Its center is therefore at
c the centroid of those bins = vsbin. 
c The number of uniform bins that they contain is csbin.
c But we invert csbin so that a multiplication is needed and zeros are 
c corrected for. 
         do k=1,nsbins
            if(csbin(k,id).eq.0)then
               write(*,*)'Fixing csbin zero',k,id,fv(k,id)
     $              ,vhbin(k,id)
               vsbin(k,id)=vhbin(k,id)
            else
               vsbin(k,id)=vsbin(k,id)/csbin(k,id)
               csbin(k,id)=1./csbin(k,id)
            endif
         enddo
      enddo
c      write(*,*)'Bincalc has chosen',nsbins,' bin placement.'
c      write(*,*)' vsbin',vsbin
c      write(*,*)(vhbin(k,1),k=0,nsbins)
c      write(*,*)' csbin',csbin
c      write(*,*)' fsv  ',fsv
c      write(*,*)' ibinmap',ibinmap
      end
c**********************************************************************
      subroutine subaccum(mdims,xpart,ifpart,iocpart,
     $     isfull,vlimit,xnewlim)
c Over the whole particle population accumulate to space-resolved bins
      real xpart(3*mdims,iocpart)
      integer ifpart(iocpart)
      integer isfull(mdims)
c Spatial limits bottom-top, dimensions
      real xnewlim(2,mdims)
c Velocity limits
      real vlimit(2,mdims)
c Subbin addressing.
      integer isind(mdims)
c Local storage with adequate dimensional length:
      real xsf(10)
c The equal bins fill the ranges xnewlim with index lengths isfull.
      do id=1,mdims
         xsf(id)=.999999/(xnewlim(2,id)-xnewlim(1,id))
      enddo
      do j=1,iocpart
c Only for filled slots
         if(ifpart(j).eq.1)then
            do id=1,mdims
c Calculate the spatial sub-bin pointer.
               isind(id)=1+int(isfull(id)*
     $              (xpart(id,j)-xnewlim(1,id))*xsf(id))
               if(isind(id).gt.isfull(id))isind(id)=isfull(id)
            enddo
            ip=ipfindex(mdims,isfull,isind)
c            if(ip.gt.1000)write(*,*)isfull,isind,ip
            call subvaccum(xpart(1,j),vlimit,ip+1)
         endif
      enddo
             
      end
c*****************************************************************
      subroutine subvaccum(xr,vlimit,ip)
c Accumulate a particle into velocity bin corresponding to position ip
      include 'meshcom.f'
      include 'ptaccom.f'
      include 'plascom.f'
      include 'griddecl.f'
      real vlimit(2,mdims)
      real xr(3*mdims)

      do id=1,mdims
c Assign velocities to bins.
         v=xr(mdims+id)
         if(v.lt.vlimit(1,id))v=vlimit(1,id)
         if(v.gt.vlimit(2,id))v=vlimit(2,id)
         ibin=nint(nptdiag*(.000005+0.99999*
     $        (v-vlimit(1,id))/(vlimit(2,id)-vlimit(1,id)))+0.5)
         if(ibin.lt.1.or.ibin.gt.nptdiag) write(*,*)id,' ibin',ibin
     $        ,nptdiag,v,vlimit(1,id),vlimit(2,id)
         if(csbin(1,1).ne.-1.)then
c Doing summed bin accumulation
            ibs=ibinmap(ibin,id)
c Test for corruption ought not eventually to be necessary.
            if(ip.le.0.or.ibs.lt.0.or.id.lt.0.or.ip.gt.nsub_tot)then
               write(*,*)'ibs,id,ip,nsub_tot',ibs,id,ip,nsub_tot
            endif
            fvx(ibs,id,ip)=fvx(ibs,id,ip)+1
         endif
      enddo
      denfvx(ip)=denfvx(ip)+1
      end
c*******************************************************************
      subroutine fvxinit(xnewlim,cellvol)
      include 'meshcom.f'
      include 'ptaccom.f'
      real xnewlim(2,mdims)

      isfull(1)=nsub_i
      isfull(2)=nsub_j
      isfull(3)=nsub_k
      do ip=1,nsub_tot
         denfvx(ip)=0.
         do id=1,mdims
            do ibs=1,nsbins
               fvx(ibs,id,ip)=0.
            enddo
         enddo
      enddo
      cellvol=1.
      do id=1,mdims
         cellvol=cellvol*(xnewlim(2,id)-xnewlim(1,id))/isfull(id)
c         write(*,*)xnewlim(2,id),xnewlim(1,id)
      enddo
c      write(*,*)'cellvol=',cellvol
      end
c**********************************************************************
      subroutine distread(xlimit,vlimit,xnewlim,name,cellvol)

      include 'meshcom.f'
      include 'ptaccom.f'
      real xlimit(2,mdims),vlimit(2,mdims),xnewlim(2,mdims)
      character*(*) name

      open(25,file=name,status='old',form='unformatted',err=101)
      read(25)nptdiagfile,mdimsfile
      read(25)(xlimit(1,j),xlimit(2,j),vlimit(1,j),vlimit(2,j),
     $     j=1,mdims)
      read(25)((xdiag(i,j),px(i,j),i=1,nptdiag),j=1,mdims)
      read(25)((vdiag(i,j),fv(i,j),i=1,nptdiag),j=1,mdims)
      read(25)(xnewlim(1,j),xnewlim(2,j),j=1,mdims)
      read(25)nsbf,(isfull(i),i=1,mdims),cellvol
      if(nsbf.ne.nsbins .or. isfull(1).ne.nsub_i .or.
     $     isfull(2).ne.nsub_j .or. isfull(3).ne.nsub_k) goto 103
      read(25)(((fvx(i,j,k),i=1,nsbins),j=1,mdims),k=1,nsub_tot)
      read(25)(denfvx(k),k=1,nsub_tot)
      read(25)((vsbin(i,j),csbin(i,j),fsv(i,j),i=1,nsbins),
     $     (vhbin(i,j),i=0,nsbins),j=1,mdims)
      read(25)((ibinmap(i,j),i=1,nptdiag),j=1,mdims)
      close(25)

      return
 101  write(*,*)'Error opening file:',name
      close(25,status='delete')
      stop
 103  write(*,*)'Mismatch of declared fv-cell dimensions',
     $     nsbins,nsub_i,nsub_j,nsub_k
      write(*,*)'               with those in file read:',nsbf,isfull
      write(*,*)'Adjust ptaccom.f to match file and recompile.'
      stop
 
      end
c**********************************************************************
      subroutine distwrite(xlimit,vlimit,xnewlim,name,cellvol)
      include 'meshcom.f'
      include 'ptaccom.f'
      real xlimit(2,mdims),vlimit(2,mdims),xnewlim(2,mdims)
      character*(*) name

      open(25,file=name,status='unknown',err=101)
      close(25,status='delete')
      open(25,file=name,status='new',form='unformatted',err=101)
      write(25)nptdiag,mdims
      write(25)(xlimit(1,j),xlimit(2,j),vlimit(1,j),vlimit(2,j),
     $     j=1,mdims)
      write(25)((xdiag(i,j),px(i,j),i=1,nptdiag),j=1,mdims)
      write(25)((vdiag(i,j),fv(i,j),i=1,nptdiag),j=1,mdims)
      write(25)(xnewlim(1,j),xnewlim(2,j),j=1,mdims)
      write(25)nsbins,(isfull(i),i=1,mdims),cellvol
      write(25)(((fvx(i,j,k),i=1,nsbins),j=1,mdims),k=1,nsub_tot)
      write(25)(denfvx(k),k=1,nsub_tot)
      write(25)((vsbin(i,j),csbin(i,j),fsv(i,j),i=1,nsbins),
     $     (vhbin(i,j),i=0,nsbins),j=1,mdims)
      write(25)((ibinmap(i,j),i=1,nptdiag),j=1,mdims)
      close(25)

      return
 101  write(*,*)'Error opening file:',name
      close(25,status='delete')

      end
c******************************************************************
      subroutine pltsubdist(icellin,jcellin,kcellin,vlimit,xnewlim
     $     ,cellvol)
      implicit none
      integer icellin,jcellin,kcellin
      include 'meshcom.f'
      include 'partcom.f'
      include 'ptaccom.f'
c Distributions in ptaccom.f
      real vlimit(2,mdims),xnewlim(2,mdims)
      real cellvol
      real fvplt(nsbins)
      integer ip,id,iv
      integer ip3index
      external ip3index
      integer lentrim
      external lentrim
      character*100 string
      integer icell,jcell,kcell
      real wicell,wjcell,wkcell
c------------------------------------------
      icell=icellin
      jcell=jcellin
      kcell=kcellin
 1    ip=ip3index(isfull,icell,jcell,kcell)+1
      call multiframe(3,1,2)
      wicell=(xnewlim(2,1)-xnewlim(1,1))/nsub_i
      wjcell=(xnewlim(2,2)-xnewlim(1,2))/nsub_j
      wkcell=(xnewlim(2,3)-xnewlim(1,3))/nsub_k
      write(string,'(''Cell'',3i3,''  x=('',f6.2,'','',f6.2,'')'//
     $     ' y=('',f6.2,'','',f6.2,'') z=('',f6.2,'','',f6.2,'')'')')
     $     icell,jcell,kcell
     $     ,xnewlim(1,1)+(icell-1)*wicell,xnewlim(1,1)+icell*wicell
     $     ,xnewlim(1,2)+(jcell-1)*wjcell,xnewlim(1,2)+jcell*wjcell
     $     ,xnewlim(1,3)+(kcell-1)*wkcell,xnewlim(1,3)+kcell*wkcell
      do id=1,3
c Do correct scaling:
         do iv=1,nsbins
            fvplt(iv)=fvx(iv,id,ip)*csbin(iv,id)
     $              *nptdiag/(vlimit(2,id)-vlimit(1,id))/cellvol
         enddo
c         call automark(vsbin(1,id),fvx(1,id,ip),nsbins,1)
         call automark(vsbin(1,id),fvplt,nsbins,1)
         if(id.eq.1)
     $        call boxtitle(string(1:lentrim(string)))
         write(string,'(a,i1,a)')'f(v!d',id,'!d)'
         if(id.eq.3)then
            call axlabels('Velocity',string(1:lentrim(string)))
         else
            call axlabels('',string(1:lentrim(string)))
         endif
         call winset(.true.)
c         call polybox(vhbin(0,id),fvx(1,id,ip),nsbins)
         call polybox(vhbin(0,id),fvplt,nsbins)
      enddo
      call eye3d(ip)
c Increment one of the dimensions using arrow keys.
      if(ip.eq.65361)then
         icell=mod(icell,nsub_i)+1
      elseif(ip.eq.65362)then
         jcell=mod(jcell,nsub_j)+1
      elseif(ip.eq.65363)then
         kcell=mod(kcell,nsub_k)+1
      else
         return
      endif
      goto 1
c--------------------------------------
      end
