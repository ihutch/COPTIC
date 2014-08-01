c********************************************************************
      subroutine partdistup(xlimit,vlimit,xnewlim,cellvol,myid,ibset
     $     ,ispecies)
c Routine to update the particle distribution diagnostics for using the
c current particle information for species ispecies.
c
c If this is called for the first time, i.e. with cellvol.le.0, then
c calculate the bins for accumulating the distribution. Also, in this
c case, if cellvol=-1 project in the direction of Bfield, else cellvol=0
c just use cartesian components.  For this first call, in a
c multiprocessor situation, we all-reduce the accumulated fine
c distribution fv, nvaccum (and fsv) on which the binning is calculated.
c Thus the bin calculation is based upon the total particles, not just
c the particles for the head node.  However, for use with multiple
c particle files by partexamine, only the first file read is used to
c determine the bin mapping. [Otherwise we would have to read all files
c a second time for the subbinning.]
c
c Otherwise if cellvol.gt.0 just accumulate to bins.
      implicit none
      include 'ndimsdecl.f'
      include 'plascom.f'
      include 'meshcom.f'
      include 'partcom.f'
      include 'ptaccom.f'
      real xlimit(2,ndimsmax),vlimit(2,ndimsmax),xnewlim(2,ndimsmax)
      real cellvol
      integer myid,ibset
      integer nvlist
      integer ispecies

c Use the first occasion to establish the accumulation range.
c Indicated by cellvolume le 0.
      if(cellvol.le.0.)then
c Decide whether to project velocity in the direction of Bfield.
         if(cellvol.eq.-1)then
            ivproj=1
         else
            ivproj=0
         endif
c This ought to be determined based on the number of samples.
c But xlimits mean that's problematic.
c            nvlist=100
c A small nvlist means we go right to the edge of actual particles.
         nvlist=5
         call vlimitdeterm(vlimit,nvlist,ivproj,Bfield,ispecies)
         if(myid.eq.0)write(*,'('' Velocity limits:'',6f8.3)')
     $        vlimit
c         write(*,*)'ispecies=',ispecies
         call minmaxreduce(ndimsmax,vlimit)
c            write(*,'('' Velocity reduced:'',6f7.3)') vlimit
c Indicate csbin not initialized and start initialization
         csbin(1,1)=-1.
         call partacinit(vlimit)

c Do the accumulation for this file up to maximum relevant slot. 
         nfvaccum=0
         call partsaccum(xlimit ,vlimit,xnewlim,nfvaccum,ispecies)
         if(myid.eq.0)write(*,*)'Accumulated',nfvaccum,' of',ioc_part
     $        ,' total',' in',xlimit
c Reduce back the data for MPI cases.
         call ptdiagreduce()
         call minmaxreduce(ndimsmax,xnewlim)
         if(myid.eq.0)write(*,'(a,i8,a,6f8.3)')'Reduced',nfvaccum
     $        ,' Xnewlim=',xnewlim
c Should do this only the first time.
         call bincalc()
         call fvxinit(xnewlim,cellvol,ibset)
      else
         call partsaccum(xlimit ,vlimit,xnewlim,nfvaccum,ispecies)
      endif
c         write(*,*)'isfull',isfull,cellvol
c         write(*,*)'calling subaccum'
      call subaccum(isuds,vlimit,xnewlim,ispecies)
      end
c****************************************************************
      subroutine partacinit(vlimit)
c Initialize uniform bins for accumulation of the particles.
      include 'ndimsdecl.f'
      include 'meshcom.f'
      include 'ptaccom.f'
      include 'plascom.f'
      real vlimit(2,ndimsmax)

c Initialization.
      do id=1,ndimsmax
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
      include 'ndimsdecl.f'
      include 'meshcom.f'
      include 'ptaccom.f'
      include 'plascom.f'
      real vlimit(2,ndimsmax)
      real xr(3*ndimsmax)

      do id=1,ndimsmax
c Assign velocities to bins.
         if(ivproj.eq.0)then
            v=xr(ndimsmax+id)
         else
c Projected.
            v=vproject(ndimsmax,id,xr,Bfield)
         endif

         if(v.lt.vlimit(1,id))v=vlimit(1,id)
         if(v.gt.vlimit(2,id))v=vlimit(2,id)
         ibin=nint(nptdiag*(.000005+0.99999*
     $        (v-vlimit(1,id))/(vlimit(2,id)-vlimit(1,id)))+0.5)
         if(ibin.lt.1.or.ibin.gt.nptdiag)then
            write(*,*)'For id ',id,' erroneous ibin',ibin,v
            write(*,*)nptdiag,vlimit(1,id),vlimit(2,id)
         endif
         fv(ibin,id)=fv(ibin,id)+1.
         if(csbin(1,1).ne.-1.)then
c Doing summed bin accumulation
            ibs=ibinmap(ibin,id)
            fsv(ibs,id)=fsv(ibs,id)+1.
         endif
c Assign positions to bins
         x=(xr(id)-xmeshstart(id))/(xmeshend(id)-xmeshstart(id))
         ibin=nint(0.500001+x*(nptdiag-.0002))
         if(ibin.lt.1 .or. ibin.gt.nptdiag)then
            write(*,*)'ibin=',ibin,id,x,xr(id)
         else
            px(ibin,id)=px(ibin,id)+1.
         endif
      enddo

      end
c**********************************************************************
      subroutine partsaccum(xlimit,vlimit,xnewlim,nfvaccum,ispecies)
      implicit none
      integer nfvaccum,ispecies
      include 'ndimsdecl.f'
      include 'partcom.f'

c Accumulate all particles in the xlimit range into velocity bins.
c If on entry xnewlim(2).eq.xnewlim(1) then adjust those limits.
c Spatial limits bottom-top, dimensions
      real xlimit(2,ndims),xnewlim(2,ndims)
c Velocity limits
      real vlimit(2,ndims)
      integer limadj,j,id
      real x

      if(xnewlim(1,1).eq.xnewlim(2,1))then
         limadj=1
      else
         limadj=0
      endif
      do j=iicparta(ispecies),iocparta(ispecies)
c Only for filled slots
         if(x_part(iflag,j).ne.0)then
            do id=1,ndims
               x=x_part(id,j)
               if(x.lt.xlimit(1,id).or.x.gt.xlimit(2,id))goto 12
               if(limadj.eq.1)then
c Adjust the xnewlim.
                  if(x.lt.xnewlim(1,id))xnewlim(1,id)=x
                  if(x.gt.xnewlim(2,id))xnewlim(2,id)=x
               endif
            enddo
            nfvaccum=nfvaccum+1
            call oneaccum(x_part(1,j),vlimit)
         endif
 12      continue
      enddo
c      write(*,*)'Accumulated',nfvaccum,' in',xlimit,' of',iocpart
c     $        ,' total'
c      if(limadj.eq.1)write(*,'(a,6f10.5)')' xnewlim=',xnewlim

             
      end
c**********************************************************************
      subroutine vlimitdeterm(vlimit,nvlist,ivproj,Bfield,ispecies)
      implicit none
      include 'ndimsdecl.f'
      include 'partcom.f'
      integer nvlist
      integer ivproj
      integer ispecies

c Determine the required vlimits for this particle distribution.
c On entry
c    x_part contains all the particles.
c    nvlist is the number of particles in the top and bottom bins.
c    ispecies is the species number to operate on => iic, ioc
c On exit
c    vlimit contains the velocity of the nvlist'th particle from the 
c           bottom and the top of the list.
c Velocity limits
      real vlimit(2,ndims)
c Velocity Sorting arrays
      integer nvlistmax
      parameter (nvlistmax=200)
      real vtlist(nvlistmax),vblist(nvlistmax)
      real Bfield(ndims)
      integer id,j
      real v,vproject
      external vproject

      if(nvlist.gt.nvlistmax)nvlist=nvlistmax
      do id=1,ndims
         do j=1,nvlist
            vblist(j)=vlimit(1,id)
            vtlist(j)=vlimit(2,id)
         enddo
         do j=iicparta(ispecies),iocparta(ispecies)
c Only for filled slots
            if(x_part(iflag,j).ne.0)then
               if(ivproj.eq.0)then
c Straight cartesian.
                  v=x_part(id+ndims,j)
               else
c Projected
                  v=vproject(ndims,id,x_part(1,j),Bfield)
               endif
c Insert the value v into its ordered place in vtlist, retaining the top
c nvlist values, and into vblist retaining the bottom-most nvlist values.
               call sorttoplimit(v,vtlist,nvlist)
               call sortbottomlimit(v,vblist,nvlist)
            endif
         enddo
c         write(*,*)vblist
c         write(*,*)vtlist
         if(vblist(nvlist).lt.vtlist(nvlist))then
c Take the vlimits to be the position of the nvlist velocity.
            vlimit(1,id)=vblist(nvlist)
            vlimit(2,id)=vtlist(nvlist)
         else
            write(*,*)'All velocities equal',vblist(nvlist),' for id',id
            vlimit(1,id)=-4.
            vlimit(2,id)=4.
            write(*,*)'vlimits',vlimit
         endif
      enddo

      end
c========================================================================
      subroutine sortbottomlimit(v,vlist,nvlist)
      real vlist(nvlist)
c Insert the value v into its ordered place in vlist, retaining the bottom
c nvlist values. 
c So on return vlist(1:nvlist) are the ordered bottommost v-values.
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
c On return vlist(1:nvlist) contain the (reversed) topmost nvlist values.
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
c Calculate the non-uniform bins for particle distribution accumlation.
c There are nsbins combined bins for each dimension. 
c These are in common subdiag:
c ibinmap(nptdiag,ndimsmax) is the map from uniform to combined bins
c pointing to an nsbins bin for each of the nptdiag uniform bins.
c vsbin(nsbins,ndimsmax) is the center velocity of the combined bins
c csbin(nsbins,ndimsmax) is the number of fine bins in each combined bin
c fsv is the sum of fv in each combined bin during the initial
c     accumulation and bin calculation.
c vhbin(0:nsbins,ndimsmax) is the histogram boundaries of the combined bins
c If nptdiag.le.nsbins, implying no compression of the data, then
c a simple identification of summed and uniform bins is used so that the
c summed bins are actually uniform. The lengths must be set in ptaccom.f

      include 'ndimsdecl.f'
      include 'meshcom.f'
      include 'ptaccom.f'

c 
      ifixed=0
      do id=1,ndimsmax
         cumfv(0,id)=0.
         do j=1,nsbins
            vsbin(j,id)=0.
            csbin(j,id)=0.
            fsv(j,id)=0.
         enddo
         dv=(vdiag(nptdiag,id)-vdiag(1,id))/(nptdiag-1)
         vhbin(0,id)=vdiag(1,id)-dv*0.5
         ib=1
         if(nsbins.lt.nptdiag)then
c Normal summed binning compression.
            do k=1,nptdiag
               cumfv(k,id)=cumfv(k-1,id)+fv(k,id)/float(nfvaccum)
c Prevent rounding overflow here rather than by increments:
               cumfv(k,id)=min(1.,cumfv(k,id))
c            write(*,*)id,k,fv(k,id),cumfv(k,id),ib
c This linear mapping does not work well.
c            ib=1+ int(cumfv(k,id)*(nsbins)*(.99999))
               bx=float(ib)/nsbins
c cubic progression.
c               cfn=1.0001*(3.*bx**2-2.*bx**3)
c quintic progression puts more bins further out.
c            cfn=1.00002*bx**3*(10.-15.*bx+6*bx**2) -.00001
c Does not help to reduce the number of zeroes in 1.0002.
c septic progression is broadest, but still good.
               cfn=(((-20.*bx+70.)*bx-84.)*bx+35.)*bx**4
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
                  ifixed=ifixed+1
c               write(*,*)'Fixing csbin zero',k,id,fv(k,id)
c     $              ,vhbin(k,id)
                  vsbin(k,id)=vhbin(k,id)
               else
                  vsbin(k,id)=vsbin(k,id)/csbin(k,id)
                  csbin(k,id)=1./csbin(k,id)
               endif
            enddo
         else
c Ignore the cumulative calculation and do uniform mapping
            if(id.eq.1)write(*,*)'****Uniform Bins for nptdiag,nsbins='
     $           ,nptdiag,nsbins
            do k=1,nptdiag
               ibinmap(k,id)=k
               vsbin(k,id)=vdiag(k,id)
               vhbin(k,id)=vhbin(0,id)+(k-1)*dv
               csbin(k,id)=1.
               fsv(k,id)=fv(k,id)
            enddo
         endif
      enddo
      if(ifixed.gt.0)then
         write(*,*)'Bincalc had to fix some zero csbins.'
      endif
c      write(*,*)'Bincalc has chosen',nsbins,' bin placement.'
c      write(*,*)' vsbin',vsbin
c      write(*,*)(vhbin(k,1),k=0,nsbins)
c      write(*,*)' csbin',csbin
c      write(*,*)' fsv  ',fsv
c      write(*,*)' ibinmap',ibinmap
      end
c**********************************************************************
      subroutine subaccum(isuds,vlimit,xnewlim,ispecies)
      implicit none
      include 'ndimsdecl.f'
      include 'partcom.f'
      integer isuds(ndims)
      integer ispecies
c Spatial limits bottom-top, dimensions
      real xnewlim(2,ndims)
c Velocity limits
      real vlimit(2,ndims)
c Subbin addressing.
      integer isind(ndims)
c Local storage
      real xsf(ndims)
c The equal bins fill the ranges xnewlim with index lengths isuds.
      integer id,j,ip,ipfindex
      external ipfindex
      do id=1,ndims
         xsf(id)=.999999/(xnewlim(2,id)-xnewlim(1,id))
      enddo
      do j=iicparta(ispecies),iocparta(ispecies)
c Only for filled slots
         if(x_part(iflag,j).ne.0)then
            do id=1,ndims
c Calculate the spatial sub-bin pointer.
               isind(id)=1+int(isuds(id)*
     $              (x_part(id,j)-xnewlim(1,id))*xsf(id))
c If we are outside the range skip this particle.
               if(isind(id).gt.isuds(id).or.isind(id).le.0)goto 1
            enddo
            ip=ipfindex(ndims,isuds,isind)
c            if(ip.gt.1000)write(*,*)isuds,isind,ip
            if(ip.lt.0)write(*,*)isind,ip
            call subvaccum(x_part(1,j),vlimit,ip+1)
         endif
 1       continue
      enddo
             
      end
c*****************************************************************
      subroutine subvaccum(xr,vlimit,ip)
c Accumulate a particle into velocity bin corresponding to position ip
      include 'ndimsdecl.f'
      include 'meshcom.f'
      include 'ptaccom.f'
      include 'plascom.f'
      real vlimit(2,ndimsmax)
      real xr(3*ndimsmax)

      do id=1,ndimsmax
c Assign velocities to bins.
         if(ivproj.eq.0)then
            v=xr(ndimsmax+id)
         else
            v=vproject(ndimsmax,id,xr,Bfield)
         endif
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
               write(*,*)'Bin address corruption ibs,id,ip,nsub_tot'
     $              ,ibs,id,ip,nsub_tot
            endif
            fvx(ibs,id,ip)=fvx(ibs,id,ip)+1
         endif
      enddo
      denfvx(ip)=denfvx(ip)+1
      end
c*******************************************************************
      subroutine fvxinit(xnewlim,cellvol,ibset)
      include 'ndimsdecl.f'
      include 'meshcom.f'
      include 'ptaccom.f'
      real xnewlim(2,ndimsmax)

      isfull(1)=nsub_i
      isfull(2)=nsub_j
      isfull(3)=nsub_k
      if(ibset.eq.0)then
         do i=1,ndimsmax
            isuds(i)=isfull(i)
         enddo
      endif
      do ip=1,nsub_tot
         denfvx(ip)=0.
         do id=1,ndimsmax
            do ibs=1,nsbins
               fvx(ibs,id,ip)=0.
            enddo
         enddo
      enddo
      cellvol=1.
      do id=1,ndimsmax
c         cellvol=cellvol*(xnewlim(2,id)-xnewlim(1,id))/isfull(id)
         cellvol=cellvol*(xnewlim(2,id)-xnewlim(1,id))/isuds(id)
c         write(*,*)xnewlim(2,id),xnewlim(1,id)
      enddo
c      write(*,*)'cellvol=',cellvol
      end
c**********************************************************************
      subroutine distread(xlimit,vlimit,xnewlim,name,cellvol)

      include 'ndimsdecl.f'
      include 'meshcom.f'
      include 'ptaccom.f'
      real xlimit(2,ndimsmax),vlimit(2,ndimsmax),xnewlim(2,ndimsmax)
      character*(*) name

      open(25,file=name,status='old',form='unformatted',err=101)
      read(25)nptdiagfile,ndimsmaxfile
      read(25)(xlimit(1,j),xlimit(2,j),vlimit(1,j),vlimit(2,j),
     $     j=1,ndimsmax)
      read(25)((xdiag(i,j),px(i,j),i=1,nptdiag),j=1,ndimsmax)
      read(25)((vdiag(i,j),fv(i,j),i=1,nptdiag),j=1,ndimsmax)
      read(25)(xnewlim(1,j),xnewlim(2,j),j=1,ndimsmax)
      read(25)nsbf,(isfull(i),i=1,ndimsmax),cellvol
      read(25)(isuds(i),i=1,ndimsmax)
      isftot=1
      do i=1,ndimsmax
         isftot=isftot*isfull(i)
      enddo
c Check if the nsbins is correct and there's enough storage.
      if(nsbf.ne.nsbins)goto 103
      if(isftot.gt.nsub_tot)goto 103
      read(25)(((fvx(i,j,k),i=1,nsbins),j=1,ndimsmax),k=1,isftot)
      read(25)(denfvx(k),k=1,isftot)
      read(25)((vsbin(i,j),csbin(i,j),fsv(i,j),i=1,nsbins),
     $     (vhbin(i,j),i=0,nsbins),j=1,ndimsmax)
      read(25)((ibinmap(i,j),i=1,nptdiag),j=1,ndimsmax)
      close(25)

      return
 101  write(*,*)'Error opening file:',name(1:lentrim(name))
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
      include 'ndimsdecl.f'
      include 'meshcom.f'
      include 'ptaccom.f'
      real xlimit(2,ndimsmax),vlimit(2,ndimsmax),xnewlim(2,ndimsmax)
      character*(*) name

      open(25,file=name,status='unknown',err=101)
      close(25,status='delete')
      open(25,file=name,status='new',form='unformatted',err=101)
      write(25)nptdiag,ndimsmax
      write(25)(xlimit(1,j),xlimit(2,j),vlimit(1,j),vlimit(2,j),
     $     j=1,ndimsmax)
      write(25)((xdiag(i,j),px(i,j),i=1,nptdiag),j=1,ndimsmax)
      write(25)((vdiag(i,j),fv(i,j),i=1,nptdiag),j=1,ndimsmax)
c      write(*,'(2f10.4)')((vdiag(i,j),fv(i,j),i=1,5),j=1,ndimsmax)
      write(25)(xnewlim(1,j),xnewlim(2,j),j=1,ndimsmax)
      write(25)nsbins,(isfull(i),i=1,ndimsmax),cellvol
      write(25)(isuds(i),i=1,ndimsmax)
      isftot=1
      do i=1,ndimsmax
         isftot=isftot*isfull(i)
      enddo
      write(25)(((fvx(i,j,k),i=1,nsbins),j=1,ndimsmax),k=1,isftot)
      write(25)(denfvx(k),k=1,isftot)
      write(25)((vsbin(i,j),csbin(i,j),fsv(i,j),i=1,nsbins),
     $     (vhbin(i,j),i=0,nsbins),j=1,ndimsmax)
      write(25)((ibinmap(i,j),i=1,nptdiag),j=1,ndimsmax)
      close(25)

      return
 101  write(*,*)'Error opening file:',name(1:lentrim(name))
      close(25,status='delete')

      end
c******************************************************************
c 3-D only code.
      subroutine pltsubdist(icellin,jcellin,kcellin,vlimit,xnewlim
     $     ,cellvol,ndfirst,ndlast)
      implicit none
      integer icellin,jcellin,kcellin
      integer ndfirst,ndlast
      include 'ndimsdecl.f'
      include 'meshcom.f'
      include 'partcom.f'
      include 'ptaccom.f'
      include 'plascom.f'
c Distributions in ptaccom.f
      real vlimit(2,ndimsmax),xnewlim(2,ndimsmax)
      real cellvol,ftot,vtot,v2tot
      real fvplt(nsbins)
      integer ip,id,iv
      integer ip3index
      external ip3index
      integer lentrim
      external lentrim
      character*100 string
      integer icell,jcell,kcell
      real wicell,wjcell,wkcell
      integer jicell,jjcell,jkcell,ii,jj,kk,ip3,idw
      integer iup,isw
      real xylim(4)
      integer ips,itrace
      logical loverplot
c      data ndfirst/3/ndlast/3/
      data ips/0/itrace/0/loverplot/.false./
      data idw/0/
      data iup/1/jicell/0/jjcell/0/jkcell/0/
c------------------------------------------
      if(ndims.ne.3)then
         write(*,*)'Skipping pltsubdist 3-D only code for ndims=',ndims
         write(*,*)'But it might work ok if cellvol and xnewlim are ok'
         return
      endif
      isw=0
      icell=icellin
      jcell=jcellin
      kcell=kcellin
      write(*,*)'Adjust cell position with arrow keys.'
      write(*,*)'Toggle integration over direction with x,y,z.'
      write(*,*)'Adjust axis ends with a'
 1    ip=ip3index(isuds,icell,jcell,kcell)+1
      call multiframe(ndlast-ndfirst+1,1,2)
      wicell=(xnewlim(2,1)-xnewlim(1,1))/isuds(1)
      wjcell=(xnewlim(2,2)-xnewlim(1,2))/isuds(2)
      wkcell=(xnewlim(2,3)-xnewlim(1,3))/isuds(3)
      if(ivproj.eq.0)then
      write(string,'(''Cell'',3i3,''  x=('',f6.2,'','',f6.2,'')'//
     $     ' y=('',f6.2,'','',f6.2,'') z=('',f6.2,'','',f6.2,'')'')')
     $     icell*(1-jicell),jcell*(1-jjcell),kcell*(1-jkcell)
     $     ,xnewlim(1,1)+(icell-1)*wicell*(1-jicell)
     $        ,xnewlim(1,1)+wicell*(icell*(1-jicell)+jicell*isuds(1))
     $     ,xnewlim(1,2)+(jcell-1)*wjcell*(1-jjcell)
     $        ,xnewlim(1,2)+wjcell*(jcell*(1-jjcell)+jjcell*isuds(2))
     $     ,xnewlim(1,3)+(kcell-1)*wkcell*(1-jkcell)
     $        ,xnewlim(1,3)+wkcell*(kcell*(1-jkcell)+jkcell*isuds(3))
      else
 51      format('Cell',3i3,'  x=(',f6.2,',',f6.2,') y=(',f6.2,',',f6.2
     $     ,') z=(',f6.2,',',f6.2,') projected:(',3f6.2,')')
      write(string,51)
     $     icell,jcell,kcell
     $     ,xnewlim(1,1)+(icell-1)*wicell,xnewlim(1,1)+icell*wicell
     $     ,xnewlim(1,2)+(jcell-1)*wjcell,xnewlim(1,2)+jcell*wjcell
     $     ,xnewlim(1,3)+(kcell-1)*wkcell,xnewlim(1,3)+kcell*wkcell
     $     ,Bfield
      endif
c      write(*,'(a,$)')'v, T in bin'
      do id=ndfirst,ndlast
c Do correct scaling:
         do iv=1,nsbins
c Alternatively we combine bins in possibly multiple
c dimensions, if jicell, jjcell, or jkcell are non-zero.
            fvplt(iv)=0.
            do ii=0,(isuds(1)-1)*jicell
               do jj=0,(isuds(2)-1)*jjcell
                  do kk=0,(isuds(3)-1)*jkcell
                     ip3=ip3index(isuds,1+ii+(1-jicell)*(icell-1),
     $                    1+jj+(1-jjcell)*(jcell-1),
     $                    1+kk+(1-jkcell)*(kcell-1))+1
            fvplt(iv)=fvplt(iv)+fvx(iv,id,ip3)*csbin(iv,id)
     $              *nptdiag/(vlimit(2,id)-vlimit(1,id))/cellvol
                  enddo
               enddo 
            enddo
         enddo
         ftot=0.
         vtot=0.
         v2tot=0.
         do iv=1,nsbins
            vtot=vtot+vsbin(iv,id)*fvx(iv,id,ip)
            v2tot=v2tot+vsbin(iv,id)**2*fvx(iv,id,ip)
            ftot=ftot+fvx(iv,id,ip)
         enddo
c         write(*,'('','',2f10.5,$)')vtot/ftot,sqrt(v2tot/ftot-(vtot/ftot
c     $        )**2)
         if(loverplot.and.ndlast-ndfirst.eq.0)then
            itrace=itrace+1
            call color(itrace)
            call polymark(vsbin(1,id),fvplt,nsbins,1)
         else
c            call automark(vsbin(1,id),fvplt,nsbins,1)
c            call autoinit(vsbin(1,id),fvplt,nsbins)
c            write(*,*)'manautoinit',isw,xylim
            call manautoinit(vsbin(1,id),fvplt,nsbins,isw,xylim(1)
     $           ,xylim(2),xylim(3),xylim(4))
            call axis()
         endif
         if(btest(idw,id-1))then
c         if(id.eq.3)then
            write(*,*)'#',icell*(1-jicell),jcell*(1-jjcell),kcell*(1
     $        -jkcell),' dimension',id
            write(*,*)nsbins
            write(*,'(2g14.6)')(vsbin(kk,id),fvplt(kk),kk=1,nsbins)
c            write(*,*)
         endif

         if(id.eq.ndfirst)call boxtitle(string(1:lentrim(string)))
         if(ivproj.eq.0)then
            write(string,'(a,i1,a)')'f(v!d',id,'!d)'
         else
            string=' '
         endif
         if(id.eq.1)then
            if(string(1:1).eq.' ')
     $           write(string,'(a)')'f(v!d!A|!@!d)'
         elseif(id.eq.2)then
            if(string(1:1).eq.' ')
     $           write(string,'(a)')'2!Ap!@v!d!A`!@!df(v!d!A`!@!d)'
         else
            if(string(1:1).eq.' ')
     $           write(string,'(a,i1,a)')'f(v!dz!d)'
         endif
         call axlabels(' ',string(1:lentrim(string)))
         call winset(.true.)
            call polymark(vsbin(1,id),fvplt,nsbins,1)
c         call polybox(vhbin(0,id),fvx(1,id,ip),nsbins)
         call polybox(vhbin(0,id),fvplt,nsbins)
      enddo
      call winset(.false.)
      call axlabels('Velocity',' ')
      call eye3d(ip)
c      write(*,*)'ip',ip
c Increment one of the dimensions using arrow keys.
      if(ip.eq.65361)then
         icell=mod(icell-1+isuds(1)+iup,isuds(1))+1
      elseif(ip.eq.65362)then
         kcell=mod(kcell,isuds(3))+1
      elseif(ip.eq.65364)then
         kcell=mod(kcell-2+isuds(3),isuds(3))+1
      elseif(ip.eq.65363)then
         jcell=mod(jcell-1+isuds(2)+iup,isuds(2))+1
      elseif(ip.eq.65505)then
c ffe1 shift_L=2^16-31=65505
         iup=-iup
c Control averaging over dimensions
      elseif(ip.eq.ichar('x'))then
         jicell=1-jicell
      elseif(ip.eq.ichar('y'))then
         jjcell=1-jjcell
      elseif(ip.eq.ichar('z'))then
         jkcell=1-jkcell
c Control overplotting
      elseif(ip.eq.ichar('o'))then
         if(loverplot)itrace=0
         loverplot=.not.loverplot
      elseif(ip.eq.ichar('p'))then
         if(ips.eq.0)then
            write(*,*)'Postscript output on'
            call pfset(3)
            ips=3
         else
            write(*,*)'Postscript output off'
            ips=0
            call prtend()
            call pfset(0)
         endif
      elseif(ip.eq.ichar('a'))then
c------------------------------
c Control axis settings.
c Get another key input
         write(*,*)'Axis control.'
     $        ,' Type axis-end switch 1-4, value, return.'
 200     call eye3d(ip)
c         write(*,*)ip
         if((ip.ge.49 .and. ip.le.52))then
c            write(*,*)char(ip)
c Set the isw accordingly.
            isw=ieor(isw,2**(ip-49))
c            isw=ibset(isw,ip-49)
            write(*,*)'isw set to',isw,'  enter end value:'
c Get the real number
            call eyereadreal(xylim(ip-48))
c            write(*,*)'eyeread return',ip-48,xylim(ip-48)
         else
            write(*,*)'Not a number in [1,4]:', char(ip)
            goto 200
         endif
c------------------------------
      else
c         write(*,*)ip
         return
      endif
      goto 1
c--------------------------------------
      end
c****************************************************************
      subroutine eyereadreal(x)
      real x
      logical ldec

 2    ldec=.false.
      x=0
      sign=1.
      decfac=1
 1    call eye3d(i)
      if((i.ge.48 .and. i.le.57))then
         if(ldec)then
            x=x+sign*(i-48.)/10**decfac
            decfac=decfac+1
         else
            x=x*10.+sign*(i-48)
         endif
         write(*,*)x
c,'  decfac=',decfac
         goto 1
      else
         if(i.eq.46)then
            ldec=.not.ldec
            goto 1
         elseif(i.eq.45)then
            x=-x
            sign=-sign
c            write(*,*)x
            goto 1
         elseif(i.ne.65293)then
c Return indicates end. Anything else starts over.           
c Backspace is 65288
            write(*,*)'Starting Read Real Over:',i,char(i)
            goto 2
         endif
      endif
      end
c****************************************************************
      subroutine fsvzero()
      include 'ndimsdecl.f'
      include 'meshcom.f'
      include 'partcom.f'
      include 'ptaccom.f'
      do id=1,ndimsmax
         do jj=1,nsbins
            fsv(jj,id)=0.
         enddo
      enddo
      end
c****************************************************************
      real function vproject(ndims,id,xr,Bfield)
c Return the id'th component of the "Projected" velocity xr(ndims+ii)
c relative to the direction cosines Bfield(ndims). 
c id=1 Parallel, id=2 Perpendicular, id=3 3rd (=id'th) component.
      implicit none
      integer ndims,id
      real xr(ndims),Bfield(ndims) 
      real v,v2
      integer ii
c Projected. Very clumsy.
      v=0
      do ii=1,ndims
         v=v+xr(ii+ndims)*Bfield(ii)
      enddo
      if(id.eq.2)then
c Perpendicular
         v2=0
         do ii=1,ndims
            v2=v2+(xr(ii+ndims)-v*Bfield(ii))**2
         enddo
         v=sqrt(v2)
      elseif(id.eq.3)then
c Z-component.
         v=xr(id+ndims)
      endif
      vproject=v
      end

c********************************************************************
c Convert 3-D indices into pointer. Only used in pltsubdist.
      function ip3index(ifull,i,j,k)
      integer ifull(3),i,j,k 
      ip3index=0
      ip3index=(k-1)+ip3index*ifull(3)
      ip3index=(j-1)+ip3index*ifull(2)
      ip3index=(i-1)+ip3index*ifull(1)
      end
