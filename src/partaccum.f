!********************************************************************
      subroutine partdistup(xlimit,vlimit,xnewlim,cellvol,myid,ibinit
     $     ,ispecies)
! Routine to update the particle distribution diagnostics for using the
! current particle information for species ispecies.
!
! If this is called for the first time, i.e. with cellvol.le.0, then
! calculate the bins for accumulating the distribution. Also, in this
! case, if cellvol=-1 project in the direction of Bfield, else cellvol=0
! just use cartesian components.  For this first call, in a
! multiprocessor situation, we all-reduce the accumulated fine
! distribution fv, nvaccum (and fsv) on which the binning is calculated.
! Thus the bin calculation is based upon the total particles, not just
! the particles for the head node.  However, for use with multiple
! particle files by partexamine, only the first file read is used to
! determine the bin mapping. [Otherwise we would have to read all files
! a second time for the subbinning.]
!
! Otherwise if cellvol.gt.0 just accumulate to bins.
      implicit none
      include 'ndimsdecl.f'
      include 'plascom.f'
      include 'meshcom.f'
      include 'partcom.f'
      include 'ptaccom.f'
      real xlimit(2,ndims),vlimit(2,ndims),xnewlim(2,ndims)
      real cellvol
      integer myid,ibinit
      integer nvlist
      integer ispecies
      integer id

! Use the first occasion to establish the accumulation range.
! Indicated by cellvolume le 0.
      if(cellvol.le.0.)then
! Decide whether to project velocity in the direction of Bfield.
         if(cellvol.eq.-1)then
            ivproj=1
         else
            ivproj=0
         endif
! This ought to be determined based on the number of samples.
! But xlimits mean that's problematic.
!            nvlist=100
! A small nvlist means we go right to the edge of actual particles.
         nvlist=5
         call vlimitdeterm(vlimit,nvlist,ivproj,Bfield,ispecies)
         call minmaxreduce(ndims,vlimit)
         do id=1,ndims
            call adjustlimit(nptdiag,vlimit(1,id))
         enddo
         if(myid.eq.0)write(*,'('' Velocity limits:'',6f8.3)')
     $        vlimit
! Indicate csbin not initialized and start initialization
         csbin(1,1)=-1.
         call partacinit(vlimit)
!         write(*,*)'partsaccum0',fsv(10,1),ifsv(10,1)

! Do the accumulation for this file up to maximum relevant slot. 
         nfvaccum=0
         call partsaccum(xlimit,vlimit,xnewlim,nfvaccum,ispecies)
         if(myid.eq.0)write(*,'(a,i8,a,i8,a,6f6.1)')' Accumulated'
     $        ,nfvaccum,' of',ioc_part,' total in',xlimit
! Reduce back the data for MPI cases.
         call ptdiagreduce()
         call minmaxreduce(ndims,xnewlim)
         if(myid.eq.0)write(*,'(a,i8,a,6f8.3)')' Reduced',nfvaccum
     $        ,' Xnewlim=',xnewlim
! Should do this only the first time.
         call bincalc()
! Clean up by reinitializing everything.         
         call fvxinit(xnewlim,cellvol,ibinit)
         call fsvzero()
         call partacinit(vlimit)
      else
         call partsaccum(xlimit,vlimit,xnewlim,nfvaccum,ispecies)
      endif
      call subaccum(isuds,vlimit,xnewlim,ispecies)
      end
!****************************************************************
      subroutine partacinit(vlimit)
! Initialize uniform bins for accumulation of the particles.
      include 'ndimsdecl.f'
      include 'meshcom.f'
      include 'ptaccom.f'
      include 'plascom.f'
      real vlimit(2,ndims)

! Initialization.
      do id=1,ndims
         do i=1,nptdiag
            fv(i,id)=0.
            ifv(i,id)=0
            px(i,id)=0.
            vdiag(i,id)=vlimit(1,id)
     $           +(vlimit(2,id)-vlimit(1,id))*(i-0.5)/nptdiag
            xdiag(i,id)=xmeshstart(id)+(i-0.5)*
     $           (xmeshend(id)-xmeshstart(id))/(nptdiag)
         enddo
!         write(*,*)'Position cell-center range',
!     $        id,xdiag(1,id),xdiag(nptdiag,id)
!         write(*,*)'Velocity cell-center range',
!     $        id,vdiag(1,id),vdiag(nptdiag,id)
      enddo
      end
!*****************************************************************
      subroutine oneaccum(xr,vlimit)
! Accumulate a particle into velocity bins.
      include 'ndimsdecl.f'
      include 'meshcom.f'
      include 'ptaccom.f'
      include 'plascom.f'
      real vlimit(2,ndims)
      real xr(3*ndims)

      do id=1,ndims
! Assign velocities to bins.
         if(ivproj.eq.0)then
            v=xr(ndims+id)
         else
! Projected.
            v=vproject(ndims,id,xr,Bfield)
         endif

         if(v.lt.vlimit(1,id))v=vlimit(1,id)
         if(v.gt.vlimit(2,id))v=vlimit(2,id)
         ibin=nint(nptdiag*(.000005+0.99999*
     $        (v-vlimit(1,id))/(vlimit(2,id)-vlimit(1,id)))+0.5)
         if(ibin.lt.1.or.ibin.gt.nptdiag)then
            write(*,*)'For id ',id,' erroneous ibin',ibin,v
            write(*,*)nptdiag,vlimit(1,id),vlimit(2,id)
         endif
         ifv(ibin,id)=ifv(ibin,id)+1
!         fv(ibin,id)=real(ifv(ibin,id)) ! Avoid real rounding loss.
         fv(ibin,id)=fv(ibin,id)+1.
         if(csbin(1,1).ne.-1.)then
! Doing summed bin accumulation
            ibs=ibinmap(ibin,id)
! Avoid floating underflow miscounting by using integers.
            ifsv(ibs,id)=ifsv(ibs,id)+1
            fsv(ibs,id)=real(ifsv(ibs,id))
!            fsv(ibs,id)=fsv(ibs,id)+1.
         endif
! Assign positions to bins
         x=(xr(id)-xmeshstart(id))/(xmeshend(id)-xmeshstart(id))
         ibin=nint(0.500001+x*(nptdiag-.0002))
         if(ibin.lt.1 .or. ibin.gt.nptdiag)then
            write(*,*)'ibin=',ibin,id,x,xr(id)
         else
            px(ibin,id)=px(ibin,id)+1.
         endif
      enddo

      end
!**********************************************************************
      subroutine partsaccum(xlimit,vlimit,xnewlim,nfvaccum,ispecies)
      implicit none
      integer nfvaccum,ispecies
      include 'ndimsdecl.f'
      include 'partcom.f'

! Accumulate all particles in the xlimit range into velocity bins.
! If on entry xnewlim(2).eq.xnewlim(1) then adjust those limits.
! Spatial limits bottom-top, dimensions
      real xlimit(2,ndims),xnewlim(2,ndims)
! Velocity limits
      real vlimit(2,ndims)
      integer limadj,j,id
      real x

      if(xnewlim(1,1).eq.xnewlim(2,1))then
         limadj=1
      else
         limadj=0
      endif
      do j=iicparta(ispecies),iocparta(ispecies)
! Only for filled slots
         if(x_part(iflag,j).ne.0)then
            do id=1,ndims
               x=x_part(id,j)
               if(x.lt.xlimit(1,id).or.x.gt.xlimit(2,id))goto 12
               if(limadj.eq.1)then
! Adjust the xnewlim.
                  if(x.lt.xnewlim(1,id))xnewlim(1,id)=x
                  if(x.gt.xnewlim(2,id))xnewlim(2,id)=x
               endif
            enddo
            nfvaccum=nfvaccum+1
            call oneaccum(x_part(1,j),vlimit)
         endif
 12      continue
      enddo

             
      end
!**********************************************************************
      subroutine vlimitdeterm(vlimit,nvlist,ivproj,Bfield,ispecies)
      implicit none
      include 'ndimsdecl.f'
      include 'partcom.f'
      integer nvlist
      integer ivproj
      integer ispecies

! Determine the required vlimits for this particle distribution.
! On entry
!    x_part contains all the particles.
!    nvlist is the number of particles in the top and bottom bins.
!    ispecies is the species number to operate on => iic, ioc
! On exit
!    vlimit contains the velocity of the nvlist'th particle from the 
!           bottom and the top of the list.
! Velocity limits
      real vlimit(2,ndims)
! Velocity Sorting arrays
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
! Only for filled slots
            if(x_part(iflag,j).ne.0)then
               if(ivproj.eq.0)then
! Straight cartesian.
                  v=x_part(id+ndims,j)
               else
! Projected
                  v=vproject(ndims,id,x_part(1,j),Bfield)
               endif
! Insert the value v into its ordered place in vtlist, retaining the top
! nvlist values, and into vblist retaining the bottom-most nvlist values.
               call sorttoplimit(v,vtlist,nvlist)
               call sortbottomlimit(v,vblist,nvlist)
            endif
         enddo
!         write(*,*)vblist
!         write(*,*)vtlist
         if(vblist(nvlist).lt.vtlist(nvlist))then
! Take the vlimits to be the position of the nvlist velocity.
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
!========================================================================
      subroutine sortbottomlimit(v,vlist,nvlist)
      real vlist(nvlist)
! Insert the value v into its ordered place in vlist, retaining the bottom
! nvlist values. 
! So on return vlist(1:nvlist) are the ordered bottommost v-values.
!      write(*,*)'bot',v,nvlist,vlist
      if(v.ge.vlist(nvlist))return
      if(v.ge.vlist(1))then
         i1=1
         i2=nvlist
! Find my position by bisection.
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
! Here i2 is the position of list value just greater than v.
!      write(*,*)'botend',i1,i2,i,vlist(i2),v
      do i=nvlist,i2+1,-1
         vlist(i)=vlist(i-1)
      enddo
      vlist(i2)=v

      end
!*********************************************************************
      subroutine sorttoplimit(v,vlist,nvlist)
      real vlist(nvlist)
! Insert the value v into its reverse-ordered place in vlist, retaining
! the top nvlist values. 
! On return vlist(1:nvlist) contain the (reversed) topmost nvlist values.
!      write(*,*)'topstart',nvlist,v,vlist
      if(v.le.vlist(nvlist))return
      if(v.le.vlist(1))then
         i1=1
         i2=nvlist
! Find my position by bisection.
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
! Here i2 is the position of list value just less than v.
!      write(*,*)'topend',i1,i2,i,vlist(i2),v
      do i=nvlist,i2+1,-1
         vlist(i)=vlist(i-1)
      enddo
      vlist(i2)=v
      end

!***********************************************************************
      subroutine bincalc()
! Calculate the non-uniform bins for particle distribution accumlation.
! There are nsbins combined bins for each dimension. 
! These are in common subdiag:
! ibinmap(nptdiag,ndimsmax) is the map from uniform to combined bins
! pointing to an nsbins bin for each of the nptdiag uniform bins.
! vsbin(nsbins,ndimsmax) is the center velocity of the combined bins
! csbin(nsbins,ndimsmax) is the number of fine bins in each combined bin
! fsv is the sum of fv in each combined bin during the initial
!     accumulation and bin calculation.
! vhbin(0:nsbins,ndimsmax) is the histogram boundaries of the combined bins
! If nptdiag.le.nsbins, implying no compression of the data, then
! a simple identification of summed and uniform bins is used so that the
! summed bins are actually uniform. The lengths must be set in ptaccom.f

      include 'ndimsdecl.f'
      include 'meshcom.f'
      include 'ptaccom.f'

! 
      ifixed=0
      do id=1,ndims
         cumfv(0,id)=0.
         do j=1,nsbins
            vsbin(j,id)=0.
            csbin(j,id)=0.
            fsv(j,id)=0.
            ifsv(j,id)=0
         enddo
         dv=(vdiag(nptdiag,id)-vdiag(1,id))/(nptdiag-1)
         vhbin(0,id)=vdiag(1,id)-dv*0.5
         ib=1
         if(nsbins.lt.nptdiag)then
! Normal summed binning compression.
            do k=1,nptdiag
               cumfv(k,id)=cumfv(k-1,id)+fv(k,id)/float(nfvaccum)
! Prevent rounding overflow here rather than by increments:
               cumfv(k,id)=min(1.,cumfv(k,id))
!            write(*,*)id,k,fv(k,id),cumfv(k,id),ib
! This linear mapping does not work well.
!            ib=1+ int(cumfv(k,id)*(nsbins)*(.99999))
               bx=float(ib)/nsbins
! cubic progression.
!               cfn=1.0001*(3.*bx**2-2.*bx**3)
! quintic progression puts more bins further out.
!            cfn=1.00002*bx**3*(10.-15.*bx+6*bx**2) -.00001
! Does not help to reduce the number of zeroes in 1.0002.
! septic progression is broadest, but still good.
               cfn=(((-20.*bx+70.)*bx-84.)*bx+35.)*bx**4
!            write(*,*)'cfn',cfn,cumfv(k,id)
               if(cumfv(k,id).gt.cfn)then
! find the histogram bin-boundary.
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
! Also accumulate this data into the summed bins
               fsv(ib,id)=fsv(ib,id)+fv(k,id)
            enddo
! Fix any summed bins that are off the end, and hence zero width.
            do ibk=ib,nsbins
               vhbin(ibk,id)=vhbin(0,id)+nptdiag*dv
            enddo
! Now cumfv is the cumulative probability distribution (to 1.0) over the
! uniform bins, and ibinmap maps those bins to nonuniform bins.  The
! uniform bin centers are in vdiag. Each non-uniform bin is composed of
! the sum of the uniform bins that map to it. Its center is therefore at
! the centroid of those bins = vsbin. 
! The number of uniform bins that they contain is csbin.
! But we invert csbin so that a multiplication is needed and zeros are 
! corrected for. 
            do k=1,nsbins
               if(csbin(k,id).eq.0)then
                  ifixed=ifixed+1
!               write(*,*)'Fixing csbin zero',k,id,fv(k,id)
!     $              ,vhbin(k,id)
                  vsbin(k,id)=vhbin(k,id)
               else
                  vsbin(k,id)=vsbin(k,id)/csbin(k,id)
                  csbin(k,id)=1./csbin(k,id)
               endif
            enddo
         else
! Ignore the cumulative calculation and do uniform mapping
            if(id.eq.1)write(*,*)'****Uniform Bins for nptdiag,nsbins='
     $           ,nptdiag,nsbins
            do k=1,nptdiag
               ibinmap(k,id)=k
               vsbin(k,id)=vdiag(k,id)
!               vhbin(k,id)=vhbin(0,id)+(k-1)*dv
               vhbin(k,id)=vhbin(0,id)+(k)*dv
               csbin(k,id)=1.
               fsv(k,id)=fv(k,id)
            enddo
         endif
      enddo
      if(ifixed.gt.0)then
         write(*,*)'Bincalc had to fix some zero csbins.'
      endif
!      write(*,*)'Bincalc has chosen',nsbins,' bin placement.'
!      write(*,*)' vsbin',vsbin
!      write(*,*)(vhbin(k,1),k=0,nsbins)
!      write(*,*)' csbin',csbin
!      write(*,*)' fsv  ',fsv
!      write(*,*)' ibinmap',ibinmap
      end
!**********************************************************************
      subroutine subaccum(isuds,vlimit,xnewlim,ispecies)
      implicit none
      include 'ndimsdecl.f'
      include 'partcom.f'
      integer isuds(ndims)
      integer ispecies
! Spatial limits bottom-top, dimensions
      real xnewlim(2,ndims)
! Velocity limits
      real vlimit(2,ndims)
! Subbin addressing.
      integer isind(ndims)
! Local storage
      real xsf(ndims)
! The equal bins fill the ranges xnewlim with index lengths isuds.
      integer id,j,ip,ipfindex
      external ipfindex
      do id=1,ndims
         xsf(id)=.999999/(xnewlim(2,id)-xnewlim(1,id))
      enddo
      do j=iicparta(ispecies),iocparta(ispecies)
! Only for filled slots
         if(x_part(iflag,j).ne.0)then
            do id=1,ndims
! Calculate the spatial sub-bin pointer, corrected 21 Aug 19.
               isind(id)=nint(0.5+isuds(id)*
     $              (x_part(id,j)-xnewlim(1,id))*xsf(id))
! If we are outside the range skip this particle.
               if(isind(id).gt.isuds(id).or.isind(id).le.0)goto 1
            enddo
            ip=ipfindex(ndims,isuds,isind)
!            if(ip.gt.1000)write(*,*)isuds,isind,ip
            if(ip.lt.0)write(*,*)isind,ip
            call subvaccum(x_part(1,j),vlimit,ip+1)
         endif
 1       continue
      enddo
             
      end
!*****************************************************************
      subroutine subvaccum(xr,vlimit,ip)
! Accumulate a particle into velocity bin corresponding to position ip
      include 'ndimsdecl.f'
      include 'meshcom.f'
      include 'ptaccom.f'
      include 'plascom.f'
      real vlimit(2,ndims)
      real xr(3*ndims)
      integer ibin3(ndims)

      do id=1,ndims
! Assign velocities to bins.
         if(ivproj.eq.0)then
            v=xr(ndims+id)
         else
            v=vproject(ndims,id,xr,Bfield)
         endif
         if(v.lt.vlimit(1,id))v=vlimit(1,id)
         if(v.gt.vlimit(2,id))v=vlimit(2,id)
         ibin=nint(nptdiag*(.000005+0.99999*
     $        (v-vlimit(1,id))/(vlimit(2,id)-vlimit(1,id)))+0.5)
         ibin3(id)=ibin
         if(ibin.lt.1.or.ibin.gt.nptdiag) write(*,*)id,' ibin',ibin
     $        ,nptdiag,v,vlimit(1,id),vlimit(2,id)
         if(csbin(1,1).ne.-1.)then
! Doing summed bin accumulation
            ibs=ibinmap(ibin,id)
! Test for corruption ought not eventually to be necessary.
            if(ip.le.0.or.ibs.lt.0.or.id.lt.0.or.ip.gt.nsub_tot)then
               write(*,*)'Bin address corruption ibs,id,ip,nsub_tot'
     $              ,ibs,id,ip,nsub_tot
            endif
            fvx(ibs,id,ip)=fvx(ibs,id,ip)+1
         endif
      enddo
      denfvx(ip)=denfvx(ip)+1
! Now we have the bin number for all three dimensions. Increment f2vx.
      if(nptdiag.eq.nsbins)then
         do id=1,ndims
            j=mod(id,ndims)+1
            k=mod(id+1,ndims)+1
            f2vx(ibin3(j),ibin3(k),id,ip)=f2vx(ibin3(j),ibin3(k),id,ip)
     $           +1
         enddo
      endif
      end
!*******************************************************************
      subroutine fvxinit(xnewlim,cellvol,ibinit)
      include 'ndimsdecl.f'
      include 'meshcom.f'
      include 'ptaccom.f'
      real xnewlim(2,ndims)

      isfull(1)=nsub_i
      isfull(2)=nsub_j
      isfull(3)=nsub_k
      if(ibinit.eq.0)then
         do i=1,ndims
            isuds(i)=isfull(i)
         enddo
      endif
      do ip=1,nsub_tot
         denfvx(ip)=0.
         do id=1,ndims
            do ibs=1,nsbins
               fvx(ibs,id,ip)=0.
               do ib2=1,nsbins
                  f2vx(ib2,ibs,id,ip)=0.
               enddo
            enddo
         enddo
      enddo
      cellvol=1.
      do id=1,ndims
         cellvol=cellvol*(xnewlim(2,id)-xnewlim(1,id))/isuds(id)
!         write(*,*)xnewlim(2,id),xnewlim(1,id)
      enddo
!      write(*,*)'cellvol=',cellvol
      end
!**********************************************************************
      subroutine distread(xlimit,vlimit,xnewlim,name,cellvol)

      include 'ndimsdecl.f'
      include 'meshcom.f'
      include 'ptaccom.f'
      real xlimit(2,ndims),vlimit(2,ndims),xnewlim(2,ndims)
      character*(*) name

      open(25,file=name,status='old',form='unformatted',err=101)
      read(25)nptdiag,ndimsfile
      read(25)(xlimit(1,j),xlimit(2,j),vlimit(1,j),vlimit(2,j),
     $     j=1,ndims)
      read(25)((xdiag(i,j),px(i,j),i=1,nptdiag),j=1,ndims)
      read(25)((vdiag(i,j),fv(i,j),i=1,nptdiag),j=1,ndims)
      read(25)(xnewlim(1,j),xnewlim(2,j),j=1,ndims)
      read(25)nsbf,(isfull(i),i=1,ndims),cellvol
      read(25)(isuds(i),i=1,ndims)
      isftot=1
      do i=1,ndims
         isftot=isftot*isfull(i)
      enddo
! Check id the nsbins is correct and there's enough storage.
      if(nsbf.ne.nsbins)goto 103
      if(isftot.gt.nsub_tot)goto 103
      read(25)(((fvx(i,j,k),i=1,nsbins),j=1,ndims),k=1,isftot)
      read(25)(denfvx(k),k=1,isftot)
      read(25)((vsbin(i,j),csbin(i,j),fsv(i,j),i=1,nsbins),
     $     (vhbin(i,j),i=0,nsbins),j=1,ndims)
      read(25)((ibinmap(i,j),i=1,nptdiag),j=1,ndims)
      if(nsbf.eq.nptdiag)then
         read(25,err=104,end=104)((((f2vx(m,i,j,k),m=1,nsbf),i=1,nsbf),j
     $        =1,ndims),k=1,isftot)
      endif
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
 104  write(*,*)'Error reading f2vx',nsbf

      end
!**********************************************************************
      subroutine distwrite(xlimit,vlimit,xnewlim,name,cellvol)
      include 'ndimsdecl.f'
      include 'meshcom.f'
      include 'ptaccom.f'
      real xlimit(2,ndims),vlimit(2,ndims),xnewlim(2,ndims)
      character*(*) name

      open(25,file=name,status='unknown',err=101)
      close(25,status='delete')
      open(25,file=name,status='new',form='unformatted',err=101)
      write(25)nptdiag,ndims
      write(25)(xlimit(1,j),xlimit(2,j),vlimit(1,j),vlimit(2,j),
     $     j=1,ndims)
      write(25)((xdiag(i,j),px(i,j),i=1,nptdiag),j=1,ndims)
      write(25)((vdiag(i,j),fv(i,j),i=1,nptdiag),j=1,ndims)
!      write(*,'(2f10.4)')((vdiag(i,j),fv(i,j),i=1,5),j=1,ndims)
      write(25)(xnewlim(1,j),xnewlim(2,j),j=1,ndims)
      write(25)nsbins,(isfull(i),i=1,ndims),cellvol
      write(25)(isuds(i),i=1,ndims)
      isftot=1
      do i=1,ndims
         isftot=isftot*isfull(i)
      enddo
      write(25)(((fvx(i,j,k),i=1,nsbins),j=1,ndims),k=1,isftot)
      write(25)(denfvx(k),k=1,isftot)
      write(25)((vsbin(i,j),csbin(i,j),fsv(i,j),i=1,nsbins),
     $     (vhbin(i,j),i=0,nsbins),j=1,ndims)
      write(25)((ibinmap(i,j),i=1,nptdiag),j=1,ndims)
      if(nptdiag.eq.nsbins)then
         write(25)((((f2vx(m,i,j,k),m=1,nsbins),i=1,nsbins),j=1,ndims),k
     $        =1,isftot)
      endif
      close(25)

      return
 101  write(*,*)'Error opening file:',name(1:lentrim(name))
      close(25,status='delete')

      end
!******************************************************************
! 3-D only code.
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
! Distributions in ptaccom.f
      real vlimit(2,ndims),xnewlim(2,ndims)
      real cellvol,ftot,vtot,v2tot
      real fvplt(nsbins)
      integer ip,id,iv,iw
      integer ip3index
      external ip3index
      integer lentrim
      external lentrim
      character*100 string,cellstring
! 2-D variables
      character pp(4,nsbins,nsbins)
      integer nl,j
      parameter (nl=10)
      real cl(nl)
      integer ieye3d
      real f2vc(nsbins,nsbins,ndims)
      include '../accis/world3.h'
      integer icell,jcell,kcell
      real wicell,wjcell,wkcell
      integer jicell,jjcell,jkcell,ii,jj,kk,ip3,idw
      integer iup,isw,iwidth
      real xylim(4)
      integer ips,itrace,ies,level,ics
      logical loverplot
!      data ndfirst/3/ndlast/3/
      data ips/0/itrace/-1/loverplot/.false./
      data idw/0/level/1/
      data iup/1/jicell/0/jjcell/0/jkcell/0/
!------------------------------------------
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
      if(ndfirst.eq.ndlast)write(*,*)'For single frame overplot hit o'
 1    ip=ip3index(isuds,icell,jcell,kcell)+1
      wicell=(xnewlim(2,1)-xnewlim(1,1))/isuds(1)
      wjcell=(xnewlim(2,2)-xnewlim(1,2))/isuds(2)
      wkcell=(xnewlim(2,3)-xnewlim(1,3))/isuds(3)
      if(ivproj.eq.0)then
      write(cellstring,'(''Cell'',3i3,''  x=('',f6.1,'','',f6.1,'')'//
     $     ' y=('',f6.1,'','',f6.1,'') z=('',f6.1,'','',f6.1,'')'')')
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
      write(cellstring,51)
     $     icell,jcell,kcell
     $     ,xnewlim(1,1)+(icell-1)*wicell,xnewlim(1,1)+icell*wicell
     $     ,xnewlim(1,2)+(jcell-1)*wjcell,xnewlim(1,2)+jcell*wjcell
     $     ,xnewlim(1,3)+(kcell-1)*wkcell,xnewlim(1,3)+kcell*wkcell
     $     ,Bfield
      endif
!      write(*,'(a,$)')'v, T in bin'
! ----------------- Start of 1-d distribution drawing.
      call multiframe(ndlast-ndfirst+1,1,2)
      do id=ndfirst,ndlast
! Do correct scaling:
         do iv=1,nsbins
! Alternatively we combine bins in possibly multiple
! dimensions, if jicell, jjcell, or jkcell are non-zero.
            fvplt(iv)=0.
            do iw=1,nsbins
               f2vc(iv,iw,id)=0.
            enddo
            do ii=0,(isuds(1)-1)*jicell
               do jj=0,(isuds(2)-1)*jjcell
                  do kk=0,(isuds(3)-1)*jkcell
                     ip3=ip3index(isuds,1+ii+(1-jicell)*(icell-1),
     $                    1+jj+(1-jjcell)*(jcell-1),
     $                    1+kk+(1-jkcell)*(kcell-1))+1
            fvplt(iv)=fvplt(iv)+fvx(iv,id,ip3)*csbin(iv,id)
     $              *nptdiag/(vlimit(2,id)-vlimit(1,id))/cellvol
            do iw=1,nsbins
! I'm not clear about the normalization here.
               f2vc(iv,iw,id)=f2vc(iv,iw,id)+f2vx(iv,iw,id,ip3)
!*csbin(iv,id)
!     $              *nsbins/(vlimit(2,id)-vlimit(1,id))/cellvol
            enddo
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
!         write(*,'('','',2f10.5,$)')vtot/ftot,sqrt(v2tot/ftot-(vtot/ftot
!     $        )**2)
!         write(*,*)'itrace',itrace,loverplot,ndlast-ndfirst,id
         if(loverplot.and.ndlast-ndfirst.eq.0.and.itrace.ne.-1)then
            itrace=itrace+1
            call color(itrace)
            call polymark(vsbin(1,id),fvplt,nsbins,1)
         else
            call manautoinit(vsbin(1,id),fvplt,nsbins,isw,xylim(1)
     $           ,xylim(2),xylim(3),xylim(4))
            call axis()
!            write(*,*)'Called manauto',id,ndfirst,loverplot
            if(id.eq.ndfirst.and..not.loverplot)call
     $           boxtitle(cellstring(1:lentrim(cellstring)))
            itrace=0
         endif
         if(btest(idw,id-1))then
!         if(id.eq.3)then
            write(*,*)'#',icell*(1-jicell),jcell*(1-jjcell),kcell*(1
     $        -jkcell),' dimension',id
            write(*,*)nsbins
            write(*,'(2g14.6)')(vsbin(kk,id),fvplt(kk),kk=1,nsbins)
         endif

         ics=18*ndfirst
         if(loverplot)call legendline(.05,1.-.05*(itrace+1),258,
     $        cellstring(ics-2:ics+14))
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
         if(itrace.eq.0)call axlabels(' ',string(1:lentrim(string)))
         call winset(.true.)
         call polymark(vsbin(1,id),fvplt,nsbins,1)
!         call polybox(vhbin(0,id),fvx(1,id,ip),nsbins)
         if(nsbins.ne.nptdiag)then
            call polybox(vhbin(0,id),fvplt,nsbins)
         else
            call polyline(vsbin(1,id),fvplt,nsbins)
         endif
      enddo
      call winset(.false.)
      if(itrace.eq.0)call axlabels('Velocity',' ')
      call eye3d(ip)
! ------------ End of 1-d distribution drawing.
      if(nsbins.eq.nptdiag)then
! Start of 2-d distribution drawing.
         write(*,*)'Space key advances v-projection direction. f-count'
         call multiframe(0,0,0)
         id=3
 98      call pltinit(0.,1.,0.,1.)
         call jdrwstr(0.1,.7,cellstring(1:lentrim(cellstring)),1.)
         call iwrite(id,iwidth,string)
         call jdrwstr(0.02,0.02,'Direction ',1.)
         call drcstr(string)
!         call drcstr(' Advance with space.')
!       Plot the surface. With axes (1). Web color 10, axis color 7.
         j=level + 256*10 + 256*256*7
         call hidweb(vsbin(1,1),vsbin(1,2),f2vc(1,1,id),nsbins,nsbins
     $        ,nsbins,j)
         call color(15)
         call charsize(0.02,0.02)
         if(id.eq.1)then
            call ax3labels('!Bv!dy!d!@','!Bv!dz!d!@','!Bf(v)!@     ')
         elseif(id.eq.2)then
            call ax3labels('!Bv!dz!d!@','!Bv!dx!d!@','!Bf(v)!@     ')
         elseif(id.eq.3)then
            call ax3labels('!Bv!dx!d!@','!Bv!dy!d!@','!Bf(v)!@     ')
         endif
         call charsize(0.0,0.0)
         level=5
         call cont3proj(f2vc(1,1,id),pp,nsbins,nsbins,nsbins,cl,0,
     $        vsbin(1,1),vsbin(1,2),1,1.)
         call hdprset(0,0.)
! How to enable interactive plot rotation. Just this call instead of pltend:
         ies=ieye3d()
         call rotatezoom(ies)
! Space key changes f-integration direction.
         if(ies.eq.ichar(' '))id=mod(id,ndims)+1
! Exit on 0,q, or an arrow key.
         if(ies.ne.0 .and. ies.ne.ichar('q') .and. ies.lt.65300 )goto 98    
      endif
! -------------- End of 2-d. 
!      write(*,*)'ip',ip
! Increment one of the dimensions using arrow keys.
      if(ip.eq.65361)then
         icell=mod(icell-1+isuds(1)+iup,isuds(1))+1
      elseif(ip.eq.65362)then
         kcell=mod(kcell,isuds(3))+1
      elseif(ip.eq.65364)then
         kcell=mod(kcell-2+isuds(3),isuds(3))+1
      elseif(ip.eq.65363)then
         jcell=mod(jcell-1+isuds(2)+iup,isuds(2))+1
      elseif(ip.eq.65505)then
! ffe1 shift_L=2^16-31=65505
         iup=-iup
! Control averaging over dimensions
      elseif(ip.eq.ichar('x'))then
         jicell=1-jicell
      elseif(ip.eq.ichar('y'))then
         jjcell=1-jjcell
      elseif(ip.eq.ichar('z'))then
         jkcell=1-jkcell
! Control overplotting
      elseif(ip.eq.ichar('o'))then
         loverplot=.not.loverplot
         if(loverplot)itrace=-1
!         write(*,*)'loverplot=',loverplot,itrace
!         if(loverplot)call prtend(' ')
      elseif(ip.eq.ichar('p'))then
         if(ips.eq.0)then
            write(*,*)'Postscript output on'
            call pfset(3)
            ips=3
         else
            write(*,*)'Postscript output off'
            ips=0
            call prtend(' ')
            call pfset(0)
         endif
      elseif(ip.eq.ichar('a'))then
!------------------------------
! Control axis settings.
! Get another key input
         write(*,*)'Axis control.'
     $        ,' Type axis-end switch 1-4, xmin,xmax,ymin,ymax.'
 200     call eye3d(ip)
!         write(*,*)ip
         if((ip.ge.49 .and. ip.le.52))then
!            write(*,*)char(ip)
! Toggle the isw accordingly.
            isw=ieor(isw,2**(ip-49))
            write(*,'(a,i3,$)')'isw set to',isw
! Get the real number
            if(isw/2**(ip-49)-(isw/2**(ip-48))*2.ne.0)then
               write(*,*)'  enter axis end value:'
               call eyereadreal(xylim(ip-48))
            else
               write(*,*)' manual value toggled off.'
            endif
!            write(*,*)'eyeread return',ip-48,xylim(ip-48)
         else
            write(*,*)'Not a number in [1,4]:', char(ip)
            goto 200
         endif
!------------------------------
      else
!         write(*,*)ip
         call prtend(' ')
         return
      endif
      goto 1
!--------------------------------------
      end
!****************************************************************
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
!,'  decfac=',decfac
         goto 1
      else
         if(i.eq.46)then
            ldec=.not.ldec
            goto 1
         elseif(i.eq.45)then
            x=-x
            sign=-sign
!            write(*,*)x
            goto 1
         elseif(i.ne.65293)then
! Return indicates end. Anything else starts over.           
! Backspace is 65288
            write(*,*)'Starting Read Real Over:',i,char(i)
            goto 2
         endif
      endif
      end
!****************************************************************
      subroutine fsvzero()
      include 'ndimsdecl.f'
      include 'meshcom.f'
      include 'partcom.f'
      include 'ptaccom.f'
      do id=1,ndims
         do jj=1,nsbins
            fsv(jj,id)=0.
            ifsv(jj,id)=0
         enddo
      enddo
      end
!****************************************************************
      real function vproject(ndims,id,xr,Bfield)
! Return the id'th component of the "Projected" velocity xr(ndims+ii)
! relative to the direction cosines Bfield(ndims). 
! id=1 Parallel, id=2 Perpendicular, id=3 3rd (=id'th) component.
      implicit none
      integer ndims,id
      real xr(ndims),Bfield(ndims) 
      real v,v2
      integer ii
! Projected. Very clumsy.
      v=0
      do ii=1,ndims
         v=v+xr(ii+ndims)*Bfield(ii)
      enddo
      if(id.eq.2)then
! Perpendicular
         v2=0
         do ii=1,ndims
            v2=v2+(xr(ii+ndims)-v*Bfield(ii))**2
         enddo
         v=sqrt(v2)
      elseif(id.eq.3)then
! Z-component.
         v=xr(id+ndims)
      endif
      vproject=v
      end

!********************************************************************
! Convert 3-D indices into pointer. Only used in pltsubdist.
      function ip3index(ifull,i,j,k)
      integer ifull(3),i,j,k 
      ip3index=0
      ip3index=(k-1)+ip3index*ifull(3)
      ip3index=(j-1)+ip3index*ifull(2)
      ip3index=(i-1)+ip3index*ifull(1)
      end
!********************************************************************
      block data ptaccomdata
      include 'ndimsdecl.f'
      include 'ptaccom.f'
      data nptdiag/nptdiagmax/
      end
!*******************************************************************
      subroutine adjustlimit(nptdiag,vlimit)
! Make the value zero lie either at the center of a v-cell or on the
! boundary between cells, so that uniform cells are symmetrically placed
! with respect to zero.
      integer nptdiag
      real vlimit(2)
      dv=(vlimit(2)-vlimit(1))/nptdiag
      svi=2.*vlimit(1)/dv
      delta=(svi-int(svi))*dv/2.
      vlimit(1)=vlimit(1)-delta
      vlimit(2)=vlimit(2)-delta
!      write(*,*)'Adjusted vlimits to',vlimit,vlimit(1)/dv
      end
