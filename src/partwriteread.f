!**********************************************************************
      subroutine datawrite(myid,partfilename,restartpath,ifull
     $     ,iuds,u,uave,qave)
! Write data in u,uave,qave to files in path restartpath, with 
! names constructed from the parameters and suitable extensions.
      integer myid,ifull(*),iuds(*)
      real u(*),uave(*),qave(*)
      character*(*) partfilename,restartpath
      include 'ndimsdecl.f'
      include 'griddecl.f'
      include 'ptchcom.f'
      character*100 localfilename

      partfilename=restartpath
      call partwrite(partfilename)
!      write(*,*)'Returned from partwrite'
!      write(*,*)partfilename(1:lentrim(partfilename)),myid
      if(myid.eq.0)then
         if(iptch_mask.ne.0)then
            localfilename=restartpath
            call namewrite(localfilename,ifull,iuds,1,uave,'.ua')
            call mditeradd(u,ndims,ifull,iuds,0,uci)
            call mditeradd(uave,ndims,ifull,iuds,0,uci)
         endif
         localfilename=restartpath
         call namewrite(localfilename,ifull,iuds,1,uave,'.pha')
         localfilename=restartpath
         call namewrite(localfilename,ifull,iuds,1,qave,'.den')
         localfilename=restartpath
         call namewrite(localfilename,ifull,iuds,1,u,'.phi')
         localfilename=restartpath
         call writefluxfile(localfilename)
      endif

      end
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Write particle data to disk.
      subroutine partwrite(name)
! File name:
      character*(*) name
! My mpi id
      include 'ndimsdecl.f'
      include 'partcom.f'
      include 'plascom.f'
      include 'rancom.f'
      include 'myidcom.f'
      character*(100) charout

      call nameconstruct(name)
! Use an extension length sufficient for the total number of processes
! but always at least 3.
      call nameappendint(name,'.',myid,
     $     max(3,int(log10(float(nprocs))+1.)))
      open(22,file=name,status='unknown',err=101)
      close(22,status='delete')
      open(22,file=name,status='new',form='unformatted',err=101)

! New MV format
      write(charout,52)debyelen,rs,phip,dt,ldiags
 52   format('MV2. debyelen,rs,phip,dt,ldiags:',4f10.4,l5)
      write(22)charout
      write(22)debyelen,rs,phip,dt,ldiags
!      write(22)ranstate
      write(22)nspecies,rhoinf,nrein,phirein,numprocs
      write(22)Bt,Bfield,caverein,chi
      write(22)iocparta,iicparta,nparta,eoverms,vpars,vperps,vds
      write(22)((x_part(j,i),j=1,idtp),i=iicparta(1)
     $        ,iocparta(nspecies))
! write out ranlux settings.
      call rluxut(ranluxstate)
      write(22)ranluxstate
      close(22)
      write(*,*)'Wrote particle data to ',name(1:lentrim(name))
      return

 101  continue
      write(*,*)'Error opening file:',name(1:lentrim(name))
      close(22,status='delete')

      end
!*****************************************************************
      subroutine partread(name,ierr)
! Return ierr bit(0) no file. bit(1) no dtprec. bit(2) no Bfield etc.
      character*(*) name
      integer ierr
      include 'ndimsdecl.f'
      include 'partcom.f'
      include 'plascom.f'
      include 'rancom.f'
      include 'myidcom.f'
      character*(100) charout

      ierr=0
      open(23,file=name,status='old',form='unformatted',err=101)
      read(23)charout
      if(charout(1:2).eq.'MV')then
! Multispecies versions:
!         if(myid.eq.0)write(*,*)'Partread MV version detected'
         read(23)debyelen,rs,phip,dt,ldiags
         if(.not.ichar(charout(3:3)).ge.2)read(23)ranstate
         read(23)nspecies,rhoinf,nrein,phirein,numprocs
         read(23)Bt,Bfield,caverein,chi
         read(23)iocparta,iicparta,nparta,eoverms,vpars,vperps,vds
         read(23,err=102,end=102)
     $        ((x_part(j,i),j=1,idtp),i=1,iocparta(nspecies))
         if(ichar(charout(3:3)).ge.2)then
! Read back ranlux settings and initialize.
            read(23)ranluxstate
            call rluxin(ranluxstate)
         endif
      else
! Older versions.
         read(23)debyelen,Ti,vds(1),rs,phip
         read(23)ranstate
         read(23)ioc_part
         if(charout(1:2).eq.'de')then
            if(myid.eq.0)write(*,*)'Version 1 detected'
            read(23)iic_part,n_part,dt,ldiags,rhoinf,nrein,
     $           phirein,numprocs,
     $           ((x_part(j,i),j=1,3*ndims),ifp,i=1,ioc_part)
            x_part(iflag,i)=ifp
         elseif(charout(1:2).eq.'V2')then
            if(myid.eq.0)write(*,*)'Version 2 detected'
            read(23)iic_part,n_part,dt,ldiags,rhoinf,nrein,
     $           phirein,numprocs,
     $           ((x_part(j,i),j=1,iflag),i=1,ioc_part)
         endif
! Extra particle data written since 30 July 2010.
         read(23,err=102,end=102)(x_part(idtp,i),i=1,ioc_part)
         read(23,err=104,end=104)eoverm,Bt,Bfield,vpar,vperp
         read(23,err=105,end=105)caverein,chi
      endif
      if(numprocs.gt.nprocs)then
! Adjust values for reduced restart processes.
         if(myid.eq.0)write(*,*)'Partread numprocs'
     $        ,' bigger than nprocs.',numprocs,nprocs,' Adjusted.'
         rhoinf=rhoinf*nprocs/numprocs
         numprocs=nprocs
         ierr=ierr+16
      endif
      goto 103
 102  write(*,*)'=========== Partread Failed reading x_part ========='
      write(*,*)'nspecies,idtp,iocparta'
     $     ,nspecies,idtp,iocparta(nspecies),i,j
      ierr=2
      return
 104  write(*,*)'=========== Partread Failed reading Bt etc. ========='
      ierr=ierr+4
      return
 105  write(*,*)'=========== Partread Failed reading caverein ========='
      ierr=ierr+8
      return
 103  close(23)
!      write(*,*)'Finished reading back particle data from '
!     $     ,name(1:lentrim(name))
!      write(*,*)'Charout=',charout(1:lentrim(charout))
! Check that the read back data is sane
      do i=1,iocparta(nspecies)
         if(x_part(iflag,i).ne.0)then 
            if(x_part(7,i).eq.0)then
               write(*,*)'Bizarre particle data read back',
     $              i,ioc_part,(x_part(j,i),j=1,9)
            endif
         endif
      enddo

! Zero the flags of higher slots.
      do i=iocparta(nspecies)+1,n_partmax
         x_part(iflag,i)=0
      enddo

      return
 101  continue
      write(*,*)'No file: ',name(1:lentrim(name))
      ierr=1
      end
!******************************************************************
      subroutine array3write(name,ifull,iuds,ied,u)
! 3-Dimensions assumed. But extra dimension ied allowed.
! File name:
      character*(*) name
      integer ifull(3),iuds(3),ied
      real u(ifull(1),ifull(2),ifull(3),ied)
      include 'ndimsdecl.f'
      include 'plascom.f'
      include 'meshcom.f'
      character*(130) charout

!      write(*,*)'ifull',ifull
      write(charout,51)debyelen,Ti,vds(1),rs,phip,ixnlength,ifull
 51   format('V4 debyelen,Ti,vd,rs,phip',5f9.4,' ixnlength',i6,' ifull'
     $     ,3i6)
      open(22,file=name,status='unknown',err=101)
      close(22,status='delete')
      open(22,file=name,status='new',form='unformatted',err=101)
      write(22)charout
      write(22)debyelen,Ti,vd,rs,phip
      write(22)ixnp,xn
      write(22)iuds
      write(22)ied
      write(22)((((u(i,j,k,l),i=1,iuds(1)),j=1,iuds(2)),k=1,iuds(3)),l=1
     $     ,ied)
      close(22)
!      write(*,'(''Wrote array data to '',a,3i4)')
!     $     name(1:lentrim(name)),iuds
      return

 101  continue
      write(*,*)'Error opening file:',name(1:lentrim(name))
      close(22,status='delete')

      end
!******************************************************************
      subroutine array3read(name,ifull,iuds,ied,u,ierr)
! 3-Dimensions assumed. Version detection and extra dimension.
! On entry, ied is the maximum allowed extra dimension.
!           ierr if .ne.0 indicates write informational messages.
! On exit, ied is the actual number of extra dimension. That is, the
! number of 3d arrays actually read.
! File name:
      character*(*) name
      integer ifull(3),iuds(3),ied
      real u(ifull(1),ifull(2),ifull(3),ied)
      include 'ndimsdecl.f'
      include 'plascom.f'
      include 'meshcom.f'
      character*(130) charout
      integer ifulr(3),istop
      data istop/0/

      open(23,file=name,status='old',form='unformatted',err=101)
      if(ierr.ne.0)write(*,*)'Opened ',name(1:lentrim(name))
      read(23)charout
! --------- Parsing the leading string to test parameter consistency.
!      write(*,'(2a)')'Charout=',charout(1:lentrim(charout))
      irst=istrstr(charout,'ifull')
      if(irst.ne.0)then
         read(charout(irst+5:),*)ifulr
         do i=1,3
            if(ifulr(i).gt.ifull(i))then
               write(*,'(a,3i5,a,3i5)')
     $              'Allocated array dimension mismatch. File:'
     $              ,ifulr,' Code:',ifull
               istop=2
               goto 102
            endif
         enddo
      endif
      irst=istrstr(charout,'ixnlength')
      if(irst.ne.0)then
! String contains ixnlength value. Get it and check it.
         read(charout(irst+9:),*)ixnlen
         if(.not.ixnlen.le.ixnlength)then
            write(*,'(a,2i5,a)')'**** ixnlength mismatch in array3read',
     $           ixnlen,ixnlength,' possible griddecl problem'
            istop=1
         endif
      endif
 102  if(istop.gt.0)then 
         read(23)debyelen
         read(23)ixnp
         read(23)iuds
         write(*,*)charout(1:lentrim(charout))
         write(*,'(a,3i5,a)')'To read this file run $ ./setdimens'
     $              ,ifulr,'  to adjust griddecl.f'
         stop '********  array3read fatal error  *********'
      endif
!-----------
!      write(*,*)charout(irst+9:)
      read(23)debyelen,Ti,vd,rs,phip
      read(23)ixnp,(xn(kk),kk=1,ixnlen)
      read(23)iuds
      if(iuds(1).gt.ifull(1).or.iuds(2).gt.ifull(2).or.iuds(3).gt
     $     .ifull(3))then 
         write(*,*)'iuds',iuds,' bigger than ifull',ifull
         stop '********  array3read fatal error  *********'
      endif
      if(charout(1:2).eq.'de')then
! First version
         if(ierr.ne.0)write(*,*)'Old version file'
         ied=0
         read(23)(((u(i,j,k,1),i=1,iuds(1)),j=1,iuds(2)),k=1,iuds(3))
      elseif(charout(1:1).eq.'V')then
         read(23)ie
         if(ierr.ne.0)write(*,*
     $        )'New version. Number of quantities in file=',ie
         if(ie.gt.ied)then
            write(*,*)'Greater than allowed number:',ied,' not all read'
            ie=ied
         else
            ied=ie
         endif
         read(23)((((u(i,j,k,l),i=1,iuds(1)),j=1,iuds(2)),k=1,iuds(3)),l
     $        =1,ied)
      endif
      close(23)
      if(ierr.ne.0)write(*,'(''Read back array data from '',a,3i4)')
     $     name(1:lentrim(name)),iuds
      ierr=0
      return

 101  continue
      write(*,*)'Failed to open file:',name(1:lentrim(name))
      ierr=1

      end
!******************************************************************
!******************************************************************
      subroutine namewrite(name,ifull,iuds,ied,u,extension)
! Construct name (extending input string name) and write data. 
      character*(*) name,extension
      integer ifull(*),iuds(*),ied
      real u(*)
      call nameconstruct(name)
      i=nbcat(name,extension)
      call array3write(name,ifull,iuds,ied,u)
      end
!******************************************************************
      subroutine nameconstruct(name)
      character*(*) name
      include 'ndimsdecl.f'
      include 'plascom.f'
      include 'colncom.f'
      include 'meshcom.f'
! Construct a filename that contains many parameters
! Using the routines in strings_names.f
! Initializing the name must done earlier.
      if(Ti.ge.1 .and. Ti.lt.10)then
         call nameappendint(name,'T',nint(Ti),1)
      else
         call nameappendexp(name,'T',Ti,1)
      endif
      if(vd.eq.0.)then
      elseif(vd.lt.9.5)then
         call nameappendint(name,'v',nint(100*vd),3)
      else
         call nameappendint(name,'v',nint(10*vd),3)
      endif
      if(phip.eq.0.)then
      elseif(phip.lt.9.5)then
         call nameappendint(name,'P',nint(100*phip),3)
      else
         call nameappendint(name,'P',nint(10*phip),3)
      endif
      if(debyelen.ge.1 .and. debyelen.lt.10.)then
         call nameappendint(name,'L',nint(debyelen),1)
      else
         call nameappendexp(name,'L',debyelen,1)
      endif
      if(ixnp(2)-ixnp(1).gt.3)call nameappendint(name,'x'
     $     ,int(xmeshend(1)) ,int(alog10(xmeshend(1))+1))      
      if(ixnp(3)-ixnp(2).gt.3)call nameappendint(name,'y'
     $     ,int(xmeshend(2)),int(alog10(xmeshend(2))+1))
      if(ixnp(4)-ixnp(3).gt.3)call nameappendint(name,'z'
     $     ,int(xmeshend(3)),int(alog10(xmeshend(3))+1))
      if(colntime.ne.0)call nameappendexp(name,'c',colntime,1)
      if(Bt.gt..01)call nameappendexp(name,'B',Bt,1)
      end
!************************************************************************
      subroutine reportprogress(nstep,nsteps,nsubc,ndropped,ierr)
! Output parameter development to stdout.
      implicit none
      integer nstep,nsteps,nsubc,ndropped,ierr
      include 'ndimsdecl.f'
      include 'partcom.f'
      integer k
      real fluxdiag
      external fluxdiag

      if(ierr.eq.0)then
         write(*,'('' '',i5.5,$)')nstep
      else
! Step number print
         if(nstep.gt.9999.or.abs(ierr).gt.999)then
            write(*,'(i5.5,i5,$)')nstep,ierr
         else
            write(*,'(i4.4,i4,$)')nstep,ierr
         endif
         
! write out flux to object 1.
         write(*,'(f6.3,''| '',$)')fluxdiag()
         if(nspecies.gt.1)then
            write(*,*)(k,nparta(k),k=1,nspecies)
         else
            if(mod(nstep,5).eq.0)write(*,*)
         endif
         if(mod(nstep,(nsteps/25+1)*5).eq.0)then
            write(*,*)'nrein   n_part  ioc_part  rhoinf     dt'
     $           ,'   passthrus'
            write(*,'(i6,i9,i9,f11.2,f9.4,i6)')
     $           nrein,n_part,ioc_part,rhoinf,dt,npassthrough
            npassthrough=0
            if(nsubc.ne.0)write(*,'(''Subcycled:'',i5,$)')nsubc
            if(ndropped.ne.0)then
! Report dropped ions because of excessive acceleration.
               write(*,'(a,i5,a,f8.3)'
     $              )' dropped-ion period-total:',ndropped
     $              ,'  per step average:'
     $              ,float(ndropped/((nsteps/25+1)*5))
               ndropped=0
            else
               write(*,*)
            endif
         endif
      endif
      end
!************************************************************************
      subroutine periodicwrite(ifull,iuds,iLs,diagsum,uave,lmyidhead
     $     ,ndiags,ndiagmax,nstep,nsteps,idistp,vlimit
     $     ,xnewlim,cellvol,ibinit,idcount,restartpath,iavesteps)
! Periodic reduction, reporting, and writing of information on the 
! state of the simulation.
      implicit none
      integer ndiags,ndiagmax,nstep,nsteps,idistp,ibinit,idcount
      integer iavesteps
      logical lmyidhead
      real cellvol
      include 'ndimsdecl.f'
      include 'partcom.f'
      integer ifull(ndims),iuds(ndims),iLs(ndims+1)
      real diagsum(ifull(1),ifull(2),ifull(3),ndiagmax+1,nspeciesmax)
      real uave(ifull(1),ifull(2),ifull(3))
      real xlimit(2,ndimsmax),vlimit(2,ndimsmax)
      real xnewlim(2,ndimsmax)
      character*(*) restartpath

! Local variables:
      character*100 diagfilename,argument
      integer ipin,idiag,ispecies
      integer lentrim
      external lentrim

      if(ndiags.gt.0)then
! Reduce the data
         do ispecies=1,nspecies
            call diagreduce(diagsum(1,1,1,1,ispecies),ndims,ifull
     $           ,iuds,iLs,ndiags)
            call diagperiod(diagsum(1,1,1,1,ispecies),ifull,iuds
     $           ,iLs,ndiags)
! Do any other processing? 
! Divide diagsums by iavsteps, to give average diagnostic density.
            do idiag=1,ndiags
               call mditermults(diagsum(1,1,1,idiag,ispecies),ndims
     $              ,ifull,iuds,0,1./iavesteps,0.)
            enddo
! Write the ave potential into the ndiags+1 slot of diagsum (by adding
! to the previously zeroed values).
            ipin=0
            call mditeradd(diagsum(1,1,1,ndiags+1,ispecies),ndims
     $           ,ifull,iuds,ipin,uave)
            if(lmyidhead)then
! If I'm the head, write it, but use shorthand mostly.
               if(iavesteps.lt.10)then
                  write(*,'(a,$)')'D'
                  if(mod(nstep,10).eq.0)write(*,*)
               else
                  if(ispecies.eq.1)write(*,'(a,i3,a,$)')'Diags',ndiags
     $                 ,' species'
                  write(*,'(a,i1,$)')' ',ispecies
                  if(ispecies.eq.nspecies)write(*,*)
               endif
               if(ispecies.eq.1)then
                  if(nsteps.gt.9999.or.nstep.gt.9999)then
                     write(argument,'(''.dia'',i5.5)')nstep
                  else
                     write(argument,'(''.dia'',i4.4)')nstep
                  endif
               else
                  if(nsteps.gt.9999.or.nstep.gt.9999)then
                     write(argument,'(''.d'',i1.1,''a'',i5.5)')
     $                 ispecies,nstep
                  else
                     write(argument,'(''.d'',i1.1,''a'',i4.4)')
     $                 ispecies,nstep
                  endif
               endif
               diagfilename=restartpath
               call namewrite(diagfilename,ifull,iuds,ndiags+1
     $              ,diagsum(1,1,1,1,ispecies),argument)
            endif
! Now reinit diagsum
            do idiag=1,ndiags+1
               call mditerset(diagsum(1,1,1,idiag,ispecies),ndims
     $              ,ifull,iuds,0,0.)
            enddo
         enddo
      endif
      if(idistp.ne.0)then
! Particle distribution diagnostics
         if(idcount.ne.0)then
! Not for the first time, print or plot.
! Reduce the data from nodes.
            call ptdiagreduce()
            if(lmyidhead)then 
               if(2*(idistp/2)-4*(idistp/4).ne.0)
     $              call pltsubdist(5,9,9,vlimit,xnewlim,cellvol,1,3)
! I think pltsubdist ought to get its first 3 arguments from isuds.
               diagfilename=' '
               call nameconstruct(diagfilename)
               if(nsteps.gt.9999)then
                  write(diagfilename(lentrim(diagfilename)+1:)
     $                 ,'(''.pex'',i5.5)')nstep
               else
                  write(diagfilename(lentrim(diagfilename)+1:)
     $                 ,'(''.pex'',i4.4)')nstep
               endif
               if(idistp-2*(idistp/2).ne.0)
     $              call distwrite(xlimit,vlimit,xnewlim,
     $              diagfilename,cellvol)
            endif
! (Re)initialize the accumulation
            ibinit=1
            call fvxinit(xnewlim,cellvol,ibinit)
            call partacinit(vlimit)
! Unless we want fsv to accumulate all the particle counts, it also
! ought to be reinitialized here. (partexamine expects this.)
            call fsvzero()
         endif
         idcount=idcount+1
      endif
      
      end
