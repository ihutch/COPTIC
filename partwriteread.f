c**********************************************************************
      subroutine datawrite(myid,partfilename,restartpath,ifull
     $     ,iuds,u,uave,qave)
c Write data in u,uave,qave to files in path restartpath, with 
c names constructed from the parameters and suitable extensions.
      integer myid,ifull(*),iuds(*)
      real u(*),uave(*),qave(*)
      character*(*) partfilename,restartpath
      include 'ndimsdecl.f'
      include 'griddecl.f'
      include 'ptchcom.f'
      character*100 localfilename

      partfilename=restartpath
      call partwrite(partfilename,myid)
c      write(*,*)'Returned from partwrite'
c      write(*,*)partfilename(1:lentrim(partfilename)),myid
      if(myid.eq.0)then
         if(iptch_copy.ne.0)then
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
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c Write particle data to disk.
      subroutine partwrite(name,myid)
c File name:
      character*(*) name
c My mpi id
      integer myid
      include 'ndimsdecl.f'
      include 'partcom.f'
      include 'plascom.f'
      include 'rancom.f'
      character*(100) charout

      call nameconstruct(name)
      call nameappendint(name,'.',myid,3)
      open(22,file=name,status='unknown',err=101)
      close(22,status='delete')
      open(22,file=name,status='new',form='unformatted',err=101)

c New MV format
      write(charout,52)debyelen,rs,phip,dt,ldiags
 52   format('MV2. debyelen,rs,phip,dt,ldiags:',4f10.4,l5)
      write(22)charout
      write(22)debyelen,rs,phip,dt,ldiags
c      write(22)ranstate
      write(22)nspecies,rhoinf,nrein,phirein,numprocs
      write(22)Bt,Bfield,caverein,chi
      write(22)iocparta,iicparta,nparta,eoverms,vpars,vperps,vds
      write(22)((x_part(j,i),j=1,idtp),i=iicparta(1)
     $        ,iocparta(nspecies))
c write out ranlux settings.
      call rluxut(ranluxstate)
      write(22)ranluxstate
      close(22)
      write(*,*)'Wrote particle data to ',name(1:lentrim(name))
      return

 101  continue
      write(*,*)'Error opening file:',name(1:lentrim(name))
      close(22,status='delete')

      end
c*****************************************************************
      subroutine partread(name,ierr)
c Return ierr bit(0) no file. bit(1) no dtprec. bit(2) no Bfield etc.
      character*(*) name
      integer ierr
      include 'ndimsdecl.f'
      include 'partcom.f'
      include 'plascom.f'
      include 'rancom.f'
      character*(100) charout

      ierr=0
      open(23,file=name,status='old',form='unformatted',err=101)
      read(23)charout
      if(charout(1:2).eq.'MV')then
c Multispecies versions:
         write(*,*)'Partread MV version detected'
         read(23)debyelen,rs,phip,dt,ldiags
         if(.not.ichar(charout(3:3)).ge.2)read(23)ranstate
         read(23)nspecies,rhoinf,nrein,phirein,numprocs
         read(23)Bt,Bfield,caverein,chi
         read(23)iocparta,iicparta,nparta,eoverms,vpars,vperps,vds
         read(23,err=102,end=102)
     $        ((x_part(j,i),j=1,idtp),i=1,iocparta(nspecies))
         if(ichar(charout(3:3)).ge.2)then
c Read back ranlux settings and initialize.
            read(23)ranluxstate
            call rluxin(ranluxstate)
         endif
      else
c Older versions.
         read(23)debyelen,Ti,vd,rs,phip
         read(23)ranstate
         read(23)ioc_part
         if(charout(1:2).eq.'de')then
            write(*,*)'Version 1 detected'
            read(23)iic_part,n_part,dt,ldiags,rhoinf,nrein,
     $           phirein,numprocs,
     $           ((x_part(j,i),j=1,3*ndims),ifp,i=1,ioc_part)
            x_part(iflag,i)=ifp
         elseif(charout(1:2).eq.'V2')then
            write(*,*)'Version 2 detected'
            read(23)iic_part,n_part,dt,ldiags,rhoinf,nrein,
     $           phirein,numprocs,
     $           ((x_part(j,i),j=1,iflag),i=1,ioc_part)
         endif
c Extra particle data written since 30 July 2010.
         read(23,err=102,end=102)(x_part(idtp,i),i=1,ioc_part)
         read(23,err=104,end=104)eoverm,Bt,Bfield,vpar,vperp
         read(23,err=105,end=105)caverein,chi
      endif
      goto 103
 102  write(*,*)'=========== No dtprec data in partfile.========='
      ierr=2
 104  write(*,*)'=========== No Bfield etc. in partfile.========='
      ierr=ierr+4
 105  write(*,*)'=========== No caverein  . in partfile.========='
      ierr=ierr+8
 103  close(23)
c      write(*,*)'Finished reading back particle data from '
c     $     ,name(1:lentrim(name))
c      write(*,*)'Charout=',charout(1:lentrim(charout))
c Check that the read back data is sane
      do i=1,iocparta(nspecies)
         if(x_part(iflag,i).ne.0)then 
            if(x_part(7,i).eq.0)then
               write(*,*)'Bizarre particle data read back',
     $              i,ioc_part,(x_part(j,i),j=1,9)
            endif
         endif
      enddo

c Zero the flags of higher slots.
      do i=iocparta(nspecies)+1,n_partmax
         x_part(iflag,i)=0
      enddo

      return
 101  continue
      write(*,*)'No file: ',name(1:lentrim(name))
      ierr=1
      end
c******************************************************************
      subroutine array3write(name,ifull,iuds,ied,u)
c 3-Dimensions assumed. But extra dimension ied allowed.
c File name:
      character*(*) name
      integer ifull(3),iuds(3),ied
      real u(ifull(1),ifull(2),ifull(3),ied)
      include 'ndimsdecl.f'
      include 'plascom.f'
      include 'meshcom.f'
      character*(130) charout

c      write(*,*)'ifull',ifull
      write(charout,51)debyelen,Ti,vd,rs,phip,ixnlength
 51   format('V3 debyelen,Ti,vd,rs,phip',5f9.4,' ixnlength',i6)
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
c      write(*,'(''Wrote array data to '',a,3i4)')
c     $     name(1:lentrim(name)),iuds
      return

 101  continue
      write(*,*)'Error opening file:',name(1:lentrim(name))
      close(22,status='delete')

      end
c******************************************************************
      subroutine array3read(name,ifull,iuds,ied,u,ierr)
c 3-Dimensions assumed. Version detection and extra dimension.
c On entry, ied is the maximum allowed extra dimension.
c           ierr if .ne.0 indicates write informational messages.
c On exit, ied is the actual number of extra dimension. That is, the
c number of 3d arrays actually read.
c File name:
      character*(*) name
      integer ifull(3),iuds(3),ied
      real u(ifull(1),ifull(2),ifull(3),ied)
      include 'ndimsdecl.f'
      include 'plascom.f'
      include 'meshcom.f'
      character*(130) charout

      open(23,file=name,status='old',form='unformatted',err=101)
      write(*,*)'Opened',name
      read(23)charout
c      write(*,'(2a)')'Charout=',charout(1:lentrim(charout))
      irst=istrstr(charout,'ixnlength')
      if(irst.ne.0)then
c String contains ixnlength value. Get it and check it.
         read(charout(irst+9:),*)ixnlen
         if(ixnlen.ne.ixnlength)then
            write(*,*)'ixnlength mismatch',
     $           ' in array3read. Written with different griddecl.',
     $           ixnlen,ixnlength
            stop 'array3read fatal'
         endif
      endif
c      write(*,*)charout(irst+9:)
      read(23)debyelen,Ti,vd,rs,phip
      read(23)ixnp,xn
      read(23)iuds
      if(charout(1:2).eq.'de')then
c First version
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
c******************************************************************
c******************************************************************
      subroutine namewrite(name,ifull,iuds,ied,u,extension)
c Construct name (extending input string name) and write data. 
      character*(*) name,extension
      integer ifull(*),iuds(*),ied
      real u(*)
      call nameconstruct(name)
      i=nbcat(name,extension)
      call array3write(name,ifull,iuds,ied,u)
      end
c******************************************************************
      subroutine nameconstruct(name)
      character*(*) name
      include 'ndimsdecl.f'
      include 'plascom.f'
      include 'colncom.f'
      include 'meshcom.f'
c Construct a filename that contains many parameters
c Using the routines in strings_names.f
c We do not now start by zeroing the string. Has to be done earlier.
c      name=' '
      call nameappendexp(name,'T',Ti,1)
      if(vd.lt.9.5)then
         call nameappendint(name,'v',nint(100*vd),3)
      else
         call nameappendint(name,'v',nint(10*vd),3)
      endif
      if(abs(phip).lt.9.5)then
         call nameappendint(name,'P',nint(abs(phip*100)),3)
      else
         call nameappendint(name,'P',nint(abs(phip*10)),3)
      endif
      call nameappendexp(name,'L',debyelen,1)
      call nameappendint(name,'z',nint(xmeshend(3)),3)
      call nameappendint(name,'x',nint(xmeshend(2)),2)
      if(colntime.ne.0)call nameappendexp(name,'c',colntime,1)
      if(Bt.gt..01)call nameappendexp(name,'B',Bt,1)
      end
c*****************************************************************
      subroutine phipset(myid)
      include 'ndimsdecl.f'
      include 'plascom.f'
      include '3dcom.f'

      if(obj_geom(oabc,1).ne.0)then
         phip=-obj_geom(oabc+2,1)/obj_geom(oabc,1)
         if(myid.eq.0)write(*,*)'Object 1 potential=',phip
      elseif(obj_geom(oradius,1).ne.0.)then
         phip=obj_geom(omag,1)*obj_geom(oradius,1)
         if(myid.eq.0)write(*,*)'Potential from point charge'
     $        ,obj_geom(omag,1),' at radius ',obj_geom(oradius,1)
     $        ,' Charge:',phip
      else
         if(myid.eq.0)write(*,*)'Potential phip not set from objects.'
      endif

      end
c******************************************************************
c************************************************************************
      subroutine reportprogress(nf_step,nsteps,nsubc,ndropped)
c Output parameter development to stdout.
      implicit none
      integer nf_step,nsteps,nsubc,ndropped
      include 'ndimsdecl.f'
      include 'partcom.f'
      integer k
      real fluxdiag
      external fluxdiag

c write out flux to object 1.
            write(*,'(f6.3,''| '',$)')fluxdiag()
            if(nspecies.gt.1)then
               write(*,*)(k,nparta(k),k=1,nspecies)
            else
               if(mod(nf_step,5).eq.0)write(*,*)
            endif
            if(mod(nf_step,(nsteps/25+1)*5).eq.0)then
               write(*,*)'nrein   n_part  ioc_part  rhoinf    dt'
     $              ,'   passthrus'
               write(*,'(i6,i9,i9,f9.2,f9.4,i6)')
     $              nrein,n_part,ioc_part,rhoinf,dt,npassthrough
               npassthrough=0
               if(nsubc.ne.0)write(*,'(''Subcycled:'',i5,$)')nsubc
               if(ndropped.ne.0)then
c Report dropped ions because of excessive acceleration.
                  write(*,'(a,i5,a,f8.3)'
     $             )' dropped-ion period-total:',ndropped
     $                 ,'  per step average:'
     $                 ,float(ndropped/((nsteps/25+1)*5))
                  ndropped=0
               else
                  write(*,*)
               endif
            endif
            end
c************************************************************************
      subroutine periodicwrite(ifull,iuds,iLs,diagsum,uave,lmyidhead
     $     ,ndiags,ndiagmax,nf_step,nsteps,idistp,vlimit
     $     ,xnewlim,cellvol,ibset,idcount)
c Periodic reduction, reporting, and writing of information on the 
c state of the simulation.
      implicit none
      integer ndiags,ndiagmax,nf_step,nsteps,idistp,ibset,idcount
      logical lmyidhead
      real cellvol
      include 'ndimsdecl.f'
      include 'partcom.f'
      integer ifull(ndims),iuds(ndims),iLs(ndims+1)
      real diagsum(ifull(1),ifull(2),ifull(3),ndiagmax+1,nspeciesmax)
      real uave(ifull(1),ifull(2),ifull(3))
      real xlimit(2,ndimsmax),vlimit(2,ndimsmax)
      real xnewlim(2,ndimsmax)

c Local variables:
      character*100 diagfilename,argument
      integer ipin,idiag,ispecies
      integer lentrim
      external lentrim

      if(ndiags.gt.0)then
c Reduce the data
         do ispecies=1,nspecies
            call diagreduce(diagsum(1,1,1,1,ispecies),ndims,ifull
     $           ,iuds,iLs,ndiags)
            call diagperiod(diagsum(1,1,1,1,ispecies),ifull,iuds
     $           ,iLs,ndiags)
c Do any other processing? Here or later?
c Write the ave potential into the ndiags+1 slot of diagsum (by adding
c to the previously zeroed values).
            ipin=0
            call mditeradd(diagsum(1,1,1,ndiags+1,ispecies),ndims
     $           ,ifull,iuds,ipin,uave)
            if(lmyidhead)then
c If I'm the head, write it.
               if(ispecies.eq.1)write(*,'(a,i3,a,$)')'Diags',ndiags
     $              ,' species'
               write(*,'(a,i1,$)')' ',ispecies
               if(ispecies.eq.nspecies)write(*,*)
               if(ispecies.eq.1)then
                  if(nsteps.gt.9999)then
                     write(argument,'(''.dia'',i5.5)')nf_step
                  else
                     write(argument,'(''.dia'',i4.4)')nf_step
                  endif
               else
                  if(nsteps.gt.9999)then
                     write(argument,'(''.d'',i1.1,''a'',i5.5)')
     $                 ispecies,nf_step
                  else
                     write(argument,'(''.d'',i1.1,''a'',i4.4)')
     $                 ispecies,nf_step
                  endif
               endif
               diagfilename=' '
               call namewrite(diagfilename,ifull,iuds,ndiags+1
     $              ,diagsum(1,1,1,1,ispecies),argument)
            endif
c Now reinit diagsum
            do idiag=1,ndiags+1
               call mditerset(diagsum(1,1,1,idiag,ispecies),ndims
     $              ,ifull,iuds,0,0.)
            enddo
         enddo
      endif
      if(idistp.ne.0)then
c Particle distribution diagnostics
         if(idcount.ne.0)then
c Not for the first time, print or plot.
c Reduce the data from nodes.
            call ptdiagreduce()
            if(lmyidhead)then 
               if(2*(idistp/2)-4*(idistp/4).ne.0)
     $              call pltsubdist(5,9,9,vlimit,xnewlim,cellvol,1,3)
               diagfilename=' '
               call nameconstruct(diagfilename)
               if(nsteps.gt.9999)then
                  write(diagfilename(lentrim(diagfilename)+1:)
     $                 ,'(''.pex'',i5.5)')nf_step
               else
                  write(diagfilename(lentrim(diagfilename)+1:)
     $                 ,'(''.pex'',i4.4)')nf_step
               endif
               if(idistp-2*(idistp/2).ne.0)
     $              call distwrite(xlimit,vlimit,xnewlim,
     $              diagfilename,cellvol)
            endif
c (Re)initialize the accumulation
            ibset=1
            call fvxinit(xnewlim,cellvol,ibset)
            call partacinit(vlimit)
c Unless we want fsv to accumulate all the particle counts, it also
c ought to be reinitialized here. (partexamine expects this.)
            call fsvzero()
         endif
         idcount=idcount+1
      endif
      
      end
c******************************************************************
