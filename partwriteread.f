c**********************************************************************
      subroutine datawrite(myid,partfilename,phifilename,ifull
     $     ,iuds,u,uave,qave)
      integer myid,ifull(*),iuds(*)
      real u(*),uave(*),qave(*)
      character*(*) partfilename,phifilename
      include 'griddecl.f'
      include 'ptchcom.f'
      character*100 localfilename

      call partwrite(partfilename,myid)
      if(myid.eq.0)then
         if(iptch_copy.ne.0)then
            call namewrite(phifilename,ifull,iuds,1,uave,'.ua')
            call mditeradd(u,ndims,ifull,iuds,0,uci)
            call mditeradd(uave,ndims,ifull,iuds,0,uci)
         endif
         call namewrite(phifilename,ifull,iuds,1,uave,'.pha')
         call namewrite(phifilename,ifull,iuds,1,qave,'.den')
         call namewrite(phifilename,ifull,iuds,1,u,'.phi')
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
      include 'partcom.f'
      include 'plascom.f'
      include 'ran1com.f'
      character*(100) charout

      call nameconstruct(name)
      call nameappendint(name,'.',myid,3)
c      write(*,*)name
      write(charout,51)debyelen,Ti,vd,rs,phip
 51   format('debyelen,Ti,vd,rs,phip:',5f10.4)

      open(22,file=name,status='unknown',err=101)
      close(22,status='delete')
      open(22,file=name,status='new',form='unformatted',err=101)

      write(22)charout
      write(22)debyelen,Ti,vd,rs,phip

      write(22)ranstate
      write(22)ioc_part
      write(22)iregion_part,n_part,dt,ldiags,rhoinf,nrein,
     $     phirein,numprocs,
     $     ((x_part(j,i),j=1,3*npdim),if_part(i),i=1,ioc_part)
      write(22)(dtprec(i),i=1,ioc_part)
      write(22)rmtoz,Bt,Bfield,vpar,vperp
      close(22)
c      write(*,*)'Wrote particle data to ',name(1:lentrim(name))
      return

 101  continue
      write(*,*)'Error opening file:',name
      close(22,status='delete')

      end
c*****************************************************************
      subroutine partread(name,ierr)
c Return ierr bit(0) no file. bit(1) no dtprec. bit(2) no Bfield etc.
      character*(*) name
      integer ierr
      include 'partcom.f'
      include 'plascom.f'
      include 'ran1com.f'
      character*(100) charout

      ierr=0
      open(23,file=name,status='old',form='unformatted',err=101)
      read(23)charout
      read(23)debyelen,Ti,vd,rs,phip
      read(23)ranstate
      read(23)ioc_part
      read(23)iregion_part,n_part,dt,ldiags,rhoinf,nrein,
     $     phirein,numprocs,
     $     ((x_part(j,i),j=1,3*npdim),if_part(i),i=1,ioc_part)
c Extra particle data written since 30 July 2010.
      read(23,end=102)(dtprec(i),i=1,ioc_part)
      read(23,end=104)rmtoz,Bt,Bfield,vpar,vperp
      goto 103
 102  write(*,*)'=========== No dtprec data in partfile.========='
      ierr=2
 104  write(*,*)'=========== No Bfield etc. in partfile.========='
      ierr=ierr+4
 103  close(23)
c      write(*,*)'Finished reading back particle data from '
c     $     ,name(1:lentrim(name))
c      write(*,*)'Charout=',charout(1:lentrim(charout))
c Check that the read back data is sane
      do i=1,ioc_part
         if(if_part(i).ne.0)then 
            if(x_part(7,i).eq.0)then
               write(*,*)'Bizarre particle data read back',
     $              i,ioc_part,(x_part(j,i),j=1,9)
            endif
         endif
      enddo

c Zero the flags of higher slots.
      do i=ioc_part+1,n_partmax
         if_part(i)=0
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
      include 'plascom.f'
      include 'meshcom.f'
      character*(100) charout

c      write(*,*)'ifull',ifull
      write(charout,51)debyelen,Ti,vd,rs,phip
 51   format('V2 debyelen,Ti,vd,rs,phip:',5f9.4)
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
      write(*,*)'Error opening file:',name
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
      include 'plascom.f'
      include 'meshcom.f'
      character*(100) charout

      open(23,file=name,status='old',form='unformatted',err=101)
      read(23)charout
c      write(*,'(2a)')'Charout=',charout(1:lentrim(charout))
      read(23)debyelen,Ti,vd,rs,phip
      read(23)ixnp,xn
      read(23)iuds
      if(charout(1:2).eq.'de')then
c First version
         if(ierr.ne.0)write(*,*)'Old version file'
         ied=0
         read(23)(((u(i,j,k,1),i=1,iuds(1)),j=1,iuds(2)),k=1,iuds(3))
      elseif(charout(1:2).eq.'V2')then
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
      write(*,*)'Error opening file:',name
      ierr=1
      end
c******************************************************************
c******************************************************************
      subroutine namewrite(name,ifull,iuds,ied,u,extension)
      character*(*) name,extension
      integer ifull(3),iuds(3),ied
      real u(ifull(1),ifull(2),ifull(3))
      call nameconstruct(name)
      i=nbcat(name,extension)
      call array3write(name,ifull,iuds,ied,u)
      end
c******************************************************************
      subroutine nameconstruct(name)
      character*(*) name
      include 'plascom.f'
      include 'colncom.f'
      include 'meshcom.f'
c Construct a filename that contains many parameters
c Using the routines in strings_names.f
      name=' '
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
      end
c Below here are the obsolete versions which can be deleted once
c we are convinced there are no bugs or needs. Done 17 Aug 2010
c******************************************************************
