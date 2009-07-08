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

      name=' '
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

      close(22)
      write(*,*)'Wrote particle data to ',name(1:lentrim(name))
      return

 101  continue
      write(*,*)'Error opening file:',name
      close(22,status='delete')

      end
c*****************************************************************
      subroutine partread(name,ierr)
      character*(*) name
      include 'partcom.f'
      include 'plascom.f'
      include 'ran1com.f'
      character*(100) charout

      open(23,file=name,status='old',form='unformatted',err=101)
      read(23)charout
      read(23)debyelen,Ti,vd,rs,phip
      read(23)ranstate
      read(23)ioc_part
      read(23)iregion_part,n_part,dt,ldiags,rhoinf,nrein,
     $     phirein,numprocs,
     $     ((x_part(j,i),j=1,3*npdim),if_part(i),i=1,ioc_part)
      close(23)
      write(*,*)'Read back particle data from ',name(1:lentrim(name))
      write(*,*)'Charout=',charout(1:lentrim(charout))
      ierr=0
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
 101  write(*,*)'Error opening file:',name
      ierr=1
      end
c******************************************************************
      subroutine phiwrite(name,ifull,iuds,u)
c 3-Dimensions assumed.
c File name:
      character*(*) name
      integer ifull(3),iuds(3)
      real u(ifull(1),ifull(2),ifull(3))
      include 'plascom.f'
      character*(100) charout

      name=' '
      call nameconstruct(name)
      i=nbcat(name,'.phi')
      write(charout,51)debyelen,Ti,vd,rs,phip
 51   format('debyelen,Ti,vd,rs,phip:',5f10.4)
      open(22,file=name,status='unknown',err=101)
      close(22,status='delete')
      open(22,file=name,status='new',form='unformatted',err=101)
      write(22)charout
      write(22)debyelen,Ti,vd,rs,phip
      write(22)iuds
      write(22)(((u(i,j,k),i=1,iuds(1)),j=1,iuds(2)),k=1,iuds(3))
      close(22)
      write(*,*)'Wrote potential data to ',name(1:lentrim(name))
      return

 101  continue
      write(*,*)'Error opening file:',name
      close(22,status='delete')

      end
c******************************************************************
      subroutine phiread(name,ifull,iuds,u,ierr)
c 3-Dimensions assumed.
c File name:
      character*(*) name
      integer ifull(3),iuds(3)
      real u(ifull(1),ifull(2),ifull(3))
c Don't update the plasma data based on this routine:
c      include 'plascom.f'
      character*(100) charout

      open(23,file=name,status='old',form='unformatted',err=101)
      read(23)charout
      read(23)debyelen,Ti,vd,rs,phip
      read(23)iuds
      read(23)(((u(i,j,k),i=1,iuds(1)),j=1,iuds(2)),k=1,iuds(3))
      close(23)
      write(*,*)'Read back potential data from ',name(1:lentrim(name))
      ierr=0
      return

 101  continue
      write(*,*)'Error opening file:',name
      ierr=1
      end
c******************************************************************
      subroutine nameconstruct(name)
      character*(*) name
      include 'plascom.f'
c Construct a filename that contains many parameters
c Using the routines in strings_names.f
      call nameappendexp(name,'T',Ti,1)
      call nameappendint(name,'v',nint(100*vd),3)
      if(rs.ge.100)then
         call nameappendint(name,'R',ifix(rs/10.),2)
      else
         call nameappendint(name,'r',ifix(rs),2)
      endif
      call nameappendint(name,'P',ifix(abs(phip)),2)
      call nameappendexp(name,'L',debyelen,1)

      end
