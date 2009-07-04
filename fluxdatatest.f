c**************************************************************
      program fluxdatatest
      include '3dcom.f'
      character*100 filename,argument

      filename='flux_T1e0v000r05P02L1e0.dat'
      do i=1,iargc()
         call getarg(i,argument)
         if(argument(1:2).eq.'-f')
     $        read(argument(3:),'(a)')filename
         if(argument(1:2).eq.'-h')goto 201
         if(argument(1:2).eq.'-?')goto 201
      enddo

      call readfluxfile(filename,ierr)

      write(*,*) '   No. quantities, ',
     $     'No. objects, No. steps'
      write(*,*)mf_quant,mf_obj,nf_step
      write(*,*) 'First 2 step addresses ',
     $     (((nf_address(i,j,k),i=1,mf_quant),' ,'
     $     ,j=1,mf_obj),k=0,2),'...'

      write(*,'(a,i3,a)') 'Position data for ',nf_posdim,' dimensions.'
      write(*,'(10f8.4)')((ff_data(nf_address(1,1,1-j)+i-1)
     $     ,i=1,nf_posno(1,1)),j=1,nf_posdim)

      write(*,*)nf_step
      do k=1,nf_step,max(nf_step/5,1)
         write(*,'(a,i3,a)') 'Step(',k,') data'
         write(*,'(10f8.4)')(ff_data(nf_address(1,1,k)+i-1)
     $     ,i=1,nf_posno(1,1))
      enddo

      call fluxave()
      call exit(1)
 201  write(*,*)'Usage: fluxdatatest [-ffilename]'

      write(*,*)'Read back data from file written unformatted as:'
      write(*,*)'nf_step,mf_quant,mf_obj'
      write(*,*)'((nf_posno(i,j),i=1,mf_quant),j=1,mf_obj)'
      write(*,*)'(((nf_address(i,j,k),i=1,mf_quant),j=1,mf_obj),',
     $     'k=0,nf_step+1)'
      write(*,*)'(ff_data(i),i=1,nf_address(1,1,nf_step+1))'

      end
