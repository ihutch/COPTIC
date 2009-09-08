c**************************************************************
      program fluxdatatest
      include '3dcom.f'
      character*100 filename,argument
      logical lplot

      data lplot/.true./
      
      fn1=0.5
      fn2=1.

      filename='T1e0v000r05P02L1e0.flx'
      do i=1,iargc()
         call getarg(i,argument)
         if(argument(1:3).eq.'-n1')
     $        read(argument(4:),'(f10.4)')fn1
         if(argument(1:3).eq.'-n2')
     $        read(argument(4:),'(f10.4)')fn2
         if(argument(1:2).eq.'-f')
     $        read(argument(3:),'(a)')filename
         if(argument(1:2).eq.'-p')lplot=.false.
         if(argument(1:2).eq.'-h')goto 201
         if(argument(1:2).eq.'-?')goto 201
         if(argument(1:1).ne.'-') read(argument(1:),'(a)')filename
      enddo

      call readfluxfile(filename,ierr)

      write(*,*)'found',ff_data(nf_address(nf_flux,1,-1)+1-1)
     $     ,nf_address(nf_flux,1,-1)

      write(*,*) '   No. quantities, ',
     $     'No. objects, No. steps'
      write(*,*)mf_quant,mf_obj,nf_step
      write(*,*) 'Posn and first 2 step addresses ',
     $     (((nf_address(i,j,k),i=1,mf_quant),' ,'
     $     ,j=1,mf_obj),k=1-nf_posdim,2),'...'

      do k=1,mf_obj
      write(*,'(a,i3,a,2i4,a,i3)') 'Position data for ',nf_posdim
     $     ,' dimensions. nf_dimlens='
     $     ,nf_dimlens(1,k,1),nf_dimlens(1,k,2),' Object',k
      write(*,'(10f8.4)')((ff_data(nf_address(nf_flux,k,1-j)+i-1)
     $     ,i=1,nf_posno(1,k)),j=1,nf_posdim)

      do kk=1,nf_step,max(nf_step/5,1)
         write(*,'(a,i3,a,f10.4,a)')'Step(',kk,') rho=',ff_rho(kk)
     $        ,'  Flux data'
         write(*,'(10f8.4)')(ff_data(nf_address(nf_flux,1,kk)+i-1)
     $     ,i=1,nf_posno(1,k))
      enddo
      
      
      n1=fn1*nf_step
      n2=fn2*nf_step
      call fluxave(n1,n2,k,lplot)
      write(*,'(''==================== End of Object'',i2,'' ->'',i3)')
     $     k,nf_geommap(k)
      enddo
      call exit(1)
 201  write(*,*)'Usage: fluxdatatest [-ffilename,-n1fff,-n2fff]'

      write(*,*)'Read back data from file written unformatted as:'
      write(*,*)'nf_step,mf_quant,mf_obj'
      write(*,*)'((nf_posno(i,j),i=1,mf_quant),j=1,mf_obj)'
      write(*,*)'(((nf_address(i,j,k),i=1,mf_quant),j=1,mf_obj),',
     $     'k=0,nf_step+1)'
      write(*,*)'(ff_data(i),i=1,nf_address(1,1,nf_step+1))'

      end
