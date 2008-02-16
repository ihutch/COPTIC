c**************************************************************
      program fluxdatatest
      include '3dcom.f'
c      call fluxdatainit()
      call readfluxfile('flux_T1e0v000r05P02L1e0.dat')

      write(*,*) 'Got data for mq,mo,nsteps=',mf_quant,mf_obj,nf_step
      write(*,*) 'First 2 step addresses ',
     $     (((nf_address(i,j,k),i=1,mf_quant),' ,'
     $     ,j=1,mf_obj),k=0,2),'...'

      write(*,*) 'Angle data'
      write(*,'(10f7.4)')(ff_data(nf_address(1,1,0)+i-1)
     $     ,i=1,nf_posno(1,1))

      do k=1,nf_step,nf_step/5
      write(*,'(a,i3,a)') 'Step(',k,') data'
      write(*,'(10f7.4)')(ff_data(nf_address(1,1,1)+i-1)
     $     ,i=1,nf_posno(1,1))
      enddo

      call fluxave()

      end
