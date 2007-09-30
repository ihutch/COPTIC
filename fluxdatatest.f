c**************************************************************
      program fluxdatatest
      include '3dcom.f'
c      call fluxdatainit()
      call readfluxfile('flux_T1e0v000r05P02L1e0.dat')

      write(*,*) 'nf_addresses for mq,mo=',mf_quant,mf_obj,' :',
     $     (((nf_address(i,j,k),i=1,mf_quant),' ,'
     $     ,j=1,mf_obj),k=0,3),'...'

      write(*,*) 'Angle data'
      write(*,'(10f7.4)')(ff_data(nf_address(1,1,0)+i-1)
     $     ,i=1,nf_posno(1,1))

      write(*,*) 'Step1 data'
      write(*,'(10f7.4)')(ff_data(nf_address(1,1,1)+i-1)
     $     ,i=1,nf_posno(1,1))

      call fluxave()

      end
