c**************************************************************
      program fluxdatatest
      include '../3dcom.f'
      include '../plascom.f'
      include '../sectcom.f'

      real plotdata(10000,5),stepdata(10000)
      character*100 filename,argument
      integer iplot,iprint
      data iplot/1/iprint/1/
      
      fn1=0.5
      fn2=1.
      rview=1.
      iosw=3

      filename='T1e0v000r05P02L1e0.flx'
      do i=1,iargc()
         call getarg(i,argument)
         if(argument(1:3).eq.'-n1')
     $        read(argument(4:),'(f10.4)')fn1
         if(argument(1:3).eq.'-n2')
     $        read(argument(4:),'(f10.4)')fn2
         if(argument(1:2).eq.'-f')
     $        read(argument(3:),'(a)')filename
         if(argument(1:2).eq.'-p')
     $        read(argument(3:),'(i5)')iplot
         if(argument(1:2).eq.'-w')
     $        read(argument(3:),'(i5)')iprint
         if(argument(1:2).eq.'-r')
     $        read(argument(3:),'(f10.4)')rview
         if(argument(1:2).eq.'-i')
     $        read(argument(3:),'(i5)')iosw
         if(argument(1:2).eq.'-h')goto 201
         if(argument(1:2).eq.'-?')goto 201
         if(argument(1:1).ne.'-') read(argument(1:),'(a)')filename
      enddo

      call readfluxfile(filename,ierr)

c      write(*,*)'found',ff_data(nf_address(nf_flux,1,-1)+1-1)
c     $     ,nf_address(nf_flux,1,-1)

      write(*,*)'   No. objects,   No. steps,    No. quantities(obj)'
      write(*,'(i12,i12,i12,20i3)')mf_obj,nf_step,
     $     (mf_quant(j),j=1,mf_obj)
      if(iprint.gt.0)write(*,*) 'Posn and first 2 step addresses ',
     $     (((nf_address(i,j,k),i=1,mf_quant(j)),' ,'
     $     ,j=1,mf_obj),k=1-nf_posdim,2),'...'


c For all the objects being flux tracked.
      do k=1,mf_obj
         if(iprint.gt.0)then
         write(*,'(a,i3,a,3i4,a,i3)') 'Position data for ',nf_posdim
     $        ,' flux-indices. nf_dimlens='
     $        ,(nf_dimlens(1,k,kd),kd=1,nf_posdim-1),' Object',k
         write(*,'(10f8.4)')((ff_data(nf_address(nf_flux,k,1-j)+i-1)
     $        ,i=1,nf_posno(1,k)),j=1,nf_posdim)
         
         do kk=1,nf_step,max(nf_step/5,1)
            write(*,'(a,i3,a,f10.4,a)')'Step(',kk,') rho=',ff_rho(kk)
     $           ,'  Flux data'
            write(*,'(10f8.2)')(ff_data(nf_address(nf_flux,k,kk)+i-1)
     $           ,i=1,nf_posno(nf_flux,k))
            if(mf_quant(k).ge.2)then
               write(*,'(''x-momentum'',i4)')nf_posno(nf_gx,k)
               write(*,'(10f8.3)')(ff_data(nf_address(nf_gx,k,kk)+i-1)
     $              ,i=1,nf_posno(nf_gx,k))
            endif
            if(mf_quant(k).ge.3)then
               write(*,'(''y-momentum'',i4)')nf_posno(nf_gy,k)
               write(*,'(10f8.3)')(ff_data(nf_address(nf_gy,k,kk)+i-1)
     $              ,i=1,nf_posno(nf_gy,k))
            endif
            if(mf_quant(k).ge.4)then
               write(*,'(''z-momentum'',i4)')nf_posno(nf_gz,k)
               write(*,'(10f8.3)')(ff_data(nf_address(nf_gz,k,kk)+i-1)
     $              ,i=1,nf_posno(nf_gz,k))
            endif
         enddo
         kk=nf_step
            write(*,'(a,i3,a,f10.4,a)')'Step(',kk,') rho=',ff_rho(kk)
     $           ,'  Flux data'
            write(*,'(10f8.2)')(ff_data(nf_address(nf_flux,k,kk)+i-1)
     $           ,i=1,nf_posno(nf_flux,k))
         endif
         plotdata(i,j)=pressforce(j,k,i)
      
         n1=fn1*nf_step
         n2=fn2*nf_step
         call fluxave(n1,n2,k,iplot,rhoinf)
      enddo

c Plots if 
      if(iplot.ne.0)then

      call pltinit(1.,float(nf_step),-100.,100.)
      call axis()
      call axlabels('step','Force')
      write(*,*)'   Field,       part,       press,'
     $     ,'       total,     charge:'
      do k=1,mf_obj
         call color(k)
         call iwrite(k,iwd,argument)
         avefield=0.
         avepart=0.
         avepress=0.
         avecharge=0.
         iavenum=0
         do i=1,nf_step
            plotdata(i,1)=fieldforce(3,k,i)*debyelen**2
            plotdata(i,2)=pressforce(3,k,i)
            plotdata(i,3)=partforce(3,k,i)               
            plotdata(i,4)=plotdata(i,1)+plotdata(i,2)+plotdata(i,3)
            plotdata(i,5)=charge_ns(k,i)               
c            write(*,*)(plotdata(i,j),j=1,3)
            stepdata(i)=i
            if(i.ge.n1 .and. i.le.n2)then
               iavenum=iavenum+1
               avefield=avefield+fieldforce(3,k,i)
               avepress=avepress+pressforce(3,k,i)
               avepart=avepart+partforce(3,k,i)
               avecharge=avecharge+charge_ns(k,i)
            endif
         enddo
         avefield=debyelen**2*avefield/float(iavenum)
         avepress=avepress/float(iavenum)
         avepart=avepart/float(iavenum)
         avecharge=avecharge/float(iavenum)
c         call autoplot(stepdata,plotdata(1,3),nf_step)
         call polyline(stepdata,plotdata(1,3),nf_step)
         call legendline(.1+.4*(k-1),.2,0,
     $        'partforce '//argument(1:iwd))
         call dashset(1)
         call polyline(stepdata,plotdata(1,1),nf_step)
         call legendline(.1+.4*(k-1),.15,0,
     $        'fieldforce '//argument(1:iwd))
         call dashset(2)
         call polyline(stepdata,plotdata(1,2),nf_step)
         call legendline(.1+.4*(k-1),.1,0,
     $        'pressforce '//argument(1:iwd))
         call dashset(4)
         call polyline(stepdata,plotdata(1,4),nf_step)
         call legendline(.1+.4*(k-1),.25,0,
     $        'total '//argument(1:iwd))
         call dashset(0)
         write(*,'(''================== End of Object'',i2,'' ->'',i3)')
     $     k,nf_geommap(k)
         write(*,'(6f12.5  )')
     $        avefield,avepart,avepress,
     $        avefield+avepart+avepress,avecharge
c     $        ,avecharge/avefield
      enddo

      call pltend()

      if(sc_ipt.ne.0)write(*,*)'Intersections=',sc_ipt
      call objplot(ndims,rview,iosw)

      endif

      call exit(1)
 201  write(*,*)'Usage: fluxdatatest [-ffilename,'//
     $     '-n1fff,-n2fff,-piii,-wiii,-rfff,-iiii]'

      write(*,*)'Read back flux data from file:'

      end
c******************************************************************

