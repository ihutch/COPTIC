c**************************************************************
      program fluxexamine
      include '../3dcom.f'
      include '../plascom.f'
      include '../sectcom.f'

      real plotdata(10000,5),stepdata(10000)
      character*100 filename,argument
      integer iplot,iprint,ifmask,idimf
      real avefield(ns_ndims),avepress(ns_ndims),avepart(ns_ndims)
      real rp
      data iplot/1/iprint/1/
      data ifmask/1023/
      data idimf/0/
      data rp/0./

      fn1=0.5
      fn2=1.
      rview=1.
      iosw=1

      filename='T1e0v000r05P02L1e0.flx'
      do i=1,iargc()
         call getarg(i,argument)
         if(argument(1:3).eq.'-n1')
     $        read(argument(4:),'(f10.4)')fn1
         if(argument(1:3).eq.'-n2')
     $        read(argument(4:),'(f10.4)')fn2
         if(argument(1:2).eq.'-f')
     $        read(argument(3:),'(i5)')idimf
c iplot is the quantity number to plot and average.
         if(argument(1:2).eq.'-p')
     $        read(argument(3:),'(i5)')iplot
         if(argument(1:2).eq.'-w')
     $        read(argument(3:),'(i5)')iprint
         if(argument(1:2).eq.'-m')
     $        read(argument(3:),'(i5)')ifmask
         if(argument(1:3).eq.'-rp')then
            read(argument(4:),*)rp
         elseif(argument(1:2).eq.'-r')then
            read(argument(3:),'(f10.4)')rview
         endif
         if(argument(1:2).eq.'-i')
     $        read(argument(3:),'(i5)')iosw
         if(argument(1:2).eq.'-h')goto 201
         if(argument(1:2).eq.'-?')goto 201
         if(argument(1:1).ne.'-') read(argument(1:),'(a)')filename
      enddo

      call readfluxfile(filename,ierr)

c      write(*,*)'found',ff_data(nf_address(nf_flux,1,-1)+1-1)
c     $     ,nf_address(nf_flux,1,-1)

      write(*,*)'   No. objects,   No. steps, dt,   No. quantities(obj)'
      write(*,'(i12,i12,f10.4,20i3)')mf_obj,nf_step,ff_dt(nf_step),
     $     (mf_quant(j),j=1,mf_obj)
c      write(*,*) 'Posn and first 2 step addresses ',
c     $     (((nf_address(i,j,k),i=1,mf_quant(j)),' ,'
c     $     ,j=1,mf_obj),k=1-nf_posdim,2),'...'

c      write(*,*)'geommap for objects',mf_obj,(nf_geommap(k),k=1,mf_obj)


c For all the objects being flux tracked.
      do k=1,mf_obj
         if(k.eq.iprint)then
            if(mf_quant(k).ge.1)then
               write(*,'(a,i3,a,3i4,a,i3)') 'Position data for '
     $              ,nf_posdim,' flux-indices. nf_dimlens='
     $              ,(nf_dimlens(1,k,kd),kd=1,nf_posdim-1),' Object',k
               write(*,'(10f8.4)')((ff_data(nf_address(nf_flux,k,1-j)+i
     $              -1),i=1,nf_posno(1,k)),j=1,nf_posdim)
            endif
            do kk=max(nf_step/5,1),nf_step,max(nf_step/5,1)
               if(mf_quant(k).ge.1)then
                  write(*,'(a,i4,a,f10.4,a)')'Step(',kk,') rho='
     $                 ,ff_rho(kk),'  Flux data'
                  write(*,'(10f8.2)')(ff_data(nf_address(nf_flux,k,kk)+i
     $                 -1),i=1,nf_posno(nf_flux,k))
               endif
               if(mf_quant(k).ge.2)then
                  write(*,'(''x-momentum'',i4)')nf_posno(nf_gx,k)
                  write(*,'(10f8.3)')(ff_data(nf_address(nf_gx,k,kk)+i
     $                 -1),i=1,nf_posno(nf_gx,k))
               endif
               if(mf_quant(k).ge.3)then
                  write(*,'(''y-momentum'',i4)')nf_posno(nf_gy,k)
                  write(*,'(10f8.3)')(ff_data(nf_address(nf_gy,k,kk)+i
     $                 -1),i=1,nf_posno(nf_gy,k))
               endif
               if(mf_quant(k).ge.4)then
                  write(*,'(''z-momentum'',i4)')nf_posno(nf_gz,k)
                  write(*,'(10f8.3)')(ff_data(nf_address(nf_gz,k,kk)+i
     $                 -1),i=1,nf_posno(nf_gz,k))
               endif
            enddo
         endif
         plotdata(i,j)=pressforce(j,k,i)
      
         n1=fn1*nf_step
         n2=fn2*nf_step
         if(mf_quant(k).ge.iplot)then
c            write(*,*)'Plotting',k,mf_quant(k),iplot
            call fluxave(n1,n2,k,iplot,rhoinf)
         endif
      enddo

c Plots if 
      if(iplot.ne.0)then
         call pltinit(1.,float(nf_step),-100.,100.)
         call axis()
         call axlabels('step','Force')
      endif
      if(rp.ne.0.)write(*,'(a,f10.4,a,f10.4)')
     $     'Radius',rp,' Potential',phip
      if(idimf.eq.0)write(*,*)'   Field,       part,       press,'
     $     ,'       total,   ave over steps',n1,n2
      nplot=0
      do k=1,mf_obj
         imk=ifmask/2**(k-1)
         imk=imk-2*(imk/2)

         do j=1,ns_ndims
            avefield(j)=0.
            avepart(j)=0.
            avepress(j)=0.
         enddo
         avecharge=0.
         iavenum=0
         do i=1,nf_step
            plotdata(i,1)=fieldforce(3,k,i)*debyelen**2
            plotdata(i,2)=pressforce(3,k,i)
            plotdata(i,3)=partforce(3,k,i)               
            plotdata(i,4)=plotdata(i,1)+plotdata(i,2)+plotdata(i,3)
            plotdata(i,5)=charge_ns(k,i)
            stepdata(i)=i
            if(i.ge.n1 .and. i.le.n2)then
               iavenum=iavenum+1
               do j=1,ns_ndims
                  avefield(j)=avefield(j)+fieldforce(j,k,i)
                  avepress(j)=avepress(j)+pressforce(j,k,i)
                  avepart(j)=avepart(j)+partforce(j,k,i)
               enddo
               avecharge=avecharge+charge_ns(k,i)
            endif
         enddo
         do j=1,ns_ndims
            avefield(j)=debyelen**2*avefield(j)/float(iavenum)
            avepress(j)=avepress(j)/float(iavenum)
            avepart(j)=avepart(j)/float(iavenum)
         enddo
         avecharge=avecharge/float(iavenum)
         if(iplot.ne.0)then
         if(imk.ne.0)then
            nplot=nplot+1
            call color(k)
            call iwrite(k,iwd,argument)
            call polyline(stepdata,plotdata(1,3),nf_step)
            call legendline(.1+.4*(nplot-1),.2,0,
     $        'partforce '//argument(1:iwd))
            call dashset(1)
            call polyline(stepdata,plotdata(1,1),nf_step)
            call legendline(.1+.4*(nplot-1),.15,0,
     $        'fieldforce '//argument(1:iwd))
            call dashset(2)
            call polyline(stepdata,plotdata(1,2),nf_step)
            call legendline(.1+.4*(nplot-1),.1,0,
     $           'pressforce '//argument(1:iwd))
            call dashset(4)
            call polyline(stepdata,plotdata(1,4),nf_step)
            call legendline(.1+.4*(nplot-1),.25,0,
     $           'total '//argument(1:iwd))
            call dashset(0)
         endif
         endif
         if(idimf.eq.0)then
            write(*,101)k,nf_geommap(k),obj_geom(oradius,nf_geommap(k))
     $           ,avecharge
 101        format('========== Object',i2,' ->'
     $           ,i3,' radius=',f7.3,' ========  Charge=',f10.4)
            do j=1,ns_ndims            
               write(*,'(5f12.5  )')
     $              avefield(j),avepart(j),avepress(j),
     $              avefield(j)+avepart(j)+avepress(j)
            enddo
         else
            write(*,*)obj_geom(oradius,nf_geommap(k))
     $           ,avefield(idimf)+avepart(idimf)+avepress(idimf)
         endif
      enddo

      if(iplot.ne.0)then
         call pltend()
         if(iplot.eq.1)call objplot(rview,iosw)
      endif

      call exit(1)
 201  write(*,*)'Usage: fluxexamine [filename '//
     $     '-n1fff -n2fff -piii -wiii -rfff -iiii]'
      write(*,*)'Read back flux data from file and display.'
      write(*,*)'-n1,-n2 step range over which to average.'
      write(*,*)'-p set quantity to average and plot'
      write(*,*)'-w set object whose data is to be written'
      write(*,*)'-m mask objects whose data is to be plotted'
      write(*,*)'-r set size of plot window'
      write(*,*)'-i set iosw for objplot'
      write(*,*)'-f<id> set dimension whose force to summarize'

      end
c******************************************************************

