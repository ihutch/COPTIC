!**************************************************************
      program fluxexamine
      include '../src/ndimsdecl.f'
      include '../src/3dcom.f'
      include '../src/plascom.f'
      include '../src/sectcom.f'
      include '../src/colncom.f'
      include '../src/vtkcom.f'

      parameter (ntr=10000)
      real plotdata(ntr,6),stepdata(ntr)
!      real traceave(ntr)
      real traceave(ntr)
      character*100 filename,argument
      integer iplot,iprint,ifmask,idimf,iomask,ivprn,ivtk,istep,irw
      real avefield(ndims),avepress(ndims),avepart(ndims)
      real avetotal(ndims),avecoln(ndims),avesq(ndims)
      real rp,yrange,cv(ndims)
      parameter (ngets=10,ngetpar=6,nfget=ngets*ngetpar)
      integer ifluxget(ngets,ngetpar)
      intrinsic ibclr,btest
      character*10 typetext(5)
      character*10 poslabel(4,5)
      integer objnumtotext(7)
      data objnumtotext/1,2,3,5,3,4,4/
      data poslabel/'cos(theta)','phi(rad)','0','Area',
     $     'x','y','z','Area','radius','theta(rad)','z','Area'
     $     ,'centroid','theta(rad)','topfrac','Area','1','2','3','Area'/
      data typetext/'Sphere ','Cuboid ','Cylinder ','SurfRev '
     $     ,'Pllelopp '/
      data iplot/1/iprint/1/ivprn/0/iquiet/1/
      data ifmask/1023/iomask/0/
      data idimf/3/irw/0/
      data rp/0./cv/ndims*0./
      data ifluxget/nfget*0/
      external getfluxdenave,getfluxden,getposcoord

      if_species=1
      fn1=0.5
      fn2=1.
      rview=1.
      iosw=1
      yrange=0.
      nbox=0
      iarg1=1
      ivtk=0
      istep=5
      vtkflag=0
      nget=0
      ipinit=0

      filename='T1e0v000r05P02L1e0.flx'
 11   continue
      do iarg=iarg1,iargc()
!         write(*,*)'iarg',iarg
         call getarg(iarg,argument)
!         write(*,*) argument
         if(argument(1:3).eq.'-n1')
     $        read(argument(4:),'(f10.4)')fn1
         if(argument(1:3).eq.'-n2')
     $        read(argument(4:),'(f10.4)')fn2
         if(argument(1:2).eq.'-f')
     $        read(argument(3:),'(i5)')idimf
         if(argument(1:2).eq.'-b')
     $        read(argument(3:),'(i5)')nbox
         if(argument(1:3).eq.'-pf')then
            call pfset(3)
         elseif(argument(1:2).eq.'-p')then
! iplot is the quantity number to plot and average.
            read(argument(3:),'(i5)')iplot
         endif
         if(argument(1:2).eq.'-q')then
! Any negative iplot value turns off initial plots.
            iplot=-1
!            ifmask=0
            iosw=-99
            iquiet=0
! Tell not to do the force plots:
            ipinit=-1
!            read(argument(3:),'(i5)')iquiet
         elseif(argument(1:2).eq.'-g')then
            nget=nget+1
            read(argument(3:),*,end=301,err=301)
     $           (ifluxget(nget,k),k=1,ngetpar)
            goto 302
 301        write(*,'(a,10i4)')'Error reading fluxget data',nget
     $           ,(ifluxget(nget,k),k=1,ngetpar)
            nget=nget-1
            stop
 302        continue
         endif
         if(argument(1:2).eq.'-w')
     $        read(argument(3:),'(i5)')iprint
         if(argument(1:2).eq.'-m')
     $        read(argument(3:),'(i5)')ifmask
         if(argument(1:4).eq.'-vtk')then
            ivtk=1
            read(argument(5:),'(i5)') istep
         elseif(argument(1:2).eq.'-v') then
            read(argument(3:),'(i5)')ivprn
         endif
         if(argument(1:3).eq.'-rp')then
            read(argument(4:),*)rp
         elseif(argument(1:3).eq.'-rw')then
            read(argument(4:),*)irw
         elseif(argument(1:2).eq.'-r')then
            read(argument(3:),*)rview
         elseif(argument(1:2).eq.'-c')then
            read(argument(3:),*,err=201,end=201)cv
         elseif(argument(1:2).eq.'-y')then
            read(argument(3:),'(f10.4)')yrange
         endif
         if(argument(1:2).eq.'-s')
     $        read(argument(3:),'(i5)')if_species
         if(argument(1:2).eq.'-i')
     $        read(argument(3:),'(i5)')iosw
         if(argument(1:3).eq.'-om')then
            read(argument(4:),'(i5)')iomask
         elseif(argument(1:2).eq.'-o')then
            read(argument(3:),'(i5)')ims
! set to mask out all objects if we are starting anew:            
!            if(iomask.eq.0)iomask=2*(2**30-1)+1
            if(iomask.eq.0)iomask=65535
            iomask=IBCLR(iomask,ims-1)
         endif

         if(argument(1:2).eq.'-h')goto 201
         if(argument(1:2).eq.'-?')goto 201
         if(argument(1:1).ne.'-')then
            read(argument(1:),'(a)')filename
            iarg1=iarg+1
            goto 10
         endif
      enddo
      iarg1=iargc()+1
 10   continue
      if(ivprn.eq.0)then
! Give informational messages.
         ierr=1
      else
         ierr=0
      endif
      call readfluxfile(filename,ierr)

! Test if species requested compatible.
      if(if_species.lt.1.or.if_species.gt.nf_species)then
         write(*,*)'Error reading file:',filename,ierr
         write(*,*)'if_species incompatible',if_species,nf_species
         stop
      else
         write(*,*)'Set iplot for species',if_species
      endif
      
!      write(*,*)'File:',filename(1:lentrim(filename)) 
      write(*,*)'debyelen,Ti,vd,rs,phip,colntime,subcycle,vneut',
     $     ',fcold,dropac,Tneut,Eneut'
      write(*,'(a,12f6.2)')':',debyelen,Ti,vd,rs,phip ,colntime,subcycle
     $     ,vneutral,fcollided,dropaccel,Tneutral,Eneutral
      if(ivprn.eq.0)then
!      write(*,*)'   No. objects,   No. steps, dt,   No. quantities(obj)'
!         write(*,'(i12,i12,f10.4,20i3)')mf_obj,nf_step,ff_dt(nf_step),
!     $     (mf_quant(j),j=1,mf_obj)
         write(*,'(a,i3,a,i5,a,f8.4,a,20i3)')'Objects=',mf_obj,' Steps='
     $        ,nf_step,' dt=',ff_dt(nf_step),' Quantities(obj)='
     $        ,(mf_quant(j),j=1,mf_obj)
!      write(*,*) 'Posn and first 2 step addresses ',
!     $     (((nf_address(i,j,k),i=1,mf_quant(j)),' ,'
!     $     ,j=1,mf_obj),k=1-nf_posdim,2),'...'
!      write(*,*)'geommap for objects',mf_obj,(nf_geommap(k),k=1,mf_obj)
      endif

! For all the objects being flux tracked.
      do k=1,mf_obj
         itype=int(obj_geom(otype,nf_geommap(k)))
         itexttype=objnumtotext(itype)
         if(k.eq.iprint.and.nf_posno(nf_flux,k).lt.200)then
!         if(k.eq.iprint)then
            if(mf_quant(k).ge.1)then
               write(*,'(a,i3,a,3i4,a,i3)') 'Position data for '
     $              ,nf_posdim,' flux-indices: x,y,z,A. nf_dimlens='
     $              ,(nf_dimlens(1,k,kd),kd=1,nf_posdim-1),' Object',k
               write(*,'(a,i2,a)')'Object type',itype
     $              ,' '//typetext(itexttype)
! Write out the position data hopefully comprehensibly.
               do j=1,nf_posdim
                  write(*,*)poslabel(j,itexttype)
                  if(abs(ff_data(nf_address(nf_flux,k,1-j))).lt.10)then
                     write(*,'(10f8.4)')(ff_data(nf_address(nf_flux,k,1
     $                    -j)+i-1),i=1,nf_posno(1,k))
                  else
                     write(*,'(10f8.1)')(ff_data(nf_address(nf_flux,k,1
     $                    -j)+i-1),i=1,nf_posno(1,k))
                  endif
               enddo
               write(*,*)'addresses',
     $              (nf_address(nf_flux,k,1-j),j=1,nf_posdim),' values'
     $              ,(ff_data(nf_address(nf_flux,k,1-j)),j=1,nf_posdim)
!               write(*,'(10f8.4)')((ff_data(nf_address(nf_flux,k,1-j)+i
!     $              -1),i=1,nf_posno(1,k)),j=1,nf_posdim)
            endif
            do kk=max(nf_step/istep,1),nf_step,max(nf_step/istep,1)
               if(mf_quant(k).ge.1)then
                  write(*,'(a,i4,a,f10.2,a)')'Step(',kk,') rho='
     $                 ,ff_rho(kk),'  Flux data'
                  write(*,'(10f8.2)')(ff_data(nf_address(nf_flux,k,kk)+i
     $                 -1),i=1,nf_posno(nf_flux,k))
               endif
               if(mf_quant(k).ge.2)then
                  write(*,'(''x-momentum'',i4)')nf_posno(nf_gx,k)
                  write(*,'(10f8.1)')(ff_data(nf_address(nf_gx,k,kk)+i
     $                 -1),i=1,nf_posno(nf_gx,k))
               endif
               if(mf_quant(k).ge.3)then
                  write(*,'(''y-momentum'',i4)')nf_posno(nf_gy,k)
                  write(*,'(10f8.1)')(ff_data(nf_address(nf_gy,k,kk)+i
     $                 -1),i=1,nf_posno(nf_gy,k))
               endif
               if(mf_quant(k).ge.4)then
                  write(*,'(''z-momentum'',i4)')nf_posno(nf_gz,k)
                  write(*,'(10f8.1)')(ff_data(nf_address(nf_gz,k,kk)+i
     $                 -1),i=1,nf_posno(nf_gz,k))
               endif
            enddo
         endif
         if(k.eq.iprint.and.mf_quant(k).ge.1.and..true.)then
            do kk=max(nf_step/istep,1),nf_step,max(nf_step/istep,1)
               if(mf_quant(k).ge.1)then
! If flag '-vtk' is true then we call vtkoutput
                  if (ivtk.eq.1)then
                     vtkflag=1
                     call vtkoutput(filename,kk,abs(iplot),iomask)
                  endif
               endif
            enddo
         endif
         n1=int(fn1*nf_step)
         n2=int(fn2*nf_step)
!         if(mf_quant(k).ge.iplot)then
!            write(*,*)'Plotting',k,mf_quant(k),iplot
         ipltp=sign(abs(iplot)+if_quant(k,if_species),iplot)
         call fluxave(n1,n2,k,ipltp,rhoinf)
!         endif
      enddo
! Write nicely the fluxdensity versus position for a sphere. 
! Object k. Face 1, indexed by i,j, [third index 0 ignored].
! We write out in row-order not column order. Which is 2nd index.
      if(irw.gt.0.)then
         if(irw.gt.mf_obj)then
            write(*,*)'Asked for non-existent object',irw,'of',mf_obj
            stop
         endif
         k=irw
         write(*,*)'Writing row-ordered fluxdensity for object',k
         imax=nf_dimlens(1,k,1)
         jmax=nf_dimlens(1,k,2)
         write(*,*)imax,jmax
         do i=1,imax
            coord=getposcoord(2,k,1,1,i,0)
            write(*,401)coord
            do j=1,jmax
               write(*,401)getfluxdenave(1,k,1,j,i,0)
            enddo
            write(*,*)
         enddo
 401     format(' ',f10.4,$)
         write(*,'(a,$)')'Coordinate:'
         do j=1,imax
            coord=getposcoord(1,k,1,j,1,0)
            write(*,401)coord
         enddo
         write(*,*)
         call exit(0)
      endif
! Section for writing out average flux to specified facets.
      write(*,*)'   dt=',ff_dt(nf_step),'   rhoinf=',rhoinf
      do i=1,nget
         if(i.eq.1)write(*,*)'get,iq,iob,fc,i1,i2,i3,fluxden, flux'
     $        ,'    fluxden/dt/rhoinf'
         fluxden=getfluxdenave(ifluxget(i,1),nf_map(ifluxget(i,2))
     $        ,ifluxget(i,3),ifluxget(i,4),ifluxget(i,5)
     $        ,ifluxget(i,6))
         flux=getfluxave(ifluxget(i,1),nf_map(ifluxget(i,2))
     $        ,ifluxget(i,3),ifluxget(i,4),ifluxget(i,5)
     $        ,ifluxget(i,6))

         write(*,'(i4,6i3,2f10.2,f10.4)')i, (ifluxget(i,k),k=1,ngetpar)
     $        ,fluxden,flux,fluxden/ff_dt(nf_step)/rhoinf
      enddo


      if(iquiet.ne.0)then
         write(*,*)'Doing npart plot',nf_step,ff_rho(1),ff_rho(nf_step)
         do i=1,nf_step
            stepdata(i)=i
            plotdata(i,1)=nf_npart(i)
!            write(*,*)i,nf_npart(i)
            plotdata(i,2)=ff_rho(i)
         enddo
         if(nf_npart(1).ne.nf_npart(nf_step))then
            call autoplot(stepdata,plotdata(1,1),nf_step)
            call axlabels('step','No. of particles')
            call scalewn(1.,stepdata(nf_step),0,2.*ff_rho(nf_step)
     $           ,.false.,.false.)
            call dashset(1)
            call color(5)
            call polyline(stepdata,plotdata(1,2),nf_step)
            call dashset(0)
            call axptset(1.,1.)
            call yaxis(first,0.)
            call axlabels(' ','rhoinf')
            call axptset(0.,0.)
            call color(15)
         else
            call autoplot(stepdata,plotdata(1,2),nf_step)
            call axlabels('step','rho-infinity')
         endif
         call pltend()
      endif
! Plots if 
      if(rp.ne.0.)write(*,'(a,f10.4,a,f10.4)')
     $     'Radius',rp,' Potential',phip
      if(ivprn.eq.0.and.ifmask.ne.0)write(*,'(2a,2i5)')
     $     '   Field,       part,       press,'
     $     ,'     collision,    total,  steps ave',n1,n2
      nplot=0
      do k=1,mf_obj
         imk=ifmask/2**(k-1)
         imk=imk-2*(imk/2)
!         write(*,*)'ifmask=',ifmask,' k=',k,' imk=',imk
         if(ifmask.eq.0)goto 50

         do j=1,ndims
            avefield(j)=0.
            avepart(j)=0.
            avepress(j)=0.
            avecoln(j)=0.
            avesq(j)=0.
         enddo
         avecharge=0.
         iavenum=0
!         write(*,'(a,$)')'Step 1 force values part/press/field:'
!         write(*,*)k,partforce(idimf,k,1),pressforce(idimf,k,1)
!     $        ,fieldforce(idimf,k,1)
         do i=1,nf_step
            plotdata(i,1)=fieldforce(idimf,k,i)*debyelen**2
            if(.not.pressforce(idimf,k,i).lt.4.e7)then
               write(*,*)'***Giant pressforce at step',i
     $              ,pressforce(idimf,k,i),'  zeroed!'
               plotdata(i,2)=0.
            else
               plotdata(i,2)=pressforce(idimf,k,i)
            endif

            if(.not.abs(partforce(idimf,k,i)).lt.4.e2)then
!               write(*,*)i-3,partforce(idimf,k,i-3)
!               write(*,*)i-2,partforce(idimf,k,i-2)
!               write(*,*)i-1,partforce(idimf,k,i-1)
               write(*,*)'***Giant partforce at step',i
     $              ,partforce(idimf,k,i),'  zeroed!'               
!               write(*,*)i+1,partforce(idimf,k,i+1)
!               stop
               plotdata(i,3)=0.
            else
               plotdata(i,3)=partforce(idimf,k,i)               
            endif
            plotdata(i,5)=charge_ns(k,i)
            plotdata(i,6)=colnforce(idimf,k,i)               
            plotdata(i,4)=plotdata(i,1)+plotdata(i,2)+plotdata(i,3)
     $           + plotdata(i,6)
            stepdata(i)=i
            if(i.ge.n1 .and. i.le.n2)then
               iavenum=iavenum+1
               do j=1,ndims
                  avefield(j)=avefield(j)+fieldforce(j,k,i)
                  avepress(j)=avepress(j)+pressforce(j,k,i)
                  avepart(j)=avepart(j)+partforce(j,k,i)
                  avecoln(j)=avecoln(j)+colnforce(j,k,i)
                  avesq(j)=avesq(j)+(fieldforce(j,k,i)+pressforce(j,k,i)
     $                 +partforce(j,k,i)+colnforce(j,k,i))**2
               enddo
               avecharge=avecharge+charge_ns(k,i)
            endif
         enddo
         do j=1,ndims
            avefield(j)=debyelen**2*avefield(j)/float(iavenum)
            avepress(j)=avepress(j)/float(iavenum)
            avepart(j)=avepart(j)/float(iavenum)
            avecoln(j)=avecoln(j)/float(iavenum)
            avesq(j)=sqrt(avesq(j)/float(iavenum))
            avetotal(j)=avefield(j)+avepart(j)+avepress(j)
         enddo
         avecharge=avecharge/float(iavenum)
         if(iplot.ne.0)then
         if(imk.ne.0)then
            nplot=nplot+1
            if(ipinit.eq.0)then
               if(yrange.eq.0.)yrange=2.*avesq(idimf)
               call pltinit(1.,float(nf_step)
     $              ,-yrange,1.5*yrange)
               call axis()
               call iwrite(idimf,iwr,argument)
               call axlabels('step','Force-'//argument)
               call winset(.true.)
               ipinit=ipinit+1
            endif
            if(ipinit.gt.0)then
!               write(*,*)'FLUXEXAMINE calling color',k
               call color(k)
               call dashset(4)
               call iwrite(k,iwd,argument)
               call boxcarave(nf_step,nbox,plotdata(1,3),traceave)
               call polyline(stepdata,traceave,nf_step)
!            call polyline(stepdata,plotdata(1,3),nf_step)
               call legendline(.1+.4*(nplot-1),.2,0,
     $              'partforce '//argument(1:iwd))
               call dashset(1)
               call boxcarave(nf_step,nbox,plotdata(1,1),traceave)
               call polyline(stepdata,traceave,nf_step)
!            call polyline(stepdata,plotdata(1,1),nf_step)
               call legendline(.1+.4*(nplot-1),.15,0,
     $              'fieldforce '//argument(1:iwd))
               call dashset(2)
               call boxcarave(nf_step,nbox,plotdata(1,2),traceave)
               call polyline(stepdata,traceave,nf_step)
               call legendline(.1+.4*(nplot-1),.1,0,
     $              'pressforce '//argument(1:iwd))
               if(avecoln(1).ne.0)then
                  call dashset(3)
                  call boxcarave(nf_step,nbox,plotdata(1,6),traceave)
                  call polyline(stepdata,traceave,nf_step)
!               call polyline(stepdata,plotdata(1,6),nf_step)
                  call legendline(.1+.4*(nplot-1),.05,0,
     $                 'collisions '//argument(1:iwd))
               endif
               call dashset(0)
               call boxcarave(nf_step,nbox,plotdata(1,4),traceave)
               call polyline(stepdata,traceave,nf_step)
               call legendline(.1+.4*(nplot-1),.25,0,
     $              'total '//argument(1:iwd))
            endif
         endif
         endif
         if(ivprn.eq.0)then
            write(*,101)k,nf_geommap(k),obj_geom(oradius,nf_geommap(k))
     $        ,obj_geom(ocenter+2,nf_geommap(k)),avecharge
 101        format('===== Object',i2,' ->'
     $        ,i3,' radius=',f7.3,' zcenter=',f7.3,' Charge='
     $        ,f10.4,' =====')
            do j=1,ndims            
               write(*,'(5f12.5  )')
     $              avefield(j),avepart(j),avepress(j),avecoln(j),
     $              avefield(j)+avepart(j)+avepress(j)+avecoln(j)
            enddo
         endif
         if(btest(ivprn,k-1))then
! v printing this object
            write(*,*)vd,avefield(idimf)+avepart(idimf)+avepress(idimf)
     $           +avecoln(idimf)
         endif

      enddo
 50   continue

!      write(*,*)'iomask=',iomask,' iosw=',iosw,' iplot=',iplot

      if(iplot.ne.0.and.vtkflag.eq.0.and.iosw.ge.-3.and.iosw.lt.3)then
         call pltend()
!         write(*,*)'abs(iplot),rview,cv,iosw,iomask',abs(iplot),rview,cv
!     $        ,iosw,iomask
         call objplot(abs(iplot),rview,cv,iosw,iomask)
      endif

! Read more arguments if there are any.
      if(iarg1.le.iargc())goto 11
      
      call exit(0)
 201  write(*,*)'Usage: fluxexamine [filename '//
     $     '-n1fff -n2fff -piii -wiii -rfff -iiii ...]'
      write(*,*)'Read back flux data from file and display.'
      write(*,*)'-n1,-n2 fractional step range over which to average.'
      write(*,*)'-p set quantity to average and plot.'
     $     ,' Default -p1. Non-positive no initial plots'
      write(*,*)'-q suppress all plots and fluxave printing.'
      write(*,*)'-w set object whose read-back data is to be written,',
     $     ' or none.'
      write(*,*)'-m mask objects whose force is to be plotted'
     $     ,' -m0 none.'
      write(*,*)'-v mask objects whose force vs v is to be printed.'
     $     ,' No other printing.'
      write(*,*)'  -v1024 masks out first 10 and stops force print.'
      write(*,*)'-r set 3-D size of plot window'
      write(*,*)'-s set species whose flux to analyse [default 1]'
      write(*,*)'-cff,ff,ff set 3-D center of plot window'
      write(*,*)'-i set iosw for objplot:'
     $     ,' Color 0 position, 1 flux, 2 flux-density. Else none.'
      write(*,*)'   If negative, no face annotation.'
      write(*,*)'-oiii add object iii to 3D objects to plot'
     $     ,' (first time masking all others).'
      write(*,*)'-omiii set full mask of 3D objects not to plot.'
     $     ,' (Use -i9 for no plot.)'
      write(*,*)'-f<id> set dimension whose force to plot'
      write(*,*)'-b<nb> set boxcar average range +-nb.'
      write(*,*)'-yfff set range of force plot'
      write(*,*)'-vtkiii outputs a vtk file every nf_step/iii steps'
      write(*,*)'-g<i,j,k,l,m,n> specify a facet to print average'
     $     ,' iquant,iobj,iface,i1,i2,i3'
      write(*,*)'-pf Write plots to files.'
!      write(*,*)'-rp... Write potential phip' 
      write(*,*)'-rwiii Row write average fluxdensity for object iii'
      write(*,*)'Turn off plots progressively: -p-1 -m0 -i9 -q'
      write(*,*)'Turn off print progressively: -w0 -v1024 -p256 ' 

      end
!********************************************************************
      subroutine fluxcheck()
      include '../src/ndimsdecl.f'
      include '../src/3dcom.f'
      include '../src/plascom.f'
      include '../src/sectcom.f'
      end
!*********************************************************************
! This subroutine calls vtkwrite to write in the common blocks our vtk data
! then it calls vtkwritescalarfacets to write vtk files for unstructred meshes
      subroutine vtkoutput(filename,kk,iplot,iomask)
      include '../src/ndimsdecl.f'
      include '../src/3dcom.f'
      include '../src/plascom.f'
      include '../src/sectcom.f'
      include '../src/colncom.f'
      include '../src/vtkcom.f'
      parameter (ntr=nvtkindmax)
      integer centering(ntr),conn(ntr),celltypes(ntr)
      character*100 filename
      integer kk,iosw,iplot,iomask
      character*4 charstep

      iosw=kk-nf_step
      call vtkwrite(iplot,iosw,iomask)
      do i=1,vtkindex
         celltypes(i)=9
         centering(i)=0
      enddo
      do j=1,4*vtkindex
         conn(j)=j-1
      enddo
      write(charstep,'(i4.4)') kk
      call vtkwritescalarfacets(4*vtkindex,vtkflx,vtkpoints,vtkindex,
     $                          celltypes,conn,centering,
     $                          0,filename(1:lentrim(filename))//
     $                           charstep//char(0),'flux'//char(0))
      end
