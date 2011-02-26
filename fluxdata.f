c Initialize the flux data, determining what we are saving and where.
c The objects whose flux is to be tracked are indicated by 
c obj_geom(ofluxtype,i). If this is zero, it is not tracked.
c The number of ibins in each (2) of the surface dimensions is indicated
c by obj_geom(ofn1[/2],i), and data space and addresses allocated.
c The uniform-spacing bin-positions are calculated and set.
      subroutine fluxdatainit(myid)
      include '3dcom.f'
      include 'sectcom.f'
c-----------------------------------------------
c Intersection diagnostic points.
      sc_ipt=0
c-----------------------------------------------
c Initialize here to avoid giant block data program.
      nf_step=0
      do i=1,nf_quant
         do j=1,nf_obj
            nf_posno(i,j)=0
            do k=1-nf_posdim,nf_maxsteps
               nf_address(i,j,k)=0
            enddo
         enddo
      enddo
      do i=1,nf_datasize
         ff_data(i)=0.
      enddo
c-------------------------------------------------
c This initialization ought not to be necessary except for partforce.
      do k=1,nf_maxsteps
         do j=1,nf_obj
            do i=1,ns_ndims
               fieldforce(i,j,k)=0.
               pressforce(i,j,k)=0.
               partforce(i,j,k)=0.
            enddo
            charge_ns(j,k)=0.
         enddo
      enddo
c-------------------------------------------------
c Initialize object number
      mf_obj=0
      do i=1,ngeomobj
         if(obj_geom(ofluxtype,i).eq.0)then
c No flux setting for this object.
         elseif(obj_geom(ofluxtype,i).ge.1
     $           .and. obj_geom(ofluxtype,i).le.nf_quant)then
c Might eventually need more interpretation of fluxtype.
            if(mf_obj.lt.nf_obj)then
               mf_obj=mf_obj+1
            else
               write(*,*)'WARNING: Flux object number exceeds nf_obj='
     $              ,nf_obj,' Omitting later objects.'
               goto 1
            endif
            mf_quant(mf_obj)= obj_geom(ofluxtype,i)
c The mapped object number != object number.
            nf_map(i)=mf_obj
            nf_geommap(mf_obj)=i
            itype=int(obj_geom(otype,i))
c Use only bottom 8 bits:
            itype=itype-256*(itype/256)
c There are nfluxes positions for each quantity.
c At present there's no way to prescribe different grids for each
c quantity in the input file. But the data structures could 
c accommodate such difference if necessary. 
            do j=1,mf_quant(mf_obj)
               if(itype.eq.1)then
c Sphere one face only n3=1
                  nfluxes=1
                  do k=1,2
                     nf_dimlens(j,mf_obj,k)=int(obj_geom(ofn1+k-1,i))
                     nfluxes=nfluxes*obj_geom(ofn1+k-1,i)
                  enddo
               elseif(itype.eq.2 .or. itype.eq.4)then
c Cuboid or parallelopiped.
c Six faces, each using two of the three array lengths. 
                  nfluxes=0
                  do kk=1,int(obj_geom(offc,i))
                     k=mod(kk-1,ns_ndims)+1
                     if(kk.le.ns_ndims)
     $                 nf_dimlens(j,mf_obj,k)=int(obj_geom(ofn1+k-1,i))
                     nf_faceind(j,mf_obj,kk)=nfluxes
                     nfluxes=nfluxes+
     $                    obj_geom(ofn1+mod(k,ns_ndims),i)
     $                    *obj_geom(ofn1+mod(k+1,ns_ndims),i)
c                     write(*,*)k,' nfluxes=',nfluxes
c     $                    ,obj_geom(ofn1+mod(k,ns_ndims),i)
                  enddo
c                  write(*,*)(nf_faceind(j,mf_obj,k),k=1,2*ns_ndims)
c                  nfluxes=2*nfluxes
               elseif(itype.eq.3 .or. itype.eq.5)then
c Cylinder specifying nr, nt, nz. 
c Three facets in the order bottom side top.
                  nfluxes=0
                  nf_faceind(j,mf_obj,1)=nfluxes
                  nf_dimlens(j,mf_obj,1)=int(obj_geom(ofn1,i))
                  nfluxes=nfluxes+
     $                 int(obj_geom(ofn1,i))*int(obj_geom(ofn2,i))
                  nf_faceind(j,mf_obj,2)=nfluxes
                  nf_dimlens(j,mf_obj,2)=int(obj_geom(ofn2,i))
                  nfluxes=nfluxes+
     $                 int(obj_geom(ofn2,i))*int(obj_geom(ofn3,i))
                  nf_faceind(j,mf_obj,3)=nfluxes
                  nf_dimlens(j,mf_obj,3)=int(obj_geom(ofn3,i))
                  nfluxes=nfluxes+
     $                 int(obj_geom(ofn1,i))*int(obj_geom(ofn2,i))
c                  write(*,*)nfluxes,(nf_faceind(j,mf_obj,k),k=1,3)
               else
                  write(*,*)'Unknown object in fluxdatainit. Quit.'
                  stop
               endif
               nf_posno(j,mf_obj)=nfluxes
               
            enddo
            if(myid.eq.0)
     $           write(*,'(a,2i2,a,i2,a,i5,a,i3,a,i3,a,i3,a,2i3)')
     $           ' Fluxinit of object',i,mf_obj
     $           ,' ftype',int(obj_geom(ofluxtype,i))
     $           ,' Total',nf_posno(1,mf_obj),' dimlens:'
     $           ,nf_dimlens(1,mf_obj,1),'x',nf_dimlens(1,mf_obj,2)
     $           ,'x',nf_dimlens(1,mf_obj,3)
     $           ,' Quantities:',mf_quant(mf_obj)
         else
            write(*,*)'==== Unknown flux type',obj_geom(ofluxtype,i)
            stop
         endif
      enddo
 1    continue
c      write(*,*)'nf_posno=',((nf_posno(j,i),j=1,2),i=1,4)
c-------------------------------------------------
c Now we create the addressing arrays etc.
      call nfaddressinit()
c After which, nf_address(i,j,k) points to the start of data for 
c quantity i, object j, step k. 
c So we can pass nf_data(nf_address(i,j,k)) as a vector start.
c-------------------------------------------------
c Initialize the position and area data for each object.
      call positioninit()

      end
c******************************************************************
      subroutine nfaddressinit()
      include '3dcom.f'
c General iteration given correct settings of nf_posno. Don't change!
c Zero nums to silence incorrect warnings.
      numdata=0
      numobj=0
      nf_address(1,1,1-nf_posdim)=1
      do k=1-nf_posdim,nf_maxsteps+2
         if(k.gt.1-nf_posdim)
     $        nf_address(1,1,k)=nf_address(1,1,k-1)+numobj
         numobj=0
         do j=1,mf_obj
            if(j.gt.1)nf_address(1,j,k)=nf_address(1,j-1,k)+numdata
            numdata=0
            do i=1,mf_quant(j)
               if(i.gt.1)nf_address(i,j,k)=
     $              nf_address(i-1,j,k)+nf_posno(i-1,j)
               numdata=numdata+nf_posno(i,j)
c               write(*,*)i,j,k,nf_posno(i,j),nf_address(i,j,k),numdata
            enddo
            numobj=numobj+numdata
         enddo
      enddo
c Check if we might overrun the datasize.
      if(nf_address(nf_quant,mf_obj,nf_maxsteps+2)+numobj
     $     .gt.nf_datasize)then
         write(*,*)'DANGER: data from',nf_quant,mf_obj,nf_maxsteps,
     $        ' would exceed nf_datasize',nf_datasize
         stop
      else
      endif
      end
c******************************************************************
      subroutine positioninit()
      include '3dcom.f'

c      real xyz(ns_ndims)
c The k=1-nf_posdim to k=0 slots exist for us to put descriptive
c information: values that provide positions to correspond to the
c fluxes.  Area is most important and is in 1-nf_posdim.
c
c Ellipses require elliptic surfaces which are a mess. For equal radii,
c the area element is just dA=2\pi a^2 dcostheta. 
      do i=1,ngeomobj
         if(obj_geom(ofluxtype,i).gt.0)then
            itype=obj_geom(otype,i)
            itype=itype-256*(itype/256)
            io=nf_map(i)
            if(itype.eq.1)then
c Sphere -----------------------
               if(obj_geom(oradius,i).ne.obj_geom(oradius+1,i) .or.
     $              obj_geom(oradius,i).ne.obj_geom(oradius+2,i))then
                  write(*,*)'Warning! Non-spherical spheroid'
     $                 ,' surface areas'
     $                 ,' not calculated correctly.'
               endif
               do j=1,mf_quant(io)
c Area of each element. They are equal only if radii are equal. 
                  ar=4.*3.1415926*obj_geom(oradius,i)**2
     $                 /(nf_dimlens(j,io,1)*nf_dimlens(j,io,2))
c                  write(*,*)'Object',i,' Quant',j,' Facet areas=',ar
                  do i2=1,nf_dimlens(j,io,2)
                     p=3.1415926*(-1.+2.*(i2-0.5)/nf_dimlens(j,io,2))
                     do i1=1,nf_dimlens(j,io,1)
                        c=-1.+2.*(i1-0.5)/nf_dimlens(j,io,1)
                        ip=i1+(i2-1)*int(nf_dimlens(j,io,1))
c Position values are cos(theta) and psi. Third not used.
c                        write(*,*)j,io,i1,i2,nf_p1,ioff,ip,c
c     $                       ,nf_address(j,io,nf_p1)+ip-1
                        ff_data(nf_address(j,io,nf_p1)+ip-1)=c
                        ff_data(nf_address(j,io,nf_p2)+ip-1)=p
                        ff_data(nf_address(j,io,nf_pa)+ip-1)=ar
                     enddo
                  enddo
               enddo
c            write(*,*)'Set ff_data',i,ioff,nf_map(i),io
c End of sphere case.
            elseif(itype.eq.2 .or. itype.eq.4)then
c Parallelopiped/Cube ----------------------------
c facet areas are face areas divided by no of facets.
               do j=1,mf_quant(io)
c Over different faces:
                  do k=1,2*ns_ndims
                     k1=mod(k-1,ns_ndims)+1
                     k2=mod(k  ,ns_ndims)+1
                     k3=mod(k+1,ns_ndims)+1
c                     write(*,'(/,a,i2,3i3,$)')'Face',k
c     $                    ,nf_faceind(j,io,k)
c     $                    ,nf_dimlens(j,io,k2),nf_dimlens(j,io,k3)
                     if(itype.eq.4)then
c Structure position of start of covariant vectors
                        i1=pp_vec+ns_ndims*(k2-1)
                        i2=pp_vec+ns_ndims*(k3-1)
c face area equals magnitude of cross product. 
c 3d assumption for convenience.
                        ar2= (obj_geom(i1+1,i)*obj_geom(i2+2,i)
     $                    -obj_geom(i1+2,i)*obj_geom(i2+1,i))**2
     $                    +(obj_geom(i1+2,i)*obj_geom(i2  ,i)
     $                    -obj_geom(i1  ,i)*obj_geom(i2+2,i))**2
     $                    +(obj_geom(i1  ,i)*obj_geom(i2+1,i)
     $                    -obj_geom(i1+1,i)*obj_geom(i2  ,i))**2
                        ar=4.*sqrt(ar2)/
     $                    (nf_dimlens(j,io,k2)*nf_dimlens(j,io,k3))
                     else
c Cube
                        ar=4.*obj_geom(oradius+k2-1,i)
     $                       *obj_geom(oradius+k3-1,i)/
     $                    (nf_dimlens(j,io,k2)*nf_dimlens(j,io,k3))
                     endif
c Store data. At the moment, just the areas.
                     do j3=1,nf_dimlens(j,io,k3)
                     do j2=1,nf_dimlens(j,io,k2)
                        ip=j2+(j3-1)*int(nf_dimlens(j,io,k2))
     $                       +nf_faceind(j,io,k)
                        ff_data(nf_address(j,io,nf_pa)+ip-1)=ar
c                        write(*,'(i4,f8.4,$)')
c     $                       j2,j3
c     $                       ,int(nf_dimlens(j,io,k2))
c     $                       ,ip,ar
                        if(itype.eq.2)then
c Cube position data x,y,z
                           xr1=-obj_geom(oradius+k1-1,i)
                           if(k.gt.ns_ndims)xr1=-xr1
                           xr2=obj_geom(oradius+k2-1,i)*
     $                          (-1.+2.*(j2-0.5)/nf_dimlens(j,io,k2))
                           xr3=obj_geom(oradius+k3-1,i)*
     $                          (-1.+2.*(j3-0.5)/nf_dimlens(j,io,k3))
c                           write(*,*)j2,j3,xr1,xr2,xr3
                           ff_data(nf_address(j,io,1-k1)+ip-1)=
     $                          obj_geom(ocenter+k1-1,i)+xr1
                           ff_data(nf_address(j,io,1-k2)+ip-1)=
     $                          obj_geom(ocenter+k2-1,i)+xr2
                           ff_data(nf_address(j,io,1-k3)+ip-1)=
     $                          obj_geom(ocenter+k3-1,i)+xr3
                        endif
                     enddo
                     enddo
                  enddo
               enddo
            elseif(itype.eq.3 .or. itype.eq.5)then
c Cylinder and Non-aligned cylinder -------------------------------
               if(itype.eq.3)then
                  ica=obj_geom(ocylaxis,i)
                  rc=obj_geom(oradius+mod(ica,ns_ndims),i)
                  zr=obj_geom(oradius+ica-1,i)
                  zc=obj_geom(ocenter+ica-1,i)
                  if(rc.ne.obj_geom(oradius+mod(ica+1,ns_ndims),i))then
                     write(*,*)'Warning! Elliptical cylinder'
     $                    ,' surface areas'
     $                    ,' not calculated correctly.'
                  endif
               else
                  rc=1.
                  zr=1.
                  zc=0.
               endif
               do j=1,mf_quant(io)
c Area of each element. Faces are in order bottom, side, top.
c Bottom and top facet areas:
                  ar=3.1415926*rc**2
     $                 /(nf_dimlens(j,io,1)*nf_dimlens(j,io,2))
c                  write(*,*)'Object',i,' Quant',j,' Facet areas=',ar
                  do i2=1,nf_dimlens(j,io,2)
                     do i1=1,nf_dimlens(j,io,1)
                        ip=i1+(i2-1)*int(nf_dimlens(j,io,1))
c Positional data:
c r (not r^2)
                        r=rc*sqrt((i1-0.5)/nf_dimlens(j,io,1))
                        ff_data(nf_address(j,io,nf_pr)+ip-1)=r
                        ff_data(nf_address(j,io,nf_pr)+ip-1
     $                       +nf_faceind(j,io,3))=r
c theta
                        t=3.1415927*
     $                       (-1.+2.*(i2-0.5)/nf_dimlens(j,io,2))
c                        write(*,*)i1,nf_dimlens(j,io,1),t
                        ff_data(nf_address(j,io,nf_pt)+ip-1)=t
                        ff_data(nf_address(j,io,nf_pt)+ip-1
     $                       +nf_faceind(j,io,3))=t
c z
                        ff_data(nf_address(j,io,nf_pz)+ip-1)=zc-zr
                        ff_data(nf_address(j,io,nf_pz)+ip-1
     $                       +nf_faceind(j,io,3))=zc+zr
c area
                        ff_data(nf_address(j,io,nf_pa)+ip-1)=ar
                        ff_data(nf_address(j,io,nf_pa)+ip-1
     $                       +nf_faceind(j,io,3))=ar
                     enddo
                  enddo
c Side 2 pi r 2 z:
                  ar=4.*3.1415926*rc*zr
     $                 /(nf_dimlens(j,io,2)*nf_dimlens(j,io,3))
c index theta,z
                  do i2=1,nf_dimlens(j,io,3)
                     do i1=1,nf_dimlens(j,io,2)
                        ip=i1+(i2-1)*int(nf_dimlens(j,io,2))
                        ff_data(nf_address(j,io,nf_pa)+ip-1
     $                       +nf_faceind(j,io,2))=ar
                        ff_data(nf_address(j,io,nf_pr)+ip-1
     $                       +nf_faceind(j,io,2))=rc
c                        write(*,'(i4,f8.4,$)')
c     $                       ,ip,ar
c theta
                        t=3.1415927*
     $                       (-1.+2.*(i1-0.5)/nf_dimlens(j,io,2))
c                        write(*,*)i1,nf_dimlens(j,io,1),t
                        ff_data(nf_address(j,io,nf_pt)+ip-1
     $                       +nf_faceind(j,io,2))=t
c z
                        z=zc+zr*(-1.+2.*(i2-0.5)/nf_dimlens(j,io,3))
                        ff_data(nf_address(j,io,nf_pz)+ip-1
     $                       +nf_faceind(j,io,2))=z
                        

                     enddo
                  enddo
               enddo
            else
c Unknown ------------------------------------------
               write(*,*)'Flux positions uncalculated for object'
     $           ,i,' type',obj_geom(otype,i)
            endif
         endif
c      write(*,*)'Initialized positional data to:',nf_address(1,1,1)-1
c      write(*,'(10f8.4)')(ff_data(k),k=1,nf_address(1,1,1)-1)
      enddo

c      write(*,*)'10 steps of initialized slots:',nf_address(1,1,6)-1
c      write(*,'(10f8.4)')(ff_data(k),k=nf_address(1,1,1)
c     $     ,nf_address(1,1,11)-1)

      end
c******************************************************************
      subroutine tallyexit(i,idiffreg,ltlyerr)
c Document the exit of this particle just happened. 
c Assign the exit to a specific object, and bin on object.
c (If it is a mapped object, decided by objsect.)      
c On entry
c        i is particle number, 
c        idiffreg is the difference between its region and active.
c On exit ltlyerr is true if an error occurred else unchanged.
c Normally, returning an error will cause this particle to be
c considered to have left the particle region, so it will be discarded.
      integer i,idiffreg
      logical ltlyerr
      include 'partcom.f'
      include '3dcom.f'

      idiff=abs(idiffreg)
      idp=idiff
c Determine (all) the objects crossed and call objsect for each.
      iobj=0
 1    if(idiff.eq.0) return
      iobj=iobj+1
      idiff=idiff/2
      if(idp.ne.idiff*2)then
         call objsect(i,iobj,ierr)
         if(ierr.ne.0)then
            ireg=insideall(npdim,x_part(1,i))
            r=0.
            r1=0.
            do id=1,3
               r=r+x_part(id,i)**2
               r1=r1+(x_part(id,i)-dt*x_part(id+3,i))**2
            enddo
            r=sqrt(r)
            r1=sqrt(r1)
            write(*,*)'Tallyexit error',ierr,i,iobj,idiffreg
            if(ierr.eq.99)write(*,*)'Unknown object type.'
            write(*,*)'xpart,r=',(x_part(k,i),k=1,6),r,ireg
            write(*,*)'xp1',(x_part(k,i)-dt*x_part(k+3,i),k=1,3),r1
            ltlyerr=.true.
            return
         endif
      endif
      idp=idp/2
      goto 1
      
      end
c******************************************************************
c****************************************************************
      subroutine objsect(j,iobj,ierr)
c Find the intersection of the last step of particle j with  
c object iobj, and update the positioned-fluxes accordingly.
c ierr is returned: 0 good. 1 no intersection. 99 unknown object.
c
      include '3dcom.f'
      include 'partcom.f'
      include 'sectcom.f'

      real x1(npdim),x2(npdim)
      integer isc
      data isc/0/

      ierr=0

c Do nothing for untracked object and report no error.
      if(nf_map(iobj).eq.0)return
      infobj=nf_map(iobj)

      itype=int(obj_geom(otype,iobj))
c Use only bottom 8 bits:
      itype=itype-256*(itype/256)

c Get the positions:
      do i=1,npdim
         x1(i)=x_part(i,j)-dt*x_part(i+3,j)
         x2(i)=x_part(i,j)
      enddo
      if(itype.eq.1)then
c Sphere intersection. Return the bin number and direction ijbin,sd.
         call spherefsect(npdim,x1,x2,iobj,ijbin,sd,fraction)
      elseif(itype.eq.2)then
c Cube intersection. Return the bin number and direction ijbin,sd.
         call cubefsect(npdim,x1,x2,iobj,ijbin,sd,fraction)
      elseif(itype.eq.3)then
c Cylinder
         call cylfsect(npdim,x1,x2,iobj,ijbin,sd,fraction)
      elseif(itype.eq.4)then
c Parallelopiped intersection.
c         write(*,*)'Calling pllel',npdim,x1,x1,iobj
         call pllelofsect(npdim,x1,x2,iobj,ijbin,sd,fraction)
c         write(*,*)'Returning from pllel',npdim,x1,x1,iobj,ijbin,sd
      elseif(itype.eq.5)then
c Non-aligned cylinder
         call cylgfsect(npdim,x1,x2,iobj,ijbin,sd,fraction)
      else
c         write(*,*)'Unknown object in objsect'
c Unknown object type.
         ierr=99
         return
      endif
      if(abs(sd).ne.1.)then
c One end seems on exactly the boundary. Don't count this intersection.
         if(fraction.ne.1.)then
c There's a worse problem. Report the error.
            write(*,*)'sd fraction,type,No',sd,fraction,itype,iobj
            ierr=1
         endif
         return
      endif
c------------------------------
c Saving the x1 and intersection for diagnostics.
      if(sc_ipt.lt.sc_npts)then
         sc_ipt=sc_ipt+1
      endif
c Count up to total intersections. Then store cyclically. 
      isc=mod(isc,sc_npts)+1
      iob_sc(isc)=iobj
      ibin_sc(isc)=ijbin
c      write(*,*)'Saving intersection',isc,iobj,ijbin
      do i=1,sc_ndims
         x_sc(i,1,isc)=x1(i)
         x_sc(i,2,isc)=x1(i)*(1.-fraction)+x2(i)*fraction
      enddo
c------------------------------
c Adding into bins. At this point all we need is sd, ijbin.
c Particle Flux.
      iaddress=ijbin+nf_address(nf_flux,infobj,nf_step)
      ff_data(iaddress)=ff_data(iaddress)+sd
c This is the way to test that one is really accessing the right bin:
c      write(*,*)'ijbin=',ijbin,nf_posno(1,infobj),infobj,sd,fraction
c      ff_data(iaddress)=ijbin
c Perhaps ought to consider velocity interpolation.
      if(mf_quant(infobj).ge.2)then
c Momentum               
         iaddress=ijbin+nf_address(nf_gx,infobj,nf_step)
         ff_data(iaddress)=ff_data(iaddress)+ sd*x_part(4,j)
      endif
      if(mf_quant(infobj).ge.3)then
         iaddress=ijbin+nf_address(nf_gy,infobj,nf_step)
         ff_data(iaddress)=ff_data(iaddress)+ sd*x_part(5,j)
      endif
      if(mf_quant(infobj).ge.4)then
         iaddress=ijbin+nf_address(nf_gz,infobj,nf_step)
         ff_data(iaddress)=ff_data(iaddress)+ sd*x_part(6,j)
      endif
      if(mf_quant(infobj).ge.5)then
c Energy
         iaddress=ijbin+nf_address(nf_heat,infobj,nf_step)
         xx=0.
         do k=1,npdim
            xx=xx+x_part(3+k,j)**2
         enddo
         ff_data(iaddress)=ff_data(iaddress)+ sd*xx
      endif
c Accumulate the particle force= momentum/time over whole object.
c Normalized to rhoinf
      do id=1,ns_ndims
         partforce(id,infobj,nf_step)=partforce(id,infobj,nf_step)
     $        +sd*x_part(ns_ndims+id,j)/dt/rhoinf
      enddo
c            write(*,'(a,3f10.4,i4,i4)')'Partforce='
c     $           ,(partforce(kk,infobj,nf_step),kk=1,3)
c     $           ,infobj,nf_step
c If the bins were different we would have to recalculate ibin. 

      end
c********************************************************************
      subroutine spherefsect(npdim,xp1,xp2,iobj,ijbin,sd,fraction)
c Given npdimensional sphere, nf_map infobj, find the intersection of
c the line joining xp1, xp2 with it.  Use the 3dcom.f data to determine
c the flux bin to which that point corresponds and return the bin index
c in ijbin (zero-based), and the direction in which the sphere was
c crossed in sd. If no intersection is found return sd=0. Return
c the fractional position between xp1 and xp2 in frac
      integer npdim
      real xp1(npdim),xp2(npdim)
      integer iobj,ijbin
      real sd,fraction
      real tiny,onemtiny
      parameter (tiny=1.e-5,onemtiny=1.-2.*tiny)
      include '3dcom.f'
c 3D here.
      parameter (nds=3)
      real x12(nds)

      if(npdim.ne.nds)stop 'Wrong dimension number in spherefsect'
      infobj=nf_map(iobj)
      fraction=1.
      ida=0
      call sphereinterp(npdim,ida,xp1,xp2,
     $     obj_geom(ocenter,iobj),obj_geom(oradius,iobj),fraction
     $     ,f2,sd,C,D)
      if(sd.eq.0 .or. fraction-1..ge.tiny .or. fraction.lt.0.)then
c This section can be triggered inappropriately if rounding causes
c fraction to be >=1 when really the point is just outside or on the 
c surface. Then we get a sd problem message.
         fraction=1.
         sd=0.
         return
      endif
c This code decides which of the nf_posno for this object
c to update corresponding to this crossing, and then update it. 
c Calculate normalized intersection coordinates.
      do i=1,npdim
         x12(i)=((1.-fraction)*xp1(i)+fraction*xp2(i)
     $        -obj_geom(ocenter+i-1,iobj))
     $        /obj_geom(oradius+i-1,iobj)
      enddo
c Bin by cos(theta)=x12(3) uniform grid in first nf_dimension. 
c ibin runs from 0 to N-1 cos = -1 to 1.
      ibin=int(nf_dimlens(nf_flux,infobj,1)*(onemtiny*x12(3)+1.)*0.5)
      psi=atan2(x12(2),x12(1))
c jbin runs from 0 to N-1 psi = -pi to pi.
      jbin=int(nf_dimlens(nf_flux,infobj,2)
     $     *(0.999999*psi/3.1415926+1.)*0.5)
      ijbin=ibin+jbin*nf_dimlens(nf_flux,infobj,1)
      if(ijbin.gt.nf_posno(1,infobj))then
         write(*,*)'ijbin error in spherefsect'
         write(*,*)infobj,ijbin,nf_posno(1,infobj),ibin,jbin
     $        ,nf_dimlens(nf_flux,infobj,1),nf_dimlens(nf_flux,infobj,2)
     $        ,x12(3),obj_geom(ocenter+2,iobj),xp1(3),xp2(3)
     $        ,fraction
      endif
      end
c*********************************************************************
      subroutine sphereinterp(npdim,ida,xp1,xp2,xc,rc,f1,f2,sd,C,D)
c Given two different npdim dimensioned vectors xp1,xp2,and a sphere
c center xc radius rc, find the intersection of the line joining x1,x2,
c with the sphere and return it as the value of the fraction f1 of
c x1->x2 to which this corresponds, chosen always positive if possible, 
c and closest to 0. The other intersection fraction in f2.
c Also return the direction of crossing in sd, and the fractional
c radial distance^2 outside the sphere of the two points in C and D. 
c (positive means inward from x1 to x2). If there is no intersection,
c return fraction=1., sd=0.  If ida is non-zero then form the radius
c only over the other dimensions.  In that case the subsurface (circle)
c is the figure whose intersection is sought.
      integer npdim,ida
      real xp1(npdim),xp2(npdim),xc(npdim),rc(npdim)
      real sd,f1,f2,C,D
c
      real x1,x2,A,B
c Prevent a singularity if the vectors coincide.
      A=1.e-25
      B=0.
      C=-1.
      D=-1.
c x1 and x2 are the coordinates in system in which sphere 
c has center 0 and radius 1.
      ni=npdim
      if(ida.ne.0)ni=npdim-1
      do ii=1,ni
         i=mod(ida+ii-1,npdim)+1
         xci=xc(i)
         rci=rc(i)
         x1=(xp1(i)-xci)/rci
         x2=(xp2(i)-xci)/rci
         A=A+(x2-x1)**2
         B=B+x1*(x2-x1)
         C=C+x1**2
         D=D+x2**2
      enddo
c Crossing direction from x1 to intersection
      sd=sign(1.,C)
      disc=B*B-A*C
      if(disc.ge.0.)then
         disc=sqrt(disc)
c (A always positive)
         if(C.gt.0. .and. B.lt.0) then
c Discrim<|B|. Can take minus sign and still get positive. 
            f1=(-B-disc)/A
            f2=(-B+disc)/A
         else
c If B positive or C negative must use plus sign. 
            f1=(-B+disc)/A
            f2=(-B-disc)/A
         endif
      else
c         write(*,*)'Sphere-crossing discriminant negative.',A,B,C
         sd=0.
         f1=1.
         f2=1.
         return
      endif
      end
c*********************************************************************
      subroutine cubefsect(npdim,xp1,xp2,iobj,ijbin,sd,fraction)
c Given a coordinate-aligned cube object iobj. Find the point of
c intersection of the line joining xp1,xp2, with it, and determine the
c ijbin to which it is therefore assigned, and the direction it is
c crossed (sd=+1 means inward from 1 to 2). 
c Intersection fractional distance from xp1 to xp2 returned.

c The cube is specified by center and radii!=0 (to faces.) which define
c two planes (\pm rc) in each coordinate. Inside corresponds to between
c these two planes, i.e. x-xc < |rc|. We define the facets of the cube
c to be the faces (planes) in the following order:
c +rc_1,+rc_2,+rc_3,-rc_1,-rc_2,-rc_3 which is the order of the
c coefficients of the adjacent vectors. 

      integer npdim,iobj,ijbin
      real xp1(npdim),xp2(npdim)
      real sd
      include '3dcom.f'
      real xn1(ns_ndims),xn2(ns_ndims)
      sd=0.
      ig1=inside_geom(npdim,xp1,iobj)
      ig2=inside_geom(npdim,xp2,iobj)
      if(ig1.eq.1 .and. ig2.eq.0)then
c xp1 inside & 2 outside
         inside=0
      elseif(ig2.eq.1 .and. ig1.eq.0)then
c xp2 inside & 1 outside
         inside=1
      else
c Both inside or outside. (Used to be neither inside). Any intersection
c will be disallowed. This means we cut off any such edges of the
c cube. But actually for coordinate-aligned cube there are no
c problematic non-normal intersections.
         fraction=1.
         return
      endif
c In direction sd
      sd=2*inside-1.
c Package from here by converting into normalized position
      do i=1,npdim
         xn1(i)=(xp1(i)-obj_geom(ocenter+i-1,iobj))
     $        /obj_geom(oradius+i-1,iobj)
         xn2(i)=(xp2(i)-obj_geom(ocenter+i-1,iobj))
     $        /obj_geom(oradius+i-1,iobj)
      enddo
c And calling the unit-cube version.
      if(inside.eq.0)then
         call cubeexplt(npdim,xn1,xn2,ijbin,iobj,fraction)
      else
         call cubeexplt(npdim,xn2,xn1,ijbin,iobj,fraction)
         fraction=1.-fraction
      endif
c Code Diagnostic:
c      do i=1,3
c         xd(i)=(1.-fraction)*xp1(i)+ fraction*xp2(i)
c         xd2(i)=(1.-fraction)*xp2(i)+ fraction*xp1(i)
c      enddo
c      write(*,'(7f8.4)')fraction,xp1,xp2,xd,xd2
      end
c*********************************************************************
      subroutine pllelofsect(npdim,xp1,xp2,iobj,ijbin,sd,fraction)
c Given a general parallelopiped object iobj. Find the point
c of intersection of the line joining xp1,xp2, with it, and determine
c the ijbin to which it is therefore assigned, and the direction it is
c crossed (sd=+1 means inward from 1 to 2).

c The object is specified by center and three vectors pqr. Each of npdim
c pairs of parallel planes consists of the points: +-p + c_q q + c_r r.
c Where p is one of the three base (covariant) vectors and qr the
c others, and c_q are real coefficients.  Inside corresponds to between
c these two planes, i.e. contravariant coefficients <1. We define the
c facets of the cube to be the faces (planes) in the following order:
c +v_1,+v_2,+v_3,-v_1,-v_2,-v_3. 
c Then within each face the facet indices are in cyclic order. But that
c is determined by the cubeexplt code.
      integer npdim,iobj,ijbin
      real xp1(npdim),xp2(npdim)
      real sd
      include '3dcom.f'
      real xn1(pp_ndims),xn2(pp_ndims)
      sd=0.

c      write(*,*)'Pllelo',npdim,xp1,xp2,iobj,pp_ndims
      ins1=0
      ins2=0
      do j=1,pp_ndims
         xn1(j)=0.
         xn2(j)=0.
c Contravariant projections.
         do i=1,npdim
c Cartesian coordinates.
            ii=(ocenter+i-1)
            xc=obj_geom(ii,iobj)
c xn1, xn2 are the contravariant coordinates with respect to the center.
            ji=(pp_contra+pp_ndims*(j-1)+i-1)
c            write(*,*)'ji',ji
            xn1(j)=xn1(j)+(xp1(i)-xc)*obj_geom(ji,iobj)
            xn2(j)=xn2(j)+(xp2(i)-xc)*obj_geom(ji,iobj)
         enddo
         if(abs(xn1(j)).ge.1.)ins1=1
         if(abs(xn2(j)).ge.1.)ins2=1
      enddo
c ins1,2 indicate inside (0) or outside (1) for each point. 
c In direction sd
      sd=2*ins1-1.
c And calling the unit-cube version.
c      write(*,*)'Calling cubeexplt',xn1,xn2
      if(ins1.eq.0 .and. ins2.eq.1)then
         call cubeexplt(npdim,xn1,xn2,ijbin,iobj,fraction)
      elseif(ins2.eq.0 .and. ins1.eq.1)then
         call cubeexplt(npdim,xn2,xn1,ijbin,iobj,fraction)
         fraction=1.-fraction
      else
         fraction=1.
         return
      endif
      end
c*******************************************************************
      subroutine pllelfrac(xp,xn,objg)
c For input point xp(ndims) return the normalized position relative to
c the parallelogram object objg in output xn(ndims)
      include '3dcom.f'
      real xp(ns_ndims),xn(ns_ndims)
      real objg(odata)
      do j=1,pp_ndims
         xn(j)=0.
c Contravariant projections.
         do i=1,ns_ndims
c Cartesian coordinates.
            ii=(ocenter+i-1)
            xc=objg(ii)
c xn1, xn2 are the contravariant coordinates with respect to the center.
            ji=(pp_contra+pp_ndims*(j-1)+i-1)
c            write(*,*)'ji',ji
            xn(j)=xn(j)+(xp(i)-xc)*objg(ji)
         enddo
      enddo
      end
c***********************************************************************
      subroutine cubeexplt(npdim,xp1,xp2,ijbin,iobj,fmin)
c For a unit cube, center 0 radii 1, find the intersection of line
c xp1,xp2 with it and return the relevant flux index ijbin.
c xp1 must be always inside the cube.
      real xp1(npdim),xp2(npdim)
      integer ijbin,iobj
      include '3dcom.f'
      integer idebug
      imin=0
      idebug=0

      fn=0.
      fmin=100.
      do i=1,2*npdim
         im=mod(i-1,npdim)+1
c First half of the i's are negative. Second half positive.
         xc1=(((i-1)/npdim)*2-1)
         xd1=(xp1(im)-xc1)
         xd2=(xp2(im)-xc1)
         if(xd1.lt.0. .neqv. xd2.lt.0)then
c Crossed this plane
            fn=xd1/(xd1-xd2)
c At fraction fn (always positive) from the inside point.
            if(fn.lt.fmin)then
               fmin=fn
               imin=i
c         write(*,*)i,xc1,' xd1,xd2=',xd1,xd2,fmin
c     $              ,xp1(im),xp2(im)
c     $        ,(1.-fmin)*xd1+fmin*xd2
c     $        ,(1.-fmin)*xp1(im)+fmin*xp2(im)
            endif
         endif
      enddo
c Now the minimum fraction is in fmin, which is the crossing point.
c imin contains the face-index of this crossing. Indexing within the
c face on equal spaced grid whose numbers have been read into
c obj_geom(ofn.,iobj):
      infobj=nf_map(iobj)
      ibstep=1
      ibin=0
      if(idebug.eq.1)write(*,'(''Position'',$)')
      do i=1,npdim-1
c The following defines the order of indexation. ofn1 is the next highest
c cyclic index following the face index. So on face 1 or 4 the other
c two indices on the face are 2,3. But on face 2,5 they are 3,1.
         k=mod(mod(imin-1,3)+1+i-1,npdim)+1
         xk=(1.-fmin)*xp1(k)+fmin*xp2(k)
c Not sure that this is the best order for the plane. Think! :
c         xcr=(1.-xk)*.5
c This has xcr run from 0 to 1 as xk goes from -1. to +1. :
         xcr=(1.+xk)*.5
         if(idebug.eq.1)write(*,'(i2,f7.2,i5,i5,'',''$)')k,xcr
     $        ,nf_dimlens(nf_flux,infobj,k)
     $        ,int(nf_dimlens(nf_flux,infobj,k)*(0.999999*xcr))
         ibin=ibin+ibstep*
     $        int(nf_dimlens(nf_flux,infobj,k)*(0.999999*xcr))
         ibstep=ibstep*nf_dimlens(nf_flux,infobj,k)
      enddo
c Now we have ibin equal to the face-position index, and ibstep equal
c to the face-position size. Add the face-offset for imin. This is
c tricky for unequal sized faces. So we need to have stored it. 
      ijbin=ibin+nf_faceind(nf_flux,infobj,imin)
      if(idebug.eq.1)write(*,*)'Ending cubefsect',ijbin
     $     ,nf_faceind(nf_flux,infobj,imin),imin
c That's it.
      end
c*********************************************************************
      subroutine cylfsect(npdim,xp1,xp2,iobj,ijbin,sdmin,fmin)
c Master routine for calling cylusect after normalization of cyl.
      integer npdim,iobj,ijbin
      real xp1(npdim),xp2(npdim)
      real sdmin
      include '3dcom.f'
      real xn1(pp_ndims),xn2(pp_ndims)

      ida=int(obj_geom(ocylaxis,iobj))
      do i=1,pp_ndims
         ii=mod(i-ida+2,pp_ndims)+1
         xn1(ii)=(xp1(i)-obj_geom(ocenter+i-1,iobj))
     $        /obj_geom(oradius+i-1,iobj)
         xn2(ii)=(xp2(i)-obj_geom(ocenter+i-1,iobj))
     $        /obj_geom(oradius+i-1,iobj)
      enddo
      call cylusect(npdim,xn1,xn2,iobj,ijbin,sdmin,fmin)
      end
c*********************************************************************
      subroutine cylusect(npdim,xp1,xp2,iobj,ijbin,sdmin,fmin)
c Find the point of intersection of the line joining xp1,xp2, with the
c UNIT cylinder, and determine the ijbin to which it is therefore
c assigned, and the direction it is crossed (sdmin=+1 means inward from
c 1 to 2). The 1-axis is where theta is measured from and the 3-axis
c is the axial direction.
c The facets of the cylinder are the end faces -xr +xr, and the curved
c side boundary. 3 altogether.  The order of faces is bottom, side, top.
      
      integer npdim,iobj,ijbin
      real xp1(npdim),xp2(npdim)
      real sdmin
      include '3dcom.f'
c 3D here.
      parameter (nds=3)
      real x12(nds)
      real fn(4),zrf(4),sdf(4)
      real xc(nds),rc(nds)
      data xc/0.,0.,0./rc/1.,1.,1./

      ida=3
c First, return if both points are beyond the same axial end.
      z1=xp1(ida)
      z2=xp2(ida)
      xd=z2-z1
      fmin=1.
      sdmin=0.
      if((z1.gt.1. .and. z2.gt.1.).or.
     $     (-z1.gt.1. .and. -z2.gt.1))return
      sds=0.
c Find the intersection (if any) with the circular surface.
      call sphereinterp(npdim,ida,xp1,xp2,
     $     xc,rc,fn(1),fn(2),sds,d1,d2)
      if(sds.ne.0)then
c Directions are both taken to be that of the closest. 
c A bit inconsistent but probably ok. 
         sdf(1)=sds
         sdf(2)=sds
         zrf(1)=(1.-fn(1))*xp1(ida)+fn(1)*xp2(ida)
         zrf(2)=(1.-fn(2))*xp1(ida)+fn(2)*xp2(ida)
      else
c No radial intersections
         zrf(1)=2.
         zrf(2)=2.
      endif
c Find the axial intersection fractions with the end planes.
      if(xd.ne.0)then
         fn(3)=(1.-z1)/xd
         sdf(3)=-1.
         if(z1.gt.1.)sdf(3)=1.
         fn(4)=(-1.-z1)/xd
         sdf(4)=-1.
         if(z1.lt.-1.)sdf(4)=1.
         zrf(3)=0.
         zrf(4)=0.
         do k=1,npdim
            if(k.ne.ida)then
               xkg1=(1.-fn(3))*xp1(k)+fn(3)*xp2(k)
               zrf(3)=zrf(3)+(xkg1)**2
               xkg2=(1.-fn(4))*xp1(k)+fn(4)*xp2(k)
               zrf(4)=zrf(4)+(xkg2)**2
            endif
         enddo
      else
c Pure radial difference. No end-intersections anywhere.
         zrf(3)=2.
         zrf(4)=2.
      endif
c Now we have 4 possible fractions fn(4). Two or none of those
c are true. Truth is determined by abs(zrf(k))<=1. Choose closest.
      fmin=10.
      kmin=0
      do k=1,4
         if(abs(zrf(k)).le.1)then
            if(fn(k).ge.0. .and. fn(k).lt.fmin)then
               kmin=k 
               fmin=fn(k)
            endif
         endif
      enddo
      if(fmin.gt.1.)then
c No crossing
         sdmin=0.
         fmin=1.
         return
      else
         sdmin=sdf(kmin)
         if(kmin.le.2)then
c radial crossing
            imin=0.
         else
c axial crossing
            imin=-1.
            zida=(1.-fmin)*z1+fmin*z2
            if(zida.gt.0.)imin=1.
         endif
      endif

c Now the minimum fraction is in fmin, which is the crossing point.
c imin contains the face-index of this crossing. -1,0, or +1.
c Calculate normalized intersection coordinates.
      do i=1,npdim
         x12(i)=(1.-fmin)*xp1(i)+fmin*xp2(i)
      enddo
c Calculate r,theta,z (normalized) relative to the ida direction as z.
      z=x12(ida)
      theta=atan2(x12(mod(ida+1,npdim)+1),x12(mod(ida,npdim)+1))
      r2=0.
      do i=1,npdim-1
         k=mod(ida+i-1,npdim)+1
         r2=r2+x12(k)**2
      enddo
c      write(*,'(a,7f7.4,3i3)')'r2,theta,z,x12,fmin,imin'
c     $     ,r2,theta,z,x12,fmin,imin
c End blocks are of size nr x nt, and the curved is nt x nz.
c 3-D only here. 
      infobj=nf_map(iobj)
      ijbin=0.
      if(imin.ne.0)then
c Ends
         if(imin.eq.1)then
c offset by (nr+nz)*nt
            ijbin=(nf_dimlens(nf_flux,infobj,1)
     $           +nf_dimlens(nf_flux,infobj,3))
     $           *nf_dimlens(nf_flux,infobj,2)
         endif
c Uniform mesh in r^2 normalized. 
         ir=int(nf_dimlens(nf_flux,infobj,1)*(0.999999*r2))
         it=int(nf_dimlens(nf_flux,infobj,2)
     $     *(theta/3.1415927+1.)*0.5)
         ijbin=ijbin+ir+it*nf_dimlens(nf_flux,infobj,1)
      else
c Side. Offset to this facet nr*nt:
         ijbin=nf_dimlens(nf_flux,infobj,1)*nf_dimlens(nf_flux,infobj,2)
         it=int(nf_dimlens(nf_flux,infobj,2)
     $     *(theta/3.1415927+1.)*0.5)
         iz=int(nf_dimlens(nf_flux,infobj,3)*(0.999999*z+1.)*0.5)
c Index in order theta,z
         ijbin=ijbin+it+iz*nf_dimlens(nf_flux,infobj,2)
      endif
c      write(*,'(6f8.4,3i3)')xp1,xp2,ir,it,iz
      end
c*********************************************************************
      subroutine cylgfsect(npdim,xp1,xp2,iobj,ijbin,sd,fraction)
c Given a general cylinder object iobj. Find the point
c of intersection of the line joining xp1,xp2, with it, and determine
c the ijbin to which it is therefore assigned, and the direction it is
c crossed (sd=+1 means inward from 1 to 2).

c The object is specified by center and three contravariant vectors.
c When the position relative to the center is dotted into the contra
c variant vector it yields the coordinate relative to the unit cylinder,
c whose third component is the axial direction. 

      integer npdim,iobj,ijbin
      real xp1(npdim),xp2(npdim)
      real sd
      include '3dcom.f'
      real xn1(pp_ndims),xn2(pp_ndims)

c j refers to transformed coordinates in which it is unit cyl
      do j=1,pp_ndims
         xn1(j)=0.
         xn2(j)=0.
c Contravariant projections.
         do i=1,npdim
c i refers to the Cartesian coordinates.
            xc=obj_geom(ocenter+i-1,iobj)
c xn1, xn2 are the contravariant coordinates with respect to the center.
            ji=(pp_contra+pp_ndims*(j-1)+i-1)
c            write(*,*)'ji',ji
            xn1(j)=xn1(j)+(xp1(i)-xc)*obj_geom(ji,iobj)
            xn2(j)=xn2(j)+(xp2(i)-xc)*obj_geom(ji,iobj)
         enddo
      enddo
c Now xn1,2 are the coordinates relative to the unit cylinder.      
      fraction=1.
c Shortcut
      z1=xn1(pp_ndims)
      z2=xn2(pp_ndims)
      if((z1.ge.1..and.z2.ge.1).or.(z1.le.-1..and.z2.le.-1.))return
c Call the unit-cylinder code.
      call cylusect(npdim,xn1,xn2,iobj,ijbin,sd,fraction)
c      write(*,*)'cylusect return',ijbin,fraction

      end
c*********************************************************************
c*********************************************************************
      subroutine timeave(nu,u,uave,ictl)
c Average a quantity u(nu) over steps with a certain decay number
c into uave.
c ictl controls the actions as follows:
c       0    do the averaging.
c       1    do nothing except set nstep=1.
c       2    do nothing except set nstave=nu
c       3    do both 1 and 2. 
c       99   do nothing except increment nstep.
c The 99 call should be done at the end of all usage of this routine
c for the present step.
      real u(nu),uave(nu)

      integer nstep,nstave
      data nstep/1/nstave/20/
c Normal call.
      if(ictl.eq.0)then
         do i=1,nu
            uave(i)=(uave(i)*(nstep-1)+u(i))/nstep
         enddo
         return
      endif
      
      if(ictl.eq.1 .or. ictl.eq.3)then
         nstep=1
      endif
      if(ictl.eq.2 .or. ictl.eq.3)then
         nstave=nu
      endif
      if(ictl.ge.99)then
         if(nstep.le.nstave) nstep=nstep+1
      endif

      end
c***********************************************************************
      subroutine fluxdiag()
c Get the total count to object 1 and convert to normalized flux.
c Assuming it's a unit sphere.
      include '3dcom.f'
c For rhoinf, dt
      include 'partcom.f'

      sum=0
      do i=1,nf_posno(nf_flux,1)
         sum=sum+ff_data(nf_address(nf_flux,1,nf_step)+i-1)
      enddo
      flux=sum/(4.*3.14159)/rhoinf/dt
      write(*,'(f6.3,''| '',$)')flux
c      write(*,'(a,f7.0,f8.3)')'Total flux',sum,
c     $     sum/(4.*3.14159)/rhoinf/dt
c      write(*,'(10f7.1)')(ff_data(nf_address(1,1,nf_step)+i-1),
c     $     i=1,nf_posno(1,1))

      end
c***********************************************************************
c Averaging the flux data for quantity abs(iquant)
c over all positions for object ifobj,
c for the steps n1 to n2.
c The positions might be described by more than one dimension, but
c that is irrelevant to the averaging.
c Plot the quantity iquant if positive (not if negative).
c Plotting does not attempt to account for the multidimensionality.
c Rhoinf is returned in rinf.
      subroutine fluxave(n1,n2,ifobj,iquant,rinf)
      integer n1,n2
c      logical lplot
      integer iquant
      include '3dcom.f'
      include 'sectcom.f'
c Use for averaging:ff_data(nf_address(iq,ifobj,nf_maxsteps+1)+i-1)
c which is ff_data(iav+i)
      real fluxofstep(nf_maxsteps),step(nf_maxsteps)
      character*30 string

      iq=abs(iquant)
c If quantity asked for is not available, do nothing.
      if(iq.gt.mf_quant(ifobj).or.iq.eq.0)then
c         write(*,*)ifobj,iq,mf_quant(ifobj)
         return
      endif
c Offset of averaging location:
c      iav=nf_address(iq,ifobj,nf_maxsteps+1)-1
      iav=nf_address(iq,ifobj,nf_step+1)-1
c Offset to average flux density location
c      iavd=nf_address(iq,ifobj,nf_maxsteps+2)-1
      iavd=nf_address(iq,ifobj,nf_step+2)-1
c Offset to area
      iaa=nf_address(iq,ifobj,nf_pa)-1

c      write(*,*)'Fluxave addresses',iav,iavd,iaa
c Check step numbers for rationality.
      if(n1.lt.1)n1=1
      if(n2.gt.nf_step)n2=nf_step
      if(n2-n1.lt.0)then
         write(*,*)'fluxave incorrect limits:',n1,n2
         return
      endif
      do i=1,nf_posno(iq,ifobj)
         ff_data(iav+i)=0.
         ff_data(iavd+i)=0.
      enddo
c Sum the data over steps.
      tot=0
      do i=1,nf_posno(iq,ifobj)
         do is=n1,n2
            ff_data(iav+i)=ff_data(iav+i)
     $           +ff_data(nf_address(iq,ifobj,is)+i-1)
         enddo
         tot=tot+ff_data(iav+i)
         ff_data(iav+i)=ff_data(iav+i)/(n2-n1+1)
c Get the area and divide to give the flux density.
         area=ff_data(iaa+i)
         ff_data(iavd+i)=ff_data(iav+i)/area
c         write(*,*)ifobj,iq,i,' Area,Flux-density',area,ff_data(iavd +i)
      enddo

c Total the flux over positions as a function of step.
      tdur=0.
      rinf=0.
      do is=1,n2
         fluxstep=0
         do i=1,nf_posno(iq,ifobj)
            fluxstep=fluxstep+ff_data(nf_address(iq,ifobj,is)+i-1)
         enddo
         step(is)=is
         fluxofstep(is)=fluxstep
         if(is.ge.n1)then
c            write(*,*)'n1,n2,is,ff_dt(is)',n1,n2,is,ff_dt(is)
            tdur=tdur+ff_dt(is)
            rinf=rinf+ff_rho(is)*ff_dt(is)
         endif
c         write(*,*)'step',is,' fluxofstep',fluxofstep(is)
      enddo
c      write(*,*)'tot,rinf,tdur,n2',tot,rinf,tdur,n1,n2
      tot=tot/tdur
      rinf=rinf/tdur

c From here on is non-general and is mostly for testing.
      write(*,'(a,i3,a,i3,a,i4,i4,a,f10.3)')' Average flux quant',iq
     $     ,', object',ifobj,', over steps',n1,n2,', per unit time:',tot
      write(*,*)'rhoinf:',rinf,' Total:',nint(tot*tdur)
     $     ,'  Average collected per step by posn:'
      write(*,'(10f8.2)')(ff_data(iav+i),i=1,nf_posno(iq,ifobj))
      fluxdensity=tot/(4.*3.14159)/rinf
      write(*,*)'Flux density*r^2, normalized to rhoinf'
     $     ,fluxdensity
c      write(*,*)'Sectcom ipt=',sc_ipt
      rmtoz=1.
      flogfac=0.5*alog(2.*3.1415926/(rmtoz*1837.))
      phifloat=alog(fluxdensity)+flogfac
      write(*,*)'Floating potential=',phifloat

      if(iquant.gt.0)then
         write(string,'(''Object '',i3,'' Quantity'',i3)')
     $        nf_geommap(ifobj),iquant
         call autoplot(step,fluxofstep,n2)
         call boxtitle(string)
         call axlabels('step','Spatially-summed flux number')
         call pltend()
         call automark(ff_data(nf_address(iq,ifobj,nf_p1))
     $        ,ff_data(iav+1),nf_posno(iq,ifobj),1)
         call boxtitle(string)
         call axlabels('First flux face variable',
     $        'Time-averaged flux number')
         call pltend()
      endif

      end
c*******************************************************************
      subroutine writefluxfile(name)
c File name:
      character*(*) name
c Common data containing the BC-object geometric information
      include '3dcom.f'
c Particle common data
      include 'partcom.f'
c Plasma common data
      include 'plascom.f'
c Intersection data
      include 'sectcom.f'

      character*(100) charout
c Zero the name first. Very Important!
c Construct a filename that contains many parameters
c Using the routines in strings_names.f
      name=' '
      call nameconstruct(name)
c     np=nbcat(name,'.flx')
      call nbcat(name,'.flx')
c      write(*,*)name
      write(charout,51)debyelen,Ti,vd,rs,phip
 51   format('debyelen,Ti,vd,rs,phip:',5f10.4)

c      write(*,*)'mf_obj=',mf_obj,nf_step,mf_quant(1)

      open(22,file=name,status='unknown',err=101)
      close(22,status='delete')
      open(22,file=name,status='new',form='unformatted',err=101)
c This write sequence must be exactly that read below.
      write(22)charout
      write(22)debyelen,Ti,vd,rs,phip
      write(22)nf_step,mf_quant,mf_obj,(nf_geommap(j),j=1,mf_obj)
      write(22)(ff_rho(k),k=1,nf_step)
      write(22)(ff_dt(k),k=1,nf_step)
      write(22)((nf_posno(i,j),(nf_dimlens(i,j,k),k=1,nf_ndims)
     $     ,(nf_faceind(i,j,k),k=1,2*nf_ndims)
     $     ,i=1,mf_quant(j)),j=1,mf_obj)
      write(22)(((nf_address(i,j,k),i=1,mf_quant(j)),j=1,mf_obj),
     $     k=1-nf_posdim,nf_step+2)
      write(22)(ff_data(i),i=1,nf_address(1,1,nf_step+2)-1)
c New force write.
      write(22)(((fieldforce(i,j,k),pressforce(i,j,k) ,partforce(i,j,k)
     $     ,charge_ns(j,k),i=1,ns_ndims),j=1,mf_obj),k=1,nf_step)
c Object data:
      write(22)ngeomobj
      write(22)((obj_geom(j,k),j=1,odata),nf_map(k),k=1,ngeomobj)
      write(22)ibool_part,ifield_mask,iptch_mask,lboundp,rjscheme
c Intersection data:
      write(22)sc_ipt
      write(22)(((x_sc(j,i,k),j=1,sc_ndims),i=1,2),iob_sc(k),
     $     ibin_sc(k),k=1,sc_ipt)

      close(22)
c      write(*,*)'Wrote flux data to ',name(1:lentrim(name))

      return
 101  continue
      write(*,*)'Error opening file:',name
      close(22,status='delete')
      end
c*****************************************************************
      subroutine readfluxfile(name,ierr)
      character*(*) name
      include '3dcom.f'
      include 'plascom.f'
      include 'sectcom.f'
      character*(100) charout

      open(23,file=name,status='old',form='unformatted',err=101)
      read(23)charout
      read(23)debyelen,Ti,vd,rs,phip
      read(23)nf_step,mf_quant,mf_obj,(nf_geommap(j),j=1,mf_obj)
      read(23)(ff_rho(k),k=1,nf_step)
      read(23)(ff_dt(k),k=1,nf_step)
      read(23)((nf_posno(i,j),(nf_dimlens(i,j,k),k=1,nf_ndims)
     $     ,(nf_faceind(i,j,k),k=1,2*nf_ndims)
     $     ,i=1,mf_quant(j)),j=1,mf_obj)
      read(23)(((nf_address(i,j,k),i=1,mf_quant(j)),j=1,mf_obj),
     $     k=1-nf_posdim,nf_step+2)
      read(23)(ff_data(i),i=1,nf_address(1,1,nf_step+2)-1)
      read(23,end=102, err=102)(((fieldforce(i,j,k),pressforce(i,j,k)
     $     ,partforce(i,j,k),charge_ns(j,k),i=1,ns_ndims),j=1,mf_obj),k
     $     =1,nf_step)
c Object data:
      read(23)ngeomobj
      read(23)((obj_geom(j,k),j=1,odata),nf_map(k),k=1,ngeomobj)
      read(23)ibool_part,ifield_mask,iptch_mask,lboundp,rjscheme
c      write(*,*)'Object data for',ngeomobj,' objects:'
c      write(*,*)((obj_geom(j,k),j=1,odata),nf_map(k),k=1,ngeomobj)
c      write(*,*)ibool_part,ifield_mask,iptch_mask,lboundp,rjscheme
c Intersection data:
      read(23)sc_ipt
      read(23)(((x_sc(j,i,k),j=1,sc_ndims),i=1,2),iob_sc(k),
     $     ibin_sc(k),k=1,sc_ipt)
      goto 103
 102  write(*,*)'Failed to read back forces. Old format?'
 103  close(23)

c Now one might have to reconstruct nf_faceind from nf_dimlens.

      write(*,*)'Read back flux data from ',name(1:lentrim(name))
      write(*,*)charout(1:lentrim(charout))
      ierr=0
      return
 101  write(*,*)'Error opening file:',name
      ierr=1
      end
c*********************************************************************
c Obsolete
c*********************************************************************
      subroutine cylfsect1(npdim,xp1,xp2,iobj,ijbin,sdmin,fmin)
c Given a coordinate-aligned cylinder object iobj. Find the point of
c intersection of the line joining xp1,xp2, with it, and determine the
c ijbin to which it is therefore assigned, and the direction it is
c crossed (sdmin=+1 means inward from 1 to 2).  

c The cylinder is specified by center and radii!=0 (to faces.)  in each
c coordinate. Plus the axial coordinate.  Inside corresponds to between
c the axial planes, i.e. x-xc < |rc|, and orthogonal radius < 1.

c The facets of the cylinder to be the end faces -xr +xr, and the curved
c side boundary. 3 altogether.  The order of faces is bottom, side, top.
      
      integer npdim,iobj,ijbin
      real xp1(npdim),xp2(npdim)
      real sdmin
      include '3dcom.f'
c 3D here.
      parameter (nds=3)
      real x12(nds)
      real fn(4),zrf(4),sdf(4)

      ida=int(obj_geom(ocylaxis,iobj))
c First, return if both points are beyond the same axial end.
      xca=obj_geom(ocenter+ida-1,iobj)
      xra=obj_geom(oradius+ida-1,iobj)
      xd1=(xp1(ida)-xca)
      xd2=(xp2(ida)-xca)
      xd=xd2-xd1
      fmin=1.
      sdmin=0.
      if((xd1.gt.xra .and. xd2.gt.xra).or.
     $     (-xd1.gt.xra .and. -xd2.gt.xra))return
      sds=0.
c Find the intersection (if any) with the circular surface.
      call sphereinterp(npdim,ida,xp1,xp2,
     $     obj_geom(ocenter,iobj),obj_geom(oradius,iobj),
     $     fn(1),fn(2),sds,d1,d2)
      if(sds.ne.0)then
c Directions are both taken to be that of the closest. 
c A bit inconsistent but probably ok. 
         sdf(1)=sds
         sdf(2)=sds
         zrf(1)=((1.-fn(1))*xp1(ida)+fn(1)*xp2(ida)
     $        -obj_geom(ocenter+ida-1,iobj))
     $        /obj_geom(oradius+ida-1,iobj)
         zrf(2)=((1.-fn(2))*xp1(ida)+fn(2)*xp2(ida)
     $        -obj_geom(ocenter+ida-1,iobj))
     $        /obj_geom(oradius+ida-1,iobj)
      else
c No radial intersections
         zrf(1)=2.
         zrf(2)=2.
      endif
c Find the axial intersection fractions with the end planes.
      if(xd.ne.0)then
         fn(3)=(xra-xd1)/(xd2-xd1)
         sdf(3)=-1.
         if(xd1.gt.xra)sdf(3)=1.
         fn(4)=(-xra-xd1)/(xd2-xd1)
         sdf(4)=-1.
         if(xd1.lt.-xra)sdf(4)=1.
         zrf(3)=0.
         zrf(4)=0.
c         z3=xp1(ida)*(1.-fn(3))+xp2(ida)*fn(3)
c         z4=xp1(ida)*(1.-fn(4))+xp2(ida)*fn(4)
c         write(*,*)'Axial crossings',fn(3),z3,fn(4),z4
         do k=1,npdim
            if(k.ne.ida)then
               xrk=obj_geom(oradius+k-1,iobj)
               xkg1=(1.-fn(3))*xp1(k)+fn(3)*xp2(k)
     $              -obj_geom(ocenter+k-1,iobj)
               zrf(3)=zrf(3)+(xkg1/xrk)**2
               xkg2=(1.-fn(4))*xp1(k)+fn(4)*xp2(k)
     $              -obj_geom(ocenter+k-1,iobj)
               zrf(4)=zrf(4)+(xkg2/xrk)**2
c               write(*,*)k,zrf(3),zrf(4),xrk,xd1,xd2,xra
            endif
         enddo
      else
c Pure radial difference. No end-intersections anywhere.
         zrf(3)=2.
         zrf(4)=2.
      endif
c Now we have 4 possible fractions fn(4). Two or none of those
c are true. Truth is determined by abs(zrf(k))<=1. Choose closest.
      fmin=10.
      kmin=0
      do k=1,4
         if(abs(zrf(k)).le.1)then
            if(fn(k).ge.0. .and. fn(k).lt.fmin)then
               kmin=k
               fmin=fn(k)
            endif
         endif
      enddo
      if(fmin.gt.1.)then
c No crossing
         sdmin=0.
         fmin=1.
         return
      else
         sdmin=sdf(kmin)
         if(kmin.le.2)then
c radial crossing
            imin=0.
         else
c axial crossing
            imin=-1.
            zida=(1.-fmin)*xp1(ida)+fmin*xp2(ida)
            if(zida.gt.xca)imin=1.
         endif
      endif

c Now the minimum fraction is in fmin, which is the crossing point.
c imin contains the face-index of this crossing. -1,0, or +1.
c Calculate normalized intersection coordinates.
      do i=1,npdim
         x12(i)=((1.-fmin)*xp1(i)+fmin*xp2(i)
     $        -obj_geom(ocenter+i-1,iobj))/obj_geom(oradius+i-1,iobj)
      enddo
c Calculate r,theta,z (normalized) relative to the ida direction as z.
      z=x12(ida)
      theta=atan2(x12(mod(ida+1,npdim)+1),x12(mod(ida,npdim)+1))
      r2=0.
      do i=1,npdim-1
         k=mod(ida+i-1,npdim)+1
         r2=r2+x12(k)**2
c         write(*,'(i2,2f8.4,'', '',$)')k,r2,x12(k)
      enddo
c      write(*,'(a,7f7.4,3i3)')'r2,theta,z,x12,fmin,imin'
c     $     ,r2,theta,z,x12,fmin,imin
c End blocks are of size nr x nt, and the curved is nt x nz.
c 3-D only here. 
      infobj=nf_map(iobj)
      ijbin=0.
      ir=9
      it=9
      iz=9
      if(imin.ne.0)then
c Ends
         if(imin.eq.1)then
c offset by (nr+nz)*nt
            ijbin=(nf_dimlens(nf_flux,infobj,1)
     $           +nf_dimlens(nf_flux,infobj,3))
     $           *nf_dimlens(nf_flux,infobj,2)
         endif
c Uniform mesh in r^2 normalized. 
         ir=int(nf_dimlens(nf_flux,infobj,1)*(0.999999*r2))
         it=int(nf_dimlens(nf_flux,infobj,2)
     $     *(theta/3.1415927+1.)*0.5)
         ijbin=ijbin+ir+it*nf_dimlens(nf_flux,infobj,1)
      else
c Side. Offset to this facet nr*nt:
         ijbin=nf_dimlens(nf_flux,infobj,1)*nf_dimlens(nf_flux,infobj,2)
         it=int(nf_dimlens(nf_flux,infobj,2)
     $     *(theta/3.1415927+1.)*0.5)
         iz=int(nf_dimlens(nf_flux,infobj,3)*(0.999999*z+1.)*0.5)
c Index in order theta,z
         ijbin=ijbin+it+iz*nf_dimlens(nf_flux,infobj,2)
      endif
c      write(*,'(6f8.4,3i3)')xp1,xp2,ir,it,iz
      end
