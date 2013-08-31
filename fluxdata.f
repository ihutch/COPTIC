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
               colnforce(i,j,k)=0.
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
            mf_quant(mf_obj)=int(obj_geom(ofluxtype,i))
            itype=int(obj_geom(otype,i))
            i2type=(itype/256)
            if(i2type.eq.2)then
c Point object. No flux accumulations.
               if(myid.eq.0)write(*,*
     $              )'Measuring force, no flux, on point object',i
     $              ,mf_obj
               mf_quant(mf_obj)=0
            endif
c The mapped object number != object number.
            nf_map(i)=mf_obj
            nf_geommap(mf_obj)=i
c Use only bottom 8 bits from now on:
            itype=itype-256*i2type
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
                     nfluxes=nfluxes*int(obj_geom(ofn1+k-1,i))
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
     $                    int(obj_geom(ofn1+mod(k,ns_ndims),i))
     $                    *int(obj_geom(ofn1+mod(k+1,ns_ndims),i))
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
            itype=int(obj_geom(otype,i))
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
                  ica=int(obj_geom(ocylaxis,i))
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
      subroutine tallyexit(i,idiffreg,ltlyerr,dtpos)
c Document the exit of this particle just happened. 
c Assign the exit to a specific object, and bin on object.
c (If it is a mapped object, decided by objsect.)      
c On entry
c        i is particle number, 
c        idiffreg is the xor between its region and previous.
c        dtpos is the time step duration that brought us here.
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
         call objsect(i,iobj,ierr,dtpos)
         if(ierr.gt.0)then
            ireg=insideall(npdim,x_part(1,i))
            r=0.
            r1=0.
            do id=1,3
               r=r+x_part(id,i)**2
               r1=r1+(x_part(id,i)-dtpos*x_part(id+3,i))**2
            enddo
            r=sqrt(r)
            r1=sqrt(r1)
            write(*,*)'Tallyexit error',ierr,i,iobj,idiffreg
            if(ierr.eq.99)write(*,*)'Unknown object type.'
            write(*,*)'xpart,r=',(x_part(k,i),k=1,6),r,ireg
            write(*,*)'xp1',(x_part(k,i)-dtpos*x_part(k+3,i),k=1,3),r1
            ltlyerr=.true.
            return
         elseif(ierr.eq.-1)then
c This Pass-through should not happen because tallyexit is not called
c unless the region is changed. But leave for error detection.
            write(*,*)'idiffreg,i,iobj=',idiffreg,i,iobj
            write(*,'(a,6f10.5)')'xp2',(x_part(k,i),k=1,6)
            write(*,'(a,6f10.5)')'xp1',(x_part(k,i)-dtpos*x_part(k+3,i)
     $           ,k=1,3)
         endif
      endif
      idp=idp/2
      goto 1
      
      end
c******************************************************************
c****************************************************************
      subroutine objsect(j,iobj,ierr,dtpos)
c Find the intersection of the last step of particle j (length dtpos)
c with object iobj, and update the positioned-fluxes accordingly.  ierr
c is returned: 0 good. 1 no intersection. 99 unknown object.
c
      include '3dcom.f'
      include 'partcom.f'
      include 'sectcom.f'

      real x1(npdim),x2(npdim)
      integer isc
      data isc/0/

      ierr=0
      ijbin2=-1

c Do nothing for untracked objects and report no error.
      if(nf_map(iobj).eq.0)return
      infobj=nf_map(iobj)
      if(mf_quant(infobj).eq.0)return

      itype=int(obj_geom(otype,iobj))
c Use only bottom 8 bits:
      itype=itype-256*(itype/256)

c Get the positions:
      do i=1,npdim
         x1(i)=x_part(i,j)-dtpos*x_part(i+3,j)
         x2(i)=x_part(i,j)
      enddo
      if(itype.eq.1)then
c Sphere intersection. Return the bin number and direction ijbin,sd.
         call spherefsect(npdim,x1,x2,iobj,ijbin,sd,fraction,ijbin2)
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
c Do the bin adding in a subroutine.
      call binadding(j,infobj,sd,ijbin)

      if(ijbin2.ne.-1)then
c This should not happen.
         write(*,*)'Pass-through j,ijbin2,infobj=',j,ijbin2,infobj
         ierr=-1
      endif
      end
c*******************************************************************
      subroutine binadding(j,infobj,sd,ijbin)
c Add particle-j data to infobj bin ijbin with crossing-direction sd.
      implicit none
      include '3dcom.f'
      include 'partcom.f'
      integer j,infobj,ijbin
      real sd
      
      integer iaddress,k,id
      real xx

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
      
      end
c*******************************************************************
c*******************************************************************
      subroutine pllelfrac(xp,xn,iobj)
c For input point xp(ndims) return the normalized position relative to
c the parallelogram object objg in output xn(ndims)
      include '3dcom.f'
      real xp(ns_ndims),xn(ns_ndims)
      integer iobj
c      real objg(odata)
      do j=1,pp_ndims
         xn(j)=0.
c Contravariant projections.
         do i=1,ns_ndims
c Cartesian coordinates.
            ii=(ocenter+i-1)
c            xc=objg(ii)
            xc=obj_geom(ii,iobj)
c xn1, xn2 are the contravariant coordinates with respect to the center.
            ji=(pp_contra+pp_ndims*(j-1)+i-1)
c            write(*,*)'ji',ji
            xn(j)=xn(j)+(xp(i)-xc)*obj_geom(ji,iobj)
         enddo
      enddo
      end
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
      real function fluxdiag()
c Get the total count to object 1 and convert to normalized flux.
c Assuming it's a unit sphere.
      include '3dcom.f'
c For rhoinf, dt
      include 'partcom.f'

      sum=0
      do i=1,nf_posno(nf_flux,1)
         sum=sum+ff_data(nf_address(nf_flux,1,nf_step)+i-1)
      enddo
      fluxdiag=sum/(4.*3.14159)/rhoinf/dt
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
      subroutine fluxave(n1in,n2in,ifobj,iquant,rinf)
      integer n1in,n2in
c      logical lplot
      integer iquant
      include '3dcom.f'
      include 'sectcom.f'
c Use for averaging:ff_data(nf_address(iq,ifobj,nf_maxsteps+1)+i-1)
c which is ff_data(iav+i)
      real fluxofstep(nf_maxsteps),step(nf_maxsteps)
      character*30 string

      n1=n1in
      n2=n2in
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
c Initialize the name before entry. Very Important!
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
      include 'colncom.f'
      character*(100) charout

c Construct a filename that contains many parameters
c Using the routines in strings_names.f
      call nameconstruct(name)
c     np=nbcat(name,'.flx')
      call nbcat(name,'.flx')
c      write(*,*)name
      write(charout,51)debyelen,Ti,vd,rs,phip
 51   format('debyelen,Ti,vd,rs,phip:',5f10.4,' Version: 3')

c      write(*,*)'mf_obj=',mf_obj,nf_step,mf_quant(1)

      open(22,file=name,status='unknown',err=101)
      close(22,status='delete')
      open(22,file=name,status='new',form='unformatted',err=101)
c This write sequence must be exactly that read below.
      write(22)charout
      write(22)debyelen,Ti,vd,rs,phip ,colntime,subcycle,vneutral
     $     ,fcollided,dropaccel,Tneutral,Eneutral
      write(22)nf_step,mf_quant,mf_obj,(nf_geommap(j),j=1,mf_obj)
c      write(*,*)'geommap',(nf_geommap(j),j=1,mf_obj)
      write(22)(ff_rho(k),k=1,nf_step)
      write(22)(ff_dt(k),k=1,nf_step)
      write(22)((nf_posno(i,j),(nf_dimlens(i,j,k),k=1,nf_ndims)
     $     ,(nf_faceind(i,j,k),k=1,2*nf_ndims)
     $     ,i=1,mf_quant(j)),j=1,mf_obj)
      write(22)(((nf_address(i,j,k),i=1,mf_quant(j)),j=1,mf_obj),
     $     k=1-nf_posdim,nf_step+2)
      ndatalen=nf_address(1,1,nf_step+2)-1
      write(22)(ff_data(i),i=1,ndatalen)
c New force write.
      write(22)(((fieldforce(i,j,k),pressforce(i,j,k) ,partforce(i,j,k)
     $     ,colnforce(i,j,k),i=1,ns_ndims)
     $     ,charge_ns(j,k),j=1,mf_obj),k=1,nf_step)
c Object data:
      write(22)ngeomobj
      write(22)((obj_geom(j,k),j=1,odata),nf_map(k),k=1,ngeomobj)
      write(22)ibool_part,ifield_mask,iptch_mask,lboundp,rjscheme
c Intersection data:
      write(22)sc_ipt
      write(22)(((x_sc(j,i,k),j=1,sc_ndims),i=1,2),iob_sc(k),
     $     ibin_sc(k),k=1,sc_ipt)
c n_part data
      write(22)(nf_npart(k),k=1,nf_step)
      close(22)
c      write(*,*)'Wrote flux data to ',name(1:lentrim(name))
      do k=1,nf_step
         do j=1,mf_obj
            do i=1,ns_ndims
               if(.not.colnforce(i,j,k).lt.1.e30)then
                  write(*,*)'Strange colnforce',i,j,k,colnforce(i,j,k)
               endif
            enddo
         enddo
      enddo

      return
 101  continue
      write(*,*)'Error opening file:',name(1:lentrim(name))
      close(22,status='delete')
      end
c*****************************************************************
      subroutine readfluxfile(name,ierr)
c On entry ierr .ne.0 indicates write informational messages.
c On exit  ierr .ne.0 indicates error.
      character*(*) name
      include '3dcom.f'
      include 'plascom.f'
      include 'sectcom.f'
      include 'colncom.f'
      character*(100) charout

      Eneutral=0.
      open(23,file=name,status='old',form='unformatted',err=101)
      read(23)charout
c Figure out the version:
      iend=lentrim(charout)
      if(charout(iend-9:iend-3).eq.'Version')then
         read(charout(iend-1:iend),*)iversion
         if(ierr.ne.0)write(*,*)'Flux file version',iversion
         ierr=0
      else
         iversion=0
      endif
      if(iversion.le.1)then
         read(23)debyelen,Ti,vd,rs,phip
      elseif(iversion.le.2)then
         read(23)debyelen,Ti,vd,rs,phip
     $        ,colntime,subcycle,vneutral,fcollided,dropaccel,Tneutral
      elseif(iversion.le.3)then
         read(23)debyelen,Ti,vd,rs,phip ,colntime,subcycle,vneutral
     $        ,fcollided,dropaccel,Tneutral,Eneutral
      endif
      read(23)nf_step,mf_quant,mf_obj,(nf_geommap(j),j=1,mf_obj)
c      write(*,*)'geommap',(nf_geommap(j),j=1,mf_obj)
      read(23)(ff_rho(k),k=1,nf_step)
      read(23)(ff_dt(k),k=1,nf_step)
      read(23)((nf_posno(i,j),(nf_dimlens(i,j,k),k=1,nf_ndims)
     $     ,(nf_faceind(i,j,k),k=1,2*nf_ndims)
     $     ,i=1,mf_quant(j)),j=1,mf_obj)
      read(23)(((nf_address(i,j,k),i=1,mf_quant(j)),j=1,mf_obj),
     $     k=1-nf_posdim,nf_step+2)
c      read(23)(ff_data(i),i=1,nf_address(1,1,nf_step+2)-1)
      ndatalen=0
      do j=1,mf_obj
         if(mf_quant(j).ge.1)then
            ndatalen=nf_address(1,j,nf_step+2)-1
            goto 201
         endif
      enddo
 201  continue
c      write(*,*)'Datalen',ndatalen
      read(23)(ff_data(i),i=1,ndatalen)
      if(iversion.eq.0)then
         write(*,*)'Old force data version',iversion
         read(23,end=102, err=102)(((fieldforce(i,j,k),pressforce(i,j,k)
     $     ,partforce(i,j,k)
     $     ,charge_ns(j,k),i=1,ns_ndims),j=1,mf_obj),k=1,nf_step)
      else
         read(23,end=102, err=102)(((fieldforce(i,j,k),pressforce(i,j,k)
     $     ,partforce(i,j,k)
     $     ,colnforce(i,j,k)
     $     ,i=1,ns_ndims),charge_ns(j,k),j=1,mf_obj),k=1,nf_step)
      endif
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
c n_part data
      if(iversion.ge.2)read(23,end=104)(nf_npart(k),k=1,nf_step)
    
      goto 103
 102  write(*,*)'Failed to read back forces. Old format? Version='
     $     ,iversion
      ierr=2
      goto 103
 104  write(*,*)'No nf_npart data. Mismatch of versions? Version='
     $     ,iversion
 103  close(23)

c Hack to fix nans when colnforce was wrong. Delete when that data is 
c obsolete.
      isc=0
      do k=1,nf_step
         do j=1,mf_obj
            do i=1,ns_ndims
               if(.not.colnforce(i,j,k).lt.1.e30)then
                  if(isc.eq.0)write(*,*)'Strange colnforce',i,j,k
     $                 ,colnforce(i,j,k)
                  colnforce(i,j,k)=0.
                  isc=isc+1
               endif
            enddo
         enddo
      enddo

c Now one might have to reconstruct nf_faceind from nf_dimlens.
      if(ierr.ne.0)write(*,*)'Read back flux data from '
     $     ,name(1:lentrim(name))
c      write(*,*)charout(1:lentrim(charout))
c      ierr=0
      return
 101  write(*,*)'Error opening file:',name(1:lentrim(name))
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
            imin=0
         else
c axial crossing
            imin=-1
            zida=(1.-fmin)*xp1(ida)+fmin*xp2(ida)
            if(zida.gt.xca)imin=1
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
      ijbin=0
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
