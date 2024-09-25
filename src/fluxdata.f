! Initialize the flux data, determining what we are saving and where.
! The objects whose flux is to be tracked are indicated by 
! obj_geom(ofluxtype,i). If this is zero, it is not tracked.
! The number of ibins in each (2) of the surface dimensions is indicated
! by obj_geom(ofn1[/2],i), and data space and addresses allocated.
! The uniform-spacing bin-positions are calculated and set.
      subroutine fluxdatainit(myid)
      include 'ndimsdecl.f'
      include '3dcom.f'
      include 'sectcom.f'
!-----------------------------------------------
! Intersection diagnostic points.
      sc_ipt=0
! Default working species 1
      if_species=1
!-----------------------------------------------
! Initialize here to avoid giant block data program.
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
!-------------------------------------------------
! This initialization ought not to be necessary except for partforce.
      do k=1,nf_maxsteps
         do j=1,nf_obj
            do i=1,ndims
               fieldforce(i,j,k)=0.
               pressforce(i,j,k)=0.
               partforce(i,j,k)=0.
               colnforce(i,j,k)=0.
            enddo
            charge_ns(j,k)=0.
         enddo
      enddo
!-------------------------------------------------
! Initialize object number
      mf_obj=0
      do i=1,ngeomobj
         if(obj_geom(ofluxtype,i).eq.0)then
! No flux setting for this object.
         elseif(obj_geom(ofluxtype,i).ge.1
     $           .and. obj_geom(ofluxtype,i).le.nf_quant)then
! Might eventually need more interpretation of fluxtype.
            if(mf_obj.lt.nf_obj)then
               mf_obj=mf_obj+1
            else
               write(*,*)'WARNING: Flux object number exceeds nf_obj='
     $              ,nf_obj,' Omitting later objects.'
               goto 1
            endif
! Here we ought perhaps to have a way to set the number of fluxes
! different for different species. For now, they are the same.
            mf_quant(mf_obj)=0
            nq=int(obj_geom(ofluxtype,i))
            if_quant(mf_obj,1)=0
! Increment number of quantities provided there is allocated address room.
            do j=1,nf_species
               nt=mf_quant(mf_obj)+nq
               if(nt.le.nf_quant)then
                  kf_quant(mf_obj,j)=nq
                  mf_quant(mf_obj)=nt
               else
                  kf_quant(mf_obj,j)=0
               endif
               if(j.gt.1)if_quant(mf_obj,j)=if_quant(mf_obj,j-1)
     $              +kf_quant(mf_obj,j)
!               write(*,'(a,7i4)')'Species,if,kf,mf_quant',j
!     $              ,if_quant(mf_obj,j),kf_quant(mf_obj,j)
!     $              ,mf_quant(mf_obj)
            enddo
            itype=int(obj_geom(otype,i))
            i2type=(itype/256)
            if(i2type.eq.2)then
! Point object. No flux accumulations.
               if(myid.eq.0)write(*,*
     $              )'Measuring force, no flux, on point object',i
     $              ,mf_obj
               mf_quant(mf_obj)=0
            endif
! The mapped object number != object number.
            nf_map(i)=mf_obj
            nf_geommap(mf_obj)=i
! Use only bottom 8 bits from now on:
            itype=itype-256*i2type
! There are nfluxes positions for each quantity.
! At present there's no way to prescribe different grids for each
! quantity in the input file. But the data structures could 
! accommodate such difference if necessary. 
            do j=1,mf_quant(mf_obj)
               if(itype.eq.1)then
! Sphere one face only n3=1
                  nfluxes=1
                  do k=1,2
                     nf_dimlens(j,mf_obj,k)=int(obj_geom(ofn1+k-1,i))
                     nfluxes=nfluxes*int(obj_geom(ofn1+k-1,i))
                  enddo
               elseif(itype.eq.2 .or. itype.eq.4)then
! Cuboid or parallelopiped.
! Six faces, each using two of the three array lengths. 
                  nfluxes=0
                  do kk=1,int(obj_geom(offc,i))
                     k=mod(kk-1,ndims)+1
                     if(kk.le.ndims)
     $                 nf_dimlens(j,mf_obj,k)=int(obj_geom(ofn1+k-1,i))
                     nf_faceind(j,mf_obj,kk)=nfluxes
                     nfluxes=nfluxes+
     $                    int(obj_geom(ofn1+mod(k,ndims),i))
     $                    *int(obj_geom(ofn1+mod(k+1,ndims),i))
!                     write(*,*)k,' nfluxes=',nfluxes
!     $                    ,obj_geom(ofn1+mod(k,ndims),i)
                  enddo
!                  write(*,*)(nf_faceind(j,mf_obj,k),k=1,2*ndims)
!                  nfluxes=2*nfluxes
               elseif(itype.eq.3 .or. itype.eq.5)then
! Cylinder specifying nr, nt, nz. 
! Three facets in the order bottom side top. Check for enough n's
                  if(obj_geom(ofn3,i).eq.0.)then
                     write(*,*)'Zero third flux array length for object'
     $                    ,j,'Probably an error in object file spec'
                     stop
                  endif
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
!                  write(*,*)nfluxes,(nf_faceind(j,mf_obj,k),k=1,3)
               elseif(itype.eq.6.or.itype.eq.7)then
! Surface of revolution has theta (ofn1) 
! and each face (up to ovlen-1) is a line segment rotated.
! onpair-1 is the number of faces. opdiv the divisions of each.
!                  write(*,*)'Surface of Revolution flux initialization'
                  nfluxes=0
                  ntheta=int(obj_geom(ofn1,i))
                  nf_dimlens(j,mf_obj,1)=ntheta
                  nr=0
                  do kk=1,int(obj_geom(onpair,i)-1)
! faceind is the offset of each face.
                     nf_faceind(j,mf_obj,kk)=nfluxes
                     nr=nr+int(obj_geom(opdiv+kk-1,i))
                     nfluxes=nfluxes+int(obj_geom(opdiv+kk-1,i))*ntheta
                  enddo
! dimlens of r is the total number of r-divisions of all facets.
                  nf_dimlens(j,mf_obj,2)=nr
!                  stop
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
!      write(*,*)'nf_posno=',((nf_posno(j,i),j=1,2),i=1,4)
!-------------------------------------------------
! Now we create the addressing arrays etc.
      call nfaddressinit()
! After which, nf_address(i,j,k) points to the start of data for 
! quantity i, object j, step k. 
! So we can pass nf_data(nf_address(i,j,k)) as a vector start.
!-------------------------------------------------
! Initialize the position and area data for each object.
      call positioninit()

      end
!******************************************************************
      subroutine nfaddressinit()
      include 'ndimsdecl.f'
      include '3dcom.f'
! General iteration given correct settings of nf_posno. Don't change!
! Zero to silence incorrect warnings.
      numdata=0
      numobj=0
      i=2
      j=2
      nf_address(1,1,1-nf_posdim)=1
      do k=1-nf_posdim,nf_nsteps+2
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
!               write(*,*)i,j,k,nf_posno(i,j),nf_address(i,j,k),numdata
            enddo
            numobj=numobj+numdata
         enddo
      enddo
      nend=nf_address(max(i-1,1),max(j-1,1),k-1)+numobj
! Check if we might overrun the datasize.
!      write(*,*)'Flux data',nend,' for',mf_obj,' objects,'
!     $        ,nf_nsteps,' steps,',numdata,' quantities*positions'
      if(nend.gt.nf_datasize)then
         write(*,*)'DANGER: flux data',nend,' for',mf_obj,' objects,'
     $        ,nf_nsteps,' steps,',numdata,' quantities*positions'
     $        ,' would overrun nf_datasize',nf_datasize
         write(*,*)'Use fewer quantities, positions, (max)steps,'
     $        ,' or increase nf_datasize'
         stop
      else
      endif
      end
!******************************************************************
      subroutine positioninit()
      include 'ndimsdecl.f'
      include '3dcom.f'

!      real xyz(ndims)
! The k=1-nf_posdim to k=0 slots exist for us to put descriptive
! information: values that provide positions to correspond to the
! fluxes.  Area is most important and is in 1-nf_posdim.
!
! Ellipses require elliptic surfaces which are a mess. For equal radii,
! the area element is just dA=2\pi a^2 dcostheta. 
      do i=1,ngeomobj
         if(obj_geom(ofluxtype,i).gt.0)then
            itype=int(obj_geom(otype,i))
            itype=itype-256*(itype/256)
            io=nf_map(i)
            if(itype.eq.1)then
! Sphere -----------------------
               if(obj_geom(oradius,i).ne.obj_geom(oradius+1,i) .or.
     $              obj_geom(oradius,i).ne.obj_geom(oradius+2,i))then
                  write(*,*)'Warning! Non-spherical spheroid'
     $                 ,' surface areas'
     $                 ,' not calculated correctly.'
               endif
               do j=1,mf_quant(io)
! Area of each element. They are equal only if radii are equal. 
                  ar=4.*3.1415926*obj_geom(oradius,i)**2
     $                 /(nf_dimlens(j,io,1)*nf_dimlens(j,io,2))
!                  write(*,*)'Object',i,' Quant',j,' Facet areas=',ar
                  do i2=1,nf_dimlens(j,io,2)
                     p=3.1415926*(-1.+2.*(i2-0.5)/nf_dimlens(j,io,2))
                     do i1=1,nf_dimlens(j,io,1)
                        c=-1.+2.*(i1-0.5)/nf_dimlens(j,io,1)
                        ip=i1+(i2-1)*int(nf_dimlens(j,io,1))
! Position values are cos(theta) and psi. Third not used.
!                        write(*,*)j,io,i1,i2,nf_p1,ioff,ip,c
!     $                       ,nf_address(j,io,nf_p1)+ip-1
                        ff_data(nf_address(j,io,nf_p1)+ip-1)=c
                        ff_data(nf_address(j,io,nf_p2)+ip-1)=p
                        ff_data(nf_address(j,io,nf_pa)+ip-1)=ar
                     enddo
                  enddo
               enddo
!            write(*,*)'Set ff_data',i,ioff,nf_map(i),io
! End of sphere case.
            elseif(itype.eq.2 .or. itype.eq.4)then
! Parallelopiped/Cube ----------------------------
! facet areas are face areas divided by no of facets.
               do j=1,mf_quant(io)
! Over different faces:
                  do k=1,2*ndims
                     k1=mod(k-1,ndims)+1
                     k2=mod(k  ,ndims)+1
                     k3=mod(k+1,ndims)+1
!                     write(*,'(/,a,i2,3i3,$)')'Face',k
!     $                    ,nf_faceind(j,io,k)
!     $                    ,nf_dimlens(j,io,k2),nf_dimlens(j,io,k3)
                     if(itype.eq.4)then
! Structure position of start of covariant vectors
                        i1=ovec+ndims*(k2-1)
                        i2=ovec+ndims*(k3-1)
! face area equals magnitude of cross product. 
! 3d assumption for convenience.
                        ar2= (obj_geom(i1+1,i)*obj_geom(i2+2,i)
     $                    -obj_geom(i1+2,i)*obj_geom(i2+1,i))**2
     $                    +(obj_geom(i1+2,i)*obj_geom(i2  ,i)
     $                    -obj_geom(i1  ,i)*obj_geom(i2+2,i))**2
     $                    +(obj_geom(i1  ,i)*obj_geom(i2+1,i)
     $                    -obj_geom(i1+1,i)*obj_geom(i2  ,i))**2
                        ar=4.*sqrt(ar2)/
     $                    (nf_dimlens(j,io,k2)*nf_dimlens(j,io,k3))
                     else
! Cube
                        ar=4.*obj_geom(oradius+k2-1,i)
     $                       *obj_geom(oradius+k3-1,i)/
     $                    (nf_dimlens(j,io,k2)*nf_dimlens(j,io,k3))
                     endif
! Store data. At the moment, just the areas.
                     do j3=1,nf_dimlens(j,io,k3)
                     do j2=1,nf_dimlens(j,io,k2)
                        ip=j2+(j3-1)*int(nf_dimlens(j,io,k2))
     $                       +nf_faceind(j,io,k)
                        ff_data(nf_address(j,io,nf_pa)+ip-1)=ar
!                        write(*,'(i4,f8.4,$)')
!     $                       j2,j3
!     $                       ,int(nf_dimlens(j,io,k2))
!     $                       ,ip,ar
                        if(itype.eq.2)then
! Cube position data x,y,z
                           xr1=-obj_geom(oradius+k1-1,i)
                           if(k.gt.ndims)xr1=-xr1
                           xr2=obj_geom(oradius+k2-1,i)*
     $                          (-1.+2.*(j2-0.5)/nf_dimlens(j,io,k2))
                           xr3=obj_geom(oradius+k3-1,i)*
     $                          (-1.+2.*(j3-0.5)/nf_dimlens(j,io,k3))
!                           write(*,*)j2,j3,xr1,xr2,xr3
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
! Cylinder and Non-aligned cylinder -------------------------------
               if(itype.eq.3)then
                  ica=int(obj_geom(ocylaxis,i))
                  rc=obj_geom(oradius+mod(ica,ndims),i)
                  zr=obj_geom(oradius+ica-1,i)
                  zc=obj_geom(ocenter+ica-1,i)
                  if(rc.ne.obj_geom(oradius+mod(ica+1,ndims),i))then
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
! Area of each element. Faces are in order bottom, side, top.
! Bottom and top facet areas:
                  ar=3.1415926*rc**2
     $                 /(nf_dimlens(j,io,1)*nf_dimlens(j,io,2))
!                  write(*,*)'Object',i,' Quant',j,' Facet areas=',ar
                  do i2=1,nf_dimlens(j,io,2)
                     do i1=1,nf_dimlens(j,io,1)
                        ip=i1+(i2-1)*int(nf_dimlens(j,io,1))
! Positional data:
! r (not r^2)
                        r=rc*sqrt((i1-0.5)/nf_dimlens(j,io,1))
                        ff_data(nf_address(j,io,nf_pr)+ip-1)=r
                        ff_data(nf_address(j,io,nf_pr)+ip-1
     $                       +nf_faceind(j,io,3))=r
! theta
                        t=3.1415927*
     $                       (-1.+2.*(i2-0.5)/nf_dimlens(j,io,2))
!                        write(*,*)i1,nf_dimlens(j,io,1),t
                        ff_data(nf_address(j,io,nf_pt)+ip-1)=t
                        ff_data(nf_address(j,io,nf_pt)+ip-1
     $                       +nf_faceind(j,io,3))=t
! z
                        ff_data(nf_address(j,io,nf_pz)+ip-1)=zc-zr
                        ff_data(nf_address(j,io,nf_pz)+ip-1
     $                       +nf_faceind(j,io,3))=zc+zr
! area
                        ff_data(nf_address(j,io,nf_pa)+ip-1)=ar
!                        write(*,*)'Area address',nf_address(j,io,nf_pa)
!     $                       +ip-1,' Area',ar
                        ff_data(nf_address(j,io,nf_pa)+ip-1
     $                       +nf_faceind(j,io,3))=ar
                     enddo
                  enddo
! Side 2 pi r 2 z:
                  ar=4.*3.1415926*rc*zr
     $                 /(nf_dimlens(j,io,2)*nf_dimlens(j,io,3))
! index theta,z
                  do i2=1,nf_dimlens(j,io,3)
                     do i1=1,nf_dimlens(j,io,2)
                        ip=i1+(i2-1)*int(nf_dimlens(j,io,2))
                        ff_data(nf_address(j,io,nf_pa)+ip-1
     $                       +nf_faceind(j,io,2))=ar
                        ff_data(nf_address(j,io,nf_pr)+ip-1
     $                       +nf_faceind(j,io,2))=rc
!                        write(*,'(i4,f8.4,$)')
!     $                       ,ip,ar
! theta
                        t=3.1415927*
     $                       (-1.+2.*(i1-0.5)/nf_dimlens(j,io,2))
!                        write(*,*)i1,nf_dimlens(j,io,1),t
                        ff_data(nf_address(j,io,nf_pt)+ip-1
     $                       +nf_faceind(j,io,2))=t
! z
                        z=zc+zr*(-1.+2.*(i2-0.5)/nf_dimlens(j,io,3))
                        ff_data(nf_address(j,io,nf_pz)+ip-1
     $                       +nf_faceind(j,io,2))=z
                     enddo
                  enddo
               enddo
            elseif(itype.eq.6 .or. itype.eq.7)then
! Surface of revolution ------------------------------
!               write(*,*)'Flux positions for Surface of Revolution',i
! Each line segment is a truncated cone with ofn1 and opdiv
! equal divisions. The areas of the subsegments are to be equal. 
! For a segment with ends (rb,zb),(rt,zt), the total area is 
! A=\pi (rt+rb)\sqrt{(rt-rb)^2+(zt-zb)^2}. Equal area means
! equal division in r^2, so when dividing into m=0,1,...,N positions
! we should take rm^2=rb^2+(rt^2-rb^2)*m/N 
! and corresponding distance fraction fm=(rm-rb)/(rt-rb).
! The total array of points lying at the end of subsegments can be
! considered to be 0,((f_mn,m=1,N_n),n=1,npair-1). The m=0 point
! of all but the first segment is the m=N_n of the previous.
! The total array of area centroids is ((f_(m-1/2)n,m=1,N),n=1,npair-1)
! These fractions are the positions we must store. We must also store
! the areas, which are just A/N_n/ntheta.
! Do in usual order even though there can't be different arrangements
! for different flux types.
               do j=1,mf_quant(io)

! The covariant vectors' length defines the scale factor in the 
! z and r directions, from object-normalized to world.
! The rscale is already stored.
                  zscale=0.
                  do kk=1,ndims
                     zscale=zscale+obj_geom(ovec+2*ndims+kk-1,i)**2
                  enddo
                  zscale=sqrt(zscale)
                  rscale=obj_geom(orscale,i)
!                  write(*,*)'rscale,zscale,=',rscale,zscale
! Contour positions:
                  i0=0
                  do i3=1,int(obj_geom(onpair,i)-1)
                  rb=obj_geom(opr+i3-1,i)*rscale
                  rt=obj_geom(opr+i3,i)*rscale
                  rdiff=rt-rb
                  zdiff=(obj_geom(opz+i3,i)-obj_geom(opz+i3-1,i))*zscale
! Areas are equal through segments. Needs to be scaled.
                  area=3.1415927*(rb+rt)*sqrt(rdiff**2+zdiff**2)
     $                 /obj_geom(opdiv+i3-1,i)
!                  write(*,*)'area,rb,rt,zdiff,rscale,zscale',area,rb,rt
!     $                 ,zdiff,rscale,zscale
                  do i2=1,int(obj_geom(opdiv+i3-1,i))
                     i0=i0+1
! Calculate contour fractional index:
! for end of subsegment
                     fr=float(i2)/obj_geom(opdiv+i3-1,i)
                     if(abs(rdiff).gt.1.e-5*rb)then
                        rm2=rb**2*(1-fr)+rt**2*fr
                        fm=(sqrt(rm2)-rb)/(rt-rb)
                     else
                        fm=fr
                     endif
                     pm=fm+i3
! for centroid of subsegment
                     fr=(i2-0.5)/obj_geom(opdiv+i3-1,i)
                     if(abs(rdiff).gt.1.e-5*rb)then
                        rm2=rb**2*(1-fr)+rt**2*fr
                        fm=(sqrt(rm2)-rb)/(rt-rb)
                     else
                        fm=fr
                     endif
                     cm=fm+i3
! Theta positions:
                     do i1=1,nf_dimlens(j,io,1)
                        ip=i1+(i0-1)*int(nf_dimlens(j,io,1))
! Centroid frac index, theta, end frac, area
                        t=3.1415927*
     $                       (-1.+2.*(i1-0.5)/nf_dimlens(j,io,1))
!                        write(*,*)'cm,t,pm,area',cm,t,pm,area
                        ff_data(nf_address(j,io,nf_p1)+ip-1)=cm
                        ff_data(nf_address(j,io,nf_p2)+ip-1)=t
                        ff_data(nf_address(j,io,nf_p3)+ip-1)=pm
                        ff_data(nf_address(j,io,nf_p4)+ip-1)=area
                     enddo
                  enddo
                  enddo
               enddo
            else
! Unknown ------------------------------------------
               write(*,*)'Flux positions uncalculated for object'
     $           ,i,' type',obj_geom(otype,i)
            endif
         endif
!      write(*,*)'Initialized positional data to:',nf_address(1,1,1)-1
!      write(*,'(10f8.4)')(ff_data(k),k=1,nf_address(1,1,1)-1)
      enddo

!      write(*,*)'10 steps of initialized slots:',nf_address(1,1,6)-1
!      write(*,'(10f8.4)')(ff_data(k),k=nf_address(1,1,1)
!     $     ,nf_address(1,1,11)-1)

      end
!******************************************************************
      subroutine tallyexit(xprior,xi,idiffreg,ltlyerr,ispecies,fmax)
! Document the exit of this particle just happened. 
! Assign the exit to a specific object, and bin on object.
! (If it is a mapped object, decided by objsect.)      
! On entry
!        xprior is particle prior position
!        xi is particle position/velocity etc, 
!        idiffreg is all the objects the step crosses.
!        fmax is the maximum considered step fraction.
! On exit ltlyerr is true if an error occurred else unchanged.
! Normally, returning an error will cause this particle to be
! considered to have left the particle region, so it will be discarded.
      integer idiffreg
      logical ltlyerr
      include 'ndimsdecl.f'
      real xi(3*ndims)
      real xprior(ndims)
      include 'partcom.f'
      include '3dcom.f'

!      if(ispecies.ge.2)write(*,*)'Tallying',ispecies

      idiff=abs(idiffreg)
      idp=idiff
! Determine (all) the objects crossed and call objsect for each.
      iobj=0
 1    if(idiff.eq.0) return
      iobj=iobj+1
      idiff=idiff/2
      if(idp.ne.idiff*2)then
         call objsect(xprior,xi,iobj,ierr,ispecies,fmax)
         if(ierr.gt.0)then
! There's a serious error. Give info and stop
            ireg=insideall(ndims,xi(1))
            r=0.
            r1=0.
            do id=1,3
               r=r+xi(id)**2
            enddo
            r=sqrt(r)
            r1=sqrt(r1)
            write(*,*)'Tallyexit error',ierr,iobj,idiffreg
            if(ierr.eq.99)write(*,*)'Unknown object type.'
            write(*,*)'xpart,r=',(xi(k),k=1,6),r,ireg
            write(*,*)'xp1',(xprior(k),k=1,3),r1
            ltlyerr=.true.
            stop
            return
         endif
      endif
      idp=idp/2
      goto 1
      
      end
!******************************************************************
!****************************************************************
      subroutine objsect(x1,xi,iobj,ierr,ispecies,fmax)
! Find the intersection of the step of particle xi from x1
! with object iobj, and update the positioned-fluxes accordingly.
! fmax is the maximum position fraction that should be counted.
! It prevents counting intersections after leaving the region.
!  ierr is returned: 0 good. 1 no intersection. 99 unknown object.
!
      include 'ndimsdecl.f'
      real xi(3*ndims)
      real x1(ndims)
      include '3dcom.f'
      include 'partcom.f'
      include 'sectcom.f'
      include 'dbgcom.f'

      real xn1(ndims),xn2(ndims)
      integer isc
      integer ijbin
      integer imin(ovlen)
      real fraction
      real fmin(ovlen)
! This does not seem to work as expected:
      equivalence (fraction,fmin(1))
      

      data isc/0/

      idebug=idbug
      ierr=0
      ijbin=0 ! silence warnings.
! Do nothing for untracked objects and report no error.
      infobj=nf_map(iobj)
      if(infobj.eq.0)return
      if(mf_quant(infobj).eq.0)return

      ijbin2=-1
      itype=int(obj_geom(otype,iobj))
! Use only bottom 8 bits:
      itype=itype-256*(itype/256)
! sd1 is 1 if x1 is outside this object, else -1.
      sd1=1.-2.*inside_geom(ndims,x1,iobj)
!      write(*,*)'sd1=',sd1

      if(itype.eq.1)then
! New Sphere code ----------------. 
! ida=0 means use all ndims dimensions (not cyl).
         ida=0
         call spheresect(ndims,ida,x1,xi,obj_geom(ocenter,iobj)
     $        ,obj_geom(oradius,iobj),fraction,f2,sd,C,D)
         if(sd.eq.0 .or. fraction-1..gt.0. .or. fraction.lt.0.)then
            fraction=1.
            sd=0.
         else
! This code decides which of the nf_posno for this object
! to update corresponding to this crossing, and then update it.
! Calculate normalized intersection coordinates.
            call ijbinsphere(iobj,fraction,x1,xi,ijbin)
            call binadding(xi,infobj,sd,ijbin,ispecies)
            if(f2.gt.0. .and. f2.lt.1.)then
               call ijbinsphere(iobj,f2,x1,xi,ijbin)
               call binadding(xi,infobj,sd,ijbin,ispecies)
            endif
         endif
      elseif(itype.eq.2)then
! Cube intersection --------------- new code:
! Normalize to unit cube
         do i=1,ndims
            xn1(i)=(x1(i)-obj_geom(ocenter+i-1,iobj))
     $           /obj_geom(oradius+i-1,iobj)
            xn2(i)=(xi(i)-obj_geom(ocenter+i-1,iobj))
     $           /obj_geom(oradius+i-1,iobj)
         enddo
         call cubeusect(xn1,xn2,nsect,fmin,imin)
!         if(nsect.gt.1)write(*,*)'nsect multiple',nsect,fmin(1),fmin(2)
! The first is the minimum fraction from 1 to 2. But really we need
! to account for all intersections.
         do i=1,nsect
! Crossings alternate direction with initial direction given by sd1
            sd=sd1*(-1)**(i-1)
            if(fmin(i).le.fmax)then
               call ijbincube(iobj,imin(i),fmin(i),xn1,xn2,ijbin,idebug)
               call binadding(xi,infobj,sd,ijbin,ispecies)
            else
!               write(*,'(a,i2,2f10.6,f5.1)')'fmin>fmax',i,fmin(i),fmax
!     $              ,sd
            endif
        enddo
      elseif(itype.eq.3)then
! Cylinder ----------------
         ida=int(obj_geom(ocylaxis,iobj))
         do i=1,ndims
            ii=mod(i-ida+2,ndims)+1
            xn1(ii)=(x1(i)-obj_geom(ocenter+i-1,iobj))
     $           /obj_geom(oradius+i-1,iobj)
            xn2(ii)=(xi(i)-obj_geom(ocenter+i-1,iobj))
     $           /obj_geom(oradius+i-1,iobj)
         enddo
         call cylusect(xn1,xn2,iobj,nsect,fmin,imin)
         do i=1,nsect
            if(fmin(i).le.fmax)then
! Need to decide sd correctly
               sd=sd1*(-1)**(i-1)
               call ijbincyl(iobj,imin(i),fmin(i),xn1,xn2,ijbin)
               call binadding(xi,infobj,sd,ijbin,ispecies)
            endif
         enddo
      elseif(itype.eq.4)then
! Parallelopiped ------------------------
         call xp2contra(iobj,x1,xi,xn1,xn2,ins1,ins2)
         call cubeusect(xn1,xn2,nsect,fmin,imin)
!         do i=1,min(nsect,1)
         do i=1,nsect
            if(fmin(i).le.fmax)then
! First crossing is inward if ins1 is 0.
!               sd=-(2.*ins1-1.)
               sd=sd1*(-1)**(i-1)
               if(inside_geom(ndims,x1,iobj).eq.0)sd=1.
               call ijbincube(iobj,imin(i),fmin(i),xn1,xn2,ijbin,idebug)
               call binadding(xi,infobj,sd,ijbin,ispecies)
            endif
        enddo
      elseif(itype.eq.5)then
! Non-aligned cylinder ---------------------------------
         call xp2contra(iobj,x1,xi,xn1,xn2,ins1,ins2)         
         call cylusect(xn1,xn2,iobj,nsect,fmin,imin)
         do i=1,nsect
            if(fmin(i).le.fmax)then
               sd=sd1*(-1)**(i-1)
               call ijbincyl(iobj,imin(i),fmin(i),xn1,xn2,ijbin)
               call binadding(xi,infobj,sd,ijbin,ispecies)
            endif
         enddo
      elseif(itype.eq.6.or.itype.eq.7)then
! Surface of revolution -------------------------------
         ijbin=-1
         call xp2contra(iobj,x1,xi,xn1,xn2,ins1,ins2)
         call srvsect(xn1,xn2,iobj,nsect,fmin,imin)
         do i=1,nsect
            if(fmin(i).le.fmax)then
               sd=sd1*(-1)**(i-1)
               call ijbinsrv(iobj,imin(i),fmin(i),x1,xi,ijbin)
               call binadding(xi,infobj,sd,ijbin,ispecies)
            endif
         enddo
      else
!         write(*,*)'Unknown object in objsect'
! Unknown object type.
         ierr=99
         return
      endif
!------------------------------
! Saving the x1 and intersection for diagnostics.
      if(sc_ipt.lt.sc_npts)then
         sc_ipt=sc_ipt+1
      endif
! Count up to total intersections. Then store cyclically. 
      isc=mod(isc,sc_npts)+1
      iob_sc(isc)=iobj
      ibin_sc(isc)=ijbin
!      write(*,*)'Saving first intersection',isc,iobj,ijbin
      do i=1,ndims
         x_sc(i,1,isc)=x1(i)
         x_sc(i,2,isc)=x1(i)*(1.-fraction)+xi(i)*fraction
      enddo
!------------------------------
      end
!*******************************************************************
      subroutine binadding(xi,infobj,sd,ijbin,ispecies)
! Add particle xi data to infobj bin ijbin with 
! Surface-crossing-direction sd.
      implicit none
      include 'ndimsdecl.f'
      real xi(3*ndims)
      include '3dcom.f'
      include 'partcom.f'
      include 'plascom.f'
      integer infobj,ispecies
! Not array in this routine.
      integer ijbin
      real sd
      
      integer iaddress,k,id,iq
      real xx

! Adding into bins. At this point all we need is sd, ijbin.
! Multispecies approach: add if_quant to the flux index.
      iq=if_quant(infobj,ispecies)
! Particle Flux.
      iaddress=ijbin+nf_address(nf_flux+iq,infobj,nf_step)
      if(iaddress.le.1)write(*,*)'BIN ERROR',iaddress,ijbin,infobj
     $     ,ispecies,mf_quant(infobj),ff_data(iaddress)
      ff_data(iaddress)=ff_data(iaddress)+sd
! This is the way to test that one is really accessing the right bin:
!      write(*,*)'ijbin=',ijbin,nf_posno(1,infobj),infobj,sd,fraction
!      ff_data(iaddress)=ijbin
      if(kf_quant(infobj,ispecies).ge.2)then
! Momentum               
         iaddress=ijbin+nf_address(nf_gx+iq,infobj,nf_step)
         if(iaddress.le.1)write(*,*)'BIN ERROR',iaddress,ijbin,infobj
     $        ,ispecies,mf_quant(infobj),ff_data(iaddress)
         ff_data(iaddress)=ff_data(iaddress)+ sd*xi(4)
      endif
      if(kf_quant(infobj,ispecies).ge.3)then
         iaddress=ijbin+nf_address(nf_gy+iq,infobj,nf_step)
         if(iaddress.le.1)write(*,*)'BIN ERROR',iaddress,ijbin,infobj
     $        ,ispecies,mf_quant(infobj),ff_data(iaddress)
         ff_data(iaddress)=ff_data(iaddress)+ sd*xi(5)
      endif
      if(kf_quant(infobj,ispecies).ge.4)then
         iaddress=ijbin+nf_address(nf_gz+iq,infobj,nf_step)
         if(iaddress.le.1)write(*,*)'BIN ERROR',iaddress,ijbin,infobj
     $        ,ispecies,mf_quant(infobj),ff_data(iaddress)
         ff_data(iaddress)=ff_data(iaddress)+ sd*xi(6)
      endif
      if(kf_quant(infobj,ispecies).ge.5)then
! Energy
         iaddress=ijbin+nf_address(nf_heat+iq,infobj,nf_step)
         if(iaddress.le.1)write(*,*)'BIN ERROR',iaddress,ijbin,infobj
     $        ,ispecies,mf_quant(infobj),ff_data(iaddress)
         xx=0.
         do k=1,ndims
            xx=xx+xi(3+k)**2
         enddo
         ff_data(iaddress)=ff_data(iaddress)+ sd*xx
      endif
! Accumulate the particle force= momentum/time over whole object.
! Normalized to rhoinf, assuming 1/eoverm is equal to the mass.
      do id=1,ndims
         partforce(id,infobj,nf_step)=partforce(id,infobj,nf_step)
     $        +sd*xi(ndims+id)/(dt*rhoinf*abs(eoverms(ispecies)))
      enddo

      end
!*******************************************************************
      subroutine srvsectplot(iobj,xn1,xn2,fmin)
      include 'ndimsdecl.f'
      include '3dcom.f'
!      include 'partcom.f'
!      include 'sectcom.f'

      real xn1(ndims),xn2(ndims)
!      integer imin(ovlen)
      real fmin(ovlen)
      integer na
      parameter (na=20)
      real ra(na),za(na),xf(ndims)
      
      call autoplot(obj_geom(opr,iobj),obj_geom(opz,iobj),
     $     int(obj_geom(onpair,iobj)))
      call axlabels('r','z')
      do i=1,int(obj_geom(onpair,iobj))-1
! The ends of the segment chosen.
         rb=obj_geom(opr+i-1,iobj)
         rt=obj_geom(opr+i,iobj)
         zb=obj_geom(opz+i-1,iobj)
         zt=obj_geom(opz+i,iobj)
         nz=int(obj_geom(opdiv+i-1,iobj))
         dr2=(rt**2-rb**2)/nz
!          ifct=int((rm**2-rb**2)/dr2)
         do j=1,nz
            rm=sqrt((j-1.)*dr2+rb**2)
            if(abs(dr2).gt.1.e-6)then
               fr=(rm-rb)/(rt-rb)
            else
               fr=(j-1.)/(nz)
            endif
            zm=zb*(1-fr)+zt*fr
            call polymark(rm,zm,1,j+1)
         enddo
      enddo
      do kk=1,2
         do i=1,na
            if(kk.eq.2)then
               fa=fmin(1)+(fmin(2)-fmin(1))*(i-1.)/(na-1.)
            else
               fa=(i-1.)/(na-1.)
            endif
            do k=1,ndims
               xf(k)=xn1(k)*(1.-fa)+xn2(k)*fa
            enddo
            ra(i)=sqrt(xf(1)**2+xf(2)**2)
            za(i)=xf(3)
         enddo
         call color(1+(kk-1)*5)
         call dashset(0)
         if(kk.eq.1)call polymark(ra,za,1,1)
         call polyline(ra,za,na)
         call color(15)
      enddo
      call pltend()
      end
!*******************************************************************
      subroutine pllelfrac(xp,xn,iobj)
! For input point xp(ndims) return the normalized position relative to
! the parallelogram object objg in output xn(ndims)
      include 'ndimsdecl.f'
      include '3dcom.f'
      real xp(ndims),xn(ndims)
      integer iobj
!      real objg(odata)
      do j=1,ndims
         xn(j)=0.
! Contravariant projections.
         do i=1,ndims
! Cartesian coordinates.
            ii=(ocenter+i-1)
!            xc=objg(ii)
            xc=obj_geom(ii,iobj)
! xn1, xn2 are the contravariant coordinates with respect to the center.
            ji=(ocontra+ndims*(j-1)+i-1)
!            write(*,*)'ji',ji
            xn(j)=xn(j)+(xp(i)-xc)*obj_geom(ji,iobj)
         enddo
      enddo
      end
!*********************************************************************
      subroutine timeave(nu,u,uave,ictl)
! Average a quantity u(nu) over steps with a certain decay number
! into uave.
! ictl controls the actions as follows:
!       0    do the averaging.
!       1    do nothing except set nstep=1.
!       2    do nothing except set nstave=nu
!       3    do both 1 and 2. 
!       99   do nothing except increment nstep.
! The 99 call should be done at the end of all usage of this routine
! for the present step.
      real u(nu),uave(nu)

      integer nstep,nstave
      data nstep/1/nstave/20/
! Normal call.
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
!***********************************************************************
      real function fluxdiag()
! Get the total count to object 1 and convert to normalized flux.
! Assuming it's a unit sphere.
      include 'ndimsdecl.f'
      include '3dcom.f'
! For rhoinf, dt
      include 'partcom.f'

      sum=0
      do i=1,nf_posno(nf_flux,1)
         sum=sum+ff_data(nf_address(nf_flux,1,nf_step)+i-1)
      enddo
      fluxdiag=sum/(4.*3.14159)/rhoinf/dt
!      write(*,'(a,f7.0,f8.3)')'Total flux',sum,
!     $     sum/(4.*3.14159)/rhoinf/dt
!      write(*,'(10f7.1)')(ff_data(nf_address(1,1,nf_step)+i-1),
!     $     i=1,nf_posno(1,1))

      end
!***********************************************************************
! Averaging the flux data for quantity abs(iquant)
! over all positions for object ifobj,
! for the steps n1 to n2.
! The positions might be described by more than one dimension, but
! that is irrelevant to the averaging.
! Plot the quantity iquant if positive (not if negative).
! Plotting does not attempt to account for the multidimensionality.
! Rhoinf is returned in rinf.
      subroutine fluxave(n1in,n2in,ifobj,iquant,rinf)
      integer n1in,n2in
!      logical lplot
      integer iquant
      include 'ndimsdecl.f'
      include '3dcom.f'
      include 'sectcom.f'
! Use for averaging:ff_data(nf_address(iq,ifobj,nf_step+1)+i-1)
! which is ff_data(iav+i)
      real fluxofstep(nf_maxsteps),step(nf_maxsteps)
      character*30 string

      ndatamax=200
      n1=n1in
      n2=n2in
      iq=abs(iquant)
! Interrogate the second byte. If non-zero then don't print.
      inoprint=iq/256
      iq=iq-256*inoprint
! If quantity asked for is not available, do nothing.
      if(iq.gt.mf_quant(ifobj).or.iq.eq.0)then
!         write(*,*)'No such quantity',ifobj,mf_quant(ifobj),iq
         return
      endif
! Offset of averaging location:
      iav=nf_address(iq,ifobj,nf_step+1)-1
! Offset to average flux density location
      iavd=nf_address(iq,ifobj,nf_step+2)-1
! Offset to area
      iaa=nf_address(iq,ifobj,nf_pa)-1

!      write(*,*)'Fluxave addresses',iav,iavd,iaa
! Check step numbers for rationality.
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
! Sum the data over steps. 
      tot=0
! Sum the absolute count over steps
      atot=0
      do i=1,nf_posno(iq,ifobj)
         do is=n1,n2
            ff_data(iav+i)=ff_data(iav+i)
     $           +ff_data(nf_address(iq,ifobj,is)+i-1)
         enddo
         tot=tot+ff_data(iav+i)
         atot=atot+abs(ff_data(iav+i))
         ff_data(iav+i)=ff_data(iav+i)/(n2-n1+1)
! Get the area and divide to give the flux density.
         area=ff_data(iaa+i)
         ff_data(iavd+i)=ff_data(iav+i)/area
!         write(*,'(4i6,a,2f10.4)')ifobj,iq,i,iaa+i,' Area,Flux-density'
!     $        ,area,ff_data(iavd +i)
      enddo

! Total the flux over positions as a function of step.
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
!            write(*,*)'n1,n2,is,ff_dt(is)',n1,n2,is,ff_dt(is)
            tdur=tdur+ff_dt(is)
            rinf=rinf+ff_rho(is)*ff_dt(is)
         endif
!         write(*,*)'step',is,' fluxofstep',fluxofstep(is)
      enddo
!      write(*,*)'tot,rinf,tdur,n2',tot,rinf,tdur,n1,n2
      tot=tot/tdur
      rinf=rinf/tdur

! From here on is mostly for testing.
      if(inoprint.eq.0)then
         write(*,'(a,i3,a,i3,a,i5,i5,a,f9.2)')' Average flux quant',iq
     $        ,', object',ifobj,', over steps',n1,n2,', per unit time:'
     $        ,tot
         write(*,'(a,f9.2,a,i7,a,i7,a)')' rhoinf:',rinf,' Total:'
     $        ,nint(tot*tdur),' Abs',nint(atot)
     $        ,'  Average collected per step by posn:'
         if(nf_posno(iq,ifobj).lt.ndatamax)then
            write(*,'(10f8.2)')(ff_data(iav+i),i=1,nf_posno(iq,ifobj))
         else
            write(*,*)'[Too much data to write out.]'
         endif
         fluxdensity=tot/(4.*3.14159)/rinf
         write(*,*)'Flux density*r^2, normalized to rhoinf'
     $        ,fluxdensity
!      write(*,*)'Sectcom ipt=',sc_ipt
         eovermlocal=1.
         flogfac=0.5*alog(2.*3.1415926*eovermlocal/1836.)
         phifloat=alog(fluxdensity)+flogfac
         write(*,*)'Floating potential=',phifloat
      endif

      if(iquant.gt.0)then
         write(string,'(''Object '',i3,'' Quantity'',i3)')
     $        nf_geommap(ifobj),iquant
         call autoplot(step,fluxofstep,n2)
         call boxtitle(string)
         call axlabels('step','Spatially-summed flux number')
         call pltend()
         call multiframe(1,2,0)
         call automark(ff_data(nf_address(iq,ifobj,nf_p1))
     $        ,ff_data(iav+1),nf_posno(iq,ifobj),1)
         call boxtitle(string)
         call axlabels('First flux face variable',
     $        'Time-averaged flux number')
         call automark(ff_data(nf_address(iq,ifobj,nf_p2))
     $        ,ff_data(iav+1),nf_posno(iq,ifobj),1)
         call boxtitle(string)
         call axlabels('Second flux face variable',
     $        '')
         call pltend()
         call multiframe(0,0,0)
      endif
      end
!*******************************************************************
      subroutine writefluxfile(name,nstep)
! Initialize the name before entry. Very Important!
! File name:
      character*(*) name
      include 'ndimsdecl.f'
! Point charge info:
      include 'griddecl.f'
      include 'ptchcom.f'
! Common data containing the BC-object geometric information
      include '3dcom.f'
! Particle common data
      include 'partcom.f'
! Plasma common data
      include 'plascom.f'
! Intersection data
      include 'sectcom.f'
      include 'colncom.f'
      character*(100) charout

! Construct a filename that contains many parameters
! Using the routines in strings_names.f
      call nameconstruct(name)
!     np=nbcat(name,'.flx')
      call nbcat(name,'.flx')
!      write(*,*)name
      write(charout,51)debyelen,Ti,vd,rs,phip
 51   format('debyelen,Ti,vd,rs,phip:',5f10.4,' Version: 6')

!      write(*,*)'mf_obj=',mf_obj,nf_step,mf_quant(1)

      open(22,file=name,status='unknown',err=101)
      close(22,status='delete')
      open(22,file=name,status='new',form='unformatted',err=101)
! This write sequence must be exactly that read below.
      write(22)charout
      write(22)debyelen,Ti,vd,rs,phip ,colntime,subcycle,vneutral
     $     ,fcollided,dropaccel,Tneutral,Eneutral
      write(22)nf_step,mf_quant,mf_obj,(nf_geommap(j),j=1,mf_obj)
!      write(*,*)'geommap',(nf_geommap(j),j=1,mf_obj)
      write(22)(ff_rho(k),k=1,nf_step)
      write(22)(ff_dt(k),k=1,nf_step)
      write(22)((nf_posno(i,j),(nf_dimlens(i,j,k),k=1,ndims)
     $     ,(nf_faceind(i,j,k),k=1,2*ndims)
     $     ,i=1,mf_quant(j)),j=1,mf_obj)
      write(22)(((nf_address(i,j,k),i=1,mf_quant(j)),j=1,mf_obj),
     $     k=1-nf_posdim,nf_step+2)
      ndatalen=nf_address(1,1,nf_step+2)-1
      write(22)(ff_data(i),i=1,ndatalen)
! Debug:
!      write(*,*)'Writefile addresses',
!     $     (nf_address(nf_flux,1,1-j),j=1,nf_posdim),' values'
!     $     ,(ff_data(nf_address(nf_flux,1,1-j)),j=1,nf_posdim)
! New force write.
      write(22)(((fieldforce(i,j,k),pressforce(i,j,k) ,partforce(i,j,k)
     $     ,colnforce(i,j,k),i=1,ndims)
     $     ,charge_ns(j,k),j=1,mf_obj),k=1,nf_step)
! Object data:
      write(22)ngeomobj
      write(22)((obj_geom(j,k),j=1,odata),nf_map(k),k=1,ngeomobj)
      write(22)ibool_part,ifield_mask,iptch_mask,lboundp,rjscheme
! Intersection data:
      write(22)sc_ipt
      write(22)(((x_sc(j,i,k),j=1,ndims),i=1,2),iob_sc(k),
     $     ibin_sc(k),k=1,sc_ipt)
! n_part data
      write(22)(nf_npart(k),k=1,nf_step)
! Species information 
      write(22)nf_species
      write(22)((if_quant(j,k),kf_quant(j,k),j=1,mf_obj),k=1,nf_species)
      write(22)nstep
      close(22)
!      write(*,*)'Wrote flux data to ',name(1:lentrim(name))
      do k=1,nf_step
         do j=1,mf_obj
            do i=1,ndims
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
!*****************************************************************
      subroutine readfluxfile(name,ierr)
! On entry ierr .ne.0 indicates write informational messages.
! On exit  ierr .lt.0 indicates error, else ierr=nstep read.
      character*(*) name
      include 'ndimsdecl.f'
! Point charge info:
      include 'griddecl.f'
      include 'ptchcom.f'
      include '3dcom.f'
      include 'partcom.f'
      include 'plascom.f'
      include 'sectcom.f'
      include 'colncom.f'
      character*(100) charout

      Eneutral=0.
      open(23,file=name,status='old',form='unformatted',err=101)
      read(23)charout
! Figure out the version:
      iend=lentrim(charout)
      if(charout(iend-9:iend-3).eq.'Version')then
         read(charout(iend-1:iend),*)iversion
         if(ierr.ne.0)write(*,*)'Flux file version',iversion
!         ierr=0
      else
         iversion=0
      endif
      if(iversion.le.1)then
         read(23)debyelen,Ti,vd,rs,phip
      elseif(iversion.le.2)then
         read(23)debyelen,Ti,vd,rs,phip
     $        ,colntime,subcycle,vneutral,fcollided,dropaccel,Tneutral
      else
         read(23)debyelen,Ti,vd,rs,phip ,colntime,subcycle,vneutral
     $        ,fcollided,dropaccel,Tneutral,Eneutral
      endif
      ndlen=odata
! Earlier version data length was shorter.
      if(iversion.lt.5)ndlen=42
      read(23)nf_step,mf_quant,mf_obj,(nf_geommap(j),j=1,mf_obj)
!      write(*,*)'geommap',(nf_geommap(j),j=1,mf_obj)
      read(23)(ff_rho(k),k=1,nf_step)
      read(23)(ff_dt(k),k=1,nf_step)
      read(23)((nf_posno(i,j),(nf_dimlens(i,j,k),k=1,ndims)
     $     ,(nf_faceind(i,j,k),k=1,2*ndims)
     $     ,i=1,mf_quant(j)),j=1,mf_obj)
      read(23)(((nf_address(i,j,k),i=1,mf_quant(j)),j=1,mf_obj),
     $     k=1-nf_posdim,nf_step+2)
!      read(23)(ff_data(i),i=1,nf_address(1,1,nf_step+2)-1)
      ndatalen=0
      do j=1,mf_obj
         if(mf_quant(j).ge.1)then
            ndatalen=nf_address(1,j,nf_step+2)-1
            goto 201
         endif
      enddo
 201  continue
!      write(*,*)'Datalen',ndatalen
      read(23)(ff_data(i),i=1,ndatalen)
! Debug:
!      write(*,*)'Readfile addresses',
!     $     (nf_address(nf_flux,1,1-j),j=1,nf_posdim),' values'
!     $     ,(ff_data(nf_address(nf_flux,1,1-j)),j=1,nf_posdim)

      if(iversion.eq.0)then
         write(*,*)'Old force data version',iversion
         read(23,end=102, err=102)(((fieldforce(i,j,k),pressforce(i,j,k)
     $     ,partforce(i,j,k)
     $     ,charge_ns(j,k),i=1,ndims),j=1,mf_obj),k=1,nf_step)
      else
         read(23,end=102, err=102)(((fieldforce(i,j,k),pressforce(i,j,k)
     $     ,partforce(i,j,k)
     $     ,colnforce(i,j,k)
     $     ,i=1,ndims),charge_ns(j,k),j=1,mf_obj),k=1,nf_step)
      endif
! Object data:
      read(23)ngeomobj
      read(23)((obj_geom(j,k),j=1,ndlen),nf_map(k),k=1,ngeomobj)
      read(23)ibool_part,ifield_mask,iptch_mask,lboundp,rjscheme
!      write(*,*)'Object data for',ngeomobj,' objects:'
!      write(*,*)((obj_geom(j,k),j=1,ndlen),nf_map(k),k=1,ngeomobj)
!      write(*,*)((obj_geom(j,k),j=1,8),nf_map(k),k=1,ngeomobj)
! Intersection data:
      read(23)sc_ipt
      read(23)(((x_sc(j,i,k),j=1,ndims),i=1,2),iob_sc(k),
     $     ibin_sc(k),k=1,sc_ipt)
! n_part data
      if(iversion.ge.2)read(23,end=104)(nf_npart(k),k=1,nf_step)
      if(iversion.ge.4)then
         read(23,end=105)nf_species
         read(23,end=105)((if_quant(j,k),kf_quant(j,k),j=1,mf_obj),k=1
     $        ,nf_species)
         if(ierr.ne.0)write(*,*)
     $        'Flux reading: nf_species=',nf_species,
     $        ' Quantity start',(if_quant(1,k),k=1,nf_species),
     $        ' Quantity count',(kf_quant(1,k),k=1,nf_species)
      endif
      if(iversion.ge.6)then
         read(23)ierr
!         write(*,*)'Read nstep from flux file',ierr,' into ierr'
      endif
      goto 103
 102  write(*,*)'Failed to read back forces. Old format? Version='
     $     ,iversion
      ierr=2
      goto 103
 104  write(*,*)'No nf_npart data. Mismatch of versions? Version='
     $     ,iversion
      goto 103
 105  write(*,*)'Incomplete nf_species data. iversion=',iversion
 103  close(23)

! Hack to fix nans when colnforce was wrong. Delete when that data is 
! obsolete.
      isc=0
      do k=1,nf_step
         do j=1,mf_obj
            do i=1,ndims
               if(.not.colnforce(i,j,k).lt.1.e30)then
                  if(isc.eq.0)write(*,*)'Strange colnforce',i,j,k
     $                 ,colnforce(i,j,k)
                  colnforce(i,j,k)=0.
                  isc=isc+1
               endif
            enddo
         enddo
      enddo

! Now one might have to reconstruct nf_faceind from nf_dimlens.
!      if(ierr.ne.0)write(*,*)'Read back flux data from '
!     $     ,name(1:lentrim(name)),' ierr=',ierr
!      ierr=0 ! This prevents the step number being returned. 
      return
 101  write(*,*)'Error opening file:',name(1:lentrim(name))
      ierr=1
      end
!*********************************************************************
!*********************************************************************
      integer function ijbinindex(infobj,iface,i1,i2,i3)
! Given the face and indexes within it, return the ijbin address
! of the flux entry.
! On entry
!      infobj is the object in flux collection index
!      iface is the face number of this object 1-6 cubes. -1,0,+1 cyl.
!      i1, i2, [i3] are the face indexes within it. 
! On exit
!      ijbinindex is the address requested. Or -1 if error.
! For spheres,   itype=1 there is just one face and two angles i1,i2
!     cubes,     itype=2,4:  6 faces, three indices i1,i2,i3
!     cylinders, itype=3,5:  3 faces, three indices ir,it,iz
! For SurfRef,   itype=6,7:  npair+-1 faces, 2 indices each: it,irz
! In reality, only two indices are required for each facet.
! However, it seems less confusing to provide cartesian or r,theta,z 
! indices for cubes or cylinders respectively.
! Thus, for x-faces or circular faces e.g., only i2,i3 matter.
      implicit none
      integer iobj,iface,i1,i2,i3,ijbin
      include 'ndimsdecl.f'
      include '3dcom.f'
      integer ii(3),i,k,ibstep,infobj,itype

      ijbinindex=-1
      ijbin=0
      if(infobj.le.0.or.infobj.gt.mf_obj)then
         write(*,*)'ijbinindex incorrect infobj',infobj
         return
      endif
      iobj=nf_geommap(infobj)
      itype=int(obj_geom(otype,iobj))
!      write(*,*)'ijbinindex',iobj,itype,(obj_geom(k,iobj),k=1,12)
      if(itype.eq.1)then
! sphere ----------------------
         if(i1.gt.nf_dimlens(nf_flux,infobj,1).or.
     $       i2.gt.nf_dimlens(nf_flux,infobj,2))then
            write(*,'(a,i2,a,2i4)')'ijbinindex sphere asked too big'
     $           ,infobj,' i1,i2=',i1,i2
            return
         endif
         ijbin=i1-1+(i2-1)*nf_dimlens(nf_flux,infobj,1)
      elseif(itype.eq.2.or.itype.eq.4)then
! cube ------------------------
         ii(1)=i1
         ii(2)=i2
         ii(3)=i3
         ibstep=1
         if(iface.le.0.or.iface.gt.6)then
            write(*,*)'Incorrect cube face',iface
            return
         endif
         do i=1,ndims-1
            k=mod(mod(iface-1,3)+1+i-1,ndims)+1
            if(ii(k).gt.nf_dimlens(nf_flux,infobj,k))then
               write(*,'(a,i2,a,i1,a,i5)')
     $              'ijbinindex cube asked too big. Object=',infobj
     $              ,' i',k,'=',ii(k)
               return
            endif
            ijbin=ijbin+ibstep*(ii(k)-1)
            ibstep=ibstep*nf_dimlens(nf_flux,infobj,k)
         enddo
         ijbin=ijbin+nf_faceind(nf_flux,infobj,iface)
!         write(*,*)'ijbin',ijbin
      elseif(itype.eq.3.or.itype.eq.5)then
! cylinder---------------------
! iface -1,0,+1 -> 1,2,3
         ijbin=nf_faceind(nf_flux,infobj,iface+2)
         if(iface.eq.0)then
! Index in order theta,z
            if(i2.gt.nf_dimlens(nf_flux,infobj,2).or.
     $           i3.gt.nf_dimlens(nf_flux,infobj,3))then
               write(*,*)'ijbinindex cyl index too big',infobj,i2,i3
            endif
            ijbin=ijbin+i2-1+(i3-1)*nf_dimlens(nf_flux,infobj,2)
         elseif(abs(iface).eq.1)then
! Index in order r,theta
            if(i1.gt.nf_dimlens(nf_flux,infobj,1).or.
     $           i2.gt.nf_dimlens(nf_flux,infobj,2))then
               write(*,*)'ijbinindex end index too big',infobj,i1,i2
            endif
            ijbin=ijbin+(i1-1)+(i2-1)*nf_dimlens(nf_flux,infobj,1)
         else
            write(*,*)'Wrong ijbinindex cylinder face number',iface
            return
         endif
         
      elseif(itype.eq.6.or.itype.eq.7)then
! SoRevolution ------------------------------
! Need test for indexes sensible. Not done.
         ijbin=(i1-1)+nf_dimlens(nf_flux,infobj,1)*(i2-1)
     $        +nf_faceind(nf_flux,infobj,iface)
      endif

      ijbinindex=ijbin

      end
!********************************************************************
      real function getfluxdenave(iquant,ifobj,iface,i1,i2,i3)
! Return the flux density of quantity iquant on object ifobj
! face iface, at facet indices i1,i2,i3
      include 'ndimsdecl.f'
      include '3dcom.f'
      ijbin=ijbinindex(ifobj,iface,i1,i2,i3)
      if(ijbin.ne.-1.)then
         iavd=nf_address(iquant,ifobj,nf_step+2)
     $        +ijbin
         getfluxdenave=ff_data(iavd)
      else
         getfluxdenave=0.
      endif

      end
!********************************************************************
      real function getfluxave(iquant,ifobj,iface,i1,i2,i3)
! Return the flux density of quantity iquant on object ifobj
! face iface, at facet indices i1,i2,i3
      include 'ndimsdecl.f'
      include '3dcom.f'
      ijbin=ijbinindex(ifobj,iface,i1,i2,i3)
      if(ijbin.ne.-1.)then
!      write(*,*)'getfluxave',iquant,ifobj,iface,i1,i2,i3,ijbin
         iavd=nf_address(iquant,ifobj,nf_step+1)
     $        +ijbin
         getfluxave=ff_data(iavd)
      else
         getfluxave=0.
      endif
      end
!********************************************************************
      real function getposcoord(ic,ifobj,iface,i1,i2,i3)
! Return the face position relative coordinate ic on object ifobj
! face iface, at facet indices i1,i2,i3
      include 'ndimsdecl.f'
      include '3dcom.f'
!      write(*,*)'Entering getposcoord',ic,ifobj,iface,i1,i2,i3
      ijbin=ijbinindex(ifobj,iface,i1,i2,i3)
      if(ijbin.ne.-1.)then
         iavd=nf_address(nf_flux,ifobj,1-ic)
     $        +ijbin
         getposcoord=ff_data(iavd)
!      write(*,*)'getposcoord',ic,ifobj,iface,i1,i2,i3,ijbin,getposcoord
      else
         getposcoord=1.e20
      endif
      end
