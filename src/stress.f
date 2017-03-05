! Routines for calculating the force on objects via integrations of the
! stress. 
! This requires a "surface representation" of the object.
! Any object is represented by a set of KM facets. Each facet has
! a position x_k(ndims) at which the stress is evaluated, and a surface element
! A_k(ndims) which provides the force vector when the scalar product with
! the stress is taken. The total force is the sum of the KM force vectors.
! The surface representation is here a linear array of 2*ndims*km length.
! But this can be considered to be an array surfobj(2,ndims,km)

      subroutine maxwellforce(ndims,km,surfobj,fieldforce,  u,cij,iLs)
! Calculate the total maxwell force on a surface from its surface
! representation which consists of km facets each of which has 
! ndims position + ndims surface coefficients. 2.ndims.km.

      integer km,ndims
      real surfobj(2*ndims*km),fieldforce(ndims)
      real u(*),cij(*)
      integer iLs(*)

! Local storage
      parameter (mdims=3)
      real field(mdims)

      do i=1,ndims
         fieldforce(i)=0.
      enddo
      do k=1,km
         koff=2*ndims*(k-1)
         call fieldatpoint(surfobj(koff+1),u,cij,iLs,field)
!         call fieldatpointtest(surfobj(koff+1),u,cij,iLs,field)
         es=0.
         do i=1,ndims
            es=es+field(i)**2
         enddo
         do i=1,ndims
            do j=1,ndims
               stress=field(i)*field(j)
               if(i.eq.j)stress=stress-es/2.
               fieldforce(i)=fieldforce(i)+surfobj(koff+ndims+j)*stress
            enddo
         enddo
      enddo

      end

!***************************************************************
! Electron pressure force calculation.
      subroutine pressureforce(ndims,km,surfobj,force,u,cij,iLs)
! Calculate the total pressure force on a surface from its surface
! representation which consists of km facets each of which has 
! ndims position + ndims surface coefficients. 2.ndims.km.
! The normalized pressure is equal to the exponential of the potential,
! which is the factor by which the electron density is reduced cf infty.

      integer km,ndims
      real surfobj(2*ndims*km),force(ndims)
      real u(*),cij(*)
      integer iLs(*)


      do i=1,ndims
         force(i)=0.
      enddo
      do k=1,km
         koff=2*ndims*(k-1)
         phi=potentialatpoint(surfobj(koff+1),u,cij,iLs)
!         phi=potentialatpointtest(surfobj(koff+1),u,cij,iLs)
         do i=1,ndims
! 2 Feb 11. The sign was previously +, which was incorrect
! for an outward directed surface normal. 
            force(i)=force(i)-surfobj(koff+ndims+i)*exp(phi)
! 5 Nov 13 trap for unphysical forces.
            if(.not.abs(force(i)).lt.1.e20)then 
               write(*,*)'k,i,koff,phi,force(i)',k,i,koff,phi,force(i)
               write(*,*)'x=',(surfobj(koff+kk),kk=1,3)
               stop
     $           '**** pressureforce error.'
            endif
         enddo
         if(phi.eq.9999)then
            r2=0.
            do i=1,ndims
               r2=r2+surfobj(koff+i)**2
            enddo
            r=sqrt(r2)
            write(*,'(a,3i4,4f10.5)')'pressureforce',k,km,koff
     $           ,(surfobj(koff+j),j=1,3),r
         endif
      enddo
      
      end


!***************************************************************

      subroutine totalcharge(ndims,km,surfobj,charge,  u,cij,iLs)
! Calculate the total charge within a surface from the surface
! representation which consists of km facets each of which has 
! ndims position + ndims surface coefficients. 2.ndims.km.
! The charge is simply the integral E.dA over the surface.

      integer km,ndims
      real surfobj(2*ndims*km),charge
      real u(*),cij(*)
      integer iLs(*)

! Local storage
      parameter (mdims=3)
      real field(mdims)

      charge=0.
      do k=1,km
         koff=2*ndims*(k-1)
         call fieldatpoint(surfobj(koff+1),u,cij,iLs,field)
!         call fieldatpointtest(surfobj(koff+1),u,cij,iLs,field)
         do j=1,ndims
            charge=charge+surfobj(koff+ndims+j)*field(j)
         enddo
      enddo

      end

!***************************************************************
! Calculate the forces for objects that we are tracking.
      subroutine calculateforces(mdims,iLs,cij,u)
      integer mdims,iLs(mdims+1)
      real u(*),cij(*)
      include 'ndimsdecl.f'
      include '3dcom.f'

      do i=1,nf_obj
         if(ns_flags(i).eq.1)then
! Spheroid. Accumulate.
            km=ns_nt*ns_np
            call maxwellforce(ndims,km,
     $           surfobj(1,1,1,i),fieldforce(1,i,nf_step),
     $           u,cij,iLs)
!            write(*,'(i3,i3,a,3f10.5)')nf_step,i,'  Maxwellforce=',
!     $           (fieldforce(k,i,nf_step),k=1,3)
            call pressureforce(ndims,km,
     $           surfobj(1,1,1,i),pressforce(1,i,nf_step),
     $           u,cij,iLs)
            call totalcharge(ndims,km,
     $           surfobj(1,1,1,i),charge_ns(i,nf_step),
     $           u,cij,iLs)
         elseif(ns_flags(i).eq.513)then
! Point charge object. Just fieldforce.
            call fieldatpoint(obj_geom(ocenter,nf_geommap(i)),
     $           u,cij,iLs,fieldforce(1,i,nf_step))
            charge_ns(i,nf_step)=obj_geom(oradius,nf_geommap(i))
     $           *obj_geom(ofn1,nf_geommap(i))*4.*3.14159
            do nd=1,ndims
               pressforce(nd,i,nf_step)=0.
               fieldforce(nd,i,nf_step)=fieldforce(nd,i,nf_step)*
     $              charge_ns(i,nf_step)
            enddo
         else
! Unknown type. Do nothing.
         endif
      enddo

      end
!***************************************************************
! Initialize the force tracking for objects
! Must be called after fluxdatainit.
      subroutine forcetrackinit()
      include 'ndimsdecl.f'
      include '3dcom.f'

      do i=1,nf_obj
         ns_flags(i)=0
      enddo

      do i=1,ngeomobj
         nfmap=nf_map(i)
         if(nfmap.ne.0)then
! we are tracking this object
            itype=int(obj_geom(otype,i))
            it=int(itype - 256*(itype/256))
!            write(*,*)'Forceinit:'
!     $           ,i,nfmap,'  type',it
            if(it.eq.1)then
! Spheroid.
               if(itype-it.ne.512)then
                  ns_flags(nfmap)=it
! Initialize the area facet mesh for a unit sphere
                  call spheremesh(ns_nt,ns_np,surfobj(1,1,1,nfmap))
! Now apply the transformation to a spheroid with possibly different
! radii in each direction:
                  r3=obj_geom(oradius,i)*obj_geom(oradius+1,i)
     $                 *obj_geom(oradius+2,i)
!               write(*,*)'Object radius cubed',i,r3
                  do k=1,ns_np
                     do j=1,ns_nt
!                     write(*,*)(surfobj(ii,j,k,nfmap),ii=1,3)
                        do id=1,ndims
                           r=obj_geom(oradius+id-1,i)*1.00001
                           surfobj(id,j,k,nfmap)=
     $                          obj_geom(ocenter+id-1,i)+
     $                          surfobj(id,j,k,nfmap)*r
                           surfobj(ndims+id,j,k,nfmap)=
     $                          surfobj(ndims+id,j,k,nfmap)*r3/r
                        enddo
                     enddo
                  enddo
               else
                  ns_flags(nfmap)=itype
               endif
            else
! An object type I don't know how to handle. Ignore.
            endif
         endif
      enddo

!      do i=1,nf_obj
!         iobj=nf_geommap(i)
!         if(iobj.ne.0)then
!         do k=1,ns_np
!            do j=1,ns_nt
!               r=0.
!               do id=1,3
!                  r=r+surfobj(id,j,k,i)**2
!               enddo
!               write(*,*)'Radius of surface:',i,j,k,r
!            enddo
!         enddo
!         endif
!      enddo


      end
!***************************************************************
! Calculate the surfobj surface integration structure for a sphere
! of unit radius, in 3-Dimensions,
! based on uniform mesh of coordinates in cos\theta and psi. 
      subroutine spheremesh(nt,np,surfobj)
! nt,np (INPUT) are the number of theta,psi cells
! surfobj on (OUTPUT) contains the structure 2*ndims*nt*np.
      parameter (ndims=3)

      integer nt,np
      real surfobj(2*ndims*nt*np)

      parameter (PI=3.141593)

! Pointer to start of data
      ipt=0
      cim=1.
      acm=0.
      dpsi=2.*PI/np
      do it=1,nt
         ci=1.-(2*it)/float(nt)
         ac=acos(ci)
         ft=0.5*(ac-acm-ci*sqrt(1.-ci**2)+cim*sqrt(1.-cim**2))
         fz=-0.5*(ci**2-cim**2)
         thi=0.5*(ac+acm)
         cci=cos(thi)
         sci=sin(thi)
         psim=-PI
         cpim=-1.
         spim=0.
         do ip=1,np
            psi=PI*(-1.+(2*ip)/float(np))
            cpi=cos(psi)
            spi=sin(psi)
            psic=0.5*(psi+psim)
            cpc=cos(psic)
            spc=sin(psic)
! Cartesian Position of this surface:
            surfobj(ipt+1)=cpc*sci
            surfobj(ipt+2)=spc*sci
            surfobj(ipt+3)=cci
! Cartesian Area vector of this surface:
            surfobj(ipt+4)=(spi-spim)*ft
            surfobj(ipt+5)=-(cpi-cpim)*ft
            surfobj(ipt+6)=dpsi*fz
!            write(*,'(i4,i4,7f8.4,6f8.4)')it,ip
!     $           ,cci,psic,ci,psi,cpi,spi
!     $           ,ft,fz
!           write(*,'(6f8.4)')
!     $           ,(surfobj(ipt+k),k=1,6)
! Pointer to start of data for next facet
            ipt=ipt+2*ndims
            cpim=cpi
            spim=spi
            psim=psi
         enddo
         cim=ci
         acm=ac

      enddo

      end

!**********************************************************************
      subroutine chargeforce(ndims,km,surfobj,fieldforce,u,iLs)
! Calculate the force on a unit point charge inside a surface 
! By averaging the electric field on the surface whose
! representation consists of km facets each of which has 
! ndims position + ndims surface coefficients. 2.ndims.km.

      integer km,ndims
      real surfobj(2*ndims*km),fieldforce(ndims)
      real u(*)
      integer iLs(*)

! Local storage
      parameter (mdims=3)
      real field(mdims)

      do i=1,ndims
         fieldforce(i)=0.
      enddo
      totarea=0.
      do k=1,km
         koff=2*ndims*(k-1)
         call fieldatpoint(surfobj(koff+1),u,cij,iLs,field)
!         call fieldsimple3atpoint(surfobj(koff+1),u,iLs,field)
         area=0
         do j=1,ndims
            area=area+surfobj(koff+ndims+j)**2
         enddo
         area=sqrt(area)
!         write(*,*)'k=',k,' area=',area,' psn',(surfobj(koff+i),i=1,3)
         totarea=totarea+area
         do i=1,ndims
            fieldforce(i)=fieldforce(i)+field(i)*area
         enddo
!         write(*,*)' field=',field
      enddo
!      write(*,*)'surfobj',surfobj
!      write(*,*)'totarea',totarea,' fieldforce',fieldforce
      do i=1,ndims
         fieldforce(i)=fieldforce(i)/totarea
      enddo
      end
!***************************************************************
! Initialize the force tracking for one object with center xc and 
! radius r. With angular mesh ns_nt, ns_np.
      subroutine forceinitone(ns_nt,ns_np,ndims,r,xc,surfobj)
      real surfobj(2*ndims,ns_nt,ns_np)
      real r,xc(ndims)

! Initialize the area facet mesh for a unit sphere
      call spheremesh(ns_nt,ns_np,surfobj(1,1,1))
! Now apply the transformation to a spheroid with radius r
      do k=1,ns_np
         do j=1,ns_nt
            do id=1,ndims
               surfobj(id,j,k)=
     $              xc(id)+surfobj(id,j,k)*r
               surfobj(ndims+id,j,k)=
     $              surfobj(ndims+id,j,k)*r**2
            enddo
         enddo
      enddo

      end
