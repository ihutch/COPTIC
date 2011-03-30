c Routines for calculating the force on objects via integrations of the
c stress. 
c This requires a "surface representation" of the object.
c Any object is represented by a set of KM facets. Each facet has
c a position x_k(ndims) at which the stress is evaluated, and a surface element
c A_k(ndims) which provides the force vector when the scalar product with
c the stress is taken. The total force is the sum of the KM force vectors.
c The surface representation is here a linear array of 2*ndims*km length.
c But this can be considered to be an array surfobj(2,ndims,km)

      subroutine maxwellforce(ndims,km,surfobj,fieldforce,  u,cij,iLs)
c Calculate the total maxwell force on a surface from its surface
c representation which consists of km facets each of which has 
c ndims position + ndims surface coefficients. 2.ndims.km.

      real surfobj(2*ndims*km),fieldforce(ndims)
      integer km,ndims
      real u(*),cij(*)
      integer iLs(*)

c Local storage
      parameter (mdims=10)
      real field(mdims)

      do i=1,ndims
         fieldforce(i)=0.
      enddo
      do k=1,km
         koff=2*ndims*(k-1)
         call fieldatpoint(surfobj(koff+1),u,cij,iLs,field)
c         call fieldatpointtest(surfobj(koff+1),u,cij,iLs,field)
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

c***************************************************************
c Electron pressure force calculation.
      subroutine pressureforce(ndims,km,surfobj,force,  u,cij,iLs)
c Calculate the total pressure force on a surface from its surface
c representation which consists of km facets each of which has 
c ndims position + ndims surface coefficients. 2.ndims.km.
c The normalized pressure is equal to the exponential of the potential,
c which is the factor by which the electron density is reduced cf infty.

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
c         phi=potentialatpointtest(surfobj(koff+1),u,cij,iLs)
         do i=1,ndims
c 2 Feb 11. The sign was previously +, which was incorrect
c for an outward directed surface normal. 
            force(i)=force(i)-surfobj(koff+ndims+i)*exp(phi)
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


c***************************************************************

      subroutine totalcharge(ndims,km,surfobj,charge,  u,cij,iLs)
c Calculate the total charge within a surface from the surface
c representation which consists of km facets each of which has 
c ndims position + ndims surface coefficients. 2.ndims.km.
c The charge is simply the integral E.dA over the surface.

      real surfobj(2*ndims*km),charge
      integer km,ndims
      real u(*),cij(*)
      integer iLs(*)

c Local storage
      parameter (mdims=10)
      real field(mdims)

      charge=0.
      do k=1,km
         koff=2*ndims*(k-1)
         call fieldatpoint(surfobj(koff+1),u,cij,iLs,field)
c         call fieldatpointtest(surfobj(koff+1),u,cij,iLs,field)
         do j=1,ndims
            charge=charge+surfobj(koff+ndims+j)*field(j)
         enddo
      enddo

      end

c***************************************************************
c Calculate the forces for objects that we are tracking.
      subroutine calculateforces(ndims,iLs,cij,u)
      integer ndims,iLs(ndims+1)
      real u(*),cij(*)
      include '3dcom.f'

      do i=1,nf_obj
         if(ns_flags(i).eq.1)then
c Spheroid. Accumulate.
            km=ns_nt*ns_np
            call maxwellforce(ndims,km,
     $           surfobj(1,1,1,i),fieldforce(1,i,nf_step),
     $           u,cij,iLs)
c            write(*,'(i3,i3,a,3f10.5)')nf_step,i,'  Maxwellforce=',
c     $           (fieldforce(k,i,nf_step),k=1,3)
            call pressureforce(ndims,km,
     $           surfobj(1,1,1,i),pressforce(1,i,nf_step),
     $           u,cij,iLs)
            call totalcharge(ndims,km,
     $           surfobj(1,1,1,i),charge_ns(i,nf_step),
     $           u,cij,iLs)
         elseif(ns_flags(i).eq.513)then
c Point charge object. Just fieldforce.
            call fieldatpoint(obj_geom(ocenter,nf_geommap(i)),
     $           u,cij,iLs,fieldforce(1,i,nf_step))
            charge_ns(i,nf_step)=obj_geom(oradius,nf_geommap(i))
     $           *obj_geom(ofn1,nf_geommap(i))*4.*3.14159
            do nd=1,ns_ndims
               pressforce(nd,i,nf_step)=0.
               fieldforce(nd,i,nf_step)=fieldforce(nd,i,nf_step)*
     $              charge_ns(i,nf_step)
            enddo
         else
c Unknown type. Do nothing.

         endif
      enddo

      end
c***************************************************************
c Initialize the force tracking for objects
c Must be called after fluxdatainit.
      subroutine forcetrackinit()
      include '3dcom.f'

      do i=1,nf_obj
         ns_flags(i)=0
      enddo

      do i=1,ngeomobj
         nfmap=nf_map(i)
         if(nfmap.ne.0)then
c we are tracking this object
            itype=obj_geom(otype,i)
            it=int(itype - 256*(itype/256))
c            write(*,*)'Forceinit:'
c     $           ,i,nfmap,'  type',it
            if(it.eq.1)then
c Spheroid.
               if(itype-it.ne.512)then
                  ns_flags(nfmap)=it
c Initialize the area facet mesh for a unit sphere
                  call spheremesh(ns_nt,ns_np,surfobj(1,1,1,nfmap))
c Now apply the transformation to a spheroid with possibly different
c radii in each direction:
                  r3=obj_geom(oradius,i)*obj_geom(oradius+1,i)
     $                 *obj_geom(oradius+2,i)
c               write(*,*)'Object radius cubed',i,r3
                  do k=1,ns_np
                     do j=1,ns_nt
c                     write(*,*)(surfobj(ii,j,k,nfmap),ii=1,3)
                        do id=1,ns_ndims
                           r=obj_geom(oradius+id-1,i)*1.00001
                           surfobj(id,j,k,nfmap)=
     $                          obj_geom(ocenter+id-1,i)+
     $                          surfobj(id,j,k,nfmap)*r
                           surfobj(ns_ndims+id,j,k,nfmap)=
     $                          surfobj(ns_ndims+id,j,k,nfmap)*r3/r
                        enddo
                     enddo
                  enddo
               else
                  ns_flags(nfmap)=itype
               endif
            else
c An object type I don't know how to handle. Ignore.
            endif
         endif
      enddo

      do i=1,nf_obj
         iobj=nf_geommap(i)
         if(iobj.ne.0)then
         do k=1,ns_np
            do j=1,ns_nt
               r=0.
               do id=1,3
                  r=r+surfobj(id,j,k,i)**2
               enddo
c               write(*,*)'Radius of surface:',i,j,k,r
            enddo
         enddo
         endif
      enddo


      end
c***************************************************************
c Calculate the surfobj surface integration structure for a sphere
c of unit radius, in 3-Dimensions,
c based on uniform mesh of coordinates in cos\theta and psi. 
      subroutine spheremesh(nt,np,surfobj)
c nt,np (INPUT) are the number of theta,psi cells
c surfobj on (OUTPUT) contains the structure 2*ndims*nt*np.
      parameter (ndims=3)

      integer nt,np
      real surfobj(2*ndims*nt*np)

      parameter (PI=3.141593)

c Pointer to start of data
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
c Cartesian Position of this surface:
            surfobj(ipt+1)=cpc*sci
            surfobj(ipt+2)=spc*sci
            surfobj(ipt+3)=cci
c Cartesian Area vector of this surface:
            surfobj(ipt+4)=(spi-spim)*ft
            surfobj(ipt+5)=-(cpi-cpim)*ft
            surfobj(ipt+6)=dpsi*fz
c            write(*,'(i4,i4,7f8.4,6f8.4)')it,ip
c     $           ,cci,psic,ci,psi,cpi,spi
c     $           ,ft,fz
c           write(*,'(6f8.4)')
c     $           ,(surfobj(ipt+k),k=1,6)
c Pointer to start of data for next facet
            ipt=ipt+2*ndims
            cpim=cpi
            spim=spi
            psim=psi
         enddo
         cim=ci
         acm=ac

      enddo

      end

