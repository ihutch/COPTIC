c Particle advancing routine
      subroutine padvnc(ndims,cij,u,iLs)

c Number of dimensions: ndims
      integer ndims
c Storage size of the mesh arrays.
c      real cij(2*ndims+1,nx,ny,nz)
c      real u(nx,ny,nz)
      real cij(*),u(*)

c Array structure vectors: (1,nx,nx*ny,nx*ny*nz)
      integer iLs(ndims+1)

c sormesh provides ixnp, xn, the mesh spacings. (+ndims_mesh)
      include 'meshcom.f'
c Alternatively they could be passed, but we'd then need parameter.
c Local storage
c      integer ixp(ndims_mesh)
      real field(ndims_mesh)
c      real xfrac(ndims_mesh)

      common /myidcom/myid
c Make this always last to use the checks.
      include 'partcom.f'

      if(ndims.ne.ndims_mesh)
     $     stop 'Padvnc incorrect ndims number of dimensions'

      ic1=2*ndims+1

      ndimsx2=2*ndims
      do i=1,npart
         if(if_part(i).ne.0)then
         iregion=insideall(ndims,x_part(1,i))
c Inline version. See also the subroutine.
c Find out where we are (if we don't already know).
c Should not be necessary if chargetomesh has been called.
c         iu=0
c         do id=1,ndims
cc Offset to start of dimension-id-position-array.
c            ioff=ixnp(id)
c            ix=interp(xn(ioff+1),ixnp(id+1)-ioff,x_part(id,i),xm)
c            x_part(ndimsx2+id,i)=xm
cc            xfrac(id)=xm-ix
cc            iu=iu+(ix-1)*iLs(id)
c         enddo

c         call partlocate(i,iLs,iu,ixp,xfrac,iregion)
c         write(*,*)(x_part(ndimsx2+kk,i)-xfrac(kk),kk=1,3)
         r2=x_part(1,i)**2+x_part(2,i)**2+x_part(3,i)**2
         r=sqrt(r2)

c Get the ndims field components at this point. 
c We only use x_part information for location.
         do idf=1,ndims
            call getfield(
     $           ndims,cij(ic1),u,iLs
c     $           ndims,cij(ic1+ic1*iu),u(1+iu),iLs
     $           ,xn(ixnp(idf)+int(x_part(ndimsx2+idf,i)))
     $           ,idf
     $           ,x_part(ndimsx2+1,i)
c     $           ,xfrac
     $           ,iregion,field(idf))
c analytic hack for testing.
c            field(idf)=-x_part(idf,i)*2.*.18/r**3
         enddo

c         write(*,'(''iu='',i6,'' field,anal='',6f9.5)')iu,
c     $        (field(k),-x_part(k,i)*2.*.18/r**3,k=1,3)

c Accelerate          
         do j=4,6
            x_part(j,i)=x_part(j,i)+field(j-3)*dt
         enddo
c Move
         do j=1,3
            x_part(j,i)=x_part(j,i)+x_part(j+3,i)*dt
         enddo          
         
         if(i.le.norbits)then
            iorbitlen(i)=iorbitlen(i)+1
            xorbit(iorbitlen(i),i)=x_part(1,i)
            yorbit(iorbitlen(i),i)=x_part(2,i)
            zorbit(iorbitlen(i),i)=x_part(3,i)
         endif

         endif
      enddo

      end
c***********************************************************************
      subroutine partlocate(i,iLs,iu,ixp,xfrac,iregion)
c Locate the particle numbered i (from common partcom) 
c in the mesh (from common sormesh).
c Return the offset of the base of its cell in iu.
c Return the integer cell-base coordinates in ixp(ndims)
c Return the fractions of cell width at which located in xfrac(ndims)
c Return the region identifier in iregion.
c Store the mesh position into common partcom (x_part).

c sormesh provides ixnp, xn, the mesh spacings. (+ndims_mesh)
      include 'meshcom.f'
      parameter (ndimsx2=ndims_mesh*2)
      integer i,iu,iregion
      integer iLs(ndims_mesh+1)
      integer ixp(ndims_mesh)
      real xfrac(ndims_mesh)

      include 'partcom.f'

      iregion=insideall(ndims_mesh,x_part(1,i))
      iu=0
      do id=1,ndims_mesh
c Offset to start of dimension-id-position array.
         ioff=ixnp(id)
c xn is the position array for each dimension arranged linearly.
c Find the index of xprime in the array xn:
         ix=interp(xn(ioff+1),ixnp(id+1)-ioff,x_part(id,i),xm)
         xfrac(id)=xm-ix
         x_part(ndimsx2+id,i)=xm
         ixp(id)=ix
c should be ix-1
         iu=iu+(ix-1)*iLs(id)
      enddo
      
      end
