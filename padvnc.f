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
      integer ixp(ndims_mesh)
      real field(ndims_mesh)
      real xfrac(ndims_mesh)

      common /myidcom/myid,nprocs
c Make this always last to use the checks.
      include 'partcom.f'

      if(ndims.ne.ndims_mesh)
     $     stop 'Padvnc incorrect ndims number of dimensions'

      ic1=2*ndims+1

      ndimsx2=2*ndims
c We ought not to need to calculate the iregion, since it should be
c known and if a particle is outside it, it would have been reinjected:
c But for now:
      iregion=insideall(ndims,x_part(1,1))
      do i=1,n_part
c If this particle slot is occupied.
         if(if_part(i).ne.0)then
            dtprec=dt
            dtpos=dt
c Find out where we are (if we don't already know).
c Should not be necessary if chargetomesh has been called.
c         call partlocate(i,iLs,iu,ixp,xfrac,iregion)
c         write(*,*)(x_part(ndimsx2+kk,i)-xfrac(kk),kk=1,3)

c Subcycle start.
 101        continue
c Use dtaccel for acceleration. May be different from dt if there was
c a reinjection (or collision).
            dtaccel=0.5*(dt+dtprec)
c Get the ndims field components at this point. 
c We only use x_part information for location.
            do idf=1,ndims
               call getfield(
     $              ndims,cij(ic1),u,iLs
     $              ,xn(ixnp(idf)+int(x_part(ndimsx2+idf,i)))
     $              ,idf
     $              ,x_part(ndimsx2+1,i)
     $              ,iregion,field(idf))
c analytic hack for testing.
c            field(idf)=-x_part(idf,i)*2.*.18/r**3
            enddo

c         write(*,'(''iu='',i6,'' field,anal='',6f9.5)')iu,
c     $        (field(k),-x_part(k,i)*2.*.18/r**3,k=1,3)

c Accelerate          
            do j=4,6
               x_part(j,i)=x_part(j,i)+field(j-3)*dtaccel
            enddo
c Move
            do j=1,3
               x_part(j,i)=x_part(j,i)+x_part(j+3,i)*dtpos
            enddo          

            inewregion=insideall(ndims,x_part(1,i))

            if(inewregion.ne.iregion) then
c We left the region. Reinject.
               call reinject(i,x_part(1,i))
c Find where we are, since we don't yet know?
c Might not be needed if we insert needed information in reinject,
c which might be less costly.
               call partlocate(i,iLs,iu,ixp,xfrac,iregion)
               dtpos=dtpos*ran0(idum)
               dtprec=0.
c Complete reinjection by advancing by random remaining.
               goto 101
            endif

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
