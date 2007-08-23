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
      include 'sormesh.f'
c Alternatively they could be passed, but we'd then need parameter.
c Local storage
      integer ixp(ndims_mesh)
      real field(ndims_mesh),xfrac(ndims_mesh)

      common /myidcom/myid
c Make this always last to use the checks.
      include 'partcom.f'

      if(ndims.ne.ndims_mesh)
     $     stop 'Padvnc incorrect ndims number of dimensions'

      ic1=2*ndims+1

      do i=1,npart

c Find out where we are (if we don't already know).
         iregion=insideall(ndims,x_part(1,i))
c         write(*,*)'iregion=',iregion
         iu=0
         do id=1,ndims
c Offset to start of idf position array.
            ioff=ixnp(id)
c xn is the position array for each dimension arranged linearly.
c Find the index of xprime in the array xn:
            ix=interp(xn(ioff+1),ixnp(id+1)-ioff,x_part(id,i),xm)
            xfrac(id)=xm-ix
            ixp(id)=ix
            iu=iu+ix*iLs(id)
         enddo

c Get the ndims field components at this point.
         do idf=1,ndims
            ioff=ixnp(idf)
c            write(*,*)'idf,iu,ioff',idf,iu,ioff,xfrac
            call getfield(
     $           ndims
     $           ,cij(ic1+ic1*iu)
     $           ,u(1+iu)
     $           ,iLs
     $           ,xn(ioff+ixp(idf)),idf
     $           ,xfrac,iregion,field(idf))
         enddo

c         write(*,*)'ixp=',ixp,' Field=',field

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

      enddo

      end
