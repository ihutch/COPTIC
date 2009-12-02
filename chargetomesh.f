      subroutine chargetomesh(psum,iLs,diags)
c particle weight sum, having structure iLs, and (possible) diagnostics.
      real psum(*),diags(*)
c mesh data, notably ndims_mesh, since we don't pass it:
      include 'meshcom.f'
      integer iLs(ndims_mesh+1)
      real xfrac(ndims_mesh)
      integer ixp(ndims_mesh)

      include 'partcom.f'
c On entry, psum ought to have been initialized to zero.

c No diagnostics for now.
      ldiags=.false. 

c For all (possibly-active) particles.
c      do i=1,n_part
      do i=1,ioc_part
         if(if_part(i).ne.0)then
            inewregion=insideall(ndims_mesh,x_part(1,i))
            x1=x_part(ndims_mesh*2+1,i)
c Alternative to partlocate:
            iu=0
            do id=1,ndims_mesh
               ix=int(x_part(ndims_mesh*2+id,i))
               iu=iu+(ix-1)*iLs(id)
            enddo
            x2=x_part(ndims_mesh*2+1,i)
            if(x1.ne.x2)write(*,*)'Mesh pos change',i,x1,x2
c Cycle through the vertices of the box we are in.
            do ii=0,2**ndims_mesh-1
               ii1=ii
               iinc=iu
               fac=1.
c Calculate the index and weight of this vertex.
               do ik=1,ndims_mesh
c There may be bit-manipulation routines to accelerate this.
c But likely not by very much since the main cost is divide+mult by 2
                  ii2=ii1/2
                  ip=ii1-2*ii2
                  ii1=ii2
                  iinc=iinc+ip*iLs(ik)
                  xm=x_part(ndims_mesh*2+ik,i)
                  ix=int(xm)
                  xf=xm-ix
                  if(ip.eq.1)then
c                     fac=fac*xfrac(ik)
                     fac=fac*xf
                  else
c                     fac=fac*(1.-xfrac(ik))
                     fac=fac*(1.-xf)
                  endif
               enddo 
c Add to the particle sum the fraction for this vertex.
c               write(*,*)'ii,iu,iinc,fac',ii,iu,iinc,fac
               psum(1+iinc)=psum(1+iinc)+fac
               if(ldiags)then
c Do something with diagnostics.                  
                  diags(1)=0.
               endif
            enddo
         endif
      enddo

      end

c********************************************************************
      subroutine psumtoq(inc,ipoint,indi,ndims,iused,
     $     psum,rho,volumes,u)
      integer ipoint,inc
      integer indi(ndims),iused(ndims)
      real psum(*),rho(*),volumes(*),u(*)
c Partcom gives us rhoinf:
      include 'partcom.f'

c Silence warnings with spurious access.
      ind=iused(1)
      ind=indi(1)
c This routine for use in mditerarg.
c But we iterate only over the inner mesh (not edges).
c Here, u=psum, v=rho, w=volumes. 
c Set the density
      ind=1+ipoint
      if(volumes(ind).gt.1.e20)then
c This is outside the region. Compensate the electron density.
         rho(ind)=faddu(u(ind),fprime)
c         rho(ind)=0.
      else
         rho(ind)=psum(ind)/(rhoinf*volumes(ind))
      endif
      inc=1
      end
