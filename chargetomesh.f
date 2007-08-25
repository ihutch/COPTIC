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

c For all (active) particles.
      do i=1,npart
         if(if_part(i).ne.0)then
            call partlocate(i,iLs,iu,ixp,xfrac,iregion)
c         if(iregion.ne.iregion) perhaps some action: reinject?

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
                  if(ip.eq.1)then
                     fac=fac*xfrac(ik)
                  else
                     fac=fac*(1.-xfrac(ik))
                  endif
               enddo 
c Add to the particle sum the fraction for this vertex.           
               psum(1+iinc)=fac
               if(ldiags)then

               endif
            enddo
         endif
      enddo

      end
