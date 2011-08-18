      subroutine chargetomesh(psum,iLs,diagsum,ndiags)
c Assign charge and other moments to the mesh accumulators.
c Particle weight sum, having structure iLs, and (possible) diagnostics.
c Which are enumerated up to ndiags in their trailing dimension.
      real psum(*)
      real diagsum(*)
c mesh data, notably ndims_mesh, since we don't pass it:
      include 'meshcom.f'
      integer iLs(ndims_mesh+1)

      include 'partcom.f'
c On entry, psum ought to have been initialized to zero.

c For all (possibly-active) particles.
      do i=1,ioc_part
         if(if_part(i).ne.0)then
            inewregion=insideall(ndims_mesh,x_part(1,i))
c            x1=x_part(ndims_mesh*2+1,i)
c Alternative to partlocate:
            iu=0
            do id=1,ndims_mesh
               ix=int(x_part(ndims_mesh*2+id,i))
               iu=iu+(ix-1)*iLs(id)
            enddo
c            x2=x_part(ndims_mesh*2+1,i)
c This test is irrelevant now partlocate is not used.
c            if(x1.ne.x2)write(*,*)'Mesh pos change',i,x1,x2
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
               if(ndiags.gt.0)then
c Do something with diagnostics. 
c Assumption here is diags1 n, diags2-4 v, diags5-7 v^2.
                  diagsum(1+iinc)=diagsum(1+iinc)+fac
                  do k=1,ndiags-1
c Six moments. 3 for v and 3 for v^2.
                     temp=x_part(ndims_mesh+1+mod(k-1,ndims_mesh),i)
                     if(k.gt.ndims_mesh)temp=temp*temp
                     diagsum(1+iinc+k*iLs(ndims_mesh+1))=
     $                    diagsum(1+iinc+k*iLs(ndims_mesh+1))+
     $                    fac*temp
                  enddo
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
         if(.not.u(ind).lt.1.e20)then
            write(*,*)'psumtoq error',indi,rho(ind),u(ind),ind
         endif
         rho(ind)=faddu(u(ind),fprime,ind)
c         rho(ind)=0.
      else
         rho(ind)=psum(ind)/(rhoinf*volumes(ind))
      endif
c      if(.not.rho(ind).lt.1.e20)then
c         write(*,*)'psumtoq rho error',ind,rho(ind),psum(ind),rhoinf
c      endif
      inc=1
      end
c********************************************************************
      subroutine quasineutral(inc,ipoint,indi,ndims,iused,
     $     q,u,volumes)
c     ,dum3,dum4)
      integer ipoint,inc
      integer indi(ndims),iused(ndims)
      real q(*),u(*),volumes(*)
      include 'plascom.f'
c      real dum3,dum4

c Silence warnings with spurious accesses. 
c (Which ought to be optimized away.)
      ind=iused(1)
      ind=indi(1)
c This routine for use in mditerarg.
c But we iterate only over the inner mesh (not edges) by virtue of call.
      ind=1+ipoint
c This test in ineffective, because q was set by the psumtoq to
c compensate electron density.
c      if(q(ind).gt.0.)then
c So instead tell directly from volumes if we are outside:
      if(volumes(ind).le.1.e20)then
         u(ind)=alog(q(ind))
      else
         u(ind)=phip
      endif
      inc=1
      end
