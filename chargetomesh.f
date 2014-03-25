      subroutine chargetomesh(psum,mdims,iLs,diagsum,ndiags)
c Assign charge and other moments to the mesh accumulators.
c Particle weight sum, having structure iLs, and (possible) diagnostics.
c Which are enumerated up to ndiags in their trailing dimension.
      real psum(*)
      real diagsum(*)
c mesh data, notably ndims, since we don't pass it:
      include 'ndimsdecl.f'
      include 'meshcom.f'
      integer iLs(ndims+1)

      include 'partcom.f'
c On entry, psum ought to have been initialized to zero.

c For all (possibly-active) particles.
      do i=1,ioc_part
         if(if_part(i).ne.0)then
            inewregion=insideall(ndims,x_part(1,i))
c            x1=x_part(ndims*2+1,i)
c Alternative to partlocate:
            iu=0
            do id=1,ndims
               ix=int(x_part(ndims*2+id,i))
               iu=iu+(ix-1)*iLs(id)
            enddo
c Cycle through the vertices of the box we are in.
            do ii=0,2**ndims-1
               ii1=ii
               iinc=iu
               fac=1.
c Calculate the index and weight of this vertex.
               do ik=1,ndims
c There may be bit-manipulation routines to accelerate this.
c But likely not by very much since the main cost is divide+mult by 2
                  ii2=ii1/2
                  ip=ii1-2*ii2
                  ii1=ii2
                  iinc=iinc+ip*iLs(ik)
                  xm=x_part(ndims*2+ik,i)
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
                     temp=x_part(ndims+1+mod(k-1,ndims),i)
                     if(k.gt.ndims)temp=temp*temp
                     diagsum(1+iinc+k*iLs(ndims+1))=
     $                    diagsum(1+iinc+k*iLs(ndims+1))+
     $                    fac*temp
                  enddo
               endif
            enddo
         endif
      enddo
c End of charge deposition.
      end
c********************************************************************
      subroutine psumperiod(psum,mdims,ifull,iaux,iLs)
c If there are periodic particles in any dimension, do the periodic
c exchange sum.
      integer mdims,ifull(mdims),iaux(mdims)
c,iuds(ndims)
      real psum(*)
      include 'ndimsdecl.f'
      integer iLs(ndims+1)

      include 'meshcom.f'
      integer iview(3,ndims),indi(ndims)
      integer istart,iend,istride
      parameter (istart=1,iend=2,istride=3)
      include 'partcom.f'
c This data statement serves to silence ftnchek. The first mditerator
c call actually initializes the iview and indi.
      data iview/ndims*0,ndims*0,ndims*0/
      data indi/ndims*0/

      do id=1,ndims
         if(ipartperiod(id).eq.4)then
c Set view to entire array (Offsets indi [0:iaux(id)-1]).
            icomplete=mditerator(ndims,iview,indi,4,iaux)
c Use the general iterator to sum periodically.
c Slice dimension id:
            iview(iend,id)=0
 101        ii=indexcontract(ndims,ifull,indi)
            ib1=1+iLs(id)+ii
            it=1+(iaux(id)-1)*iLs(id)+ii
            it1=1+(iaux(id)-2)*iLs(id)+ii
            ib=1+ii
c Face+1=Face+1 + Other
            psum(ib1)=psum(ib1)+psum(it)
c Other-1=Other-1 + Face
            psum(it1)=psum(it1)+psum(ib)
            if(mditerator(ndims,iview,indi,0,iaux).eq.0)goto 101
         endif
      enddo
      end
c********************************************************************
      subroutine diagperiod(diagsum,mdims,ifull,iaux,iLs,ndiags)
c If there are periodic particles in any dimension, do the periodic
c exchange sum.
      integer mdims,ifull(mdims),iaux(mdims)
      real diagsum(*)
      integer iLs(mdims+1)
      include 'ndimsdecl.f'
      include 'meshcom.f'
      integer iview(3,ndims),indi(ndims)
      integer istart,iend,istride
      parameter (istart=1,iend=2,istride=3)
      include 'partcom.f'
c This data statement serves to silence ftnchek. The first mditerator
c call actually initializes the iview and indi.
      data iview/ndims*0,ndims*0,ndims*0/
      data indi/ndims*0/

      do id=1,ndims
         if(ipartperiod(id).eq.4)then
c Set view to entire array (Offsets indi [0:iaux(id)-1]).
            icomplete=mditerator(ndims,iview,indi,4,iaux)
c Use the general iterator to sum periodically.
c Slice dimension id:
            iview(iend,id)=0
 101        ii=indexcontract(ndims,ifull,indi)
            ib1=1+iLs(id)+ii
            it=1+(iaux(id)-1)*iLs(id)+ii
            it1=1+(iaux(id)-2)*iLs(id)+ii
            ib=1+ii
c Diagnostic particles must be transferred, not just added.
            do k=0,ndiags-1
               kshift=k*iLs(ndims+1)
               diagsum(ib1+kshift)= diagsum(ib1+kshift)+ diagsum(it
     $              +kshift)
               diagsum(it+kshift)=0.
               diagsum(it1+kshift)= diagsum(it1+kshift)+ diagsum(ib
     $              +kshift)
               diagsum(ib+kshift)=0.
            enddo
            if(mditerator(ndims,iview,indi,0,iaux).eq.0)goto 101
         endif
      enddo
      end

c********************************************************************
      subroutine psumtoq(inc,ipoint,indi,ndims,iLs,iused,
     $     psum,rho,volumes,u)
      integer ipoint,inc
      integer indi(ndims),iused(ndims)
      real psum(*),rho(*),volumes(*),u(*)
      real faddu
      external faddu
c Partcom gives us rhoinf:
      include 'partcom.f'

c Silence warnings with spurious access.
      ind=iLs
      ind=iused(1)
      ind=indi(1)
c This routine for use in mditerarg.
c But we iterate only over the inner mesh (not edges).
c Here, u=psum, v=rho, w=volumes. 
c Set the density
      ind=1+ipoint
      if(volumes(ind).ge.1.e20)then
c This is outside the region. 
         if(volumes(ind).ge.1.e30)then
c And all point-charge regions.
c Compensate the electron density.
            if(.not.u(ind).lt.1.e20)then
               write(*,*)'psumtoq error',indi,rho(ind),u(ind),ind
            endif
            rho(ind)=faddu(u(ind),fprime,ind)
         else
c Outside the particle region, but inside a point-charge region.
c Don't compensate for electron density. Just set rhoi=0.
c            rho(ind)=0.
c That does not seem correct. I still ought to null out the electron
c density. Do that directly, without the additional rhoci part:
            fprime=exp(u(ind))
            rho(ind)=fprime
         endif
      else
c Standard case.
         rho(ind)=psum(ind)/(rhoinf*volumes(ind))
      endif
c      if(.not.rho(ind).lt.1.e20)then
c         write(*,*)'psumtoq rho error',ind,rho(ind),psum(ind),rhoinf
c      endif
      inc=1
      end
c********************************************************************
      subroutine quasineutral(inc,ipoint,indi,ndims,iLs,iused,
     $     q,u,volumes)
c     ,dum3,dum4)
      integer ipoint,inc
      integer indi(ndims),iused(ndims)
      real q(*),u(*),volumes(*)
      include 'plascom.f'
c      real dum3,dum4

c Silence warnings with spurious accesses. 
c (Which ought to be optimized away.)
      ind=iLs
      ind=iused(1)
      ind=indi(1)
c This routine for use in mditerarg.
c But we iterate only over the inner mesh (not edges) by virtue of call.
      ind=1+ipoint
c This test is ineffective, because q was set by the psumtoq to
c compensate electron density.
c      if(q(ind).gt.0.)then
c So instead tell directly from volumes if we are outside:
      if(volumes(ind).le.1.e20)then
c There still might be too small a density, so set a floor for it.
         if(q(ind).lt.1.e-6)q(ind)=1.e-6
         u(ind)=alog(q(ind))
      else
         u(ind)=phip
      endif
      if(abs(u(ind)).gt.20.)write(*,*)'Large phi',ind,u(ind),phip
      inc=1
      end
