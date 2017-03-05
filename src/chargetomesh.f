!********************************************************************
      subroutine chargetomesh(psum,iLs,diagsum,ndiags)
! Assign charge and other moments to the mesh accumulators.
! Particle weight sum, having structure iLs, and (possible) diagnostics.
! Which are enumerated up to ndiags in their trailing dimension.
      real psum(*)
      real diagsum(*)
! mesh data, notably ndims, since we don't pass it:
      include 'ndimsdecl.f'
      integer iLs(ndims+1)
      include 'plascom.f'
      include 'partcom.f'
! On entry, psum ought to have been initialized to zero.
! Needed for boltzamp.
      include 'meshcom.f'
      include 'ptchcom.f'

! For all (possibly-active) particles.
      do ispecies=1,nspecies
         echarge=numratioa(ispecies)*sign(1.,eoverms(ispecies))
         if(ispecies.eq.2)echarge=echarge*(1.-boltzamp)
         do i=iicparta(ispecies),iocparta(ispecies)
            if(x_part(iflag,i).ne.0)then
               call achargetomesh(i,psum,iLs,diagsum,ndiags,echarge
     $              ,x_part(1,i))
            endif
         enddo
      enddo
! End of charge deposition.
      end
!********************************************************************
      subroutine achargetomesh(i,psum,iLs,diagsum,ndiags,echarge,xprior)
! Assign charge and other moments to the mesh accumulators,
! for a single particle i of echarge value, which only affects psum
! Particle weight sum, having structure iLs.
! Also (possible) diagnostics diagsum which are enumerated up to ndiags 
! in their trailing dimension.
      real psum(*)
      real diagsum(*)
! mesh data, notably ndims, since we don't pass it:
      include 'ndimsdecl.f'
      integer iLs(ndims+1)
      include 'partcom.f'
      real xprior(2*ndims)

      if(x_part(iflag,i).ne.0)then
         inewregion=insideall(ndims,x_part(1,i))
! Alternative to partlocate:
         iu=0
         do id=1,ndims
            ix=int(x_part(ndims*2+id,i))
            iu=iu+(ix-1)*iLs(id)
         enddo
! Cycle through the vertices of the box we are in.
         do ii=0,2**ndims-1
            ii1=ii
            iinc=iu
! Calculate the index and weight of this vertex.
            fac=1.
            do ik=1,ndims
! There may be bit-manipulation routines to accelerate this.
! But likely not by very much since the main cost is divide+mult by 2
               ii2=ii1/2
               ip=ii1-2*ii2
               ii1=ii2
               iinc=iinc+ip*iLs(ik)
               xm=x_part(ndims*2+ik,i)
               ix=int(xm)
               xf=xm-ix
               if(ip.eq.1)then
                  fac=fac*xf
               else
                  fac=fac*(1.-xf)
               endif
            enddo 
! Add to the particle sum the fraction for this vertex.
            psum(1+iinc)=psum(1+iinc)+fac*echarge
            if(ndiags.gt.0)then
! Do something with diagnostics. 
! Assumption here is diags1 n, diags2-4 v, diags5-7 v^2.
               diagsum(1+iinc)=diagsum(1+iinc)+fac
               do k=1,ndiags-1
! Six moments. 3 for v and 3 for v^2. Old uncentered velocity:
!                  temp=x_part(ndims+1+mod(k-1,ndims),i)
! For diagnostics use velocity advanced by half its prior step. 
                  temp=1.5*x_part(ndims+1+mod(k-1,ndims),i)
     $                 -0.5*xprior(ndims+1+mod(k-1,ndims))
                  if(k.gt.ndims)temp=temp*temp
                  kindex=1+iinc+k*iLs(ndims+1)
                  diagsum(kindex)=diagsum(kindex)+fac*temp
               enddo
            endif
         enddo
      endif

      end
!********************************************************************
      subroutine psumperiod(psum,ifull,iaux,iLs)
! If there are periodic particles in any dimension, do the periodic
! exchange sum.
! This is now obsolete because diagperiod is used everywhere despite
! being slightly less efficient.
      real psum(*)
      include 'ndimsdecl.f'
      integer ifull(ndims),iaux(ndims)
      integer iLs(ndims+1)

!      include 'meshcom.f'
      integer iview(3,ndims),indi(ndims)
      integer istart,iend,istride
      parameter (istart=1,iend=2,istride=3)
      include 'partcom.f'
! This data statement serves to silence ftnchek. The first mditerator
! call actually initializes the iview and indi.
      data iview/ndims*0,ndims*0,ndims*0/
      data indi/ndims*0/

      do id=1,ndims
         if(ipartperiod(id).eq.4)then
! Set view to entire array (Offsets indi [0:iaux(id)-1]).
            icomplete=mditerator(ndims,iview,indi,4,iaux)
! Use the general iterator to sum periodically.
! Slice dimension id:
            iview(iend,id)=0
 101        ii=indexcontract(ndims,ifull,indi)
            ib1=1+iLs(id)+ii
            it=1+(iaux(id)-1)*iLs(id)+ii
            it1=1+(iaux(id)-2)*iLs(id)+ii
            ib=1+ii
! Face+1=Face+1 + Other
            psum(ib1)=psum(ib1)+psum(it)
! Other-1=Other-1 + Face
            psum(it1)=psum(it1)+psum(ib)
            if(mditerator(ndims,iview,indi,0,iaux).eq.0)goto 101
         endif
      enddo
      end
!********************************************************************
      subroutine diagperiod(diagsum,ifull,iaux,iLs,ndiags)
! If there are periodic particles in any dimension, do the periodic
! deposit transfer. Also any half-node boundary transfers.
      real diagsum(*)
      include 'ndimsdecl.f'
      integer ifull(ndims),iaux(ndims)
      integer iLs(ndims+1)
!      include 'meshcom.f'
      integer iview(3,ndims),indi(ndims)
      integer istart,iend,istride
      parameter (istart=1,iend=2,istride=3)
      include 'partcom.f'
! This data statement serves to silence ftnchek. The first mditerator
! call actually initializes the iview and indi.
      data iview/ndims*0,ndims*0,ndims*0/
      data indi/ndims*0/

      do id=1,ndims
         if(ipartperiod(id).eq.4)then
!            write(*,*)'Diagperiod',id
! Set view to entire array (Offsets indi [0:iaux(id)-1]).
            icomplete=mditerator(ndims,iview,indi,4,iaux)
! Use the general iterator to transfer contributions periodically.
! Slice dimension id:
            iview(iend,id)=0
 101        ii=indexcontract(ndims,ifull,indi)
            ib1=1+iLs(id)+ii
            it=1+(iaux(id)-1)*iLs(id)+ii
            it1=1+(iaux(id)-2)*iLs(id)+ii
            ib=1+ii
! Diagnostic particles must be transferred, not just added.
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
         else
            if(ipartperiod(id)/64-(ipartperiod(id)/128)*2.eq.1)then
! Lower half-node position boundary. Transfer bottom to bottom+1.
               icomplete=mditerator(ndims,iview,indi,4,iaux)
               iview(iend,id)=0
 102           ii=indexcontract(ndims,ifull,indi)
               ib=1+ii
               ib1=1+iLs(id)+ii
               do k=0,ndiags-1
                  kshift=k*iLs(ndims+1)
                  diagsum(ib1+kshift)= diagsum(ib1+kshift)
     $                 + diagsum(ib+kshift)
                  diagsum(ib+kshift)=0.
               enddo
               if(mditerator(ndims,iview,indi,0,iaux).eq.0)goto 102
            endif
            if(ipartperiod(id)/128-(ipartperiod(id)/256)*2.eq.1)then
! Upper half-node position boundary. Transfer top to top-1.
               icomplete=mditerator(ndims,iview,indi,4,iaux)
               iview(iend,id)=0
 103           ii=indexcontract(ndims,ifull,indi)
               it=1+(iaux(id)-1)*iLs(id)+ii
               it1=1+(iaux(id)-2)*iLs(id)+ii
               do k=0,ndiags-1
                  kshift=k*iLs(ndims+1)
                  diagsum(it1+kshift)= diagsum(it1+kshift)
     $                 + diagsum(it+kshift)
                  diagsum(it+kshift)=0.
               enddo
               if(mditerator(ndims,iview,indi,0,iaux).eq.0)goto 103
            endif
         endif
      enddo
      end

!********************************************************************
! Obsolete version that does no background subtraction.
      subroutine psumtoqnominus(inc,ipoint,indi,ndims,iLs,iused,
     $     psum,rho,volumes,u,rhoinf)
      integer ipoint,inc
      integer indi(ndims),iused(ndims)
      real psum(*),rho(*),volumes(*),u(*),rhoinf

! Silence warnings with spurious access.
      ind=iLs
      ind=iused(1)
      ind=indi(1)
! This routine for use in mditerarg.
! But we iterate only over the inner mesh (not edges).
! Here, t=psum, u=rho, v=volumes, w=u, x=rhoinf 
! Set the density
      ind=1+ipoint
      if(volumes(ind).ge.1.e20)then
! This is outside the region. 
!         if(volumes(ind).ge.1.e30)then
! And all point-charge regions.
         rho(ind)=0.
      else
! Standard case. Use total charge density sum.
         rho(ind)=psum(ind)/(abs(rhoinf)*volumes(ind))
      endif
!      endif
      inc=1
      end
!********************************************************************
      subroutine quasineutral(inc,ipoint,indi,mdims,iLs,iused,
     $     q,u,volumes,uc)
      integer ipoint,inc
      integer indi(mdims),iused(mdims)
      real q(*),u(*),volumes(*),uc(*)
      include 'ndimsdecl.f'
      include 'plascom.f'
!      real dum3,dum4

! Silence warnings with spurious accesses. 
! (Which ought to be optimized away.)
      ind=iLs
      ind=iused(1)
      ind=indi(1)
! This routine for use in mditerarg.
! But we iterate only over the inner mesh (not edges) by virtue of call.
      ind=1+ipoint
! Tell directly from volumes if we are outside:
      if(volumes(ind).le.1.e20)then
! There still might be too small a density, so set a floor for it.
         if(q(ind).lt.1.e-6)q(ind)=1.e-6
         u(ind)=alog(q(ind))-uc(ind)
      else
         u(ind)=phip
      endif
      if(abs(u(ind)).gt.20.)write(*,*)'Large phi',ind,u(ind),phip
      inc=1
      end
!********************************************************************
      subroutine psumtoq(inc,ipoint,indi,ndims,iLs,iused,
     $     psum,rho,volumes,bckgd,rhoinf)
! This version subtracts uniform background in the particle region.
! u(*) is removed.
      integer ipoint,inc
      integer indi(ndims),iused(ndims)
      real psum(*),rho(*),volumes(*),rhoinf

! Silence warnings with spurious access.
      ind=iLs
      ind=iused(1)
      ind=indi(1)
! This routine for use in mditerarg.
! But we iterate only over the inner mesh (not edges).
! Here, t=psum, u=rho, v=volumes, w=u, x=rhoinf 
! Set the density
      ind=1+ipoint
      if(volumes(ind).ge.1.e20)then
! This is outside the region. 
!         if(volumes(ind).ge.1.e30)then
! And all point-charge regions.
         rho(ind)=0.
      else
! Standard case. Use total charge density sum.
         rho(ind)=psum(ind)/(abs(rhoinf)*volumes(ind))-bckgd
      endif
!      endif
      inc=1
      end
