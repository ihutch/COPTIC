!**********************************************************************
! Convert from previously set imeshstep, xmeshpos to ixnp.
! There must be at least 2 non-zero and monotonically increasing
! values starting imeshstep.
! Return values of iuds(ndims) per constructed mesh.
! Also return the value of rs equal to half the largest mesh box side
! length.
! Version constructs xmstart(id) to xmend(id) according to particle domain.
      subroutine meshconstruct(mdims,iuds,ifull,ipartperiod,rs)
      integer iuds(mdims),ifull(mdims),ipartperiod(mdims)
      real rs
      integer iof
      include 'ndimsdecl.f'
      include 'meshcom.f'

      iof=0
      rs=0.
      do id=1,ndims
! Pointer to start of vector.
         ixnp(id)=iof
         xn(iof+1)=xmeshpos(id,1)
         iof=iof+1
         do iblk=1,nspec_mesh-1
            msteps=imeshstep(id,iblk+1)-imeshstep(id,iblk)
            if(msteps.le.0) goto 11
! Mesh data.
            do i=1,msteps
               xn(i+iof)=((msteps-i)*xmeshpos(id,iblk)+
     $              i*xmeshpos(id,iblk+1))/msteps
            enddo
            iof=iof+msteps
         enddo
         iblk=nspec_mesh-1
 11      continue
         if(iblk.le.1) write(*,*)'Too few mesh steps. Dimension',id
         if(imeshstep(id,iblk).gt.ifull(id))then
            write(*,'(a,i1,a,i4,a,3i5)')'Meshconstruct ERROR: Meshpos('
     $           ,id,')=',imeshstep(id,iblk),'  too large for ifull='
     $           ,(ifull(k),k=1,ndims)
            stop
            endif
! Set iuds according to specified mesh
         iuds(id)=imeshstep(id,iblk)
!         write(*,'(a,i3,10f8.3)')
!     $        ' Meshspec',id,(xmeshpos(id,kk),kk=1,iblk)
         if(ipartperiod(id).eq.4.or.ipartperiod(id).eq.5)then
! Half-cell-position mesh ends (but id+1 mesh is not yet set).
            xmeshstart(id)=(xn(ixnp(id)+1)+xn(ixnp(id)+2))*.5
            xmeshend(id)=(xn(ixnp(id)+iuds(id))
     $           +xn(ixnp(id)+iuds(id)-1))*.5
         else
! Whole-cell-position ends
            xmeshstart(id)=xmeshpos(id,1)
            xmeshend(id)=xmeshpos(id,iblk)
! Unless higher byte settings used.
            if(ipartperiod(id)/64-(ipartperiod(id)/128)*2.eq.1)then
               xmeshstart(id)=(xn(ixnp(id)+1)+xn(ixnp(id)+2))*.5
            endif
            if(ipartperiod(id)/128-(ipartperiod(id)/256)*2.eq.1)then
               xmeshend(id)=(xn(ixnp(id)+iuds(id))
     $              +xn(ixnp(id)+iuds(id)-1))*.5
            endif
         endif
         rsi=0.5*abs(xmeshend(id)-xmeshstart(id))
         if(rsi.gt.rs)then 
            rs=rsi
         endif
      enddo
      ixnp(ndims+1)=iof
!      write(*,*)'Meshconstructed',xmeshstart,xmeshend
!      write(*,*)(ixnp(k),k=1,ndims+1)
!      write(*,*)(xn(k),k=1,iof)
      do id=1,ndims
! Now construct the interpolation index arrays from xmeshstart to xmeshend
         ioff=ixnp(id)
         do ipi=1,ipilen
            y=xmeshstart(id)+(ipi-1.)/(ipilen-.99999)
     $           *(xmeshend(id)-xmeshstart(id))
! Recover with ipi=(y-xmeshstart)/(xmeshend-xmeshstart)*(ipilen-1.00001)+1
            ix=interp(xn(ioff+1),ixnp(id+1)-ioff,y,x)
            if(ix.le.0)write(*,*)'Meshconstruct interp error',y,x,ix
            iposindex(ipi,id)=ix
         enddo
      enddo
!      write(*,'(20i3)')iposindex

      end
!**************************************************************
      subroutine meshshift()
      include 'ndimsdecl.f'
      include 'meshcom.f'
      include 'objcom.f'
      real shift
      parameter (shift=3.e-5)
!      write(*,'(a)')'Shifted mesh' 
      do id=1,ndims
         do j=1,nspec_mesh
            xmeshpos(id,j)=(xmeshpos(id,j)*(1.-shift)+.1*shift)
!            write(*,'(f8.4,$)')xmeshpos(id,j)
         enddo
      enddo
! Reinitialize the intersection counter.
      oi_cij=0
      end
!***************************************************************
