! Optional test main.
!      call testdiststat
!      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Obtain the statistical properties: total, centroid, and variance,
! of a scalar quantity treated as a (signed) weight, but nonzero only
! when ge qmin, over a multidimensional mesh specified in ixnp,xn.
      subroutine diststat(ndims,ifull,iused,quantity,qmin
     $     ,ixnp,xn,tot,cent,var)
      integer ndims,ifull(ndims),iused(ndims),ixnp(*)
      real quantity(*),xn(*)
      real tot,cent(ndims),var(ndims)
      integer mditer
      external mditer

! Initialize cumulators.
      tot=0.
      do id=1,ndims
         cent(id)=0.
         var(id)=0.
      enddo

! Iterate over the whole array
      icomplete=mditer(ndims,ifull,iused,index)
 1    if(quantity(index).ge.qmin)then
         tot=tot+quantity(index)
         do id=1,ndims      ! iused(id) becomes the indi during mditer.
            cent(id)=cent(id)+xn(1+ixnp(id)+iused(id))*quantity(index)
            var(id)=var(id)+xn(1+ixnp(id)+iused(id))**2*quantity(index)
         enddo
      endif
!      write(*,'(i6,3i4,3f8.3)')index,iused,xn(1+ixfirst+iused)
      if(mditer(ndims,ifull,iused,index).eq.0)goto 1

! Finalize:
      do id=1,ndims
         cent(id)=cent(id)/tot
         var(id)=var(id)/tot-cent(id)**2
      enddo

      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testdiststat
      integer ndims,nx,ny,nz
      parameter (ndims=3,nx=20,ny=30,nz=40,nxy=nx+ny,nxyz=nx+ny+nz)
      integer ifull(ndims),iused(ndims),ixnp(ndims+1)
      real xn(nx+ny+nz)
      real quantity(nx,ny,nz)
      real tot,cent(ndims),var(ndims)
      data ixnp/0,nx,nxy,nxyz/
      data xmin/-8./xmax/8./ymin/-10./ymax/10./zmin/-12./zmax/12./
      data ifull/nx,ny,nz/iused/nx,ny,nz/
      data sx/2./sy/2./sz/3./
      data cx/0./cy/1./cz/0./

      do i=1,nx
         xn(i)=xmin+(xmax-xmin)*(i-1.)/(nx-1.)
      enddo
      do i=1,ny
         xn(nx+i)=ymin+(ymax-ymin)*(i-1.)/(ny-1.)
      enddo
      do i=1,nz
         xn(nx+ny+i)=zmin+(zmax-zmin)*(i-1.)/(nz-1.)
      enddo

      do ix=1,nx
         x=xn(ix)
         do iy=1,ny
            y=xn(nx+iy)
            do iz=1,nz
               z=xn(nx+ny+iz)
               quantity(ix,iy,iz)=
     $              exp((-((x-cx)/sx)**2
     $              -((y-cy)/sy)**2-((z-cz)/sz)**2)/2)
            enddo
         enddo
      enddo
      write(*,*)ifull,iused
      write(*,'(10f8.3)')xn
      
      call diststat(ndims,ifull,iused,quantity,0.
     $     ,ixnp,xn,tot,cent,var)

      write(*,'(a,f8.2,a,3f8.4,a,3f8.4)')'tot=',tot,' cent=',cent
     $     ,' sqrtvar=',sqrt(var)
      end
! Notes:
! The hole implemented in coptic is \propto \exp[-(r/holerad)^2), which 
! has SD=sqrt(var)=holerad/sqrt(2), because \exp[-(r/s)^2/2] has SD=s.
! The SD of \sech^4(x/4) is 2.2714 from analytic expressions of the 
! integrals obtained from Wolfram.        
