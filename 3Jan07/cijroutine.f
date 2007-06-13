c************************************************************
c Specific routine for this problem. Others general.
      subroutine potlsect(id,ipm,ndims,indi,fraction,potential,dp)
c In dimension id, direction ipm, 
c from mesh point at indi(ndims) (zero-based indices, C-style),
c find any intersection of the mesh leg from this point to its neighbor
c with a potential surface. Return the fraction of the leg at which
c the intersection occurs (1 if no intersection), the potential
c of the intersection (irrelevant if fraction=1), and the +ve distance
c in computational units to the intersection (or neighbour) in dp.
      integer id,ipm,ndims
      integer indi(ndims)
      real fraction,potential,dp
      include 'sormesh.f'

      real xc(3),rc
      logical lfirst
      data lfirst/.true./
      save

      if(lfirst) then
c Initialize a sphere: center, radius, potential
         xc(1)=.5
         xc(2)=.5
         xc(3)=.5
         rc=.2
         spotl=2.
         lfirst=.false.
      endif

      call circlesect(id,ipm,ndims,indi,rc,xc,fraction,dp)
      potential=spotl
      if(fraction.ne.1)return
      
c Address of mesh point. 
      ix=indi(id)+ixnp(id)+1
c No intersection.
c Standard return.
      dp=abs(xn(ix+ipm)-xn(ix))
c      write(*,'(4i4,2f8.4)')id,ipm,indi(id),ixnp(id),xn(ix),dp
      potential=0.
      fraction=1.
      end
c****************************************************************
      subroutine circlesect(id,ipm,ndims,indi,rc,xc,fraction,dp)
      integer id,ipm,ndims
      integer indi(ndims)
      real fraction,dp
      include 'sormesh.f'

      real xc(3),rc
c      logical lfirst
c      data lfirst/.true./
      save
c The above save should not make any difference but it does.
      A=0
      B=0
      C=-rc**2
      D=-rc**2
      do i=1,ndims
         ix1=indi(i)+ixnp(i)+1
         x1=xn(ix1)-xc(i)
         x2=x1
         if(i.eq.id)then
            ix2=ix1+ipm
            x2=xn(ix2)-xc(i)
            A=A+(x2-x1)**2
            B=B+x1*(x2-x1)
            dp=abs(x2-x1)
         endif
         C=C+x1**2
         D=D+x2**2
      enddo
      fraction=1.
      if(D.ne.0. .and. D*C.le.0.)then
         if(B.ge.0. .and. A*C.le.0) then
            fraction=(-B+sqrt(B*B-A*C))/A
         elseif(B.lt.0. .and. A*C.ge.0)then
            fraction=(-B-sqrt(B*B-A*C))/A
         endif
c That should exhaust the possibilities.
         dp=fraction*dp
c         write(*,'(4i4,3f8.4)')id,ipm,indi(id),ixnp(id),xn(ix),
c     $        fraction,dp
      endif
      end

c****************************************************************
c Combined cijroutine and cijobj.
c Routine for setting cij, which is used by the general mditerate.
      subroutine cijroutine(inc,ipoint,indm1,ndims,iused,cij)
c This routine sets the object pointer in cij if there is an object
c crossing next to the the point at indi(ndims) in the 
c dimension id (plus or minus), adjusts the cij values,
c and sets the object data into the place pointed to in obj. 
      integer mdims
      parameter (mdims=10)
c Effective index in dimension, c-style (zero based)
      integer indm1(mdims)
c Shifted to correct for face omission.
      integer indi(mdims)
      integer ndims
      integer iused(ndims)
      real cij(*)

c      subroutine cijobj(id,ndims,indi,ipoint,cij,iset)
      include 'objcom.f'
      real dpm(2)
      real tiny
      parameter (tiny=1.e-15)

      do id=1,ndims
c  We are using shifted arrays; indi is here the true index.
         indi(id)=indm1(id)+1
      enddo
c ipoint here is the offset (zero-based pointer)
c The cij address is to data 2*ndims+1 long
      icb=2*ndims+1
c Object pointer defaults to zero.
      cij(icb*ipoint+ndims*2+1)=0.
c Track whether obj is set for this point.
      iset=0
c Iterate over dimensions.
      do id=1,ndims
         ifound=0
c Find the intersection, if any, of this mesh branch with the objects.
c Return the fraction of the mesh distance of the intersection, and 
c the value of the potential to be attached to it. Fraction=1. is the 
c default case, where no intersection occurs. However, we need to know
c the fractions for opposite directions before we can calculate the dds.
         icb2=2*icb
c For each direction in this dimension,
c Determine whether this is a boundary point: adjacent a fraction ne 1.
         do i=1,2
            iobj=4*(id-1)+2*i-1
            ipm=1-2*(i-1)
            call potlsect(id,ipm,ndims,indi,fraction,potential,dpm(i))
            if(fraction.ne.1.)then
               ifound=ifound+1
               if(iset.eq.0)then
c Start object data for this point if not already started.
                  objindex=objindex+1
                  if(objindex.gt.lobjmax) then
                     write(*,*)'cijobj: objindex overflow',objindex
                     stop 
                  endif
c Initialize the object data. 1s for frac,potl, 0 for diag,potterm.
                  do j=1,nobj-2
                     obj(j,objindex)=1.
                  enddo
                  obj(nobj-1,objindex)=0.
                  obj(nobj,objindex)=0.
                  cij(ipoint*icb+ndims*2+1)=objindex
                  iset=iset+1
               endif
c The object data consists of data enumerated as
c ndims*(2=forward/backward from point)*(2=fraction,potential)
c + diagonal + potential terms.
               obj(iobj,objindex)=fraction
               obj(iobj+1,objindex)=potential
            endif
         enddo
c Set the coefficients for each direction.
         do i=1,2
c Calculate the differentials avoiding divide by zero
            im=mod(i,2)+1
            dp=dpm(i)+tiny
            dm=dpm(im)
            deff=dp
            coef=(2./(deff*(dp+dm)))
            if(ifound.gt.0)then
c This is a boundary point.
               iobj=4*(id-1)+2*i-1
               if(obj(iobj,objindex).ne.1)then
c                  write(*,*)'Object intersection'
c We intersected an object in this direction. Adjust diagonal term
                  obj(icb2-1,objindex)=obj(icb2-1,objindex)
     $                 + coef
c Adjust potential sum
                  obj(icb2,objindex)=obj(icb2,objindex)
     $                 + coef*obj(iobj+1,objindex)
                  cij(icb*ipoint+2*(id-1)+i)=0.
               else
c We did not intersect an object in this direction.
                  cij(icb*ipoint+2*(id-1)+i)=coef
               endif
            else
c Standard non-boundary setting
               cij(icb*ipoint+2*(id-1)+i)=coef            
            endif
         enddo
c End of dimension iteration.
      enddo
c     Return increment of 1
      inc=1
      end
c********************************************************************
c Mostly for testing.
c Assumed 3-D routine, plots representation of the cij/obj data.
      subroutine  cijplot(ndims,ifull,iuds,cij)
      integer ndims
      parameter (mdims=10)
      integer ifull(mdims),iuds(ndims)
      real cij(ndims*2+1,ifull(1),ifull(2),ifull(3))
      include 'objcom.f'
      include 'sormesh.f'
      real xx(3),xt(3)

 51   continue
c     call pfset(3)
      call setcube(.2,.2,.2,.5,.4)
      call pltinit(0.,1.,0.,1.)
      call geteye(x2,y2,z2)
      call pltinit(0.,1.,0.,1.)
      call scale3(0.,1.,0.,1.,0.,1.)
      call trn32(0.,0.,0.,x2,y2,z2,1)
      xb=0
      yb=0
      zb=0
      if(x2.ge.0)xb=1
      if(y2.ge.0)yb=1
      if(z2.ge.0)zb=1
      icorner= (2*zb-1)*( (1 +3*yb) + (1 - 2*yb)*xb )

      call cubed(icorner)
      call axproj(icorner)

      do i=1,ifull(1)
         do j=1,ifull(2)
            do k=1,ifull(3)
               iobj=cij(2*ndims+1,i,j,k)
               if(iobj.ne.0)then
c Draw joins from the point to the fraction distance.
                  xx(1)=xn(ixnp(1)+i)
                  xx(2)=xn(ixnp(2)+j)
                  xx(3)=xn(ixnp(3)+k)
                  call vec3w(xx(1),xx(2),xx(3),0)
                  do id=1,ndims
                     do ipm=1,2
                        ispm=1-2*(ipm-1)
                        no=(id-1)*4+(ipm-1)*2+1
                        call color((no+1)/2)
                        frac=obj(no,iobj)
                        if(frac.ne.1.)then
                           xt(1)=xx(1)
                           xt(2)=xx(2)
                           xt(3)=xx(3)
                           if(id.eq.1)then
                              xt(1)=xx(1)+frac*
     $                             (xn(ixnp(1)+i+ispm)-xx(1))
                           elseif(id.eq.2)then
                              xt(2)=xx(2)+frac*
     $                             (xn(ixnp(2)+j+ispm)-xx(2))
                           else
                              xt(3)=xx(3)+frac*
     $                             (xn(ixnp(3)+k+ispm)-xx(3))
                           endif
c                           write(*,*)no,id,frac,ipm,xt
                           call vec3w(xt(1),xt(2),xt(3),1)
                           call vec3w(xx(1),xx(2),xx(3),0)
                        endif
                     enddo
                  enddo
               endif
            enddo
         enddo
      enddo

      if(ieye3d().ne.0)goto 51
c      call pltend()

      end
