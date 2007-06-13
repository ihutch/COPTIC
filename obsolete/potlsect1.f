c************************************************************
c Obsolete Specific routine for this problem.
      subroutine potlsect1(id,ipm,ndims,indi,fraction,conditions,dp)
c In dimension id, direction ipm, 
c from mesh point at indi(ndims) (zero-based indices, C-style),
c find any intersection of the mesh leg from this point to its neighbor
c with a bounding surface. Return the "fraction" of the leg at which
c the intersection occurs (1 if no intersection), the "conditions" at
c the intersection (irrelevant if fraction=1), and the +ve length
c in computational units of the full leg in dp.
      integer id,ipm,ndims
      integer indi(ndims)
      real fraction,dp
      real conditions(3)
c Equivalence: c1,c2,c3==a,b,c
c The boundary conditions are in the form c1\psi + c2\psi^\prime + c3
c =0.  Where ci=conditions(i).  Thus a fixed potential is (e.g.) c1=1
c c3=-\psi.  A fixed gradient is (e.g.) c2=1 c3=-\psi^\prime. 
c
c In multiple dimensions, if the condition is placed on the gradient
c normal to the surface, in direction n, then for axis direction i
c (assuming an orthogonal set) the value of c2 should become c2_i =
c B/(n.e_i) where e_i is the unit vector in the axis direction.
c
c [If n.e_i=0 this is a pathological intersection.  The mesh leg is
c tangential to the surface, and a normal-gradient BC is irrelevant.
c Therefore, a fraction of 1 should be returned.]
c
c If the surface conditions include non-zero c2, then it is not possible
c to apply the same BC to both sides, without a discontinuity arising.
c Continuity alone should be applied
c on the other (inactive) side of the surface.  A negative value of c2
c is used to indicate that just the appropriate continuity condition is
c to be applied, but c1,c2,c3, should all be returned with the magnitudes
c that are applied to the active side of the surface, with consistently
c reversed signs. (e.g. a=1, b=1, c=-2 -> a=-1, b=-1, c=2)
c
c A fraction of 1 causes all the bounding conditions to be ignored.

      include 'sormesh.f'

      real xc(3),xx(3),xd(3),rc
      logical lfirst
      data lfirst/.true./
      save

c Default no intersection.
      fraction=1
c Some settings that we don't want to do over and over.
      if(lfirst) then
c Initialize a sphere: center, radius, potential
         xc(1)=.5
         xc(2)=.5
         xc(3)=.5
          spotl=2.
         lfirst=.false.
      endif
c----------------------------------------------------------
c Inner constant potential circle.
      rc=.2
c      rc=.45
      call circlesect(id,ipm,ndims,indi,rc,xc,fraction,dp)
      if(fraction.ne.1.)then
         conditions(1)=1.
         conditions(2)=0.
         conditions(3)=-spotl
         return
      endif
c----------------------------------------------------------
c Second concentric sphere with zero potential.
      rd=.43
c Comment the call to disable the section.
c      call circlesect(id,ipm,ndims,indi,rd,xc,fraction,dp)
      if(fraction.ne.1)then
         conditions(1)=1.
         conditions(2)=0.
         conditions(3)=0.
         return
      endif
c
c----------------------------------------------------------
c Concentric sphere with mixed boundary condition.
      rd=.48
c Comment the call to disable the section.
      call circlesect(id,ipm,ndims,indi,rd,xc,fraction,dp)
      if(fraction.ne.1.)then
c Calculate the projection:
c Coordinates of crossing
         do i=1,ndims
c Address of mesh point. 
            ix=indi(i)+ixnp(i)+1
            xx(i)=xn(ix)
            if(i.eq.id)xx(i)=xx(i)+fraction*(xn(ix+ipm)-xn(ix))
c Component of outward vector.
            xd(i)=xx(i)-xc(i)
         enddo
         projection=ipm*(xd(id)/rd)
         if(projection.eq.0.)then
            fraction=1.
         else
c Continuity works for both directions.
            conditions(1)=sign(1./rd,projection)
            conditions(2)=1./projection
            conditions(3)=0.
c         write(*,*)(indi(i),xd(i),i=1,ndims),id,projection,
c     $           conditions(2)
c Zero outside, rather than continuity alternative
            if(.true. .and. projection.lt.0.)then
               conditions(1)=1.
               conditions(2)=0.
               conditions(3)=0.
            endif
            return
         endif
      endif
c----------------------------------------------------------
c Default return for no intersection.
c Address of mesh point. 
      ix=indi(id)+ixnp(id)+1
      dp=abs(xn(ix+ipm)-xn(ix))
      end
