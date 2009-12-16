c Interpolations.
      function boxinterp(ndm1,f,flags,d)
      real f(ndm1,ndm1),d(ndm1)
      integer flags(ndm1)
      boxinterp=0.
      if(ndm1.eq.1)then
         boxinterp=box1interp(f,flags,d)
      elseif(ndm1.eq.2)then
         boxinterp=box2interp(f,flags,d)
      else
         write(*,*)'Boxinterp Error. Unknown dimensionality:',ndm1
      endif
      end
c****************************************************************
      function box2interp(f,flags,d)
c 
c Given values of a function on up to four points adjacent on a 
c 2-D cartesian mesh: f00,f01,f10,f11; and given the fractional
c distances: d1,d2 in the two dimensions (between 0 and 1),
c interpolate the value at d1,d2 using box interpolation. 
c 
      real f(2,2)
      real d(2)
c If one or more of the f-values is absent, indicated by zero flags,
c then construct fall-back interpolation that minimizes curvature.
      integer flags(2,2)

c      write(*,*)'Box2interp flags:',flags
      if(flags(1,1).ne.0)then
c      if(.true.)then
c  f00 is not absent:
         if(flags(2,2).eq.0)then
            if(flags(1,2).ne.0)then
               if(flags(2,1).ne.0)then 
                  f(2,2)=f(2,1)+f(1,2)-f(1,1)
               else
                  f(2,1)=f(1,1)
                  f(2,2)=f(1,2)
               endif
            else
               if(flags(2,1).ne.0)then 
                  f(1,2)=f(1,1)
                  f(2,2)=f(2,1)
               else
                  f(1,2)=f(1,1)
                  f(2,1)=f(1,1)
                  f(2,2)=f(1,1)
               endif
            endif
         else
            if(flags(2,1).ne.0)then 
               if(flags(1,2).ne.0)then
c All present
               else
                  f(1,2)=f(1,1)+f(2,2)-f(2,1)
               endif
            else
               if(flags(1,2).ne.0)then
                  f(2,1)=f(1,1)+f(2,2)-f(1,2)
               else
                  f(2,1)=(f(1,1)+f(2,2))/2.
                  f(1,2)=(f(1,1)+f(2,2))/2.
               endif
            endif
         endif
      else
c f00 is absent
c         write(*,*)'f00 absent track',flags
         if(flags(1,2).ne.0)then
c f01 present
            if(flags(2,1).eq.0)then
               if(flags(2,2).eq.0)then
c all except 01 absent
                  f(1,1)=f(1,2)
                  f(2,1)=f(1,2)
                  f(2,2)=f(1,2)
               else
c 01 and 11 only present
                  f(1,1)=f(1,2)
                  f(2,1)=f(2,2)                  
               endif
            else
c Both 01 and 10 present
               if(flags(2,2).eq.0)then
                  f(1,1)=(f(1,2)+f(2,1))/2.
                  f(2,2)=(f(1,2)+f(2,1))/2.
               else
                  f(1,1)=f(1,2)+f(2,1)-f(2,2)
               endif
            endif
         else
c f01 absent
            if(flags(2,1).eq.0)then
               if(flags(2,2).eq.0)then
                  write(*,*)'**** Box no vertices!'
                  stop
               else
c all except 11 absent
                  f(1,1)=f(2,2)
                  f(1,2)=f(2,2)
                  f(2,1)=f(2,2)
               endif
            else
c f01 absent but f10 present
               if(flags(2,2).eq.0)then
c all except 10 absent
                  f(1,1)=f(2,1)
                  f(1,2)=f(2,1)
                  f(2,2)=f(2,1)
               else
                  f(1,1)=f(2,1)
                  f(1,2)=f(2,2)
               endif
            endif
         endif
      endif
c Now we have filled in all the values appropriately. 
c Interpolate:

      box2interp=f(1,1)+d(1)*(f(2,1)-f(1,1))+d(2)*(f(1,2)-f(1,1))
     $     +d(1)*d(2)*(f(2,2)-f(2,1)-f(1,2)+f(1,1))

      end
c************************************************************
      function box1interp(f,flags,d)
c 
c Given up to two f-values in a 1-D array, and a fraction d
c between them, interpolate. f00 is always present.
      real f(2)
      integer flags(2)
      real d

      if(flags(2).ne.0)then
         box1interp=(1.-d)*f(1)+d*f(2)
      else
         box1interp=f(1)
      endif
      end
c********************************************************************
c Given a monotonic function Q(x)
c on a 1-D grid x=1..nq, solve Q(x)=y for x with interpolation.
c If return is 0, then the y is outside Q's range or other error.
c The function returns the integer part of x.
      function interp(Q,nq,y,x)
      real Q(nq)
      integer nq
      real y,x
c
      integer iqr,iql,iqx
      real Qx,Qr,Ql

c      write(*,*)'nq=',nq
      interp=0
      Ql=Q(1)
      Qr=Q(nq)
      iql=1
      iqr=nq
      if((y-Ql)*(y-Qr).gt.0.) then
c Value is outside the range.
         x=0
         return
      endif
 200  if(iqr-iql.eq.1)goto 210
      iqx=(iqr+iql)/2
      Qx=Q(iqx)
c      write(*,*)y,Ql,Qx,Qr,iql,iqr
      if((Qx-y)*(Qr-y).le.0.) then
         Ql=Qx
         iql=iqx
      else
         Qr=Qx
         iqr=iqx
      endif
      goto 200
 210  continue
c Now iql and iqr, Ql and Qr bracket Q
      if(Qr-Ql.ne.0.)then
         xpart=(y-Ql)/(Qr-Ql)
         x=xpart+iql
         interp=iql
      else
         x=iql
         write(*,*)'****** Error!: interp coincident points'
      endif
      end
c**********************************************************************
c Convert index to multidimensional indices.
      subroutine indexexpand(ndims,ifull,index,ix)
c On entry index is the zero-based pointer to the position in the 
c array whose full dimensions are ifull(ndims). 
c On exit ix contains the corresponding (ndims)-dimensional indices.
      integer ndims
      integer ifull(ndims),ix(ndims)
      integer index

      ind=index
      do i=1,ndims-1
         ind2=ind/ifull(i)
         ix(i)=ind-ind2*ifull(i)+1
         ind=ind2
      enddo
      ix(ndims)=ind+1
      if(ind.gt.ifull(3)) write(*,*)'indexexpand index too big',index
      end
c********************************************************************

c*************************************************************
      function gradinterp(um,u0,up,id,icp,xm,dx0,dx1)
c Interpolate the gradient of potential,
c      real um,u0,up
c From its passed origin to adjacent points, which for this dimension
c      integer id
c The pointer to object data is
c      integer icp
c The fractional distance xm for interpolation (-1<xm<1)
c relative to the central point in the dimension id is
c      real xm
c The step sizes to adjacent points are
c      real dx0,dx1
c The object boundary data is in objcom
c      include 'objcom.f'
c The mesh vectors are in sormesh may not be needed.
c      include 'meshcom.f'
c
c The interpolation is in the form:
c u' = 2(x+dxf0/2)/(dxf0+dxf1) (up-u0)/dxd1 +
c      2(dxf1/2-x)/(dxf0+dxf1) (u0-um)/dxd0
c where u0 is the value of u at point.
c up and um are the effective forward and backward values of u
c except that they are replaced
c for object boundaries by -C/A. 
c dxp1, dxp0 are the distances to the adjacent points, or if an
c object boundary intervenes, the boundary. dxp1=dx1*fraction ...
c The boundary condition is Au + Bu' + C =0, where u' denotes the
c gradient away from the central point in this condition.
c dxf1= dxp1 (dxp1+2B/A)/(dxp1+B/A) defines the gradient control
c position (and similarly for dxf0), so that the gradient is
c controlled at dxf1/2 (=dxp1/2 if B=0, or dxp1 if A=0) 
c dxd1 = dxp1 + B/A, dxd0=dxp0 + B/A, are the divisors.
c Fractions, which give dxp1, must not be exactly zero.
      real um,u0,up
      integer id
      integer icp
      real xm
      real dx0,dx1

      include 'objcom.f'

      if(xm.lt.0)then
         x=xm*dx0
      else
         x=xm*dx1
      endif

      if(icp.ne.0)then
c There is object data. Use object boundary interpolations.
c Addresses of this dimension among objects:
         icd1=2*(id-1)*ndata_sor+1
         icd0=icd1+ndata_sor
         fraction=dob_sor(icd0,icp)
         if(fraction.lt.1. .and. fraction.ge.0.)then
            dxp0=fraction*dx0
            boa=dob_sor(icd0+1,icp)
            coa=dob_sor(icd0+2,icp)
            um=-coa
            dxf0=dxp0*(dxp0+2*boa)/(dxp0+boa)
            dxd0=dxp0+boa
         else
            dxp0=dx0
            dxf0=dx0
            dxd0=dx0
         endif
         fraction=dob_sor(icd1,icp)
         if(fraction.lt.1. .and. fraction.ge.0.)then
            dxp1=fraction*dx1
            boa=dob_sor(icd1+1,icp)
            coa=dob_sor(icd1+2,icp)
            up=-coa
            dxf1=dxp1*(dxp1+2.*boa)/(dxp1+boa)
            dxd1=dxp1+boa
         else
            dxp1=dx1
            dxf1=dx1
            dxd1=dx1
         endif
      else
c         write(*,*)'x,dx0,dx1,um,u0,up',x,dx0,dx1,um,u0,up
c Non-boundary interpolation.
         gradinterp= (2.*x+dx0)/(dx0+dx1) * (up-u0)/dx1
     $        +(dx1-2.*x)/(dx0+dx1) * (u0-um)/dx0
         return
      endif
c      write(*,'(a,/,i4,10f6.2)')
c     $     'icp,   x,   dxf0, dxf1, dxd0, dxd1, dx0,  dx1, um, u0, up',
c     $     icp,x,dxf0,dxf1,dxd0,dxd1,dx0,dx1,um,u0,up
c      write(*,*)'fraction,boa,coa',fraction,boa,coa
c General interpolation
      gradinterp= (2.*x+dxf0)/(dxf0+dxf1) * (up-u0)/dxd1
     $     +(dxf1-2.*x)/(dxf0+dxf1) * (u0-um)/dxd0

      end

c*******************************************************************
c This routine cannot use fortran-bounds-checking because it needs
c to access array elements earlier than that passed.
c
c Returns the gradient that would occur in the region iregion at the
c position relative to the passed arrays u,cij, given by node fraction
c xf in the gradient (idf) dimension.  If position is in region iregion,
c this is simple. If not, the adjacent point in the direction of the
c gradient is examined, and if it is in region iregion, extrapolation
c from that point to xf is used. If not, ix is returned as 99 indicating
c the interpolation has failed.

      subroutine gradlocalregion(cij,u,idf,icinc,iuinc,xn,
     $     xf,uprime,iregion,ix,xm)

c The cij, and u arrays, start at the element corresponding to the
c nearest node in the gradient direction.
c
c Here cij includes the offset to the object pointer so that in effect
c cij is simply the array of pointers to object data but is real.  
c
c Thus the calls generally pass: cij(nd2+1,ium2(1),ium2(2),ium2(3))
c ,u(ium2(1),ium2(2),ium2(3)). But in fact this routine should work for
c general number of dimensions.  The increments to adjacent cij,u-values
c are icinc, iuinc.
c 
c xn is the array of positions in the idf dimension, relative to point.
c (i.e. we pass xn(ixnp(idf)+ix)) needed for gradient evaluation.

      real cij(*)
      real u(*)
      include 'objcom.f'
c Not here      include 'meshcom.f'
c The direction in which we are interpolating.
      integer idf
c The increments to adjacent values in the dimension idf of cij, and u
      integer icinc,iuinc
c The position array in the idf dimension.
      real xn(*)
c The position fraction INPUT
      real xf

c The value of gradient: OUTPUT
      real uprime
c ix and xm are returned the index and fraction of interpolation.

c xn is the position array for each dimension arranged linearly.
      ix=1
      xm=xf
c Pointer to object data,
c      icp0=cij(ix)
      icp0=int(cij(1))
c If we are in wrong region, try to correct:
      if(icp0.ne.0)then
c This section is only 10% of routine cost.
         if(idob_sor(iregion_sor,icp0).ne.iregion)then
c         write(*,*)'Incorrect region',iregion,idob_sor(iregion_sor,icp0)
c     $        ,icp0,ix,xm
            jpm=1
            if(xm.lt.0.)jpm=-1
            ix=ix+jpm
            icp1=int(cij(1+(ix-1)*icinc))
c Old and seemingly incorrect version:
c            if(icp1.eq.0 .or. idob_sor(iregion_sor,icp1).ne.iregion)then
            if(icp1.ne.0.and.idob_sor(iregion_sor,icp1).ne.iregion)then
               ix=99
               uprime=0.
               return
            endif
            xm=xm-jpm
            icp0=icp1
c         write(*,*)'Base Position Adjusted ix',ix,jpm,xm,iregion
         endif
      endif
c Distance forward and backward along idf-dimension to adjacent
      dx1=xn(ix+1)-xn(ix)
      dx0=xn(ix)-xn(ix-1)

c Values of u at the points to be interpolated.
c      u0=u(1+(ix-1)*iuinc)
c      up=u(1+ix    *iuinc)
c      um=u(1+(ix-2)*iuinc)
      ixiu=1+(ix-1)*iuinc
      u0=u(ixiu)
      up=u(ixiu+iuinc)
      um=u(ixiu-iuinc)
c Do interpolation using extrapolation information pointed to by icp0.
      if(icp0.eq.0)then
c Short-cut 1: an ordinary point don't call the full routine
         if(abs(dx1-dx0).lt.1.e-6*dx0)then
c Short-cut 2: for uniform mesh:
c (saves ~25% of routine for uniform mesh even including test):
c A uprime=0. test reduces the time in this routine by about 25%.
c So this evaluation is about 25%.
            uprime= ((2.*xm+1.) * (up-u0)
     $           +(1.-2.*xm) * (u0-um))/(dx0+dx1)
         else
            if(xm.lt.0)then
               x=xm*dx0
            else
               x=xm*dx1
            endif
c         uprime= (2.*x+dx0)/(dx0+dx1) * (up-u0)/dx1
c     $        +(dx1-2.*x)/(dx0+dx1) * (u0-um)/dx0
            uprime= ((2.*x+dx0) * (up-u0)/dx1
     $           +(dx1-2.*x) * (u0-um)/dx0)/(dx0+dx1)
         endif
      else
         uprime=gradinterp(um,u0,up,idf,icp0,xm,dx0,dx1)
      endif

      end

c***************************************************************
cXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c****************************************************************
c The routines after this are unused.
c****************************************************************
cXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c*************************************************************
c Routine for interpolating gradient with correct node assignment.
c Returns the gradient corresponding to the region the assigned point
c is actually in.
      subroutine gradinterpcorrect(cij,u,idf,icinc,iuinc,
     $     xprime,uprime,ix,xm)

c The cij, and u arrays, starting at the 1st element in the idf 
c direction and the appropriate element in the other dimensions.
c (Which for cij includes the offset to the object pointer).
c Thus the calls generally pass: cij(nd2+1,ium2(1),ium2(2),ium2(3))
c ,u(ium2(1),ium2(2),ium2(3)), with ium2(idf)=1. But in fact this
c routine should work for general number of dimensions.
      real cij(*)
      real u(*)
      include 'objcom.f'
      include 'meshcom.f'
c The direction in which we are interpolating.
      integer idf
c The increments to adjacent values in the direction idf of cij, and u
      integer icinc,iuinc
c The position in the interpolated direction: INPUT
      real xprime
c The value of gradient: OUTPUT
      real uprime
c ix and xm are returned the index and fraction of interpolation.
c (Not much use except for diagnostics) 

c Dimension along which we are interpolating: idf
c Offset to start of idf position array.
      ioff=ixnp(idf)
c xn is the position array for each dimension arranged linearly.
c Find the index of xprime in the array xn:
      ix=interp(xn(ioff+1),ixnp(idf+1)-ioff,xprime,xm)
      ix=nint(xm)
      xm=xm-ix
 701  continue
c Pointer to object data,
c      icp0=cij(ndims*2+1,ium(1),ium(2),ium(3))
      icp0=int(cij(1+(ix-1)*icinc))
c               write(*,*)'i,ix,xm,icp0',i,ix,xm,icp0
c Distance forward and backward along idf-dimension to adjacent
      dx1=xn(ix+ioff+1)-xn(ix+ioff)
      dx0=xn(ix+ioff)-xn(ix+ioff-1)
c Here we use fraction if icp0!=0 to decide if we have chosen the
c correct mesh node to extrapolate from. If the distance of the point
c from the mesh node is greater than fraction, switch to the further
c node.
      if(icp0.ne.0)then
         if(xm.ge.0.)then
            jpm=1
         else
            jpm=-1
         endif
         iobj=ndata_sor*(2*(idf-1)+(1-jpm)/2)+1
         fraction=dob_sor(iobj,icp0)
c                  write(*,*)'Boundary:',icp0,iobj,jpm,fraction
c                  if(.false. .and. fraction.ge.0. .and.
         if(fraction.ge.0. .and.
     $        fraction.lt.1. .and. abs(xm).gt.fraction)then
            write(*,*)'Adjusting ix',ix,jpm,xm,fraction,
     $           ix+jpm
            ix=ix+jpm
            xm=xm-jpm
            goto 701
         endif
      endif
c Values of u at the points to be interpolated.
      u0=u(1+(ix-1)*iuinc)
      up=u(1+ix    *iuinc)
      um=u(1+(ix-2)*iuinc)
c Do interpolation.
      uprime=gradinterp(um,u0,up,idf,icp0,xm,dx0,dx1)

      end
c*******************************************************************
c Routine for interpolating gradient with correct node assignment.
c Returns the gradient that would occur in the region iregion at
c position xprime. If xprime is in region iregion, this is simple. If
c not, the adjacent point in the direction of the gradient is examined,
c and if it is in region iregion, extrapolation from that point to
c xprime is used. If not, ix is returned as negative indicating the
c interpolation has failed.
      subroutine gradinterpregion(cij,u,idf,icinc,iuinc,
     $     xprime,uprime,iregion,ix,xm)

c The cij, and u arrays, start at the 1st element in the idf 
c direction and the appropriate element in the other dimensions.
c (Which for cij includes the offset to the object pointer so that
c in effect cij is simply the array of pointers to object data).
c Thus the calls generally pass: cij(nd2+1,ium2(1),ium2(2),ium2(3))
c ,u(ium2(1),ium2(2),ium2(3)), with ium2(idf)=1. But in fact this
c routine should work for general number of dimensions.
c The increments to adjacent cij,u-values are icinc, iuinc.
      real cij(*)
      real u(*)
      include 'objcom.f'
      include 'meshcom.f'
c The direction in which we are interpolating.
      integer idf
c The increments to adjacent values in the direction idf of cij, and u
      integer icinc,iuinc
c The position in the interpolated direction: INPUT
      real xprime
c The value of gradient: OUTPUT
      real uprime
c ix and xm are returned the index and fraction of interpolation.

c Dimension along which we are interpolating: idf
c Offset to start of idf position array.
      ioff=ixnp(idf)
c xn is the position array for each dimension arranged linearly.
c Find the index of xprime in the array xn:
      ix=interp(xn(ioff+1),ixnp(idf+1)-ioff,xprime,xm)
      ix=nint(xm)
      xm=xm-ix
c Pointer to object data,
      icp0=int(cij(1+(ix-1)*icinc))
c               write(*,*)'ix,xm,icp0',ix,xm,icp0
c If we are in wrong region, try to correct:
      if(icp0.ne.0 .and. idob_sor(iregion_sor,icp0).ne.iregion)then
c         write(*,*)'Incorrect region',iregion,idob_sor(iregion_sor,icp0)
c     $        ,icp0,ix,xm
         jpm=1
         if(xm.lt.0.)jpm=-1
         ix=ix+jpm
         icp1=int(cij(1+(ix-1)*icinc))
         if(icp1.eq.0 .or. idob_sor(iregion_sor,icp1).ne.iregion)then
            ix=-ix
            uprime=0.
            return
         endif
         write(*,*)'Region Adjusted ix',ix,jpm,xm,iregion
         xm=xm-jpm
         icp0=icp1
      endif
c Distance forward and backward along idf-dimension to adjacent
      dx1=xn(ix+ioff+1)-xn(ix+ioff)
      dx0=xn(ix+ioff)-xn(ix+ioff-1)

c Values of u at the points to be interpolated.
      u0=u(1+(ix-1)*iuinc)
      up=u(1+ix    *iuinc)
      um=u(1+(ix-2)*iuinc)
c Do interpolation.
      uprime=gradinterp(um,u0,up,idf,icp0,xm,dx0,dx1)

      end
