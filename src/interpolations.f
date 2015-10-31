c Interpolations.
      function boxinterp(ndm1,f,flags,d)
      real f(ndm1,ndm1),d(ndm1)
      integer flags(ndm1)
      boxinterp=0.
      if(ndm1.eq.1)then
         boxinterp=box1interp(f,flags,d(1))
      elseif(ndm1.eq.2)then
         boxinterp=box2interp(f,flags,d(1))
      else
         write(*,*)'Boxinterp Error. Unknown dimensionality:',ndm1
      endif
      end
c****************************************************************
      function box2interpnew(f,d,iw,weights,ierr)
c 
c Given values of a function on up to four points adjacent on a 
c 2-D cartesian mesh: f00,f01,f10,f11; and given the fractional
c distances: d1,d2 in the two dimensions (between 0 and 1),
c interpolate the value at d1,d2 using box interpolation  
c weighted by weights if iw.ne.0
      real f(2,2)
      real d(2)
      integer iw
      real weights(2,2)
c Shortcut for unweighted case. 
      if(iw.eq.0)then
         box2interpnew=f(1,1)+d(1)*(f(2,1)-f(1,1))+d(2)*(f(1,2)-f(1,1))
     $     +d(1)*d(2)*(f(2,2)-f(2,1)-f(1,2)+f(1,1))
         return
      endif
      tw=0.
      fa=0.
      do i=1,2
c         whx=i+(1-2*i)*d(1)
         whx=2-i+(2*i-3)*d(1)
         do j=1,2
            wh=(2-j+(2*j-3)*d(2))*whx*weights(i,j)
            tw=tw+wh
            fa=fa+wh*f(i,j)
         enddo
      enddo
      if(tw.eq.0.)then
         write(*,*)'boxinterp zero weight everywhere error'
         write(*,*)'weights',weights,'  ds',d
         ierr=1
c Instead of stopping, just let the divide by zero throw a field 
c corruption error which will give us better diagnostics of where
c the error occurred. It usually happens because we are outside mesh.
c         stop
      endif
      box2interpnew=fa/tw
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
      function interp(Q,nq,y,x)
c Given a monotonic function Q(x)
c on a 1-D grid x=1..nq, solve Q(x)=y for x with interpolation.
c If successful, the function returns the integer part of x.
c If return is 0, then the y is outside Q's range or other error
c and then x=0 indicates y<>Q(1); x=nq+1 indicates y><Q(nq); x=1 other.
c We want interp=x=1 if y=Q(1), interp=x=nq if y=Q(nq).
      integer nq
      real Q(nq)
      real y,x
c
      integer iqr,iql,iqx
      real Qx,Qr,Ql

c      write(*,*)'nq=',nq
      interp=0
      Ql=Q(1)
      Qr=Q(nq)
      QS=Qr-Ql
c Circumlocution to catch y=NAN. [Why is this le? lt breaks geomSoR]
      if(.not.((y-Ql)*(y-Qr).le.0.)) then
c Value is outside the range.
         if((y-Ql)*(Qr-Ql).le.0.)then
            x=0
         elseif((y-Qr)*(Qr-Ql).ge.0.)then
            x=nq+1
         else
            x=1
         endif
         return
      endif
      iql=1
      iqr=nq
 200  if(iqr-iql.eq.1)goto 210
      iqx=(iqr+iql)/2
      Qx=Q(iqx)
c      write(*,*)y,Ql,Qx,Qr,iql,iqr
c      if((Qx-y)*(Qr-y).le.0.) then
      if((Qx-y)*QS.le.0.)then
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
     $        ,iql,iqx,iqr,Ql,Qr,Qx,y,Q(1),Q(nq)
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
      if(ind.gt.ifull(ndims)) write(*,*)'indexexpand index too big'
     $     ,index,' ndims=',ndims,' ifull=',ifull
      end
c********************************************************************
c Convert indices into pointer
      function ipindex(ndims,iLs,ix)
c ix is an ndims-dimensional tuple of indices which is converted into
c a zero-based pointer using the structure vector iLs
      integer iLs(ndims),ix(ndims)
      ipindex=0
      do id=1,ndims
         ipindex=ipindex+(ix(id)-1)*iLs(id)
      enddo
      end
c********************************************************************
c Convert indices into pointer using ifull
      function ipfindex(ndims,ifull,ix)
c ix is an ndims-dimensional tuple of indices which is converted into
c a zero-based pointer using the full dimensions ifull
      integer ifull(ndims),ix(ndims)
      ipfindex=0
      do id=ndims,1,-1
         ipfindex=(ix(id)-1)+ipfindex*ifull(id)
      enddo
      end
c===================================================================
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

      include 'ndimsdecl.f'
      include 'objcom.f'

      if(xm.lt.0)then
         x=xm*dx0
      else
         x=xm*dx1
      endif

      if(icp.ne.0)then
c There is object data. Use object boundary interpolations.
c Addresses of this dimension among objects:
         icd1=2*(id-1)*ndata_cij+1
         icd0=icd1+ndata_cij
         fraction=dob_cij(icd0,icp)
         if(fraction.lt.1. .and. fraction.ge.0.)then
            dxp0=fraction*dx0
            boa=dob_cij(icd0+1,icp)
            coa=dob_cij(icd0+2,icp)
            um=-coa
            dxf0=dxp0*(dxp0+2*boa)/(dxp0+boa)
            dxd0=dxp0+boa
         else
            dxp0=dx0
            dxf0=dx0
            dxd0=dx0
         endif
         fraction=dob_cij(icd1,icp)
         if(fraction.lt.1. .and. fraction.ge.0.)then
            dxp1=fraction*dx1
            boa=dob_cij(icd1+1,icp)
            coa=dob_cij(icd1+2,icp)
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
c General interpolation
      gradinterp= (2.*x+dxf0)/(dxf0+dxf1) * (up-u0)/dxd1
     $     +(dxf1-2.*x)/(dxf0+dxf1) * (u0-um)/dxd0
      if(.not.abs(gradinterp).ge.0)then
      write(*,'(a,/,i4,10f6.2,i3)')
     $     'icp,   x,   dxf0, dxf1, dxd0, dxd1, dx0,  dx1, um, u0, up ',
     $     icp,x,dxf0,dxf1,dxd0,dxd1,dx0,dx1,um,u0,up,id
      write(*,*)'fraction,boa,coa',fraction,boa,coa
      endif
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
      include 'ndimsdecl.f'
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
      icp0=int(cij(1))
c      jpm=0
      if(icp0.eq.0)then
c Short-cut 1: an ordinary point don't call the full routine
c Distance forward and backward along idf-dimension to adjacent
         dx1=xn(2)-xn(1)
c Avoid warnings about argument ranges.
         dx0=xn(1)-xn(ix-1)
c Values of u at the points to be interpolated.
c These array accesses seem about 30% of the costs of this routine.
         u0=u(1)
         up=u(1+iuinc)
         um=u(1-iuinc)
         if(abs(dx1-dx0).lt.1.e-6*dx0)then
c Short-cut 2: for uniform mesh:
c (saves ~25% of routine for uniform mesh even including test):
c A uprime=0. test reduces the time in this routine by about 25%.
c So this evaluation is about 25%.
            uprime= ((2.*xm+1.) * (up-u0)
     $           +(1.-2.*xm) * (u0-um))/(dx0+dx1)
c This makes hardly any difference.
c            uprime= 2.*((xm+0.5) * (up-u0)
c     $           +(0.5-xm) * (u0-um))/(dx0+dx1)
c This seems if anything marginally slower:
c            uprime= (2.*xm*(up-u0-u0+um) + up-um)/(dx0+dx1)
c This saves about 2 sec out of 42 (but changes phi by rounding):
c            uprime=((xm+0.5)*(up-u0)+(0.5-xm)*(u0-um))/dx1
         else
            if(xm.lt.0)then
               x=xm*dx0
            else
               x=xm*dx1
            endif
            uprime= ((2.*x+dx0) * (up-u0)/dx1
     $           +(dx1-2.*x) * (u0-um)/dx0)/(dx0+dx1)
         endif
      else
c Do interpolation using extrapolation information pointed to by icp0.
c This section is only 10% of routine cost.
c If we are in wrong region, try to correct by choosing as the 
c center of interpolation the other node adjacent to point and changing
c the xfraction accordingly:
         if(idob_cij(iregion_cij,icp0).ne.iregion)then
c         write(*,*)'Incorrect region',iregion,idob_cij(iregion_cij,icp0)
c     $        ,icp0,ix,xm
            jpm=1
            if(xm.lt.0.)jpm=-1
            ix=ix+jpm
            icp1=int(cij(1+(ix-1)*icinc))
c Old and seemingly incorrect version:
c            if(icp1.eq.0 .or. idob_cij(iregion_cij,icp1).ne.iregion)then
            if(icp1.ne.0.and.idob_cij(iregion_cij,icp1).ne.iregion)then
               ix=99
               uprime=0.
               return
            endif
            xm=xm-jpm
            icp0=icp1
c         write(*,*)'Base Position Adjusted ix',ix,jpm,xm,iregion
         endif
c Distance forward and backward along idf-dimension to adjacent
         dx1=xn(ix+1)-xn(ix)
         dx0=xn(ix)-xn(ix-1)
c Values of u at the points to be interpolated.
         ixiu=1+(ix-1)*iuinc
         u0=u(ixiu)
         up=u(ixiu+iuinc)
         um=u(ixiu-iuinc)
         uprime=gradinterp(um,u0,up,idf,icp0,xm,dx0,dx1)
      endif

      if(.not.abs(uprime).lt.1.e20)then
         write(*,*)'Getlocalregion excessive uprime',icp0,xm,ix,ixiu
     $        ,iregion
         write(*,*) 'dx0,1',dx0,dx1,' uprime',uprime,up,u0,um,u(1),iuinc
      endif

      end

c***************************************************************
cXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c****************************************************************
