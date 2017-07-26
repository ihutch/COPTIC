!******************************************************************
! Get the field value in the direction idf appropriately interpolated
! from nearby points, for a specified position.
!
! The interpolation is in the form of box-interpolation in the directions
! perpendicular to the field. Thus the position can be considered to be
! given by integer mesh indices plus real fractions of the next mesh node:
! x(ix)<=x<=x(ix+1), xf=x-x(ix). We allow passing of the whole position
! x because we do remaindering internally. 
!
! It used to be the case that:
! The arrays cij, u, and xn may be passed with local origins at the base
! of the box in question. Then we don't have to know ix, only xf for each
! dimension. This approach defeats "-ffortran-bounds-check"ing, though.
! Now, however it is assumed that the whole array is passed so that one
! can avoid overrunning the array by using the iuinc values. 
! 
! So now: the whole position x is passed, and the base:
! cij(2*ndims+1), and u(1) arrays, not the local origin, should be passed.
!
! In the gradient direction, we associate the point with the nearest
! neighbor and do interpolations from that, but the routine called knows
! how to correct for possible region crossing. In the orthogonal
! directions, we use the two nearest, i.e. the box.
!
! The field region of the point is known and passed. 
! Value will be rubbish if xff<1. or xff>nmesh because arrays will be
! overrun.

      subroutine getfield(cij,u,iuinc,xn,idf,xff,iregion,field)

      include 'ndimsdecl.f'
! Pointers and potential array with origin at the box corner.
      real cij(*),u(*)
! Increments (i.e. structure vector) of u, in each dimension.
      integer iuinc(ndims+1)
! Position array in the direction idf with origin at box corner.
      real xn(*)
! Direction (dimension) of field component: integer idf
! Fractional position to interpolate to for each dimension from the
! passed origin within cij,u. 
      real xff(ndims)
! Region code of particle
      integer iregion

! Local vector storage
      parameter (mdims=3)
      integer idn(mdims)
      real xf(mdims)

      include 'objcom.f'
! for external field need 3dcom.f
      include '3dcom.f'

! We DONT include sormesh, because xn is passed
      parameter (ipwr2nd=2**(ndims-1))
      integer iflags(ipwr2nd)
      real weights(ipwr2nd)
      real f(ipwr2nd),d(ndims-1)
      integer ii1

! Allow the passing of the real position, not just fraction.
! This is the case if ix>=1. For fractions, ix=0. 
! Calculate offset and remainder the fractions.
      iux=0
      do ii=1,ndims
         ix=int(xff(ii))
         xf(ii)=xff(ii)-ix
         if(ix.gt.1)iux=iux+iuinc(ii)*(ix-1)
      enddo
! Now iux is the index of the box lower origin. 
!      write(*,*)'iux=',iux
! Start assuming there is no non-unity weighting.
      iw=0
! xn index for passing full position
      ixn0=int(xff(idf))
! but correct it if we are passing just fractions.
      if(ixn0.lt.1)ixn0=1
! Leading dimension of cij:
      ic1=2*ndims+1
! Correct the index in field direction if fraction .gt.0.5:
      if(xf(idf).ge.0.5)then
         iu0=iuinc(idf)+iux
         ixn0=ixn0+1
         xfidf=xf(idf)-1.
      else
         iu0=0+iux
         xfidf=xf(idf)
      endif
! Now iu0 is the chosen (closest in field direction) box origin 
! General-Dimensional version without extrapolation.
      igood=0
      do ii=1,(ndims-1)
         idii=mod(idf+ii-1,ndims)+1
         d(ii)=xf(idii)
         weights(ii)=1.
         idn(ii)=idii
! Attempts at speeding up don't do much.
      enddo
      do ii=1,2**(ndims-1)
         ii1=ii-1
         iinc=iu0
! The cost of calculating iinc is non-negligible 
! but this can't be improved upon.
         do ik=1,ndims-1
            if(btest(ii1,ik-1))iinc=iinc+iuinc(idn(ik))
!            if(mod(ii1/ik,2).ne.0)iinc=iinc+iuinc(idn(ik))
         enddo
! Now iinc is the index of u including the offset from the chosen origin 
! to address the multilinear position under consideration. 
! Pass arrays with local origin.
         if(
     $        cij(1+ic1*iinc).eq.0
     $        .and.abs((xn(ixn0+1)-xn(ixn0))-(xn(ixn0)-xn(ixn0-1)))
     $              .lt.1.e-6*(xn(ixn0)-xn(ixn0-1)))then
! This is a shortcut to improve performance by not calling
! gradlocalregion when it's not necessary. It saves 25% of the cost
! of getfield on uniform mesh. >10% of total particle time.
            dx1=xn(ixn0+1)-xn(ixn0)
            dx0=xn(ixn0)-xn(ixn0-1)
            u0=u(1+iinc)
            up=u(1+iinc+iuinc(idf))
            um=u(1+iinc-iuinc(idf))
            f(ii)= ((2.*xfidf+1.) * (up-u0)
     $           +(1.-2.*xfidf) * (u0-um))/(dx0+dx1)
            ix=1
         else
            call gradlocalregion(
     $           cij(1+ic1*iinc),u(1+iinc)
     $           ,idf,ic1*iuinc(idf),iuinc(idf),xn(ixn0)
     $           ,xfidf,f(ii),iregion,ix,xm)
         endif
         if(.not.abs(f(ii)).lt.1.e20)then
            write(*,*)'Corrupt gradlocalregion'
     $           ,idf,ii,f
            write(*,*)'iinc,ic1,iuinc(idf),cij(1+ic1*iinc),u(1+iinc)'
     $           ,iinc,ic1,iuinc(idf),cij(1+ic1*iinc),u(1+iinc)
            iad=1
            do ki=1,ndims
               iad=iad+(int(xff(ki))-1)*iuinc(ki)
            enddo
            write(*,*)'xff=',xff,iad
         endif
         if(ix.ge.99)then
!            write(*,*)'Getfield no-value',ii,idf,iinc,xff
! This gradient request has failed (probably) because the lattice leg 
! has both end-points outside the region. Quite possibly there is a point 
! on the other side of one of them that is in the region, from which we
! should extrapolate. We can use gradlocalregion for that (I think) by
! adjusting the xfidf value and the base.
!            if(.false.)then
            weights(ii)=0.
            f(ii)=0.
            if(xfidf.ge.0.)then
               iincm=iinc-iuinc(idf)
               iincp=iinc+2*iuinc(idf)
               if(iincm.lt.0)then
!                  write(*,*)'1st iinc,iincm,iincp',iinc,iincm,iincp
                  icptm=0
               else
                  icptm=int(cij(1+ic1*(iinc-iuinc(idf))))
               endif
               icptp=int(cij(1+ic1*(iinc+2*iuinc(idf))))
               if((iincm.ge.0).and.(icptm.ne.0).and.
     $              (idob_cij(iregion_cij,icptm).eq.iregion))then
! xf positive, node0 (relative to 1) in region, look at node 0.
                  call gradlocalregion(cij(1+ic1*(iinc-iuinc(idf))),
     $                 u(1+iinc-iuinc(idf)) ,idf,ic1*iuinc(idf)
     $                 ,iuinc(idf),xn(ixn0-1) ,xfidf+1,f(ii),iregion,ix
     $                 ,xm)
!         if(abs(f(ii)).gt.1.e20)write(*,*)'Corrupt gradlocalregion2'
!     $           ,idf,ii,f

                  if(ix.ne.99)then
                     iw=iw+1
                     weights(ii)=weights(ii)+(1.-xfidf)
!                     write(*,*)ii,xfidf,' weights=',weights(ii)
                  endif
! The alternative is to use only the first if it works by uncommenting
! the following line and commenting the two after it.
!               elseif(icptp.ne.0.and.idob_cij(iregion_cij,icptp)
               endif
               if(iincp.le.iuinc(ndims+1)
     $              .and.icptp.ne.0.and.
     $              idob_cij(iregion_cij,icptp).eq.iregion)then
! xf positive, node3 (relative to 1) in region, look at node 3.
                  call gradlocalregion(cij(1+ic1*(iinc+2*iuinc(idf))),
     $                 u(1+iinc+2*iuinc(idf)) ,idf,ic1*iuinc(idf)
     $                 ,iuinc(idf),xn(ixn0+2) ,xfidf-2,fii,iregion,ix
     $                 ,xm)
!         if(abs(fii).gt.1.e20)write(*,*)'Corrupt gradlocalregion3'
!     $           ,idf,ii,f

                  if(ix.ne.99)then
                     if(f(ii).ne.0.)then
                        f(ii)=(1.-xfidf)*f(ii)+xfidf*fii
                     else
                        f(ii)=fii
                     endif
                     iw=iw+1
                     weights(ii)=weights(ii)+ xfidf
!                     write(*,*)'2nd ',xfidf,' weights=',weights(ii)
                  endif
               endif
            else
               iincm=iinc-2*iuinc(idf)
               iincp=iinc+iuinc(idf)
               if(iincm.le.0)then
!                  write(*,*)'idf, iinc,iincm,iincp',idf,iinc,iincm,iincp
                  icptm=0
               else
                  icptm=int(cij(1+ic1*(iinc-2*iuinc(idf))))
               endif
               icptp=int(cij(1+ic1*(iinc+iuinc(idf))))
               if((iincp.le.iuinc(ndims+1))
     $              .and.(icptp.ne.0).and.(idob_cij(iregion_cij,icptp)
     $              .eq.iregion))then
! xf negative, node2 (relative to 1) in region, look at node 2.
                  call gradlocalregion(cij(1+ic1*(iinc+iuinc(idf))),
     $                 u(1+iinc+iuinc(idf)) ,idf,ic1*iuinc(idf)
     $                 ,iuinc(idf),xn(ixn0+1) ,xfidf-1,f(ii),iregion,ix
     $                 ,xm)
!         if(abs(f(ii)).gt.1.e20)write(*,*)'Corrupt gradlocalregion4'
!     $           ,idf,ii,f

                  if(ix.ne.99)then
                     iw=iw+1
                     weights(ii)=(1.+xfidf)
!                     write(*,*)ii,xfidf,' weights=',weights(ii)
                  endif
               endif
               if((iincm.ge.0).and.(icptm.ne.0).and.
     $              (idob_cij(iregion_cij,icptm).eq.iregion))then
! xf negative, node-1 (relative to 1) in region, look at node -1.
                  call gradlocalregion(cij(1+ic1*(iinc-2*iuinc(idf))),
     $                 u(1+iinc-2*iuinc(idf)) ,idf,ic1*iuinc(idf)
     $                 ,iuinc(idf),xn(ixn0-2) ,xfidf+2,fii,iregion,ix
     $                 ,xm)
!         if(abs(fii).gt.1.e20)then
!            write(*,*)'Corrupt gradlocalregion5'
!     $           ,idf,ii,iinc,xfidf+2,ix,xm,fii
!            write(*,*)'icptm,icptp,xff',icptm,icptp,xff
!            write(*,*)iinc-2*iuinc(idf)
!         endif

                  if(ix.ne.99)then
                     if(f(ii).ne.0.)then
!                        write(*,*)'xfidf',xfidf
                        f(ii)=(1.+xfidf)*f(ii)-xfidf*fii
                     else
                        f(ii)=fii
                     endif
                     iw=iw+1
                     weights(ii)=weights(ii)-xfidf
!                     write(*,*)ii,xfidf,' weights=',weights(ii)
                  endif
               endif
            endif
!            endif
            if(ix.eq.99)then
               iflags(ii)=0
               weights(ii)=0.
               iw=iw+1
            else
               iflags(ii)=1
!               iw=iw+1
!               weights(ii)=1.
               igood=igood+1
            endif
         else
            iflags(ii)=1
            weights(ii)=1.
            igood=igood+1
! Case corresponding to extrapolation (not over extrapolation)
!            if((xm-xfidf).ne.0)then
!               weights(ii)=1.
!            endif
         endif

! Debugging code:
         if(ix.eq.98)then
            write(*,*)'ic1,iinc,idf,iregion,xfidf'
     $           ,ic1,iinc,idf,iregion,xfidf
            write(*,*)(dob_cij(k,int(cij(1+ic1*iinc))),k=1,18)
     $           ,cij(1+ic1*iinc)
         endif
!         if(abs(f(ii)).gt.1.e20)write(*,*)'Enddo corrupt',f(ii),ii
      enddo
      
      if(igood.gt.0)then
!         if(iflags(1).eq.0)write(*,*)'Zero iflags(1) error'
! Field is minus the potential gradient.
         ierr=0
         field=-box2interpnew(f,d,iw,weights,ierr)
         if(ierr.ne.0)write(*,*)xff,f
         if(.not.abs(field).lt.1.e20)then
            write(*,*)'box2interpnew field corruption in getfield'
            write(*,*)field,'  f=',f,'  d=',d,'  iw=',iw,'  weights='
     $           ,weights,'   igood=',igood,'   ix=',ix
         endif
      else
         write(*,'(''Getfield No good vertices. Region'',i3'//
     $        ','' Direction'',i2,'' Fracs'',3f8.4)')
     $        iregion,idf,xff
! This flags a problem to the calling routine. padvnc stops.
!         field=1.e13
! This instead gives a field of zero as a complete hack fallback but it
! allows the calculation to continue.
         field=0.
! Here we should look around for a point that really is in the region,
! since the whole box is not, and use that as the base node with 
! fractions greater than 1. However, this pathological case only 
! occurs in concave regions with specific relation to mesh. Probably
! it can be avoided by mesh adjustment.
! The best direction to look is probably the direction of minimum 
! lattice face distance from point. I need a suitable test case before
! I can really develop this code. 

      endif

! External field
      if(lextfield)field=field+extfield(idf)
      
      end

!********************************************************************
      subroutine getsimple3field(ndims,u,iuinc,xn,idf,xf,field)
! Do interpolation on the field (gradient of u) in a
! simple minded way. 
! u is passed with appropriate node selected as its base index,
! but if fraction .ge.0.5 in field direction, this is corrected.
! iuinc is the structure vector of the u array.
! xn is the compact position vector, local to the array position.
! idf is the direction (axis) in which to get the field.
! xf is the position fraction in 3d
! On output the field is in field.
      integer ndims
      real u(*)
      integer iuinc(3)
      real xf(3)
      real xn(*)
      real field
! Local variables
      integer iflags(2,2)
      real f(2,2)
      real d(2)
! zero circumlocution:
      data izer0/0/
! silence warnings
      data d/0.,0./iflags/0,0,0,0/f/0.,0.,0.,0./

! Silence warning about ndims unused
      ixn0=ndims
! Correct the index in field direction if fraction .gt.0.5:
      if(xf(idf).ge.0.5)then
         iu0=1+iuinc(idf)
         ixn0=1
         xfidf=xf(idf)-1.
      else
         ixn0=0
         iu0=1
         xfidf=xf(idf)
      endif

!      write(*,*)'get3 in',ndims,iuinc
! Values of u at the points to be interpolated.
!      u0=u(1+(ix-1)*iuinc)
!      up=u(1+ix    *iuinc)
!      um=u(1+(ix-2)*iuinc)
      dxp=xn(2+ixn0)-xn(1+ixn0)
! Circumlocution to prevent spurious bounds warning 
      dxm=xn(1+ixn0)-xn(izer0+ixn0)
      do id1=1,3
         do id2=id1+1,3
            if(id1.ne.idf .and. id2.ne.idf)then
               if(xf(id1).gt.1)write(*,*
     $              )'Getsimple3field fraction error',xf(id1)
               d(1)=xf(id1)
               d(2)=xf(id2)
               do ip1=0,1
               do ip2=0,1
                  u0=u(iu0+ip1*iuinc(id1)+ip2*iuinc(id2))
                  up=u(iu0+ip1*iuinc(id1)+ip2*iuinc(id2)+iuinc(idf))
                  um=u(iu0+ip1*iuinc(id1)+ip2*iuinc(id2)-iuinc(idf))
                  ugradp=(up-u0)/dxp
                  ugradm=(u0-um)/dxm
                  ug=(0.5-xfidf)*ugradm + (0.5+xfidf)*ugradp
                  f(ip1+1,ip2+1)=ug
                  iflags(ip1+1,ip2+1)=1
               enddo
               enddo
!               write(*,*)'get3',id1,id2,idf,dxp,dxm,ug
!,xf
            endif
         enddo
      enddo

      field=-box2interp(f,iflags,d)

!      write(*,'(a,10f10.4)')'get3field=',field,f
      end
!*******************************************************************
      real function getpotential(u,cij,iuinc,xff,iregion,itype)
! Pointers and potential array. Cannot here be offset [important!] 
      real cij(*),u(*)
! Increments (i.e. structure vector) of u, in each dimension.
      integer iuinc(*)
! Fractional mesh position to interpolate to for each dimension 
! accounting for offset origin.
      real xff(*)
! Region code of point
      integer iregion
! Type of interpolation
      integer itype
! If itype=0, no use of cij,iregion: skip to multilin interpolation.
! If itype=1, use "missing" vertices anyway.
! If itype=2, call fillinlin to fill in missing vertices
! If itype=3, just return the value at nearest vertex.
      
      include 'ndimsdecl.f'
      include 'objcom.f'
! Local vector storage
!      integer idn(ndims)
      real xf(ndims),uval(2**ndims)
      integer ival(2**ndims)
      integer ierr
      data ierr/0/
      
! Do not allow passing just fraction. For fractions, ix=0. 
! Calculate offset and remainder the fractions.
      iux=0
      do ii=1,ndims
         ix=int(xff(ii))
         xf(ii)=xff(ii)-ix
         if(ix.eq.0)then
            write(*,*)ndims,ii,xf
            stop 'getpotential ERROR. Passed fraction < 1' 
         else
            iux=iux+iuinc(ii)*(ix-1)
         endif
      enddo
! xff must be the total mesh position, xf is just the fraction.
! Base pointer within u.
      icb=iux+1
! Leading dimension of cij
      ic1=(2*ndims+1)
      xsqmin=1.e30
      imin=-1
      imissing=0
      do i=1,2**ndims
         i1=i-1
         inc=0
! We are calculating inc. For each box vertex use a logical offset that
! amounts to the binary representation of its sequence number, i1.
! Also calculate the square distance from point to this vertex.
         xsq=0
         do ik=1,ndims
            i2=i1/2
            if(i1-2*i2.ne.0)then 
               inc=inc+iuinc(ik)
               xsq=xsq+(1.-xf(ik))**2
            else
               xsq=xsq+xf(ik)**2 
            endif
            i1=i2
         enddo
! Now inc is the offset within u of the ith box vertex.
! Make iup the pointer to the u element.
         iup=icb+inc
         if(itype.gt.0)then
! Make icp the pointer within cij to the object pointer element. 
! Accidental expression. Think of it as (iup-1)*ic1 + 2*ndims+1.
            icp=iup*ic1
! Get that object pointer.
            ico=int(cij(icp))
            if(ico.ne.0)then
! This is an interface point
               if(idob_cij(iregion_cij,ico).ne.iregion)then
!               if(.false.)then
! This vertex is outside the region. Flag it missing
!               write(*,'(i4,$)')idob_cij(iregion_cij,ico)
                  imissing=imissing+1
                  ival(i)=0
! But store the value anyway
                  uval(i)=u(iup)
                  goto 101
               endif
            endif
         endif
! Store this valid vertex
         uval(i)=u(iup)
         ival(i)=1
! Update the minimum valid point
         if(xsq.lt.xsqmin)then
            xsqmin=xsq
            imin=i
         endif
 101     continue
      enddo
! We exit either with all vertices and uvals valid, when imissing=0,
! or there are missing vertices.
! imin refers to the closest valid vertex and that uval is valid. 
! (So we could simple mindedly use its value). 
! If imin=-1 then there was no valid vertex.
      if(imissing.eq.0)then
! Multilinear interpolate.
         getpotential=smultilinearinterp(ndims,uval,xf)      
      elseif(imin.ne.-1)then
! Fall back.
         if(itype.eq.1)then
! No correction for region just multilin anyway.
            getpotential=smultilinearinterp(ndims,uval,xf)
         elseif(itype.eq.2)then
            call fillinlin(uval,ival
     $           ,u,cij,iuinc,xff,iregion)
            getpotential=smultilinearinterp(ndims,uval,xf)
         else
! Nearest value
            getpotential=uval(imin)
         endif
      else
! No information. We ought to look around further perhaps.
         getpotential=9999
         write(*,*)'getpotential Error. no valid vertex',iregion,imin
      endif
      if(abs(getpotential).gt.20)then
         write(*,*)'Getpotential',getpotential,(xff(kk),kk=1,ndims)
     $        ,itype,imin,imissing
         ierr=ierr+1
      endif
      if(ierr.gt.50)stop
      end
!*******************************************************************
      real function smultilinearinterp(ndims,uval,xf)
! Multilinear interpolation in ndims dimensions.
! The values uval are in binary coded sequence. 
! The position fractions in xf are in order of dimension. 
      real uval(2**ndims),xf(ndims)

! Slight operational inefficiency justified by clarity and no need for
! internal storage arrays. ndims*2**ndims operations. 
      total=0.
      do i=1,2**ndims
         i1=i-1
! For each box vertex use a logical offset that
! amounts to the binary representation of its sequence number, i1.
! Weight according to whether this bit is one or zero.
         thisval=uval(i)
         do ik=1,ndims
            i2=i1/2
            if(i1-2*i2.ne.0)then 
               thisval=xf(ik)*thisval
            else
               thisval=(1.-xf(ik))*thisval
            endif
            i1=i2
         enddo
         total=total+thisval
      enddo

      smultilinearinterp=total
      end
!*****************************************************************
! Alternative. Fill in with values obtained by fitting a linear 
! expansion to the existing points. ival tells whether a vertex
! has already been filled, uval returns the values, filled in
! if necessary.
      subroutine fillinlin(uval,ival
     $     ,u,cij,iuinc,xff,iregion)
! Fill in values with the average of the others.
! Need object data 
      include 'ndimsdecl.f'
      include 'objcom.f'
! Need mesh data for xn.
      include 'meshcom.f'
! contains ndims
      real uval(2**ndims)
      integer ival(2**ndims)
! Passed derivatives for extrapolation.
      real u(*),cij(*),xff(ndims)
      integer iregion,iuinc(ndims+1)

      real grad(ndims),centroid(ndims)
! The order of the values is a binary bit representation:
! id=0 alternates every 1 = 2**0
! id=1 alternates every 2 = 2**1
! id=2 alternates every 4, etc.
! So if we are starting in the middle, dimension id+i alternates every
! 2**(mod(id+i,ndims)). 
! Within modulo 2**ndims, an increment by
! 2**(id+1) accesses the alternate id point with invariant subspace. 

!      write(*,'(a,i4,3f9.4)')'Fillinlin',iregion,xff
!      write(*,'(8i9)')ival
!      write(*,'(8f9.4)')uval
      ngoodtot=0
      aveval=0.
      do id=0,ndims-1
         idp1=id+1
         centroid(idp1)=0
         ncgood=0
         i2id=2**id
         ngood=0
         difftot=0.
         do i=0,2**(ndims-1)-1
! Over orthogonal subspace...
            itop=i/i2id
            ibot=i - itop*i2id
! Index within total array of subspace by i
            ik1=itop*(2*i2id) + ibot
! Adjacent points along id-direction.
            ik2=mod(ik1+i2id,2**ndims)
            iv1=ival(ik1+1)
            iv2=ival(ik2+1)
            uv1=uval(ik1+1)
            uv2=uval(ik2+1)
            if(iv1.eq.1)then
! zero contribution to centroid because zero position.
               ncgood=ncgood+1
               aveval=aveval+uv1
               ngoodtot=ngoodtot+1
            endif
            if(iv2.eq.1)then
               centroid(idp1)=centroid(idp1)+1
               ncgood=ncgood+1
               aveval=aveval+uv2
               ngoodtot=ngoodtot+1
               if(iv1.eq.1)then
                  difftot=difftot+uv2-uv1
                  ngood=ngood+1
               endif
            endif
         enddo
         if(.false.)then
!         if(ngood.ne.0.)then
! Actually using the field at point for gradient seems better in most
! situations where a point is missing. In which case never use this:
            grad(idp1)=difftot/ngood
         else
! Null direction. We need a better calculation.
! get field in idp1 direction at the point. 
            call getfield(cij(2*ndims+1),u,iuinc
     $           ,xn(ixnp(idp1)+1),idp1
     $           ,xff,iregion,field)
!            write(*,*)'Fillin field value',field,' direction',idp1
!     $           ,' region',iregion
            dx=xn(ixnp(idp1)+2)-xn(ixnp(idp1)+1)
            grad(idp1)=-field*dx
         endif
! ncgood can be zero only if there are no good points at all. 
         if(ncgood.eq.0)stop 'Fillin error. Zero good points'
         centroid(idp1)=centroid(idp1)/ncgood
      enddo
      aveval=aveval/ngoodtot
! Completed the linearized fit to give centroid, grad, aveval.
!      write(*,'(7f9.4)')centroid,grad,aveval
! Fill in the values
      do i=1,2**ndims
         if(ival(i).ne.1)then
! Fill in
            thisval=aveval
            i1=i-1
            do ik=1,ndims
               i2=i1/2
               if(i1-2*i2.ne.0)then 
                  thisval=thisval+grad(ik)*(1.-centroid(ik))
               else
                  thisval=thisval+grad(ik)*(-centroid(ik))
               endif
               i1=i2
            enddo
            uval(i)=thisval
         endif
      enddo

!      write(*,'(8f9.4)')uval

      end
!*****************************************************************
! Get the field at a specified position x, without having to pass lots
! of arguments. Return it as a vector in field.
      subroutine fieldatpoint(x,u,cij,iLs,field)
! u is potential, cij is stencil array, iLs is the structure of u,cij.
! All the other parameters must be obtained from commons. 
! Include the mesh xn, and ndims.
      include 'ndimsdecl.f' 
      include 'meshcom.f'
      real x(ndims),u(*),cij(*),field(ndims)
      integer iLs(ndims+1)
      external imaskregion
      real xff(ndims)

! Determine the region
!      iregion=insideall(ndims,x)
      iregion=insidemask(ndims,x)
! Locate the mesh full fractional postion of x.
      do id=1,ndims
! Offset to start of dimension-id-position array.
         ioff=ixnp(id)
! xn is the position array for each dimension arranged linearly.
! Find the index of xprime in the array xn:
         ix=interp(xn(ioff+1),ixnp(id+1)-ioff,x(id),xm)
         xff(id)=xm
      enddo

! Get each component of the field.
      do id=1,ndims
         call getfield(cij(2*ndims+1),u,iLs
     $        ,xn(ixnp(id)+1)
     $        ,id,xff,imaskregion(iregion),field(id))
      enddo

      end
!******************************************************************
! Get the potential at a specified postion x.
      real function potentialatpoint(x,u,cij,iLs)
! u is potential, cij is stencil array, iLs is the structure of u,cij.
! All the other parameters must be obtained from commons. 
! Include the mesh xn, and ndims. 
      include 'ndimsdecl.f'
      include 'meshcom.f'
      real x(ndims),u(*),cij(*)
      integer iLs(ndims+1)

      real xff(ndims)

! Determine the region
!      iregion=insideall(ndims,x)
      iregion=insidemask(ndims,x)
! Locate the mesh full fractional postion of x.
      do id=1,ndims
! Offset to start of dimension-id-position array.
         ioff=ixnp(id)
! xn is the position array for each dimension arranged linearly.
! Find the index of xprime in the array xn:
         ix=interp(xn(ioff+1),ixnp(id+1)-ioff,x(id),xm)
         xff(id)=xm
      enddo
! This call using fillinlin gives problems for quasineutral cases.
! Presumably because there are discontinuities in phi at object edge.
!      potentialatpoint=getpotential(u,cij,iLs,xff,iregion,2)
      potentialatpoint=getpotential(u,cij,iLs,xff,iregion,1)
!      if(potentialatpoint.eq.999)then
!         xs=sqrt(x(1)**2+x(2)**2+x(3)**2)
!         write(*,'(7f10.5)')x,xs,xff
!      endif

      end
!*****************************************************************
! Get the field at a specified position x, without having to pass lots
! of arguments. Return it as a vector in field.
! This simple version assumes no object boundary in the vicinity.
      subroutine fieldsimple3atpoint(x,u,iLs,field)
! u is potential, iLs is the structure of u.
! All the other parameters must be obtained from commons. 
! Include the mesh xn, and ndims. 
      include 'ndimsdecl.f'
      include 'meshcom.f'
      real x(ndims),u(*),field(ndims)
      integer iLs(ndims+1)
      real xff(ndims)
      integer kxf(ndims)

      koff=1
! Locate the mesh full fractional postion of x.
      do id=1,ndims
! Offset to start of dimension-id-position array.
         ioff=ixnp(id)
! xn is the position array for each dimension arranged linearly.
! Find the index of x in the array xn:
         ix=interp(xn(ioff+1),ixnp(id+1)-ioff,x(id),xm)
         if(ix.le.0)write(*,*)'fs3 interp error',x(id),xm
         kxf(id)=int(xm)
         koff=koff+(kxf(id)-1)*iLs(id)
         xff(id)=xm-kxf(id)
      enddo
      

! Get each component of the field.
      do id=1,ndims
         call getsimple3field(ndims,u(koff)
     $        ,iLs,xn(ixnp(id)+kxf(id)),id,xff,field(id))
      enddo
!      write(*,*)kxf,koff,xff,field

      end
!*************************************************************************
