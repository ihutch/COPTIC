c*******************************************************************
c Get the field value in the direction idf appropriately interpolated
c from nearby points, for a specified position.
c
c The interpolation is in the form of box-interpolation in the directions
c perpendicular to the field. Thus the position can be considered to be
c given by integer mesh indices plus real fractions of the next mesh node:
c x(ix)<=x<=x(ix+1), xf=x-x(ix). We allow passing of the whole position
c x because we do remaindering internally. 
c
c The arrays cij, u, and xn may be passed with local origins at the base
c of the box in question. Then we don't have to know ix, only xf for each
c dimension. This approach defeats "-ffortran-bounds-check"ing, though.
c Alternatively, if the whole position x is passed, then the base:
c cij(2*ndims+1), and u(1) arrays, not the local origin, should be passed.
c
c In the gradient direction, we associate the point with the nearest
c neighbor and do interpolations from that, but the routine called knows
c how to correct for possible region crossing. In the orthogonal
c directions, we use the two nearest, i.e. the box.
c
c The field region of the point is known and passed.

      subroutine getfield(ndims,cij,u,iuinc,xn,idf
     $     ,xff,iregion,field)

c Pointers and potential array with origin at the box corner.
      real cij(*),u(*)
c Increments (i.e. structure vector) of u, in each dimension.
      integer iuinc(ndims)
c Position array in the direction idf with origin at box corner.
      real xn(*)
c Direction (dimension) of field component:
      integer idf
c Fractional position to interpolate to for each dimension from the
c passed origin within cij,u. 
      real xff(ndims)
c Region code of particle
      integer iregion

c Local vector storage
      parameter (mdims=10)
      integer idn(mdims)
      real xf(mdims)

      include 'objcom.f'

c We DONT include sormesh, because xn is passed
      parameter (ipwr2nd=2**(ndims_sor-1))
      integer iflags(ipwr2nd)
      real f(ipwr2nd),d(ndims_sor-1)
      integer ii1,ii2

c Allow the passing of the real position, not just fraction.
c This is the case if ix>=1. For fractions, ix=0. 
c Calculate offset and remainder the fractions.
      iux=0
      do ii=1,ndims
         ix=int(xff(ii))
         xf(ii)=xff(ii)-ix
         if(ix.gt.1)iux=iux+iuinc(ii)*(ix-1)
      enddo

c      write(*,*)'iux=',iux

c xn index for passing full position
      ixn0=xff(idf)
c but correct it if we are passing just fractions.
      if(ixn0.lt.1)ixn0=1
c Leading dimension of cij:
      ic1=2*ndims+1
c Correct the index in field direction if fraction .gt.0.5:
      if(xf(idf).ge.0.5)then
         iu0=iuinc(idf)+iux
         ixn0=ixn0+1
         xfidf=xf(idf)-1.
      else
         iu0=0+iux
         xfidf=xf(idf)
      endif

c General-Dimensional version without extrapolation.
      igood=0
      iimax=1
      do ii=1,(ndims-1)
         idii=mod(idf+ii-1,ndims)+1
         d(ii)=xf(idii)
         idn(ii)=idii
c Attempts at speeding up don't do much.
c         iimax=iimax*2
      enddo
      do ii=1,2**(ndims-1)
         ii1=ii-1
         iinc=iu0
c This is likely to be most costly. We are just calculating iinc.
         do ik=1,ndims-1
c This break saves unnecessary iterations, about 15%.
            if(ii1.eq.0)goto 41
            ii2=ii1/2
            if(ii1-2*ii2.ne.0)iinc=iinc+iuinc(idn(ik))
            ii1=ii2
c This is slightly quicker but less than 25%. And not general-dimensional.
c            iinc=iinc+ipa(ii,ik)*iuinc(idn(ik))
         enddo
 41      continue
c Suppose we know there's only 3 dimensions total we could replace with
c         iinc=iinc+ipa(ii,1)*iuinc(idn(1))+ipa(ii,2)*iuinc(idn(2))
c Which saves about 25% of the original time this routine. 
c Pass arrays with local origin.
         call gradlocalregion(
     $        cij(1+ic1*iinc),u(1+iinc)
     $        ,idf,ic1*iuinc(idf),iuinc(idf),xn(ixn0)
     $        ,xfidf,f(ii),iregion,ix,xm)
         
         if(ix.ge.99)then
c           write(*,*)'Getfield no-value',ii,idf,xff
            iflags(ii)=0
         else
            iflags(ii)=1
            igood=igood+1
         endif            
      enddo
      
      if(igood.gt.0)then
c         if(iflags(1).eq.0)write(*,*)'Zero iflags(1) error'
c Field is minus the potential gradient.
         field=-boxinterp(ndims-1,f,iflags,d)
      else
         write(*,'(''Getfield No good vertices. Region'',i3'//
     $        ','' Direction'',i2,'' Fractions'',3f8.4)')
     $        iregion,idf,xff

c Here we should look around for a point that really is in the region,
c since the whole box is not, and use that as the base node with 
c fractions greater than 1. However, this pathological case only 
c occurs in concave regions with specific relation to mesh. Probably
c it can be avoided by mesh adjustment.
c The best direction to look is probably the direction of minimum 
c lattice face distance from point. I need a suitable test case before
c I can really develop this code. 

      endif

      
      end

c********************************************************************
      subroutine getsimple3field(ndims,u,iuinc,xn,idf,xf,field)
c Do interpolation on the field (gradient of u) in a totally
c simple minded way. u is passed with appropriate node selected.
c but if fraction .ge.0.5 in field direction, this is corrected.
      real u(*)
      integer iuinc(3)
      real xf(3)
      real xn(*)
      real f(2,2)
      integer iflags(2,2)
      real d(2)
c zero circumlocution:
      data izer0/0/

c Correct the index in field direction if fraction .gt.0.5:
      if(xf(idf).ge.0.5)then
         iu0=1+iuinc(idf)
         ixn0=1
         xfidf=xf(idf)-1.
      else
         ixn0=0
         iu0=1
         xfidf=xf(idf)
      endif

c      write(*,*)'get3 in',ndims,iuinc
c Values of u at the points to be interpolated.
c      u0=u(1+(ix-1)*iuinc)
c      up=u(1+ix    *iuinc)
c      um=u(1+(ix-2)*iuinc)
      dxp=xn(2+ixn0)-xn(1+ixn0)
c Circumlocution to prevent spurious bounds warning 
      dxm=xn(1+ixn0)-xn(izer0+ixn0)
      do id1=1,3
         do id2=id1+1,3
            if(id1.ne.idf .and. id2.ne.idf)then
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
c               write(*,*)'get3',id1,id2,idf,dxp,dxm,ug
c,xf
            endif
         enddo
      enddo

      field=-box2interp(f,iflags,d)

      end
c*******************************************************************
      real function getpotential(u,cij,iuinc,xff,iregion,itype)
c Pointers and potential array. Cannot here be offset [important!] 
      real cij(*),u(*)
c Increments (i.e. structure vector) of u, in each dimension.
      integer iuinc(*)
c Fractional mesh position to interpolate to for each dimension 
c accounting for offset origin.
      real xff(*)
c Region code of particle
      integer iregion
c Type of interpolation
      integer itype

c objcom, which is needed, contain ndims_... so we don't pass it.
      include 'objcom.f'
      parameter (ndims=ndims_sor)
c Local vector storage
      integer idn(ndims)
      real xf(ndims),uval(2**ndims)
      integer ival(2**ndims)
      
c Do not allow passing just fraction. For fractions, ix=0. 
c Calculate offset and remainder the fractions.
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
c xff must be the absolute mesh position, xf is just the fraction.
c Base pointer within u.
      icb=iux+1
c Leading dimension of cij
      ic1=(2*ndims+1)
      xsqmin=1.e30
      imin=-1
      imissing=0
      do i=1,2**ndims
         i1=i-1
         inc=0
c We are calculating inc. For each box vertex use a logical offset that
c amounts to the binary representation of its sequence number, i1.
c Also calculate the square distance from point to this vertex.
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
c Now inc is the offset within u of the ith box vertex.
c Make iup the pointer to the u element.
         iup=icb+inc
c Make icp the pointer within cij to the object pointer element. 
c Accidental expression. Think of it as (iup-1)*ic1 + 2*ndims+1.
         icp=iup*ic1
c Get that object pointer.
         ico=cij(icp)
         if(ico.ne.0)then
c This is an interface point
            if(idob_sor(iregion_sor,ico).ne.iregion)then
c This vertex is outside the region. Flag it missing
               imissing=imissing+1
               ival(i)=0
c But store the value anyway
               uval(i)=u(iup)
               goto 101
            endif
         endif
c Store this valid vertex
         uval(i)=u(iup)
         ival(i)=1
c Update the minimum valid point
         if(xsq.lt.xsqmin)then
            xsqmin=xsq
            imin=i
         endif
 101     continue
      enddo
c We exit either with all vertices and uvals valid, when imissing=0,
c or there are missing vertices.
c imin refers to the closest valid vertex and that uval is valid. 
c (So we could simple mindedly use its value). 
c If imin=-1 then there was no valid vertex.
      if(imissing.eq.0)then
c Multilinear interpolate.
         getpotential=smultilinearinterp(ndims,uval,xf)      
      elseif(imin.ne.-1)then
c Fall back.
         if(itype.eq.1)then
c No correction for region just multilin anyway.
            getpotential=smultilinearinterp(ndims,uval,xf)
         elseif(itype.eq.2)then
c            call fillin(ndims,uval,ival)
c            write(*,'(8i8)')ival
            call fillinlin(uval,ival,imin
     $           ,u,cij,iuinc,xff,iregion)
            getpotential=smultilinearinterp(ndims,uval,xf)
         else
c Nearest value
            getpotential=uval(imin)
         endif
      else
c No information. We ought to look around further perhaps.
         getpotential=9999
         write(*,*)'getpotential Error. no valid vertex'
      endif
      end
c*******************************************************************
      real function smultilinearinterp(ndims,uval,xf)
c Multilinear interpolation in ndims dimensions.
c The values uval are in binary coded sequence. 
c The position fractions in xf are in order of dimension. 
      real uval(2**ndims),xf(ndims)

c Slight operational inefficiency justified by clarity and no need for
c internal storage arrays. ndims*2**ndims operations. 
      total=0.
      do i=1,2**ndims
         i1=i-1
c For each box vertex use a logical offset that
c amounts to the binary representation of its sequence number, i1.
c Weight according to whether this bit is one or zero.
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
c*******************************************************************
      subroutine fillin(ndims,uval,ival)
c Fill in values with the average of the others.
      real uval(2**ndims),ival(2**ndims)

      uave=0
      utot=0
      itot=0
      do ik=1,2**ndims
         itot=itot+ival(ik)
         if(ival(ik).eq.1)utot=utot+uval(ik)
      enddo
      if(itot.ne.0) then
         uave=utot/itot
      else
      endif
      do ik=1,2**ndims
         if(ival(ik).eq.0)uval(ik)=uave
      enddo
      end
c*****************************************************************
c Alternative. Fill in with values obtained by fitting a linear 
c expansion to the existing points.
      subroutine fillinlin(uval,ival,imin
     $     ,u,cij,iuinc,xff,iregion)
c Fill in values with the average of the others.
c Need object data 
      include 'objcom.f'
c Need mesh data for xn.
      include 'meshcom.f'
      parameter (ndims=ndims_sor)
c contains ndims_sor
      real uval(2**ndims)
      integer ival(2**ndims)
c Passed derivatives for extrapolation.
      real u(*),cij(*),iuinc(ndims+1),xff(ndims)
      integer iregion

      real grad(ndims),centroid(ndims)
c The order of the values is a binary bit representation:
c id=0 alternates every 1 = 2**0
c id=1 alternates every 2 = 2**1
c id=2 alternates every 4, etc.
c So if we are starting in the middle, dimension id+i alternates every
c 2**(mod(id+i,ndims)). 
c Within modulo 2**ndims, an increment by
c 2**(id+1) accesses the alternate id point with invariant subspace. 

c      write(*,'(a,i4,3f9.4)')'Fillinlin',iregion,xff
c      write(*,'(8i9)')ival
c      write(*,'(8f9.4)')uval
      ngoodtot=0
      aveval=0.
      do id=0,ndims-1
         idp1=id+1
         centroid(idp1)=0
         ncgood=0
         i2id=2**(id)
         ngood=0
         difftot=0.
         do i=0,2**(ndims-1)-1
c Over orthogonal subspace...
            itop=i/i2id
            ibot=i - itop*i2id
c Index within total array of subspace by i
            ik1=itop*(2*i2id) + ibot
c Adjacent points along id-direction.
            ik2=mod(ik1+i2id,2**ndims)
            iv1=ival(ik1+1)
            iv2=ival(ik2+1)
            uv1=uval(ik1+1)
            uv2=uval(ik2+1)
            if(iv1.eq.1)then
c zero contribution to centroid because zero position.
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
         if(ngood.ne.0.)then
            grad(idp1)=difftot/ngood
         else
c Null direction. We need a better calculation.
c            grad(idp1)=0.
c get field in idp1 direction at the point. 
            call getfield(ndims,cij(2*ndims+1),u,iuinc
     $           ,xn(ixnp(idp1)+1),idp1
     $           ,xff,iregion,field)
c            write(*,*)'Fillin field value',field,' direction',idp1
c     $           ,' region',iregion
            dx=xn(ixnp(idp1)+2)-xn(ixnp(idp1)+1)
            grad(idp1)=-field*dx
         endif
c ncgood can be zero only if there are no good points at all. 
         if(ncgood.eq.0)stop 'Fillin error. Zero good points'
         centroid(idp1)=centroid(idp1)/ncgood
      enddo
      aveval=aveval/ngoodtot
c Completed the linearized fit to give centroid, grad, aveval.
c      write(*,'(7f9.4)')centroid,grad,aveval
c Fill in the values
      do i=1,2**ndims
         if(ival(i).ne.1)then
c Fill in
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

c      write(*,'(8f9.4)')uval

      end
