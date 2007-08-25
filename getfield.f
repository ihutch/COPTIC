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
      parameter (pwr2nd=2**(ndims_sor-1))
      integer iflags(pwr2nd)
      real f(pwr2nd),d(ndims_sor-1)

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

c Leading dimension of cij:
      ic1=2*ndims+1
c Correct the index in field direction if fraction .gt.0.5:
      if(xf(idf).ge.0.5)then
         iu0=iuinc(idf)+iux
         xfidf=xf(idf)-1.
      else
         iu0=0+iux
         xfidf=xf(idf)
      endif

c General-Dimensional version without extrapolation.
      icount=0
      igood=0
      do ii=1,(ndims-1)
         idii=mod(idf+ii-1,ndims)+1
         d(ii)=xf(idii)
         idn(ii)=idii
      enddo
      do ii=0,2**(ndims-1)-1
         ii1=ii
         iinc=iu0
         do ik=1,ndims-1
            ii2=ii1/2
            ip=ii1-2*ii2
            ii1=ii2
            iinc=iinc+ip*iuinc(idn(ik))
         enddo
         icount=icount+1
c Pass arrays with local origin.
         call gradlocalregion(
     $        cij(1+ic1*iinc),u(1+iinc)
     $        ,idf,ic1*iuinc(idf),iuinc(idf),xn
     $        ,xfidf,f(icount),iregion,ix,xm)
         
         if(ix.ge.99)then
c           write(*,*)'Getfield no-value',icount,ip1,ip2,id1,id2
            iflags(icount)=0
         else
            iflags(icount)=1
            igood=igood+1
         endif            
      enddo
      
      if(igood.gt.0)then

c         if(iflags(1).eq.0)write(*,*)'Zero iflags(1) error'
         field=boxinterp(ndims-1,f,iflags,d)
      else
         write(*,*)'getfield error: No good vertices'
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
         xfidf=xf(idf)-1.
      else
         iu0=1
         xfidf=xf(idf)
      endif

c      write(*,*)'get3 in',ndims,iuinc
c Values of u at the points to be interpolated.
c      u0=u(1+(ix-1)*iuinc)
c      up=u(1+ix    *iuinc)
c      um=u(1+(ix-2)*iuinc)
      dxp=xn(2)-xn(1)
c Circumlocution to prevent spurious bounds warning 
      dxm=xn(1)-xn(izer0)
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

      field=box2interp(f,iflags,d)

      end
