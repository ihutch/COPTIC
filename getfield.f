c Get the field value in the direction idf appropriately interpolated
c from nearby points, for a specified position.
c
c The interpolation is in the form of box-interpolation in the directions
c perpendicular to the field. Thus the position can be considered to be
c given by integer mesh indices plus real fractions of the next mesh node:
c x(ix)<=x<=x(ix+1), xf=x-x(ix).
c The arrays cij, u, and xn are passed with local origins at the base
c of the box in question. Then we don't have to know ix, only xf for each
c dimension. This approach defeats "-ffortran-bounds-check"ing, though.
c
c In the gradient direction, we associate the point with the nearest
c neighbor and do interpolations from that, but the routine called knows
c how to correct for possible region crossing. In the orthogonal
c directions, we use the two nearest, i.e. the box.
c
c The field region of the point is known and passed.

      subroutine getfield(ndims,cij,u,icinc,iuinc,xn,idf
     $     ,xf,iregion,field)

c Pointers and potential array with origin at the box corner.
      real cij(*),u(*)
c Increments of cij,u, in each dimension.
      integer icinc(ndims),iuinc(ndims)
c Position array in the direction idf with origin at box corner.
      real xn(*)
c Direction (dimension) of field component:
      integer idf
c Fractional position to interpolate to for each dimension from the
c passed origin within cij,u. 
      real xf(ndims)
c Region code of particle
      integer iregion

      include 'objcom.f'

c We DONT include sormesh, because xn is passed
      parameter (pwr2nd=2**(ndims_sor-1))
      integer iflags(pwr2nd)
      real f(pwr2nd),d(ndims_sor-1)
      logical lextrapolate
      integer iorder(4),jp1(4),jp2(4),jm1(4),jm2(4)
      data iorder/1,2,4,3/
      data jp1,jp2/1,1,0,0,0,1,0,1/
      data jm1,jm2/0,0,1,1,1,0,1,0/
      data lextrapolate/.true./
c      data lextrapolate/.false./

c This is a 3-D only version for now. But needs to be made into a
c general-dimension version.      
      do id1=1,3
         do id2=id1+1,3
            if(id1.ne.idf .and. id2.ne.idf)then
               d(1)=xf(id1)
               d(2)=xf(id2)
               icount=0
               igood=0
               do ip2=0,1
               do ip1=0,1
                  icount=icount+1
c Pass arrays with local origin.
                  call gradlocalregion(
     $              cij(1+ip1*icinc(id1)+ip2*icinc(id2))
     $              ,u(1+ip1*iuinc(id1)+ip2*iuinc(id2))
     $              ,idf,icinc(idf),iuinc(idf),xn
     $              ,xf(idf),f(icount),iregion,ix,xm)

                  if(ix.ge.99)then
c           write(*,*)'Getfield no-value',icount,ip1,ip2,id1,id2
                     iflags(icount)=0
                  else
                     iflags(icount)=1
                     igood=igood+1
                  endif            
               enddo
               enddo
c Extrapolation section
               if(lextrapolate.and.igood.ne.icount)then
c                  write(*,*)'***Extrapolation: iflags=',iflags
                  icount=0
                  do ip2=0,1
                  do ip1=0,1
                     icount=icount+1
                     if(iflags(icount).eq.0)then
c Attempt extrapolation. Examine neighbors.
                        icp=iorder(mod(iorder(icount),4)+1)
                        icm=iorder(mod(iorder(icount)+2,4)+1)
                        ifp=iflags((icp))
                        ifm=iflags((icm))
c                        write(*,*)'icp,icm,ifp,ifm',icp,icm,ifp,ifm
                        if(ifp.eq.1.and.ifm.ne.1)then
c Only ifp.
                           ip1p=2*jp1((icount))-ip1
                           ip2p=2*jp2((icount))-ip2
                           icl=(icp)
                        elseif(ifp.ne.1.and.ifm.eq.1)then
c Only ifm.
                           ip1p=2*jm1((icount))-ip1
                           ip2p=2*jm2((icount))-ip2
                           icl=(icm)
                        else
c Skip
                           goto 100
                        endif
                        call gradlocalregion(
     $                       cij(1+ip1p*icinc(id1)+ip2p*icinc(id2))
     $                       ,u(1+ip1p*iuinc(id1)+ip2p*iuinc(id2))
     $                       ,idf,icinc(idf),iuinc(idf),xn
     $                       ,xf(idf),fextrap,iregion,ix,xm)
                        if(ix.ne.99)then
c Extrapolated value present. Use it.
                           write(*,*)'Extrapolated:',icount,icl,
     $                          ' to',ip1,ip2,
     $                          jp1((icount)),jp2((icount)),
     $                          jm1((icount)),jm2((icount)),
     $                          ' from',ip1p,ip2p
                           f(icount)=2*f((icl))- fextrap
                           iflags(icount)=2
                        else
                           write(*,*)'Extrapolation failed',icount,icl,
     $                          ' to',ip1,ip2,
     $                          jp1((icount)),jp2((icount)),
     $                          jm1((icount)),jm2((icount)),
     $                          ' from',ip1p,ip2p
                        endif
 100                    continue
                     endif
                  enddo
                  enddo
               endif
c End of extrapolation section.
            endif
         enddo
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
      real u(*)
      integer iuinc(3)
      real xf(3)
      real xn(*)
      real f(2,2)
      integer iflags(2,2)
      real d(2)


c      write(*,*)'get3 in',ndims,iuinc
c Values of u at the points to be interpolated.
c      u0=u(1+(ix-1)*iuinc)
c      up=u(1+ix    *iuinc)
c      um=u(1+(ix-2)*iuinc)
      dxp=xn(2)-xn(1)
      dxm=xn(1)-xn(0)
      do id1=1,3
         do id2=id1+1,3
            if(id1.ne.idf .and. id2.ne.idf)then
               d(1)=xf(id1)
               d(2)=xf(id2)
               do ip1=0,1
               do ip2=0,1
                  u0=u(1+ip1*iuinc(id1)+ip2*iuinc(id2))
                  up=u(1+ip1*iuinc(id1)+ip2*iuinc(id2)+iuinc(idf))
                  um=u(1+ip1*iuinc(id1)+ip2*iuinc(id2)-iuinc(idf))
                  ugradp=(up-u0)/dxp
                  ugradm=(u0-um)/dxm
                  ug=(0.5-xf(idf))*ugradm + (0.5+xf(idf))*ugradp
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
