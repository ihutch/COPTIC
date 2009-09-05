c Mostly for testing.
c Assumed 3-D routine, plots representation of the cij/obj data.
      subroutine  cijplot(ndims,ifull,iuds,cij,rs,iosw)
      integer ndims,iosw
      parameter (mdims=10)
      integer ifull(mdims),iuds(ndims)
      real cij(ndims*2+1,ifull(1),ifull(2),ifull(3))
      include 'objcom.f'
      include 'meshcom.f'
      include 'partcom.f'
      real xx(3),xt(3)
      integer irx(5),iry(5),ipx(5),ipy(5),ijk(3)
      integer idelta(3,3)
      character*40 mystring
      data irx/1,0,-1,0,1/
      data iry/0,1,0,-1,0/
      data ipx/0,0,1,0,0/
      data ipy/0,0,0,1,0/
      data idelta/1,0,0,0,1,0,0,0,1/

c If iosw > 0, then wireframe the objects corresponding to its bits.
c If iosw < 0, plot the sticks as well as the wireframes.
c If iosw = 0, plot only sticks for everything.
      if(iosw.gt.0)then
         istick=0
         iwire=1
         mysw=iosw*2
      elseif(iosw.lt.0)then
         istick=1
         iwire=0
         mysw=-iosw*2
      else
         istick=2
         iwire=0
         mysw=65535
      endif

c      istick=1
c      iwire=0
      irotating=0

      call pltinit(0.,1.,0.,1.)
      call setcube(.2,.2,.2,.5,.4)
 51   continue
      call geteye(x2,y2,z2)
      call pltinit(0.,1.,0.,1.)
      call scale3(-rs,rs,-rs,rs,-rs,rs)
      call trn32(0.,0.,0.,x2,y2,z2,1)
      xb=0
      yb=0
      zb=0
      if(x2.ge.0)xb=1
      if(y2.ge.0)yb=1
      if(z2.ge.0)zb=1
      icorner= int((2*zb-1)*( (1 +3*yb) + (1 - 2*yb)*xb ))

c      call cubed(icorner)
      call axproj(icorner)

      if(istick.eq.0) goto 52
c Stick drawing.
      do i=1,iuds(1)
         do j=1,iuds(2)
            do k=1,iuds(3)
               iobj=int(cij(2*ndims+1,i,j,k))
               if(iobj.ne.0)then
c Draw joins from the point to the fraction distance.
                  xx(1)=xn(ixnp(1)+i)
                  xx(2)=xn(ixnp(2)+j)
                  xx(3)=xn(ixnp(3)+k)
                  call vec3w(xx(1),xx(2),xx(3),0)
                  do id=1,ndims
                     do ipm=1,2
                        ispm=1-2*(ipm-1)
                        no=ndata_sor*(2*(id-1)+(ipm-1))+1
                        call color((no+2)/ndata_sor)
                        frac=dob_sor(no,iobj)
                        if(frac.lt.1. .and. frac.gt.0.)then
c This a true intersection.
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
      if(iwire.eq.0 .and. ieye3d().ne.0)goto 51
      if(istick.eq.1)iwire=1
c      goto 51
 52   continue
c     Wireframe drawing.
      write(mystring,'(a)')'Object Mask:'
      call iwrite(mysw/2,ilen,mystring(13:))
      call jdrwstr(.1,.1,mystring,1.)
      do i=1,iuds(1)
         do j=1,iuds(2)
            do k=1,iuds(3)
               iobj=int(cij(2*ndims+1,i,j,k))
               iobjcode=idob_sor(iinter_sor,iobj)
               if(iobj.ne.0.and.iobjcode.ge.0)then
               isob=int(mysw/2**iobjcode-(mysw/2**(iobjcode+1))*2)
c               write(*,*)iobj,iobjcode,isob,mysw
               if(isob.ne.0)then
c Origin point
                  xx(1)=xn(ixnp(1)+i)
                  xx(2)=xn(ixnp(2)+j)
                  xx(3)=xn(ixnp(3)+k)
                  ijk(1)=i
                  ijk(2)=j
                  ijk(3)=k
c This assumes 3-D.
c Draw joins between common fractions in plane normal to id.
                  do id=1,ndims
                     i1=mod(id,3)+1
                     i2=mod(id+1,3)+1
                     ipen=0
                     do ipm=1,5       
                        no=ndata_sor*(
     $                      2*( (i1-1)*mod(ipm,2)+(i2-1)*mod(ipm+1,2))
     $                       +ipx(ipm)+ipy(ipm))+1
c                        write(*,*)'iinter',idob_sor(iinter_sor,iobj)
                        call color(iobjcode)
                        frac=dob_sor(no,iobj)
                        if(frac.lt.1. .and. frac.ge.0.)then
c This a true intersection. Draw to it.
                           xt(id)=xx(id)
                           xt(i1)=xx(i1)+frac*
     $                          (xn(ixnp(i1)+ijk(i1)+irx(ipm))-xx(i1))
                           xt(i2)=xx(i2)+frac*
     $                          (xn(ixnp(i2)+ijk(i2)+iry(ipm))-xx(i2))
c                           write(*,*)no,id,i1,i2,frac,ipm,xt
                           call vec3w(xt(1),xt(2),xt(3),ipen)
                           ipen=1
                        else
                           ipen=0
                        endif
                     enddo
                  enddo
c Draw joins between adjacent planes in the positive direction
                  do id=1,ndims
                     iobj2=int(cij(2*ndims+1,
     $                    i+idelta(1,id),
     $                    j+idelta(2,id),
     $                    k+idelta(3,id)))
                     if(iobj2.ne.0)then
                     i1=mod(id,3)+1
                     i2=mod(id+1,3)+1
                     ipen=0
                     do ipm=1,5       
                        no=ndata_sor*(
     $                      2*( (i1-1)*mod(ipm,2)+(i2-1)*mod(ipm+1,2))
     $                       +ipx(ipm)+ipy(ipm))+1
                        call color(iobjcode)
                        frac=dob_sor(no,iobj)
                        frac2=dob_sor(no,iobj2)
                        if(frac.lt.1. .and. frac.ge.0.
     $                       .and. frac2.lt.1. .and. frac2.ge.0.)then
c Double intersections. Connect them.
                           xt(id)=xx(id)
                           xt(i1)=xx(i1)+frac*
     $                          (xn(ixnp(i1)+ijk(i1)+irx(ipm))-xx(i1))
                           xt(i2)=xx(i2)+frac*
     $                          (xn(ixnp(i2)+ijk(i2)+iry(ipm))-xx(i2))
c                           write(*,*)no,id,i1,i2,frac,ipm,xt
                           call vec3w(xt(1),xt(2),xt(3),0)
                           xt(id)=xn(ixnp(id)+ijk(id)+1)
                           xt(i1)=xx(i1)+frac2*
     $                          (xn(ixnp(i1)+ijk(i1)+irx(ipm))-xx(i1))
                           xt(i2)=xx(i2)+frac2*
     $                          (xn(ixnp(i2)+ijk(i2)+iry(ipm))-xx(i2))
c                           write(*,*)no,id,i1,i2,frac,ipm,xt
                           call vec3w(xt(1),xt(2),xt(3),1)
                        endif
                     enddo
                     endif
                  enddo
               endif
               endif
            enddo
         enddo
      enddo
c Plot orbits if there are any.
c      write(*,*)'norbits,length=',norbits,iorbitlen(1)
      do kk=1,norbits
         call color(kk)
         call poly3line(xorbit(1,kk),yorbit(1,kk),zorbit(1,kk),
     $        iorbitlen(kk))
         call poly3mark(xorbit(1,kk),yorbit(1,kk),zorbit(1,kk),
     $        iorbitlen(kk),1)
      enddo
c User interface:
      call eye3d(isw)
      if(isw.eq.0)goto 55
      if(isw.eq.ichar('r').or.isw.eq.ichar('e'))irotating=1
      if(irotating.gt.0)then
c Get back current eye position xe1 etc.
         call asleep(5000)
         call trn32(xe,ye,ze,x2,y2,z2,-1)
         dtheta=.02
         if(isw.eq.ichar('r'))dtheta=-.02
         cs=cos(dtheta)
         sn=sin(dtheta)
c         write(*,*)xe,ye,ze,x2,y2,z2,cs,sn
         xex=x2-xe
         yex=y2-ye
         x2=xe+cs*xex-sn*yex
         y2=ye+sn*xex+cs*yex
         call puteye(x2,y2,z2)
c         call trn32(xe,ye,ze,x2,y2,z2,1)
         irotating=irotating-1
         goto 51
      endif
      goto 51
 55   continue
      if(istick.gt.0) then
         istick=0
         goto 51
      endif


      end
