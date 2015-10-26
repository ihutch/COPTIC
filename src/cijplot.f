c Mostly for testing.
c Assumed 3-D routine, plots representation of the cij/obj data.
      subroutine  cijplot(ifull,iuds,cij,rs,iosw)
      include 'ndimsdecl.f'
      integer iosw
      integer ifull(ndimsmax),iuds(ndimsmax)
      real cij(ndims*2+1,ifull(1),ifull(2),ifull(3))
      include 'objcom.f'
      include 'meshcom.f'
      real xx(3),xt(3),xsize(3)
      integer irx(5),iry(5),ipx(5),ipy(5),ijk(3)
      integer idelta(3,7)
c      character*40 mystring
      integer iprinting
      data irx/1,0,-1,0,1/
      data iry/0,1,0,-1,0/
      data ipx/0,0,1,0,0/
      data ipy/0,0,0,1,0/
      data idelta/1,0,0,0,1,0,0,0,1, 1,1,1,  0,1,1,1,0,1,1,1,0/
      data iprinting/0/

      if(ndims.ne.3)then
         write(*,*)'cijplot called with ndims=',ndims
         write(*,*)'Skipping this only-3D routine.'
         return
      endif

      ixud=iuds(1)
      iyud=iuds(2)
      izud=iuds(3)
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
      write(*,'(a,i3)')'Object Mask:',mysw/2
      xsm=0.
      do id=1,3
         xsize(id)=xmeshend(id)-xmeshstart(id)
         if(abs(xsize(id)).gt.xsm)xsm=abs(xsize(id))
      enddo
      xs=xsm/.25
      call pltinit(0.,1.,0.,1.)
      call setcube(xsize(1)/xs,xsize(2)/xs,xsize(3)/xs,.5,.4)
c      call setcube(.2,.2,.2,.5,.4)
 51   continue
      if(iprinting.ne.0)call pfset(3)
      call geteye(x2,y2,z2)
      call pltinit(0.,1.,0.,1.)
      call scale3(xmeshstart(1),xmeshend(1),xmeshstart(2),xmeshend(2),
     $     xmeshstart(3),xmeshend(3))
      if(rs.ne.0)call scale3(-rs,rs,-rs,rs,-rs,rs)
      call trn32(0.,0.,0.,x2,y2,z2,1)
      xb=0
      yb=0
      zb=0
      if(x2.ge.0)xb=1
      if(y2.ge.0)yb=1
      if(z2.ge.0)zb=1
      if(irotating.eq.0)then
c This is slightly different.
c         icorner= int((2*zb-1)*( (1 +3*yb) + (1 - 2*yb)*xb ))
         icorner=igetcorner()
         call ax3labels('x','y','z')
      else
         irotating=irotating-1
      endif
c      call cubed(icorner)
      call axproj(icorner)

      if(istick.eq.0) goto 52
c Stick drawing.
      do i=1,ixud
         do j=1,iyud
            do k=1,izud
               iobj=int(cij(2*ndims+1,i,j,k))
               if(iobj.ne.0)then
c Draw joins from the point to the fraction distance.
                  xx(1)=xn(ixnp(1)+i)
                  xx(2)=xn(ixnp(2)+j)
                  xx(3)=xn(ixnp(3)+k)
c                  call vec3w(xx(1),xx(2),xx(3),0)
                  do id=1,ndims
                     do ipm=1,2
                        ispm=1-2*(ipm-1)
                        no=ndata_cij*(2*(id-1)+(ipm-1))+1
                        frac=dob_cij(no,iobj)
                        if(frac.lt.1. .and. frac.gt.0.)then
c This a true intersection.
                           call vec3w(xx(1),xx(2),xx(3),0)
                           call color((no+2)/ndata_cij)
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
c                           call vec3w(xx(1),xx(2),xx(3),0)
                        endif
                     enddo
                  enddo
               endif
            enddo
         enddo
      enddo
c      if(iwire.eq.0 .and. ieye3d().ne.0)goto 51
      if(iwire.eq.0)goto 56
      if(istick.eq.1)iwire=1
c     Wireframe drawing.
 52   continue
c      write(mystring,'(a)')'Object Mask:'
c      call iwrite(mysw/2,ilen,mystring(13:))
c      call jdrwstr(.1,.1,mystring,1.)
      do i=1,ixud
         do j=1,iyud
            do k=1,izud
               iobj=int(cij(2*ndims+1,i,j,k))
               iobjcode=idob_cij(iinter_cij,iobj)
               if(iobj.ne.0.and.iobjcode.ge.0)then
               isob=int(mysw/2**iobjcode-(mysw/2**(iobjcode+1))*2)
               if(iobjcode.ne.igetcolor().and.iobjcode.ne.0)then
c                  write(*,*)'Calling color',iobj,iobjcode,isob,mysw
                  call color(iobjcode)
               endif
               if(isob.ne.0)then
c Origin point
                  xx(1)=xn(ixnp(1)+i)
                  xx(2)=xn(ixnp(2)+j)
                  xx(3)=xn(ixnp(3)+k)
                  ijk(1)=i
                  ijk(2)=j
                  ijk(3)=k
c This assumes 3-D.
c Draw joins between fractions common to a single node
c in plane normal to id.
                  do id=1,ndims
                     i1=mod(id,3)+1
                     i2=mod(id+1,3)+1
                     ipen=0
                     do ipm=1,5       
                        no=ndata_cij*(
     $                      2*( (i1-1)*mod(ipm,2)+(i2-1)*mod(ipm+1,2))
     $                       +ipx(ipm)+ipy(ipm))+1
c                        write(*,*)'iinter',idob_cij(iinter_cij,iobj)
                        frac=dob_cij(no,iobj)
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
c Draw joins between intersections of adjacent parallel lattice legs
                  do id=1,ndims
c                     id=mod(idd-1,3)+1
                     iobj2=int(cij(2*ndims+1,
     $                    i+idelta(1,id),
     $                    j+idelta(2,id),
     $                    k+idelta(3,id)))
                     if(iobj2.ne.0)then
                     i1=mod(id,3)+1
                     i2=mod(id+1,3)+1
                     ipen=0
                     do ipm=1,5       
                        no=ndata_cij*(
     $                      2*( (i1-1)*mod(ipm,2)+(i2-1)*mod(ipm+1,2))
     $                       +ipx(ipm)+ipy(ipm))+1
                        frac=dob_cij(no,iobj)
                        frac2=dob_cij(no,iobj2)
                        if(frac.lt.1. .and. frac.ge.0.
     $                       .and. frac2.lt.1. .and. frac2.ge.0.)then
c Double intersections. Connect them.
c                        call color(iobjcode)
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
c Plot orbits if there are any. (Dummy in solver.)
      call orbit3plot()
c User interface:
      iprinting=0
 56   call eye3d(isw)
      if(isw.eq.0.or.isw.eq.ichar('q'))goto 55
      if(isw.eq.ichar('r').or.isw.eq.ichar('e'))irotating=1
      if(isw.eq.ichar('p'))iprinting=mod(iprinting+1,2)
      if(isw.eq.ichar('h'))then
         write(*,*)'r,e: rotate;  p:print;  1,2,3: halve region;  ',
     $        'x,y,z: shrink region;  ','a,b,c: grow region'
      endif
      if(isw.eq.ichar('1'))then
         if(ixud.eq.iuds(1))then
            ixud=ixud/2
         else
            ixud=iuds(1)
         endif
      elseif(isw.eq.ichar('2'))then
         if(iyud.eq.iuds(2))then
            iyud=iyud/2
         else
            iyud=iuds(2)
         endif
      elseif(isw.eq.ichar('3'))then
         if(izud.eq.iuds(3))then
            izud=izud/2
         else
            izud=iuds(3)
         endif
      elseif(isw.eq.ichar('a').and.ixud.lt.iuds(1))then
         ixud=ixud+1
      elseif(isw.eq.ichar('b').and.iyud.lt.iuds(2))then
         iyud=iyud+1
      elseif(isw.eq.ichar('c').and.izud.lt.iuds(3))then
         izud=izud+1
      elseif(isw.eq.ichar('x').and.ixud.gt.1)then
         ixud=ixud-1
      elseif(isw.eq.ichar('y').and.iyud.gt.1)then
         iyud=iyud-1
      elseif(isw.eq.ichar('z').and.izud.gt.1)then
         izud=izud-1
      endif
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
c         irotating=irotating-1
         goto 51
      endif
      goto 51
 55   continue
      if(iwire.eq.0)then
         iwire=1
         goto 51
      endif
      if(istick.gt.0) then
         istick=0
         goto 51
      endif


      end
