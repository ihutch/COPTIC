c*********************************************************************
      subroutine sphereplot(objg)
c Plot the sphere objg on an already set up 3D scene.
      include '3dcom.f'
      real objg(odata)
      real xe(ns_ndims)
      parameter (nadef=20,ncosdef=20,pi=3.141593)
      parameter (ncorn=5)
      real rface(ncorn,ns_ndims)
      integer wp(ncorn),wc(ncorn)
      data wp/1,0,0,1,1/wc/0,0,1,1,0/

      ish1=objg(ofn1)
      ish2=objg(ofn2)
c Use position arrays compatible with flux array but not too coarse.
      if(ish1.gt.nadef)then 
         ncos=ish1
      elseif(ish1.gt.0)then 
         ncos=ish1*nint(float(ncosdef)/ish1)
      else
         ncos=ncosdef
         ish1=ncos
      endif
      if(ish2.gt.nadef)then
         nangle=ish2
      elseif(ish2.gt.0)then
         nangle=ish2*nint(float(nadef)/ish2)
      else
         nangle=nadef
         ish2=ncos
      endif
c      write(*,*)ish1,ish2,nangle,ncos
c Discover perspective, eye position in world coords:
      call trn32(x,y,z,xe(1),xe(2),xe(3),-1)
      call nxyz2wxyz(xe(1),xe(2),xe(3),xe(1),xe(2),xe(3))
c We are going to draw a set of cells at azimuthal angles. 
c Cell centers are at 2pi*(i-0.5)/nangle  i=1,nangles
c The eye is at angle
      psie=atan2(xe(2),xe(1))
      psio=mod((psie+pi),2.*pi)
      fpo=(nangle*(psio/(2.*pi)+0.5))
      ipo=int(fpo)
      ipp=nint(fpo)-ipo
c      write(*,*)nangle,ipo,psie/pi,psio/pi
c We should start at angle furthest away, and plot this eye angle last.
c This involves jumping from side to side. 
      do ia=1,nangle
         isign=2*mod(ia+ipp,2)-1
         i=mod(isign*(ia/2)+ipo+nangle,nangle)+1
         p1=2.*pi*(i-1)/nangle-pi
         p2=2.*pi*(i  )/nangle-pi
         itc=int(ish2*(p1/(2.*pi)+0.500001))
         sp1=sin(p1)
         sp2=sin(p2)
         cp1=cos(p1)
         cp2=cos(p2)
c         write(*,'(2i3,2f8.4)')ia,i,cp1,sp1
         do j=1,ncos
c Cosines are -1.+ 2.*(i-0.5)/ncos, i=1,ncos
c Maybe also needs to be done in order. But for now just do simply.
            c1=(-1.+2.*(j-1)/ncos)
            s1=sqrt(1.-c1**2)
            c2=(-1.+2.*(j  )/ncos)
            s2=sqrt(1.-c2**2)
c z=rcos. x=r sin(t)cos(p), y= rsin(t)sin(p)
            do k=1,ncorn
               rface(k,1)=objg(oradius)*(wp(k)*cp1+(1.-wp(k))*cp2)
     $              *(wc(k)*s1+(1.-wc(k))*s2)+objg(ocenter)
               rface(k,2)=objg(oradius+1)*(wp(k)*sp1+(1.-wp(k))*sp2)
     $              *(wc(k)*s1+(1.-wc(k))*s2)+objg(ocenter+1)
               rface(k,3)=objg(oradius+2)*(wc(k)*c1+(1.-wc(k))*c2)
     $              +objg(ocenter+2)
            enddo
c Now rface contains the positions of the face corners.
c Plot it
c            write(*,*)i,j,c1,c2
c Here we might want to access the fluxdata and color according to
c its value. For now color according to the fluxdata position. 
            ifn1=nint(ish1*(c1+1.)/2.+.5000)
c            call color(mod(itc,7)+8*mod(ifn1,2)+1)
            call color(mod(itc+mod(ifn1,2),7)+7*mod(ifn1,2)+1)
            call poly3line(rface(1,1),rface(1,2),rface(1,3),ncorn)
            call pathfill
            call color(15)
         enddo
      enddo

      end
c*********************************************************************
      subroutine cubeplot(objg)
c Plot faces of cube object on an already set up 3D scene.
      include '3dcom.f'
      real objg(odata)
      parameter (ncorn=5)
      real rface(ncorn,ns_ndims),rfc(ns_ndims)
      real xe(ns_ndims)
      integer iov(ns_ndims)
      data iov/-1,-1,-1/
      integer iof(ncorn,ns_ndims-1)
      data iof/-1,1,1,-1,-1,   -1,-1,1,1,-1/

      do i=1,ns_ndims
         if(objg(ofn1+i-1).eq.0.)objg(ofn1+i-1)=1.
      enddo
c Decide the order of face drawing.
c Get the point and eye position. Make into world units.
      call trn32(x,y,z,xe(1),xe(2),xe(3),-1)
      call nxyz2wxyz(xe(1),xe(2),xe(3),xe(1),xe(2),xe(3))
c Transform eye position into fractional cube position.      
c The starting fixed point is opposite signs from returned xn.
      do iv=1,ns_ndims
         iov(iv)=-sign(1.,xe(iv)-objg(ocenter+iv-1))
      enddo
c Now ordered.
c      write(*,*)(objg(k),k=1,4*ns_ndims),objg(ocenter),objg(oradius)
      do is=1,2*ns_ndims
         iv=mod(is-1,ns_ndims)+1
         do id=1,ns_ndims
            rfc(id)=objg(ocenter+id-1)
            if(id.eq.iv)rfc(id)=rfc(id)+iov(iv)*objg(oradius+id-1)
         enddo
         iov(iv)=-iov(iv)
c         write(*,*)'Face center',rfc
         i1=mod(1+is-2,ns_ndims)+1
         i2=mod(2+is-2,ns_ndims)+1
         i3=mod(3+is-2,ns_ndims)+1
         do k2=1,objg(ofn1+i2-1)
            f2=(k2-.5)*2.-objg(ofn1+i2-1)
            do k3=1,objg(ofn1+i3-1)
               f3=(k3-.5)*2.-objg(ofn1+i3-1)
               do ic=1,ncorn
c Set the corner offsets for this face, is.
                  rface(ic,i1)=rfc(i1)
                  rface(ic,i2)=rfc(i2)+(f2+iof(ic,1))
     $                 *objg(oradius+i2-1)/objg(ofn1+i2-1)
                  rface(ic,i3)=rfc(i3)+(f3+iof(ic,2))
     $                 *objg(oradius+i3-1)/objg(ofn1+i3-1)
               enddo
               call color(is+mod(k2+k3,2)*8)
               call poly3line(rface(1,1),rface(1,2),rface(1,3),ncorn)         
               call pathfill()
               call color(15)
            enddo
         enddo
      enddo
      end
c*********************************************************************
      subroutine cylplot(objg)
c Plot divided faces of cylinder object objg.
      include '3dcom.f'
      parameter (nadef=20,nzdef=5,pi=3.141593)
      parameter (ncorn=5)
      real rface(ncorn,ns_ndims)
      integer wp(ncorn),wc(ncorn)
      data wp/1,0,0,1,1/wc/0,0,1,1,0/
      real xe(ns_ndims)
      real objg(odata)


c Use position arrays compatible with flux array but not too coarse.
      if(objg(ofn2).gt.nadef)then
         nangle=objg(ofn2)
      elseif(objg(ofn2).gt.0)then
         nangle=objg(ofn2)*nint(float(nadef)/objg(ofn2))
      else
         nangle=nadef
      endif
      if(objg(ofn3).gt.0.)then 
         nz=objg(ofn3)
      else
         nz=nzdef
      endif
      nr=1
      if(objg(ofn1).gt.0.)nr=objg(ofn1)

      ia=objg(ocylaxis)
c Discover perspective, eye position in world coords:
      call trn32(x,y,z,xe(1),xe(2),xe(3),-1)
      call nxyz2wxyz(xe(1),xe(2),xe(3),xe(1),xe(2),xe(3))
c Flux accumulation is (zero based):
c         it=int(nf_dimlens(nf_flux,infobj,2)
c     $     *(0.999999*theta/3.1415926+1.)*0.5)
      thetae=atan2(xe(mod(ia+1,3)+1),xe(mod(ia,3)+1))
      thetao=mod((thetae+pi),2.*pi)
      fto=(nangle*(thetao/(2.*pi)+0.5))
      ito=int(fto)
      itp=nint(fto)-ito
      ie=sign(1.,xe(ia))
c      write(*,*)(objg(k),k=1,odata)
c Draw curved surface.
      do it=1,nangle
         isign=2*mod(it+itp,2)-1
         i=mod(isign*(it/2)+ito+nangle,nangle)+1
         t1=2.*pi*(i-1)/nangle-pi
         t2=2.*pi*(i  )/nangle-pi
         itc=int(objg(ofn2)*(t1/(2.*pi)+0.500001))
         do j=1,nz
            z1=(-1.+(j-1)*2./nz)*objg(oradius+ia-1)
            z2=(-1.+(j  )*2./nz)*objg(oradius+ia-1)
            do k=1,ncorn
               ids=mod(1+ia-2,ns_ndims)+1
               rface(k,ids)=(wc(k)*z1+(1.-wc(k))*z2)+objg(ocenter+ids-1)
               ids=mod(2+ia-2,ns_ndims)+1
               rface(k,ids)=objg(oradius+ids-1)
     $              *cos(wp(k)*t1+(1.-wp(k))*t2)+objg(ocenter+ids-1)
               ids=mod(3+ia-2,ns_ndims)+1
               rface(k,ids)=objg(oradius+ids-1)
     $              *sin(wp(k)*t1+(1.-wp(k))*t2)+objg(ocenter+ids-1)
            enddo
c Now rface contains the coordinates of the face corners.            
c            call color(mod(itc,7)+8*mod(j,2)+1)
            call color(mod(itc+mod(j,2),7)+7*mod(j,2)+1)
            call poly3line(rface(1,1),rface(1,2),rface(1,3),ncorn)
            call pathfill
            call color(15)
         enddo
      enddo
c Draw visible end surface, equal spacing in r^2.
      do i=1,nangle
         t1=2.*pi*(i-1)/nangle-pi
         t2=2.*pi*(i  )/nangle-pi
         itc=int(objg(ofn2)*(t1/(2.*pi)+0.50001))
         do j=1,nr
            r1=sqrt(float(j-1)/nr)
            r2=sqrt(float(j  )/nr)
            do k=1,ncorn
               rface(k,ia)=ie*objg(oradius+ia-1)+objg(ocenter+ia-1)
               ids=mod(2+ia-2,ns_ndims)+1
               rface(k,ids)=objg(oradius+ids-1)*(wc(k)*r1+(1.-wc(k))*r2)
     $              *cos(wp(k)*t1+(1.-wp(k))*t2)+objg(ocenter+ids-1)
               ids=mod(3+ia-2,ns_ndims)+1
               rface(k,ids)=objg(oradius+ids-1)*(wc(k)*r1+(1.-wc(k))*r2)
     $              *sin(wp(k)*t1+(1.-wp(k))*t2)+objg(ocenter+ids-1)
            enddo
            call color(mod(itc+mod(j,2),7)+7*mod(j,2)+1)
            call poly3line(rface(1,1),rface(1,2),rface(1,3),ncorn)
            call pathfill
            call color(15)
         enddo
      enddo

      end
c********************************************************************
      subroutine pllelplot(objg)
c Plot divided faces of parallelopiped object objg.
      include '3dcom.f'
      real objg(odata)
      real xe(ns_ndims),xn(ns_ndims)
      integer iov(ns_ndims)
       parameter (ncorn=5)
      real rface(ncorn,pp_ndims),rfc(pp_ndims)
      real xi(pp_ndims)
      integer iof(ncorn,pp_ndims-1)
      data iof/-1,1,1,-1,-1,   -1,-1,1,1,-1/
c Decide the order of face drawing. Furthest to nearest.
c Get the point and eye position. Make into world units.
      call trn32(x,y,z,xe(1),xe(2),xe(3),-1)
      call nxyz2wxyz(xe(1),xe(2),xe(3),xe(1),xe(2),xe(3))
c Transform eye position into fractional cube position.
      call pllelfrac(xe,xn,objg)
c The starting fixed point is opposite signs from returned xn.
c So the first three center vectors have values equal to minus
c the sign of xn times the three pp_vectors. 
      do iv=1,ns_ndims
         iov(iv)=-sign(1.,xn(iv))
c Fix zero flux meshes:
         if(objg(ofn1+iv-1).eq.0.)objg(ofn1+iv-1)=1.
      enddo
      
      do is=1,2*pp_ndims
         iv=mod(is-1,pp_ndims)+1
         do id=1,pp_ndims
            rfc(id)=objg(pp_orig+id-1)
     $           +iov(iv)*objg(pp_vec+(iv-1)*pp_ndims+id-1)
         enddo
         iov(iv)=-iov(iv)
         i1=mod(1+is-2,pp_ndims)+1
         i2=mod(2+is-2,pp_ndims)+1
         i3=mod(3+is-2,pp_ndims)+1
         do k2=1,objg(ofn1+i2-1)
            f2=(k2-.5)*2.-objg(ofn1+i2-1)
            do k3=1,objg(ofn1+i3-1)
               f3=(k3-.5)*2.-objg(ofn1+i3-1)
               do ic=1,ncorn
                  xi(i1)=0
                  xi(i2)=(f2+iof(ic,1))/objg(ofn1+i2-1)
                  xi(i3)=(f3+iof(ic,2))/objg(ofn1+i3-1)
                  do id=1,pp_ndims
                     rface(ic,id)=rfc(id)
                     do ivv=1,pp_ndims
                        rface(ic,id)=rface(ic,id)
     $                       +xi(ivv)*objg(pp_vec+(ivv-1)*pp_ndims+id-1)
                     enddo
                  enddo
               enddo
c Here rface(ic,id) contains the id'th coordinate of ic'th corner.
               call color(is+mod(k2+k3,2)*8)
               call poly3line(rface(1,1),rface(1,2),rface(1,3),ncorn)
               call pathfill()
               call color(15)
            enddo
         enddo
      enddo
      end
c*******************************************************************
c Plot edges/faces of objects. Window size rs.
      subroutine objplot(ndims,rs,iosw)
      integer ndims,iosw
      include '3dcom.f'
      include 'sectcom.f'
      integer index(ngeomobjmax)
      real zta(ngeomobjmax)
      data iprinting/0/
      irotating=0
      call pltinit(0.,1.,0.,1.)
      call setcube(.2,.2,.2,.5,.4)
 51   continue
      if(iprinting.ne.0)call pfset(3)
      call geteye(x2,y2,z2)
      call pltinit(0.,1.,0.,1.)
      call scale3(-rs,rs,-rs,rs,-rs,rs)
      call trn32(0.,0.,0.,x2,y2,z2,1)
      if(irotating.eq.0)then
         icorner=igetcorner()
         call ax3labels('x','y','z')
      else
         irotating=irotating-1
      endif
c      call cubed(icorner)
      call axproj(icorner)

c Decide the order in which to draw objects, based on the position of
c their centers. 
      do i=1,ngeomobj
         index(i)=i
c Get the position in view coordinates.
         call trn32(obj_geom(ocenter,i),obj_geom(ocenter+1,i),
     $        obj_geom(ocenter+2,i),xt,yt,zta(i),3)
      enddo
      call zsort(ngeomobj,zta,index)
c      write(*,*)(index(k),zta(k),k=1,ngeomobj)

c Do drawing in order
      do ik=1,ngeomobj
         iobj=index(ik)
         if(obj_geom(otype,iobj).eq.1.)then
            call sphereplot(obj_geom(1,iobj))
         elseif(obj_geom(otype,iobj).eq.2.)then
            call cubeplot(obj_geom(1,iobj))
         elseif(obj_geom(otype,iobj).eq.3.)then
            call cylplot(obj_geom(1,iobj))
         elseif(obj_geom(otype,iobj).eq.4.)then
            call pllelplot(obj_geom(1,iobj))
         endif
      enddo
      if(.true.)then
c Plot intersections
      call charsize(.005,.005)
c      write(*,*)'Intersections:'
      do i=1,sc_ipt
         call vec3w(x_sc(1,1,i),x_sc(2,1,i),x_sc(3,1,i),0)
         call vec3w(x_sc(1,2,i),x_sc(2,2,i),x_sc(3,2,i),1)
         call drcstr('!A1!@')
c         write(*,*)i,(x_sc(k,1,i),k=1,3)
      enddo
      call charsize(0.,0.)
      endif

c User interface:
      iprinting=0
 56   call eye3d(isw)
      if(isw.eq.0.or.isw.eq.ichar('q'))goto 55
      if(isw.eq.ichar('r').or.isw.eq.ichar('e'))irotating=1
      if(isw.eq.ichar('p'))iprinting=mod(iprinting+1,2)
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
         goto 51
      endif
      goto 51
 55   continue
      end
c*********************************************************************
      subroutine zsort(ngeomobj,zta,index)
      real zta(ngeomobj)
      integer index(ngeomobj)
      do j=2,ngeomobj
         a1=zta(j)
         a2=index(j)
         do i=j-1,1,-1
            if(zta(i).ge.a1)goto 10
            zta(i+1)=zta(i)
            index(i+1)=index(i)
         enddo
         i=0
 10      zta(i+1)=a1
         index(i+1)=a2
      enddo
      end
