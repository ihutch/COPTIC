!*********************************************************************
      subroutine sphereplot(iq,objg,iobj,ioswin,fmin,fmax)
! Plot the sphere objg, number iobj, on an already set up 3D scene.
! iosw determines the nature of the plot
! 0: Color code according to position.
! 1:            according to average flux already in nf_step+1
! 2:            according to average flux-density already in nf_step+2

      include 'ndimsdecl.f'
      include '3dcom.f'
      real objg(odata)
      real fmin,fmax
      real xe(ndimsmax)
      parameter (nadef=20,ncosdef=20,pi=3.141593)
      parameter (ncorn=5)
      real rface(ncorn,ndimsmax)
      integer wp(ncorn),wc(ncorn)
!      character*20 string
      logical lfw
      data wp/1,0,0,1,1/wc/0,0,1,1,0/

      ish1=int(objg(ofn1))
      ish2=int(objg(ofn2))
      ism1=1
      ism2=1
      ifobj=nf_map(iobj)
      iosw=ioswin
      if(ifobj.eq.0)iosw=0
! No flux for this object
      iav=nf_address(iq+if_quant(ifobj,if_species)
     $     ,ifobj,nf_step+abs(iosw))
! Find max and min of flux
      if(fmax.eq.fmin)then
         call minmax(ff_data(iav),nf_posno(iq,ifobj),fmin,fmax)
         fmax=1.03*fmax
         if(fmin.gt.0.)fmin=0.
      endif
! Use position arrays compatible with flux array but not too coarse.
      if(ish1.gt.ncosdef)then 
         ncos=ish1
      elseif(ish1.gt.0)then
         ism1=nint(float(ncosdef)/ish1)
         ncos=ish1*ism1
      else
         ncos=ncosdef
         ish1=ncos
      endif
      if(ish2.gt.nadef)then
         nangle=ish2
      elseif(ish2.gt.0)then
         ism2=nint(float(nadef)/ish2)
         nangle=ish2*ism2
      else
         nangle=nadef
         ish2=ncos
      endif
!      write(*,*)ish1,ish2,nangle,ncos
! Discover perspective, eye position in world coords:
      call trn32(x,y,z,xe(1),xe(2),xe(3),-1)
      call nxyz2wxyz(xe(1),xe(2),xe(3),xe(1),xe(2),xe(3))
! We are going to draw a set of cells at azimuthal angles. 
! Cell centers are at 2pi*(i-0.5)/nangle  i=1,nangles
! The eye is at angle
      psie=atan2(xe(2),xe(1))
      psio=mod((psie+pi),2.*pi)
      fpo=(nangle*(psio/(2.*pi)+0.5))
      ipo=int(fpo)
      ipp=nint(fpo)-ipo
!      write(*,*)nangle,ipo,psie/pi,psio/pi
! We should start at angle furthest away, and plot this eye angle last.
! This involves jumping from side to side. 
      do ia=1,nangle
         isign=2*mod(ia+ipp,2)-1
         i=mod(isign*(ia/2)+ipo+nangle,nangle)+1
         p1=2.*pi*(i-1)/nangle-pi
         p2=2.*pi* i   /nangle-pi
         itc=int(ish2*(p1/(2.*pi)+0.500001))
         sp1=sin(p1)
         sp2=sin(p2)
         cp1=cos(p1)
         cp2=cos(p2)
!         write(*,'(2i3,2f8.4)')ia,i,cp1,sp1
         do j=1,ncos
! Cosines are -1.+ 2.*(i-0.5)/ncos, i=1,ncos
! Maybe also needs to be done in order. But for now just do simply.
            c1=(-1.+2.*(j-1.)/ncos)
            s1=sqrt(1.-c1**2)
            c2=(-1.+2.* j   /ncos)
            s2=sqrt(1.-c2**2)
! z=rcos. x=r sin(t)cos(p), y= rsin(t)sin(p)
!            write(*,*)c1,c2,ncos
            do k=1,ncorn
               rface(k,1)=objg(oradius)*(wp(k)*cp1+(1.-wp(k))*cp2)
     $              *(wc(k)*s1+(1.-wc(k))*s2)+objg(ocenter)
               rface(k,2)=objg(oradius+1)*(wp(k)*sp1+(1.-wp(k))*sp2)
     $              *(wc(k)*s1+(1.-wc(k))*s2)+objg(ocenter+1)
               rface(k,3)=objg(oradius+2)*(wc(k)*c1+(1.-wc(k))*c2)
     $              +objg(ocenter+2)
            enddo
!            write(*,*)rface
!            write(*,*)ish1,ism1,ism2
! Now rface contains the positions of the face corners.
            k2=nint(ish1*(c1+1.)/2.+.50001)
            lfw=(mod(i+ism2,ism2).eq.ism2/2 .and.
     $           mod(j+ism1,ism1).eq.ism1/2)
!            write(*,*)'Calling facecolor',iosw
            call facecolor(iosw,1,k2,itc+1,iobj,iav,rface,fmin
     $              ,fmax,1,lfw,isign)               
         enddo
      enddo
      end
!*********************************************************************
      subroutine cubeplot(iq,objg,iobj,ioswin,fmin,fmax)
! Plot faces of cube object on an already set up 3D scene.
! iosw determines the nature of the plot
! 0: Color code according to position.
! 1:            according to average flux already in nf_step+1
! 2:            according to average flux-density already in nf_step+2

      include 'ndimsdecl.f'
      include '3dcom.f'
      real objg(odata)
!      character*20 string
      integer iov(ndimsmax)
      real xe(ndimsmax),objn1(0:ndimsmax-1)
      parameter (ncorn=5)
      real rface(ncorn,ndimsmax),rfc(ndimsmax)
      integer iof(ncorn,ndimsmax-1)
      data iof/-1,1,1,-1,-1,   -1,-1,1,1,-1/

      ifobj=nf_map(iobj)
      iosw=ioswin
      if(ifobj.eq.0)iosw=0

      iav=nf_address(iq+if_quant(ifobj,if_species)
     $     ,ifobj,nf_step+abs(iosw))
      call minmax(ff_data(iav),nf_posno(1,ifobj),fmin,fmax)
      if(fmin.gt.0.)fmin=0.

! Decide the order of face drawing.
! Get the point and eye position. Make into world units.
      call trn32(x,y,z,xe(1),xe(2),xe(3),-1)
      call nxyz2wxyz(xe(1),xe(2),xe(3),xe(1),xe(2),xe(3))
! Transform eye position into fractional cube position.      
! The starting fixed point is opposite signs from returned xn.
      do iv=1,ndims
! The divisions in different face directions:
         objn1(iv-1)=objg(ofn1+iv-1)
! Fix zero flux meshes:
         if(objn1(iv-1).eq.0.)objn1(iv-1)=1.
! The face sign to start with (furthest away):
         iov(iv)=int(-sign(1.,xe(iv)-objg(ocenter+iv-1)))
      enddo
! Now ordered.
!      write(*,*)(objg(k),k=1,4*ndims),objg(ocenter),objg(oradius)
      do is=1,2*ndims
! Face normal direction:
         iv=mod(is-1,ndims)+1
! Face index iv+ 0 or ndims depending on sign iov:
         imin=iv+ndims*(1+iov(iv))/2
         do id=1,ndims
! Set the center position of the face rfc(1:3)
            rfc(id)=objg(ocenter+id-1)
            if(id.eq.iv)rfc(id)=rfc(id)+iov(iv)*objg(oradius+id-1)
         enddo
         iov(iv)=-iov(iv)
!         write(*,*)'Face center',rfc
! The indices of objn1 (divisions) for this face
         i1=mod(1+is-2,ndims)+1
         i2=mod(2+is-2,ndims)+1
         i3=mod(3+is-2,ndims)+1
         do k2=1,int(objn1(i2-1))
!  fs run from 1-N to N-1 as ks run from 1 to N
            f2=2*k2-1.-objn1(i2-1)
            do k3=1,int(objn1(i3-1))
               f3=2*k3-1.-objn1(i3-1)
               do ic=1,ncorn
! Set the corner offsets for this face, is.
                  rface(ic,i1)=rfc(i1)
                  rface(ic,i2)=rfc(i2)+(f2+iof(ic,1))
     $                 *objg(oradius+i2-1)/objn1(i2-1)
                  rface(ic,i3)=rfc(i3)+(f3+iof(ic,2))
     $                 *objg(oradius+i3-1)/objn1(i3-1)
               enddo
               call facecolor(iosw,imin,k2,k3,iobj,iav,rface,fmin,fmax
     $              ,i2,.true.,-1)
           enddo
         enddo
      enddo
      end
!*********************************************************************
      subroutine cylplot(iq,objg,iobj,ioswin,fmin,fmax)
! Plot divided faces of cylinder object objg.
! iosw determines the nature of the plot
! 0: Color code according to position.
! 1:            according to average flux already in nf_step+1
! 2:            according to average flux-density already in nf_step+2

      include 'ndimsdecl.f'
      include '3dcom.f'
      parameter (nadef=20,nzdef=5,pi=3.141593)
      parameter (ncorn=5)
      real rface(ncorn,ndimsmax)
      real xe(ndimsmax)
      real objg(odata)
      logical lfw
      integer wp(ncorn),wc(ncorn)
      data wp/1,0,0,1,1/wc/0,0,1,1,0/

      ifobj=nf_map(iobj)
      iosw=ioswin
      if(ifobj.eq.0)iosw=0
      iav=nf_address(iq+if_quant(ifobj,if_species)
     $     ,ifobj,nf_step+abs(iosw))
      call minmax(ff_data(iav),nf_posno(1,ifobj),fmin,fmax)
      if(fmin.gt.0.)fmin=0.

      ism1=1
      ism2=1
! Use position arrays compatible with flux array but not too coarse.
      if(objg(ofn2).gt.nadef)then
         nangle=int(objg(ofn2))
      elseif(objg(ofn2).gt.0)then
         nangle=int(objg(ofn2))*nint(float(nadef)/objg(ofn2))
         ism2=nint(float(nadef)/objg(ofn2))
      else
         nangle=nadef
      endif
      if(objg(ofn3).gt.0.)then 
         nz=int(objg(ofn3))
      else
         nz=nzdef
      endif
      nr=1
      if(objg(ofn1).gt.0.)nr=int(objg(ofn1))
      ia=int(objg(ocylaxis))
! Discover perspective, eye position in world coords:
      call trn32(x,y,z,xe(1),xe(2),xe(3),-1)
      call nxyz2wxyz(xe(1),xe(2),xe(3),xe(1),xe(2),xe(3))
! Flux accumulation is (zero based):
!         it=int(nf_dimlens(nf_flux,infobj,2)
!     $     *(0.999999*theta/3.1415926+1.)*0.5)
      thetae=atan2(xe(mod(ia+1,3)+1),xe(mod(ia,3)+1))
      thetao=mod((thetae+pi),2.*pi)
      fto=(nangle*(thetao/(2.*pi)+0.5))
      ito=int(fto)
      itp=nint(fto)-ito
! int(sign(1.,xe(ia))) is the sign of the end of the cylinder closest to
! eye.  so reversing ie draws the invisible end face.
      ie=-int(sign(1.,xe(ia)))

! Draw INvisible end surface, equal spacing in r^2.
! This is necessary only for vtkfile writing. For accis plotting,
! the invisible face is always overdrawn by the others.
! The eye is at angle running from minus-pi to plus-pi.
      psie=atan2(xe(2),xe(1))
      do i=1,nangle
         t1=2.*pi*(i-1)/nangle-pi
         t2=2.*pi* i   /nangle-pi
         if(abs(mod(2.*pi+t2-psie,2.*pi)-pi).lt.pi/2.)then
            isign=-1
         else
            isign=1
         endif
         itc=int(objg(ofn2)*(t1/(2.*pi)+0.50001))
         do j=1,nr
            r1=sqrt(float(j-1)/nr)
            r2=sqrt(float(j  )/nr)
            do k=1,ncorn
               rface(k,ia)=ie*objg(oradius+ia-1)+objg(ocenter+ia-1)
               ids=mod(2+ia-2,ndims)+1
               rface(k,ids)=objg(oradius+ids-1)*(wc(k)*r1+(1.-wc(k))*r2)
     $              *cos(wp(k)*t1+(1.-wp(k))*t2)+objg(ocenter+ids-1)
               ids=mod(3+ia-2,ndims)+1
               rface(k,ids)=objg(oradius+ids-1)*(wc(k)*r1+(1.-wc(k))*r2)
     $              *sin(wp(k)*t1+(1.-wp(k))*t2)+objg(ocenter+ids-1)
            enddo
            lfw=(mod(i-1+ism2,ism2).eq.ism2/2)
            call facecolor(iosw,2+ie,j,itc+1,iobj,iav,rface,fmin
     $           ,fmax,1,lfw,isign) 
         enddo
      enddo

! Reset to visible end face
      ie=int(sign(1.,xe(ia)))
! Draw curved surface.
      do it=1,nangle
         isign=2*mod(it+itp,2)-1
         i=mod(isign*(it/2)+ito+nangle,nangle)+1
         t1=2.*pi*(i-1)/nangle-pi
         t2=2.*pi* i   /nangle-pi
         itc=int(objg(ofn2)*(t1/(2.*pi)+0.500001))
         do j=1,nz
            z1=(-1.+(j-1)*2./nz)*objg(oradius+ia-1)
            z2=(-1.+ j   *2./nz)*objg(oradius+ia-1)
            do k=1,ncorn
               ids=mod(1+ia-2,ndims)+1
               rface(k,ids)=(wc(k)*z1+(1.-wc(k))*z2)+objg(ocenter+ids-1)
               ids=mod(2+ia-2,ndims)+1
               rface(k,ids)=objg(oradius+ids-1)
     $              *cos(wp(k)*t1+(1.-wp(k))*t2)+objg(ocenter+ids-1)
               ids=mod(3+ia-2,ndims)+1
               rface(k,ids)=objg(oradius+ids-1)
     $              *sin(wp(k)*t1+(1.-wp(k))*t2)+objg(ocenter+ids-1)
            enddo
! Now rface contains the coordinates of the face corners. 
               lfw=(mod(i-1+ism2,ism2).eq.ism2/2)
               call facecolor(iosw,2,itc+1,j,iobj,iav,rface,fmin
     $              ,fmax,2,lfw,isign)               
         enddo
      enddo

! Draw visible end surface, equal spacing in r^2.
! The eye is at angle running from minus-pi to plus-pi.
      psie=atan2(xe(2),xe(1))
      do i=1,nangle
         t1=2.*pi*(i-1)/nangle-pi
         t2=2.*pi* i   /nangle-pi
         if(abs(mod(2.*pi+t2-psie,2.*pi)-pi).lt.pi/2.)then
            isign=-1
         else
            isign=1
         endif
         itc=int(objg(ofn2)*(t1/(2.*pi)+0.50001))
         do j=1,nr
            r1=sqrt(float(j-1)/nr)
            r2=sqrt(float(j  )/nr)
            do k=1,ncorn
               rface(k,ia)=ie*objg(oradius+ia-1)+objg(ocenter+ia-1)
               ids=mod(2+ia-2,ndims)+1
               rface(k,ids)=objg(oradius+ids-1)*(wc(k)*r1+(1.-wc(k))*r2)
     $              *cos(wp(k)*t1+(1.-wp(k))*t2)+objg(ocenter+ids-1)
               ids=mod(3+ia-2,ndims)+1
               rface(k,ids)=objg(oradius+ids-1)*(wc(k)*r1+(1.-wc(k))*r2)
     $              *sin(wp(k)*t1+(1.-wp(k))*t2)+objg(ocenter+ids-1)
            enddo
            lfw=(mod(i-1+ism2,ism2).eq.ism2/2)
            call facecolor(iosw,2+ie,j,itc+1,iobj,iav,rface,fmin
     $           ,fmax,1,lfw,isign) 
         enddo
      enddo

      end
!*********************************************************************
      subroutine cylgplot(iq,objg,iobj,ioswin,fmin,fmax)
! Plot divided faces of non-aligned cylinder object objg.
! iosw determines the nature of the plot
! 0: Color code according to position.
! 1:            according to average flux already in nf_step+1
! 2:            according to average flux-density already in nf_step+2

      include 'ndimsdecl.f'
      include '3dcom.f'
      parameter (nadef=20,nzdef=5,pi=3.141593)
      parameter (ncorn=5)
      real rface(ncorn,ndimsmax)
      real xe(ndimsmax),xcontra(ndimsmax)
      real objg(odata)
      logical lfw
      integer wp(ncorn),wc(ncorn)
      data wp/1,0,0,1,1/wc/0,0,1,1,0/
      ifobj=nf_map(iobj)
      iosw=ioswin
      if(ifobj.eq.0)iosw=0
      iav=nf_address(iq+if_quant(ifobj,if_species)
     $     ,ifobj,nf_step+abs(iosw))
      call minmax(ff_data(iav),nf_posno(1,ifobj),fmin,fmax)
      if(fmin.gt.0.)fmin=0.

      ism1=1
      ism2=1
! Use position arrays compatible with flux array but not too coarse.
      if(objg(ofn2).gt.nadef)then
         nangle=int(objg(ofn2))
      elseif(objg(ofn2).gt.0)then
         nangle=int(objg(ofn2))*nint(float(nadef)/objg(ofn2))
         ism2=nint(float(nadef)/objg(ofn2))
      else
         nangle=nadef
      endif
      if(objg(ofn3).gt.0.)then 
         nz=int(objg(ofn3))
      else
         nz=nzdef
      endif
      nr=1
      if(objg(ofn1).gt.0.)nr=int(objg(ofn1))

      ia=3
! Discover perspective, eye position in world coords:
      call trn32(x,y,z,xe(1),xe(2),xe(3),-1)
      call nxyz2wxyz(xe(1),xe(2),xe(3),xe(1),xe(2),xe(3))
! Transform to unit cylinder coordinates.
      call world3contra(ndims,xe,xe,iobj)
! Flux accumulation is (zero based):
      thetae=atan2(xe(mod(ia+1,3)+1),xe(mod(ia,3)+1))
      thetao=mod((thetae+pi),2.*pi)
      fto=(nangle*(thetao/(2.*pi)+0.5))
      ito=int(fto)
      itp=nint(fto)-ito

! Set to invisible end face and draw for vtkwriting.
      ie=-int(sign(1.,xe(ia)))
      psie=atan2(xe(2),xe(1))
      do i=1,nangle
         t1=2.*pi*(i-1)/nangle-pi
         t2=2.*pi* i   /nangle-pi
         if(abs(mod(2.*pi+t2-psie,2.*pi)-pi).lt.pi/2.)then
            isign=-1
         else
            isign=1
         endif
         itc=int(objg(ofn2)*(t1/(2.*pi)+0.50001))
         do j=1,nr
            r1=sqrt(float(j-1)/nr)
            r2=sqrt(float(j  )/nr)
            do k=1,ncorn
               xcontra(ia)=ie
               ids=mod(2+ia-2,ndims)+1
               xcontra(ids)=(wc(k)*r1+(1.-wc(k))*r2)
     $              *cos(wp(k)*t1+(1.-wp(k))*t2)
               ids=mod(3+ia-2,ndims)+1
               xcontra(ids)=(wc(k)*r1+(1.-wc(k))*r2)
     $              *sin(wp(k)*t1+(1.-wp(k))*t2)
! Transform to world
               call contra3world(ndims,xcontra,xcontra,iobj)
               do m=1,ndims
                  rface(k,m)=xcontra(m)
               enddo
            enddo
            lfw=(mod(i-1+ism2,ism2).eq.ism2/2)
            call facecolor(iosw,2+ie,j,itc+1,iobj,iav,rface,fmin
     $           ,fmax,1,lfw,isign) 
         enddo
      enddo
! Set back to visible face.
      ie=int(sign(1.,xe(ia)))

! Draw curved surface.
      do it=1,nangle
         isign=2*mod(it+itp,2)-1
         i=mod(isign*(it/2)+ito+nangle,nangle)+1
         t1=2.*pi*(i-1)/nangle-pi
         t2=2.*pi* i   /nangle-pi
         itc=int(objg(ofn2)*(t1/(2.*pi)+0.500001))
         do j=1,nz
            z1=(-1.+(j-1)*2./nz)
            z2=(-1.+ j   *2./nz)
            do k=1,ncorn
               ids=mod(1+ia-2,ndims)+1
               xcontra(ids)=(wc(k)*z1+(1.-wc(k))*z2)
               ids=mod(2+ia-2,ndims)+1
               xcontra(ids)=cos(wp(k)*t1+(1.-wp(k))*t2)
               ids=mod(3+ia-2,ndims)+1
               xcontra(ids)=sin(wp(k)*t1+(1.-wp(k))*t2)
! Transform back from unit-cyl to world coordinates
               call contra3world(ndims,xcontra,xcontra,iobj)
               do m=1,ndims
                  rface(k,m)=xcontra(m)
               enddo
            enddo
            lfw=(mod(i-1+ism2,ism2).eq.ism2/2)
            call facecolor(iosw,2,itc+1,j,iobj,iav,rface,fmin
     $           ,fmax,2,lfw,isign)               
         enddo
      enddo

! Draw the visible end face.
! The eye is at angle running from minus-pi to plus-pi.
      psie=atan2(xe(2),xe(1))
      do i=1,nangle
         t1=2.*pi*(i-1)/nangle-pi
         t2=2.*pi* i   /nangle-pi
         if(abs(mod(2.*pi+t2-psie,2.*pi)-pi).lt.pi/2.)then
            isign=-1
         else
            isign=1
         endif
         itc=int(objg(ofn2)*(t1/(2.*pi)+0.50001))
         do j=1,nr
            r1=sqrt(float(j-1)/nr)
            r2=sqrt(float(j  )/nr)
            do k=1,ncorn
               xcontra(ia)=ie
               ids=mod(2+ia-2,ndims)+1
               xcontra(ids)=(wc(k)*r1+(1.-wc(k))*r2)
     $              *cos(wp(k)*t1+(1.-wp(k))*t2)
               ids=mod(3+ia-2,ndims)+1
               xcontra(ids)=(wc(k)*r1+(1.-wc(k))*r2)
     $              *sin(wp(k)*t1+(1.-wp(k))*t2)
! Transform to world
               call contra3world(ndims,xcontra,xcontra,iobj)
               do m=1,ndims
                  rface(k,m)=xcontra(m)
               enddo
            enddo
            lfw=(mod(i-1+ism2,ism2).eq.ism2/2)
            call facecolor(iosw,2+ie,j,itc+1,iobj,iav,rface,fmin
     $           ,fmax,1,lfw,isign) 
         enddo
      enddo

      end
!********************************************************************
!*********************************************************************
      subroutine srvplot(iq,objg,iobj,ioswin,fmin,fmax)
! Plot divided faces of a convex surface of revolution object objg
! iosw determines the nature of the plot. Currently unused.
! 0: Color code according to position.
! 1:            according to average flux already in nf_step+1
! 2:            according to average flux-density already in nf_step+2

      include 'ndimsdecl.f'
      include '3dcom.f'
      parameter (nadef=20,nzdef=5,pi=3.141593)
      parameter (ncorn=5)
      real rface(ncorn,ndimsmax)
      real xe(ndimsmax),xcontra(ndimsmax)
      real objg(odata)
      logical lfw
      integer wp(ncorn),wc(ncorn)
      data wp/1,0,0,1,1/wc/0,0,1,1,0/
      ifobj=nf_map(iobj)
      iosw=ioswin
      if(ifobj.eq.0)iosw=0
      iav=nf_address(iq+if_quant(ifobj,if_species)
     $     ,ifobj,nf_step+abs(iosw))
      call minmax(ff_data(iav),nf_posno(1,ifobj),fmin,fmax)
      if(fmin.gt.0.)fmin=0.

      ism1=1
      ism2=1
! Use angle position arrays compatible with flux array but not too coarse.
! Unlike cylinders, theta is the first dimension. objg(ofn1) is div-No.
      if(objg(ofn1).gt.nadef)then
         nangle=int(objg(ofn1))
      elseif(objg(ofn1).gt.0)then
         nangle=int(objg(ofn1))*nint(float(nadef)/objg(ofn1))
         ism2=nint(float(nadef)/objg(ofn1))
      else
         nangle=nadef
      endif
! The axial coordinate
      ia=3
! Discover perspective, eye position in world coords:
      call trn32(x,y,z,xe(1),xe(2),xe(3),-1)
      call nxyz2wxyz(xe(1),xe(2),xe(3),xe(1),xe(2),xe(3))
! Transform to unit cylinder coordinates.
      call world3contra(ndims,xe,xe,iobj)
! The azimuthal (yx) angle of the eye
      thetae=atan2(xe(mod(ia+1,3)+1),xe(mod(ia,3)+1))
! The opposite theta angle
      thetao=mod((thetae+pi),2.*pi)
! Therefore the angles/indices to start the coloring from:
! Flux accumulation is (zero based):
      fto=(nangle*(thetao/(2.*pi)+0.5))
      ito=int(fto)
      itp=nint(fto)-ito

! Draw in axial order from furthest to nearest.
      if(xe(ia).gt.0.)then
         irz1=0
         irz2=int(objg(onpair))-2
         irzd=1
      else
         irz1=int(objg(onpair))-2
         irz2=0
         irzd=-1
      endif
!      write(*,*)'onpair',objg(onpair),irz1,irz2
! Do over line segments: faces.
      do irz=irz1,irz2,irzd
         rb=objg(opr+irz)
         rt=objg(opr+irz+1)
         zb=objg(opz+irz)
         zt=objg(opz+irz+1)
         nz=int(objg(opdiv+irz))
!         write(*,*)'irz,zb,zt,rb,rt',irz,zb,zt,rb,rt
! Draw single curved surface. Generalization of cylinder round face.
      do it=1,nangle
         isign=2*mod(it+itp,2)-1
         i=mod(isign*(it/2)+ito+nangle,nangle)+1
         t1=2.*pi*(i-1)/nangle-pi
         t2=2.*pi* i   /nangle-pi
         itc=int(objg(ofn1)*(t1/(2.*pi)+0.500001))
! Do over axial facets equally spaced in r^2 running from rb^2 to rt^2
         r2=rb
         r22=rb**2
         z2=zb
         dr2=(rt**2-rb**2)/nz
         dzdr=(zt-zb)/(rt-rb)
         dz=(zt-zb)/nz
         do j=1,nz
            r21=r22
            r22=(r21+dr2)
            r1=r2
            r2=sqrt(r22)
            z1=z2
            if(abs(dr2).ne.0.)then
               z2=zb+ (r2-rb)*dzdr
            else
               z2=zb+ j*dz
            endif
!            write(*,*)irz,j,r1,r2,r21,r22
! Find the positions of the corners
            do k=1,ncorn
               ids=mod(1+ia-2,ndims)+1
               xcontra(ids)=wc(k)*z1+(1-wc(k))*z2
               r=wc(k)*r1+(1-wc(k))*r2
               ids=mod(2+ia-2,ndims)+1
               xcontra(ids)=r*cos(wp(k)*t1+(1-wp(k))*t2)
               ids=mod(3+ia-2,ndims)+1
               xcontra(ids)=r*sin(wp(k)*t1+(1-wp(k))*t2)
! Transform back from unit-cyl to world coordinates
               call contra3world(ndims,xcontra,xcontra,iobj)
               do m=1,ndims
                  rface(k,m)=xcontra(m)
               enddo
            enddo
            lfw=(mod(it-1+ism2,ism2).eq.ism2/2)
!            write(*,*)'Facecolor',irz,it,i,j,itc+1,isign
            call facecolor(iosw,irz,itc+1,j,iobj,iav,rface
!     $           ,fmin,fmax,2,lfw,isign)               
! probably:
     $           ,fmin,fmax,1,lfw,isign)               
         enddo
      enddo
      enddo

      end
!*********************************************************************
      subroutine srvgplot(iq,objg,iobj,ioswin,fmin,fmax)
! Plot divided faces of a general surface of revolution object objg
! iosw determines the nature of the plot
! 0: Color code according to position.
! 1:            according to average flux already in nf_step+1
! 2:            according to average flux-density already in nf_step+2

      include 'ndimsdecl.f'
      include '3dcom.f'
      parameter (nadef=20,nzdef=5,pi=3.141593)
      parameter (ncorn=5)
      real rface(ncorn,ndimsmax)
      real xe(ndimsmax),xcontra(ndimsmax)
      real objg(odata)
      real rp(ovlen),zp(ovlen)
      integer iorder(ovlen)
      logical lfw
      integer wp(ncorn),wc(ncorn)
      data wp/1,0,0,1,1/wc/0,0,1,1,0/

      ifobj=nf_map(iobj)
      iosw=ioswin
      if(ifobj.eq.0)iosw=0
      iav=nf_address(iq+if_quant(ifobj,if_species)
     $     ,ifobj,nf_step+abs(iosw))
      call minmax(ff_data(iav),nf_posno(1,ifobj),fmin,fmax)
      if(fmin.gt.0.)fmin=0.

      ism1=1
      ism2=1
! Use angle position arrays compatible with flux array but not too coarse.
! Unlike cylinders, theta is the first dimension. objg(ofn1) is div-No.
      if(objg(ofn1).gt.nadef)then
         nangle=int(objg(ofn1))
      elseif(objg(ofn1).gt.0)then
         nangle=int(objg(ofn1))*nint(float(nadef)/objg(ofn1))
         ism2=nint(float(nadef)/objg(ofn1))
      else
         nangle=nadef
      endif
! The axial coordinate
      ia=3
! Discover perspective, eye position in world coords:
      call trn32(x,y,z,xe(1),xe(2),xe(3),-1)
      call nxyz2wxyz(xe(1),xe(2),xe(3),xe(1),xe(2),xe(3))
! Transform to unit coordinates.
      call world3contra(ndims,xe,xe,iobj)
! The azimuthal (yx) angle of the eye
      thetae=atan2(xe(mod(ia+1,3)+1),xe(mod(ia,3)+1))
! The opposite theta angle
      thetao=mod((thetae+pi),2.*pi)
! Therefore the angles/indices to start the coloring from:
! Flux accumulation is (zero based):
      fto=(nangle*(thetao/(2.*pi)+0.5))
      ito=int(fto)
      itp=nint(fto)-ito

! Draw single curved surface. Generalization of cylinder round face.
      do it=1,nangle
         isign=2*mod(it+itp,2)-1
         i=mod(isign*(it/2)+ito+nangle,nangle)+1
         t1=2.*pi*(i-1)/nangle-pi
         t2=2.*pi* i   /nangle-pi
         itc=int(objg(ofn1)*(t1/(2.*pi)+0.500001))

! Here's where we decide the order of drawing for rz facets.  
! Project eye onto chosen theta-plane.
         ta=(t1+t2)/2.
         re=(cos(ta)*xe(mod(ia,3)+1)+sin(ta)*xe(mod(ia+1,3)+1))
         ze=xe(ia)
! Subtract it from the contour vertexes, putting them into the rp/zp.
         np=int(objg(onpair))
         do irz=1,np
            rp(irz)=objg(opr+irz-1)-re
            zp(irz)=objg(opz+irz-1)-ze
         enddo
! Sort the order of faces into iorder, so that face=iorder(i)
         call faceorder(np,rp,zp,iorder)
!         write(*,*)np,' iorder',(iorder(irz),irz=1,np-1)
         do i=1,int(objg(onpair)-1)
            irz=iorder(i)
            rb=objg(opr+irz-1)
            rt=objg(opr+irz)
            zb=objg(opz+irz-1)
            zt=objg(opz+irz)
            nz=int(objg(opdiv+irz-1))
!         write(*,*)'irz,zb,zt,rb,rt',irz,zb,zt,rb,rt
! Do over axial facets equally spaced in r^2 running from rb^2 to rt^2
            r2=rb
            r22=rb**2
            z2=zb
            dr2=(rt**2-rb**2)/nz
            dzdr=(zt-zb)/(rt-rb)
            dz=(zt-zb)/nz
            do j=1,nz
               r21=r22
               r22=max((r21+dr2),0.)
               r1=r2
               r2=sqrt(r22)
               z1=z2
               if(abs(dr2).ne.0.)then
                  z2=zb+ (r2-rb)*dzdr
               else
                  z2=zb+ j*dz
               endif
!            write(*,*)irz,j,r1,r2,r21,r22
! Find the positions of the corners
               do k=1,ncorn
                  ids=mod(1+ia-2,ndims)+1
                  xcontra(ids)=wc(k)*z1+(1-wc(k))*z2
                  r=wc(k)*r1+(1-wc(k))*r2
                  ids=mod(2+ia-2,ndims)+1
                  xcontra(ids)=r*cos(wp(k)*t1+(1-wp(k))*t2)
                  ids=mod(3+ia-2,ndims)+1
                  xcontra(ids)=r*sin(wp(k)*t1+(1-wp(k))*t2)
! Transform back from unit-cyl to world coordinates
                  call contra3world(ndims,xcontra,xcontra,iobj)
                  do m=1,ndims
                     rface(k,m)=xcontra(m)
                  enddo
               enddo
               lfw=(mod(it-1+ism2,ism2).eq.ism2/2)
               call facecolor(iosw,irz,itc+1,j,iobj,iav,rface
     $              ,fmin,fmax,1,lfw,isign)               
            enddo
         enddo
      enddo

      end
!********************************************************************
      subroutine pllelplot(iq,objg,iobj,ioswin,fmin,fmax)
! Plot divided faces of parallelopiped object objg.
! iosw determines the nature of the plot
! 0: Color code according to position.
! 1:            according to average flux already in nf_step+1
! 2:            according to average flux-density already in nf_step+2
      include 'ndimsdecl.f'
      include '3dcom.f'
      real objg(odata)
!      character*20 string
      integer iov(ndimsmax)
      real xe(ndimsmax),xn(ndimsmax),objn1(0:ndimsmax-1)
      real xi(ndimsmax)
      parameter (ncorn=5)
      real rface(ncorn,ndimsmax),rfc(ndimsmax)
      integer iof(ncorn,ndimsmax-1)
      data iof/-1,1,1,-1,-1,   -1,-1,1,1,-1/


      ifobj=nf_map(iobj)
      iosw=ioswin
      if(ifobj.eq.0)iosw=0

      iav=nf_address(iq+if_quant(ifobj,if_species)
     $     ,ifobj,nf_step+abs(iosw))
      call minmax(ff_data(iav),nf_posno(1,ifobj),fmin,fmax)
      if(fmin.gt.0.)fmin=0.

! Decide the order of face drawing. Furthest to nearest.
! Get the point and eye position. Make into world units.
      call trn32(x,y,z,xe(1),xe(2),xe(3),-1)
      call nxyz2wxyz(xe(1),xe(2),xe(3),xe(1),xe(2),xe(3))
! Transform eye position into fractional cube position.
      call pllelfrac(xe,xn,iobj)
! The starting fixed point is opposite signs from returned xn.
! So the first three center vectors have values equal to minus
! the sign of xn times the three ovectors. 
      do iv=1,ndims
! Fix zero flux meshes:
         objn1(iv-1)=objg(ofn1+iv-1)
         if(objn1(iv-1).eq.0.)objn1(iv-1)=1.
         iov(iv)=int(-sign(1.,xn(iv)))
! Now if iov(iv) is negative that refers to the first ns_nbins bins.
      enddo
      
      do is=1,2*ndims
! Face Dimension index:
         iv=mod(is-1,ndims)+1
! Face index:
         imin=iv+ndims*(1+iov(iv))/2
         do id=1,ndims
            rfc(id)=objg(ocenter+id-1)
     $           +iov(iv)*objg(ovec+(iv-1)*ndims+id-1)
         enddo
         iov(iv)=-iov(iv)
         i1=mod(1+is-2,ndims)+1
         i2=mod(2+is-2,ndims)+1
         i3=mod(3+is-2,ndims)+1
         do k2=1,int(objn1(i2-1))
!  fs run from 1-N to N-1 as ks run from 1 to N
            f2=2*k2-1.-objn1(i2-1)
            do k3=1,int(objn1(i3-1))
               f3=2*k3-1.-objn1(i3-1)
               do ic=1,ncorn
                  xi(i1)=0
! xi's run from (1-N)+-1 to (N-1)+-1 /N
                  xi(i2)=(f2+iof(ic,1))/objn1(i2-1)
                  xi(i3)=(f3+iof(ic,2))/objn1(i3-1)
                  do id=1,ndims
                     rface(ic,id)=rfc(id)
                     do ivv=1,ndims
                        rface(ic,id)=rface(ic,id)
     $                       +xi(ivv)*objg(ovec+(ivv-1)*ndims+id-1)
                     enddo
                  enddo
               enddo
! Here rface(ic,id) contains the id'th coordinate of ic'th corner.
               call facecolor(iosw,imin,k2,k3,iobj,iav,rface,fmin,fmax
     $              ,i2,.true.,-1)
            enddo
         enddo
      enddo
      end
!*******************************************************************
      subroutine facecolor(iosw,imin,k2,k3,iobj,iav,rface,fmin,fmax
     $     ,i2,lfw,isign)
! Color the facet of the object iobj corresponding to imin,k2,k3 (is,i2)
! If iosw=0 with a color simply to delineate it.
! If iosw.ne.0 with a color corresponding to the data at ff_data(iav+offsets)
! If lfw=.true. and iosw gt 0 then annotate the facet.
! Arguments:
! imin is the face index.
! k2,k3 are the indexes of the facet within the face.
! iobj is the object. 
! iav is the address of the start of the face within ff_data.
! rface contains the 5 corners of the facet.
! fmin and fmax are the limits of the flux range contoured.
! i2 is the dimension-index of the first facet index
! lfw determines whether to write flux on this face
! isign determines the direction of such writing.

      integer iosw,imin,k2,k3,iobj,iav,i2,incorn,idim
      real fmin,fmax
      logical lfw
      include 'ndimsdecl.f'
      include '3dcom.f'
      include 'vtkcom.f'
      parameter (ncorn=5)
      real rface(ncorn,ndimsmax)
      character*20 string

      if(iosw.ne.0.or.vtkflag.eq.1)then
         ifobj=nf_map(iobj)
! Coloring by flux
         ijbin=(k2-1)+nf_dimlens(nf_flux,ifobj,i2)*(k3-1)
     $        +nf_faceind(nf_flux,ifobj,imin)
         iadd=ijbin+iav
         ff=ff_data(iadd)
! This part of the code is called when we want vtk files
! And the coloring and plotting part will be skipped
         if(vtkflag.eq.1)then
            do incorn=1,ncorn-1
              do idim=1,ndims
                 vtkpind=12*vtkindex+(incorn-1)*3+idim
                 if(vtkpind.gt.nvtkindmax)stop 'Over-ran vtpoints'
                 vtkpoints(12*vtkindex+(incorn-1)*3+idim)
     $                =rface(incorn,idim)
              enddo
           enddo
           if(vtkindex.gt.nvtkindmax)stop 'Over-ran vtkindex'
! This is the command that corrupts (even though it does not overrun).
           vtkindex=vtkindex+1
           vtkflx(vtkindex)=ff
           return
         endif
         icolor=int(240*(ff-fmin)/(fmax-fmin))+1
         call gradcolor(icolor)
      else
! Coloring just by position
         itype=int(obj_geom(otype,iobj))-256*(int(obj_geom(otype,iobj))
     $        /256)
         if(itype.eq.1)then
! Treatment for sphere.
            icolor=mod(k3-1+mod(k2,2),7)+7*mod(k2,2)+1
         elseif(itype.eq.3.or.itype.eq.5)then
! Cylinder
            icolor=mod(k2-1+mod(k3,2),7)+7*mod(k3,2)+1
         elseif(itype.eq.6.or.itype.eq.7)then
! Surface of Revolution
            icolor=mod(k3,2)+mod(imin-1+mod(k2,2),7)+7*mod(imin,2)+1
!            write(*,*)imin,k2,k3,' icolor=',icolor
!            write(*,*)(rface(ik,3),ik=1,ncorn)
         else
            icolor=(imin+mod(k2+k3,2)*8)
         endif
         call color(icolor)
      endif
      call poly3line(rface(1,1),rface(1,2),rface(1,3),ncorn)
      call pathfill()
      call color(15)
      if(iosw.gt.0.and. lfw)then
         call iwrite(ijbin,iw,string)
         string(iw+1:iw+1)=' '
         call fwrite(ff,iw,2,string(iw+2:))
         call vec3w((rface(4,1)+rface(2,1))/2.
     $        ,(rface(4,2)+rface(2,2))/2.
     $        ,(rface(4,3)+rface(2,3))/2.,0)
         call charsize(.01,.01)
         call color(15)
!         call drcstr(string)
         if(isign.lt.0)then
            call jdrcstr(string,1.)
         else
            call jdrcstr(string,-1.)
         endif
         call color(15)
         call charsize(.0,.0)
      endif
      end
!*********************************************************************
      subroutine zsort(ngeomobj,zta,index)
! Sort the values zta of length ngeomobj returning sorted index.
      real zta(ngeomobj)
      integer index(ngeomobj)
      do j=2,ngeomobj
         a1=zta(j)
         ia2=index(j)
         do i=j-1,1,-1
            if(zta(i).ge.a1)goto 10
            zta(i+1)=zta(i)
            index(i+1)=index(i)
         enddo
         i=0
 10      zta(i+1)=a1
         index(i+1)=ia2
      enddo
      end
!*******************************************************************
! Plot edges/faces of objects.
      subroutine objplot(iq,rv,cv,ioswin,iomask)
! iosw determines the nature of the plot
! 0: Color code according to position.
! 1:            according to average flux already in nf_step+1
! 2:            according to average flux density in nf_step+2
! Byte 2: 256 plot intersections, 0 don't plot intersections.
! iomask is a mask where non-zero bits mask _out_ objects.
! iq references the quantity, 1 flux, 2-4 momentum etc, to be plotted.
! For species N, iq should be added to (N-1)*kf_quant(nf_obj,ispecies).
! rv gives the Window size, cv the center of the view.
      integer iq,iosw,iomask
      real rv
      include 'ndimsdecl.f'
      include '3dcom.f'
      include 'sectcom.f'
      include 'vtkcom.f'
      real cv(ndimsmax)
      integer index(ngeomobjmax)
      real zta(ngeomobjmax)
      character*100 string
      data iprinting/0/

      ipint=ioswin/256
      iosw=ioswin- ipint*256

!      write(*,*)iosw,ipint,'iosw,iprint'
      irotating=0
      call pltinit(0.,1.,0.,1.)
      call setcube(.2,.2,.2,.5,.4)
      call geteye(x2,y2,z2)
      call trn32(0.,0.,0.,x2,y2,z2,1)
! Color gradient.
      call blueredgreenwhite()
 51   continue
      if(iprinting.ne.0)call pfset(3)
      call pltinit(0.,1.,0.,1.)
      call scale3(cv(1)-rv,cv(1)+rv,cv(2)-rv,cv(2)+rv,cv(3)-rv,cv(3)+rv)
      if(irotating.eq.0)then
         icorner=igetcorner()
         call ax3labels('x','y','z')
      else
         irotating=irotating-1
      endif
!      call cubed(icorner)
      call axproj(icorner)

! Decide the order in which to draw objects, based on the position of
! their centers. 
      do i=1,ngeomobj
         index(i)=i
! Get the position in view coordinates.
         call trn32(obj_geom(ocenter,i),obj_geom(ocenter+1,i),
     $        obj_geom(ocenter+2,i),xt,yt,zta(i),3)
      enddo
      call zsort(ngeomobj,zta,index)
!      write(*,*)(index(k),zta(k),k=1,ngeomobj)

      fmin=0.
      fmax=0.
! Do drawing in order
      do ik=1,ngeomobj
         iobj=index(ik)
         iobjmask=ibits(iomask,iobj-1,1)
         itype=int(obj_geom(otype,iobj))-256*(int(obj_geom(otype,iobj))
     $        /256)
!         write(*,*)'objplotting',ik,iobj,itype,iobjmask
!     $        ,iq,mf_quant(nf_map(iobj))
         if(iobjmask.ne.1)then
!     $        .and. 0.lt.iq.and.iq.le.mf_quant(nf_map(iobj)))then
            if(itype.eq.1.)then
               call sphereplot(iq,obj_geom(1,iobj),iobj,iosw,fmin,fmax)
            elseif(itype.eq.2.)then
               call cubeplot(iq,obj_geom(1,iobj),iobj,iosw,fmin,fmax)
            elseif(itype.eq.3.)then
               call cylplot(iq,obj_geom(1,iobj),iobj,iosw,fmin,fmax)
            elseif(itype.eq.4.)then
               call pllelplot(iq,obj_geom(1,iobj),iobj,iosw,fmin,fmax)
            elseif(itype.eq.5)then
               call cylgplot(iq,obj_geom(1,iobj),iobj,iosw,fmin,fmax)
            elseif(itype.eq.6.or.itype.eq.7)then
! srvplot works for monotonic surface of revolution.
! srvgplot should work for both.
               call srvgplot(iq,obj_geom(1,iobj),iobj,iosw,fmin,fmax)
            endif
         endif
      enddo
! This legend was removed from sphereplot to enable vtkwriting to
! use that routine without additional complexity.
      if(iosw.ne.0)then
         call gradlegend(fmin,fmax,-.35,0.,-.35,.7,-.1,.false.)
         if(abs(iosw).eq.0)then
            string='Position '
         elseif(abs(iosw).eq.1)then
            string='Flux '
         elseif(abs(iosw).eq.2)then
            string='Flux-density '
         endif
!         call iwrite(iosw,iwd,string(lentrim(string)+1:))
         call termchar(string)
         call jdrwstr(.05,.6,string,1.)
      endif
      if(ipint.eq.1)then
! Plot intersections
      call charsize(.008,.008)
!      write(*,*)'Intersections:'
      do i=1,sc_ipt
!      do i=1,300
         call vec3w(x_sc(1,1,i),x_sc(2,1,i),x_sc(3,1,i),0)
         call vec3w(x_sc(1,2,i),x_sc(2,2,i),x_sc(3,2,i),1)
         call iwrite(ibin_sc(i),iw,string)
         call drcstr('!A1!@'//string)
!         write(*,*)i,(x_sc(k,1,i),k=1,3)
      enddo
      call charsize(0.,0.)
      endif

! User interface:
      iprinting=0
      call accisflush()
      call prtend()
      call eye3d(isw)
      call rotatezoom(isw)
      if(isw.eq.ichar('p'))iprinting=mod(iprinting+1,2)
      if(isw.ne.0.and.isw.ne.ichar('q'))goto 51
      end
!*************************************************************
! This subroutine is called to write data in vtkcom common blocks
      subroutine vtkwrite(iq,ioswin,iomask)
      integer iq,iosw,iomask
      include 'ndimsdecl.f'
      include '3dcom.f'
      include 'sectcom.f'
      include 'vtkcom.f'
      integer index(ngeomobjmax)
      real zta(ngeomobjmax)
      iosw=ioswin
! Initializing vtkindex
      vtkindex=0
! Decide the order in which to draw objects, based on the position of
! their centers.
      do i=1,ngeomobj
         index(i)=i
! Get the position in view coordinates.
         call trn32(obj_geom(ocenter,i),obj_geom(ocenter+1,i),
     $        obj_geom(ocenter+2,i),xt,yt,zta(i),3)
      enddo
      call zsort(ngeomobj,zta,index)
!      write(*,*)(index(k),zta(k),k=1,ngeomobj)

      fmin=0.
      fmax=0.
! Do drawing in order
      do ik=1,ngeomobj
         iobj=index(ik)
         iobjmask=ibits(iomask,iobj-1,1)
         itype=int(obj_geom(otype,iobj))-256*(int(obj_geom(otype,iobj))
     $        /256)
!         write(*,*)'objplotting',ik,iobj,itype,iobjmask
         if(iobjmask.ne.1 .and. 0.lt.iq.and.iq.le.mf_quant(iobj))then
!            if(.false.)then
            if(itype.eq.1.)then
               call sphereplot(iq,obj_geom(1,iobj),iobj,iosw,fmin,fmax)
            elseif(itype.eq.2.)then
               call cubeplot(iq,obj_geom(1,iobj),iobj,iosw,fmin,fmax)
            elseif(itype.eq.3.)then
               call cylplot(iq,obj_geom(1,iobj),iobj,iosw,fmin,fmax)
            elseif(itype.eq.4.)then
               call pllelplot(iq,obj_geom(1,iobj),iobj,iosw,fmin,fmax)
            elseif(itype.eq.5)then
               call cylgplot(iq,obj_geom(1,iobj),iobj,iosw,fmin,fmax)
            elseif(itype.eq.6.or.itype.eq.7)then
               call srvgplot(iq,obj_geom(1,iobj),iobj,iosw,fmin,fmax)
            endif
         endif
      enddo
      end
!*************************************************************************
! Given a piecewise linear closed contour in 2-D, a polygon,
! with vertices x(n), y(n), last vertex = first vertex, n-1 sides,
! order the faces of the polygon from farthest to nearest on the basis
! of their view from the origin. A face (side) 1 is farther than face 2
! if some rays from from 1 (to origin) encounter face 2. 
! Faces are "equal" distance if neither 1 nor two is farther, for example
! if they are disjoint and no rays from either encounter the other. 
! **********************************************************************
! This version uses selection sorting to sort the integers referring to
! the face number so that iorder(1:n-1) runs through them in order from
! furthest to nearest. n is the number of vertices, n-1 the number of sides.
! The faces then should be colored in the order iorder(1:n-1).
      subroutine faceorder(n,x,y,iorder)
      integer n
      real x(n),y(n)
      integer iorder(n)
      integer imin

      do i=1,n
         iorder(i)=i
      enddo
! j is the face we are testing against the others.
      do j=1,n-2
         imin=j
! iorder(i) is the face we are testing against.
         do i=j+1,n-1
! If face iorder(i) is further than face imin, choose it the minimum.
            itf=ifacefarther(x(iorder(i)),y(iorder(i))
     $           ,x(iorder(imin)),y(iorder(imin)))
!            write(*,*)j,i,imin,' iorder(i),iorder(imin),itf',iorder(i)
!     $           ,iorder(imin),itf
            if(itf.eq.1)then
! We have found a surface that is further 
!               write(*,*)'Surface',i,' is further than',imin
               imin=i
            endif
         enddo
         if(imin.ne.j)then
! Swap the slots j and imin.
            itemp=iorder(j)
            iorder(j)=iorder(imin)
            iorder(imin)=itemp
         endif
      enddo
      end
!***********************************************************************
      integer function ifacefarther(x1,y1,x2,y2)
! Determine whether face (line segment) 1 is farther than face 2 from
! the origin, if so, return 1, if closer -1, else 0.
! Farther means there is some ray from origin to a point on 1 that 
! intersects 2 first.
      real x1(2),y1(2),x2(2),y2(2)
      real pi,twopi
      parameter (pi=3.1415926,twopi=2*3.1415926)
      ifacefarther=0
      theta11=atan2(y1(1),x1(1))
      theta12=atan2(y1(2),x1(2))
      theta21=atan2(y2(1),x2(1))
      theta22=atan2(y2(2),x2(2))
      if(abs(theta12-theta11).gt.pi.or.abs(theta22-theta21).gt.pi)then
! Adjust cut of theta not to be crossed.
         theta11=mod(theta11+twopi,twopi)
         theta12=mod(theta12+twopi,twopi)
         theta21=mod(theta21+twopi,twopi)
         theta22=mod(theta22+twopi,twopi)
      endif

!      write(*,*)'thetas',theta11,theta12,theta21,theta22
      if((theta21-theta11)*(theta12-theta21).gt.0)then
! Use 21 [21 lies between 11 and 12] or 22 if 21 is on the 1-side.
         dx=(x1(2)-x1(1))
         dy=(y1(2)-y1(1))
         side0=dx*(-y1(1))-dy*(-x1(1))
         side=dx*(y2(1)-y1(1))-dy*(x2(1)-x1(1))
         if(side.eq.0)side=dx*(y2(2)-y1(1))-dy*(x2(2)-x1(1))
      elseif((theta22-theta11)*(theta12-theta22).gt.0)then
! Use 22 [lies between 11 and 12] or else 21 if side=0
         dx=(x1(2)-x1(1))
         dy=(y1(2)-y1(1))
         side0=dx*(-y1(1))-dy*(-x1(1))
         side=dx*(y2(2)-y1(1))-dy*(x2(2)-x1(1))
         if(side.eq.0)side=dx*(y2(1)-y1(1))-dy*(x2(1)-x1(1))
      elseif((theta11-theta21)*(theta22-theta11).gt.0)then
! Use 11 [lies between 21 and 22 ] or 12
         dx=(x2(2)-x2(1))
         dy=(y2(2)-y2(1))
         side0=dx*(-y2(1))-dy*(-x2(1))
         side=dx*(y1(1)-y2(1))-dy*(x1(1)-x2(1))
         if(side.eq.0)side=dx*(y1(2)-y2(1))-dy*(x1(2)-x2(1))
         side=-side
      elseif((theta12-theta21)*(theta22-theta12).gt.0)then
! Use 12 [lies between 21 and 22 ] or 11
         dx=(x2(2)-x2(1))
         dy=(y2(2)-y2(1))
         side0=dx*(-y2(1))-dy*(-x2(1))
         side=dx*(y1(2)-y2(1))-dy*(x1(2)-x2(1))
         if(side.eq.0)side=dx*(y1(1)-y2(1))-dy*(x1(1)-x2(1))
         side=-side
      else
! Must be disjoint
         return
      endif
!      write(*,*)'dx,dy,side0,side,x1,y1',dx,dy,side0,side,x1,y1
      if(side*side0.le.0)then
         ifacefarther=-1
      else
         ifacefarther=1
      endif
      end
!********************************************************************
