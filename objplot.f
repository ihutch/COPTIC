c*********************************************************************
      subroutine sphereplot(objg,iobj,ioswin)
c Plot the sphere objg, number iobj, on an already set up 3D scene.
c iosw determines the nature of the plot
c 0: Color code according to position.
c 1:            according to average flux already in nf_maxstep+1
      include '3dcom.f'
      real objg(odata)
      real xe(ns_ndims)
      parameter (nadef=20,ncosdef=20,pi=3.141593)
      parameter (ncorn=5)
      real rface(ncorn,ns_ndims)
      integer wp(ncorn),wc(ncorn)
c      character*20 string
      logical lfw
      data wp/1,0,0,1,1/wc/0,0,1,1,0/

      ish1=objg(ofn1)
      ish2=objg(ofn2)
      ism1=1
      ism2=1
      ifobj=nf_map(iobj)
      iosw=ioswin
      if(ifobj.eq.0)iosw=0
c No flux for this object
c      write(*,*)'sphereplot',iobj,ifobj,nf_maxsteps,iosw
c     $     ,nf_posno(1,ifobj)
      iav=nf_address(1,ifobj,nf_maxsteps+1)
c Find max and min of flux
      call minmax(ff_data(iav),nf_posno(1,ifobj),fmin,fmax)
      if(fmin.gt.0.)fmin=0.
c Use position arrays compatible with flux array but not too coarse.
      if(ish1.gt.nadef)then 
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
c            write(*,*)c1,c2,ncos
            do k=1,ncorn
               rface(k,1)=objg(oradius)*(wp(k)*cp1+(1.-wp(k))*cp2)
     $              *(wc(k)*s1+(1.-wc(k))*s2)+objg(ocenter)
               rface(k,2)=objg(oradius+1)*(wp(k)*sp1+(1.-wp(k))*sp2)
     $              *(wc(k)*s1+(1.-wc(k))*s2)+objg(ocenter+1)
               rface(k,3)=objg(oradius+2)*(wc(k)*c1+(1.-wc(k))*c2)
     $              +objg(ocenter+2)
            enddo
c            write(*,*)rface
c            write(*,*)ish1,ism1,ism2
c Now rface contains the positions of the face corners.
            k2=nint(ish1*(c1+1.)/2.+.5000)
            lfw=(mod(i+ism2,ism2).eq.ism2/2 .and.
     $           mod(j+ism1,ism1).eq.ism1/2)
c            write(*,*)'Calling facecolor'
            call facecolor(iosw,1,k2,itc+1,iobj,iav,rface,fmin
     $              ,fmax,itc+1,1,lfw)               
         enddo
      enddo
c This legend ought really to account for every object. But at the 
c moment only spheres do it. 
      if(iosw.ne.0)call gradlegend(fmin,fmax,-.4,0.,-.4,.7,-.1,.false.)

      end
c*********************************************************************
      subroutine cubeplot(objg,iobj,ioswin)
c Plot faces of cube object on an already set up 3D scene.
c iosw determines the nature of the plot
c 0: Color code according to position.
c 1:            according to average flux already in nf_maxstep+1
      include '3dcom.f'
      real objg(odata)
c      character*20 string
      integer iov(ns_ndims)
      real xe(ns_ndims)
      parameter (ncorn=5)
      real rface(ncorn,ns_ndims),rfc(ns_ndims)
      integer iof(ncorn,ns_ndims-1)
      data iof/-1,1,1,-1,-1,   -1,-1,1,1,-1/

      ifobj=nf_map(iobj)
      iosw=ioswin
      if(ifobj.eq.0)iosw=0

      iav=nf_address(1,ifobj,nf_maxsteps+1)
      call minmax(ff_data(iav),nf_posno(1,ifobj),fmin,fmax)
      if(fmin.gt.0.)fmin=0.

c Decide the order of face drawing.
c Get the point and eye position. Make into world units.
      call trn32(x,y,z,xe(1),xe(2),xe(3),-1)
      call nxyz2wxyz(xe(1),xe(2),xe(3),xe(1),xe(2),xe(3))
c Transform eye position into fractional cube position.      
c The starting fixed point is opposite signs from returned xn.
      do iv=1,ns_ndims
c Fix zero flux meshes:
         if(objg(ofn1+iv-1).eq.0.)objg(ofn1+iv-1)=1.
         iov(iv)=-sign(1.,xe(iv)-objg(ocenter+iv-1))
      enddo
c Now ordered.
c      write(*,*)(objg(k),k=1,4*ns_ndims),objg(ocenter),objg(oradius)
      do is=1,2*ns_ndims
         iv=mod(is-1,ns_ndims)+1
c Face index:
         imin=iv+ns_ndims*(1+iov(iv))/2
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
c  fs run from 1-N to N-1 as ks run from 1 to N
            f2=2*k2-1.-objg(ofn1+i2-1)
            do k3=1,objg(ofn1+i3-1)
               f3=2*k3-1.-objg(ofn1+i3-1)
               do ic=1,ncorn
c Set the corner offsets for this face, is.
                  rface(ic,i1)=rfc(i1)
                  rface(ic,i2)=rfc(i2)+(f2+iof(ic,1))
     $                 *objg(oradius+i2-1)/objg(ofn1+i2-1)
                  rface(ic,i3)=rfc(i3)+(f3+iof(ic,2))
     $                 *objg(oradius+i3-1)/objg(ofn1+i3-1)
               enddo
               call facecolor(iosw,imin,k2,k3,iobj,iav,rface,fmin,fmax
     $              ,is,i2,.true.)
           enddo
         enddo
      enddo
      end
c*********************************************************************
      subroutine cylplot(objg,iobj,ioswin)
c Plot divided faces of cylinder object objg.
c iosw determines the nature of the plot
c 0: Color code according to position.
c 1:            according to average flux already in nf_maxstep+1
      include '3dcom.f'
      parameter (nadef=20,nzdef=5,pi=3.141593)
      parameter (ncorn=5)
      real rface(ncorn,ns_ndims)
      real xe(ns_ndims)
      real objg(odata)
      logical lfw
      integer wp(ncorn),wc(ncorn)
      data wp/1,0,0,1,1/wc/0,0,1,1,0/

      ifobj=nf_map(iobj)
      iosw=ioswin
      if(ifobj.eq.0)iosw=0
      iav=nf_address(1,ifobj,nf_maxsteps+1)
      call minmax(ff_data(iav),nf_posno(1,ifobj),fmin,fmax)
      if(fmin.gt.0.)fmin=0.

      ism1=0
      ism2=0
c Use position arrays compatible with flux array but not too coarse.
      if(objg(ofn2).gt.nadef)then
         nangle=objg(ofn2)
      elseif(objg(ofn2).gt.0)then
         nangle=objg(ofn2)*nint(float(nadef)/objg(ofn2))
         ism2=nint(float(nadef)/objg(ofn2))
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
               lfw=(mod(i-1+ism2,ism2).eq.ism2/2)
               call facecolor(iosw,2,itc+1,j,iobj,iav,rface,fmin
     $              ,fmax,itc+1,2,lfw)               
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
               lfw=(mod(i-1+ism2,ism2).eq.ism2/2)
               call facecolor(iosw,2+ie,itc+1,j,iobj,iav,rface,fmin
     $              ,fmax,itc+1,2,lfw)               
         enddo
      enddo

      end
c********************************************************************
      subroutine pllelplot(objg,iobj,ioswin)
c Plot divided faces of parallelopiped object objg.
c iosw determines the nature of the plot
c 0: Color code according to position.
c 1:            according to average flux already in nf_maxstep+1
      include '3dcom.f'
      real objg(odata)
c      character*20 string
      integer iov(ns_ndims)
      real xe(ns_ndims),xn(ns_ndims)
      real xi(pp_ndims)
      parameter (ncorn=5)
      real rface(ncorn,pp_ndims),rfc(pp_ndims)
      integer iof(ncorn,pp_ndims-1)
      data iof/-1,1,1,-1,-1,   -1,-1,1,1,-1/


      ifobj=nf_map(iobj)
      iosw=ioswin
      if(ifobj.eq.0)iosw=0

      iav=nf_address(1,ifobj,nf_maxsteps+1)
      call minmax(ff_data(iav),nf_posno(1,ifobj),fmin,fmax)
      if(fmin.gt.0.)fmin=0.

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
c Fix zero flux meshes:
         if(objg(ofn1+iv-1).eq.0.)objg(ofn1+iv-1)=1.
         iov(iv)=-sign(1.,xn(iv))
c Now if iov(iv) is negative that refers to the first ns_nbins bins.
      enddo
      
      do is=1,2*pp_ndims
c Face Dimension index:
         iv=mod(is-1,pp_ndims)+1
c Face index:
         imin=iv+ns_ndims*(1+iov(iv))/2
         do id=1,pp_ndims
            rfc(id)=objg(pp_orig+id-1)
     $           +iov(iv)*objg(pp_vec+(iv-1)*pp_ndims+id-1)
         enddo
         iov(iv)=-iov(iv)
         i1=mod(1+is-2,pp_ndims)+1
         i2=mod(2+is-2,pp_ndims)+1
         i3=mod(3+is-2,pp_ndims)+1
         do k2=1,objg(ofn1+i2-1)
c  fs run from 1-N to N-1 as ks run from 1 to N
            f2=2*k2-1.-objg(ofn1+i2-1)
            do k3=1,objg(ofn1+i3-1)
               f3=2*k3-1.-objg(ofn1+i3-1)
               do ic=1,ncorn
                  xi(i1)=0
c xi's run from (1-N)+-1 to (N-1)+-1 /N
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
               call facecolor(iosw,imin,k2,k3,iobj,iav,rface,fmin,fmax
     $              ,is,i2,.true.)
            enddo
         enddo
      enddo
      end
c*******************************************************************
      subroutine facecolor(iosw,imin,k2,k3,iobj,iav,rface,fmin,fmax,is
     $     ,i2,lfw)
c Color the face of the object iobj corresponding to imin,k2,k3 (is,i2)
c If iosw=0 with a color simply to delineate it.
c If iosw=1 with a color corresponding to the average flux.
c If lfw=.true. then annotate the face.
      integer iosw,imin,k2,k3,iobj,iav,is,i2
      real fmin,fmax
      logical lfw
      include '3dcom.f'
      parameter (ncorn=5)
      real rface(ncorn,pp_ndims)
      character*20 string

      if(iosw.eq.1)then
         ifobj=nf_map(iobj)
c Coloring by flux
         ijbin=(k2-1)+nf_dimlens(nf_flux,ifobj,i2)*(k3-1)
     $        +nf_faceind(nf_flux,ifobj,imin)
         iadd=ijbin+iav
         ff=ff_data(iadd)
         icolor=int(240*(ff-fmin)/(fmax-fmin))+1
c         write(*,*)k2,k3,is,i2,imin,ifobj,nf_faceind(nf_flux,ifobj,imin)
c     $        ,ijbin,ff,icolor
         call gradcolor(icolor)
      else
c Coloring just by position
         itype=obj_geom(otype,iobj)-256*(int(obj_geom(otype,iobj))/256)
         if(itype.eq.1)then
c Treatment for spheres.
            icolor=mod(k3-1+mod(k2,2),7)+7*mod(k2,2)+1
c            write(*,*)k2,k3,icolor
         elseif(itype.eq.3)then
c Cylinders
            icolor=mod(k2-1+mod(k3,2),7)+7*mod(k3,2)+1
         else
            write(*,*)'Hacked facecolor'
            icolor=(is+mod(k2+k3,2)*8)
         endif
         call color(icolor)
      endif
      call poly3line(rface(1,1),rface(1,2),rface(1,3),ncorn)
      call pathfill()
      call color(15)
      if(iosw.eq.1.and. lfw)then
         call iwrite(ijbin,iw,string)
         string(iw+1:iw+1)=' '
         call fwrite(ff,iw,2,string(iw+2:iw+2))
         call vec3w((rface(4,1)+rface(2,1))/2.
     $        ,(rface(4,2)+rface(2,2))/2.
     $        ,(rface(4,3)+rface(2,3))/2.,0)
         call charsize(.01,.01)
         call color(15)
         call drcstr(string)
         call color(15)
      endif

      end
c*******************************************************************
c Plot edges/faces of objects. Window size rs.
      subroutine objplot(ndims,rs,ioswin)
c iosw determines the nature of the plot
c 0: Color code according to position.
c 1:            according to average flux already in nf_maxstep+1
c Bit 2: 0 plot intersections [0,1], 1 don't plot intersections [2,3].
      integer ndims,iosw
      include '3dcom.f'
      include 'sectcom.f'
      integer index(ngeomobjmax)
      real zta(ngeomobjmax)
      character*10 string
      data iprinting/0/

      iosw=ioswin-(ioswin/2)*2
      ipint=ioswin-(ioswin/4)*4-iosw

      irotating=0
      call pltinit(0.,1.,0.,1.)
      call setcube(.2,.2,.2,.5,.4)
      call geteye(x2,y2,z2)
      call trn32(0.,0.,0.,x2,y2,z2,1)
c Color gradient.
      call blueredgreenwhite()
 51   continue
      call pltinit(0.,1.,0.,1.)
      call scale3(-rs,rs,-rs,rs,-rs,rs)
      if(iprinting.ne.0)call pfset(3)
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
         itype=obj_geom(otype,iobj)-256*(int(obj_geom(otype,iobj))/256)
c         write(*,*)'objplotting',ik,iobj,itype
         if(itype.eq.1.)then
            call sphereplot(obj_geom(1,iobj),iobj,iosw)
         elseif(itype.eq.2.)then
            call cubeplot(obj_geom(1,iobj),iobj,iosw)
         elseif(itype.eq.3.)then
            call cylplot(obj_geom(1,iobj),iobj,iosw)
         elseif(itype.eq.4.)then
            call pllelplot(obj_geom(1,iobj),iobj,iosw)
         endif
      enddo
      if(ipint.eq.0)then
c Plot intersections
      call charsize(.008,.008)
c      write(*,*)'Intersections:'
      do i=1,sc_ipt
c      do i=1,300
         call vec3w(x_sc(1,1,i),x_sc(2,1,i),x_sc(3,1,i),0)
         call vec3w(x_sc(1,2,i),x_sc(2,2,i),x_sc(3,2,i),1)
         call iwrite(ibin_sc(i),iw,string)
         call drcstr('!A1!@'//string)
c         write(*,*)i,(x_sc(k,1,i),k=1,3)
      enddo
      call charsize(0.,0.)
      endif

c User interface:
      iprinting=0
 56   call eye3d(isw)
      call rotatezoom(isw)
      if(isw.eq.ichar('p'))iprinting=mod(iprinting+1,2)
      if(isw.ne.0.and.isw.ne.ichar('q'))goto 51
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
